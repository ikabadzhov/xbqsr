#include <cassert>
#include "CompressedInputStream.h"

#include <array>
#include <cstdint>
#include "../../zlib.h"
#include <stdexcept>
#include <string>
#include <tbb/parallel_for.h>

// ---- GZIP/BGZF header probes ------------------------------------------------

bool CompressedInputStream::isValidFile(PlainInputStream& stream) {
    // Read the first 18 bytes (enough for GZIP + BGZF extra header)
    std::string hdr(18, '\0');
    (void)stream.read(hdr, 0, 18);
    bool ok = isValidBlockHeader(hdr);
    stream.seekg(0); // reset to beginning
    return ok;
}

CompressedInputStream::CompressedInputStream(const std::string &input): path(input), fd_(-1),
    comp_off_(0),
    decomp_off_(0),
    stash_pos_(0){
    Options opts;
    stash_cap_ = opts.stash_cap_bytes;
    seen_eof_block_ = false;
    fd_ = ::open(path.c_str(), O_RDONLY);
    if (fd_ < 0) throw std::runtime_error("open() failed for input BGZF");

    // zlib stream: enable gzip wrapper (16 + MAX_WBITS)
    std::memset(&zstrm_, 0, sizeof(zstrm_));
    if (inflateInit2(&zstrm_, 16 + MAX_WBITS) != Z_OK) {
        ::close(fd_); fd_ = -1;
        throw std::runtime_error("zlib inflateInit2 failed");
    }
#ifdef USE_LIBDEFLATE
        decomp_ = libdeflate_alloc_decompressor();
        if (!decomp_) {
            inflateEnd(&zstrm_);
            ::close(fd_); fd_ = -1;
            throw std::runtime_error("libdeflate_alloc_decompressor failed");
        }
#endif
    //if (opts.stash_reserve_bytes)
    stash_.reserve(opts.stash_reserve_bytes);

    // temp buffers (grown on demand)
    //tmp_decomp_.resize(opts.stash_cap_bytes / BGZF_BLOCK_SIZE, std::vector<uint8_t>(BGZF_BLOCK_SIZE)); // BGZF blocks decompress to ≤ 64 KiB
    //tmp_comp_.resize(opts.stash_cap_bytes / BGZF_BLOCK_SIZE, std::vector<uint8_t>(BGZF_BLOCK_SIZE)); // BGZF blocks compress to ≤ 64 KiB
    //tmp_decomp_.resize(opts.stash_cap_bytes); // BGZF blocks decompress to ≤ 64 KiB
    //tmp_comp_.resize(opts.stash_cap_bytes);
    comp_pool.resize(opts.stash_cap_bytes / BGZF_BLOCK_SIZE, std::array<uint8_t, BGZF_BLOCK_SIZE>{});
}

static inline unsigned U8(char c) { return static_cast<unsigned char>(c); }

bool CompressedInputStream::isValidBlockHeader(std::string& b) {
    if (b.size() < 18) return false;

    // GZIP magic + method
    if (U8(b[0]) != 0x1f || U8(b[1]) != 0x8b || U8(b[2]) != 0x08) return false;

    // FLG has FEXTRA bit set?
    const unsigned flg = U8(b[3]);
    if ((flg & 0x04) == 0) return false;

    // XLEN (LE) at bytes 10..11 should be 6 for BGZF
    const unsigned xlen = U8(b[10]) | (U8(b[11]) << 8);
    if (xlen < 6) return false;

    // BGZF subfield "BC" at 12..13, SLEN = 2 at 14..15
    if (U8(b[12]) != 'B' || U8(b[13]) != 'C') return false;
    const unsigned slen = U8(b[14]) | (U8(b[15]) << 8);
    if (slen != 2) return false;

    // If you want, you could read BSIZE at 16..17 here
    return true;
}

// ---- Core I/O ---------------------------------------------------------------
struct BlockPlan {
    uint64_t off;         // compressed file offset
    uint32_t comp_size;   // BSIZE + 1
    uint32_t decomp_size; // ISIZE (0 => EOF marker)
    size_t   dest_off;    // offset (relative to stash tail) where this block will write
};

ssize_t CompressedInputStream::read(std::string& buffer, ssize_t offset, ssize_t length)
{
    // ---- validate args ----
    if (length < 0 || offset < 0)
        throw std::out_of_range("CompressedInputStream::read: negative offset/length");
    const size_t uoff = static_cast<size_t>(offset);
    const size_t ulen = static_cast<size_t>(length);
    if (uoff > buffer.size() || uoff + ulen > buffer.size())
        throw std::out_of_range("CompressedInputStream::read: buffer too small for requested span");

    // Safety: comp_pool must have been resized (not reserved) in ctor
    if (comp_pool.empty())
        throw std::logic_error("CompressedInputStream: comp_pool not initialized (size==0)");

    const size_t MAX_BATCH_BLOCKS = comp_pool.size();                   // capacity of our scratch
    const size_t MAX_BATCH_BYTES  = MAX_BATCH_BLOCKS * BGZF_BLOCK_SIZE; // decompressed upper-bound per batch

    size_t want  = ulen;
    size_t total = 0;

    while (want > 0) {
        // ---- 1) serve from stash if any available ----
        const size_t avail = (stash_pos_ < stash_.size()) ? (stash_.size() - stash_pos_) : 0;
        if (avail) {
            const size_t n = std::min(want, avail);
            std::memcpy(&buffer[uoff + total], stash_.data() + stash_pos_, n);
            stash_pos_ += n;
            total      += n;
            want       -= n;
            continue;
        }

        // ---- 2) need to produce more bytes ----
        compact_stash_();                 // drop consumed prefix, keep capacity
        if (seen_eof_block_) break;       // nothing more to produce

        const size_t stash_room = (stash_cap_ > stash_.size()) ? (stash_cap_ - stash_.size()) : 0;
        if (stash_room == 0) break;

        // cover at least one block, at most free stash, and at most one batch capacity
        const size_t need_bytes = std::min(
            std::max(want, static_cast<size_t>(BGZF_BLOCK_SIZE)),
            std::min(stash_room, MAX_BATCH_BYTES)
        );

        // -------- PLAN (do NOT advance comp_off_ here) --------
        struct BlockPlan {
            uint64_t off;         // compressed file offset
            uint32_t comp_size;   // BSIZE + 1
            uint32_t decomp_size; // ISIZE (0 => EOF marker)
            size_t   dest_off;    // offset within this batch's stash tail
        };

        std::vector<BlockPlan> plan;
        plan.reserve(std::min(MAX_BATCH_BLOCKS, (need_bytes + BGZF_BLOCK_SIZE - 1) / BGZF_BLOCK_SIZE));

        size_t   batch_bytes = 0;  // sum of ISIZE for non-EOF blocks
        uint64_t off_plan    = comp_off_;
        bool     saw_eof     = false;
        size_t   planned_cnt = 0;

        while (!saw_eof && batch_bytes < need_bytes && planned_cnt < MAX_BATCH_BLOCKS) {
            BGZFBlockInfo info = probe_bgzf_block_size_(off_plan);
            if (info.comp_size == 0) { saw_eof = true; break; }

            const uint32_t comp_sz   = static_cast<uint32_t>(info.comp_size);
            const uint32_t decomp_sz = static_cast<uint32_t>(info.decomp_size);

            if (decomp_sz == 0) { // EOF marker: no output, but must consume later
                plan.push_back(BlockPlan{ off_plan, comp_sz, 0u, batch_bytes });
                off_plan += comp_sz;
                saw_eof = true;
                break;
            }

            // capacity guards BEFORE including the block
            if (batch_bytes + decomp_sz > stash_room)   break; // stash cap
            if (batch_bytes + decomp_sz > MAX_BATCH_BYTES) break; // batch cap

            plan.push_back(BlockPlan{ off_plan, comp_sz, decomp_sz, batch_bytes });
            off_plan    += comp_sz;
            batch_bytes += decomp_sz;
            ++planned_cnt;
        }

        if (plan.empty()) {
            // Nothing we can plan (either EOF or no room)
            seen_eof_block_ = seen_eof_block_ || saw_eof;
            break;
        }

        // Reserve exactly the decompressed space we will fill this round
        const size_t old_size = stash_.size();
        if (old_size + batch_bytes > stash_cap_) {
            // With guards above this should never happen
            throw std::runtime_error("Planning error: batch would exceed stash capacity");
        }
        stash_.resize(old_size + batch_bytes);
        uint8_t* const stash_base = stash_.data() + old_size;

        // -------- READ COMPRESSED (serial) --------
        if (plan.size() > MAX_BATCH_BLOCKS)
            throw std::runtime_error("plan exceeds comp_pool capacity");
        for (size_t i = 0; i < plan.size(); ++i) {
            const auto& j = plan[i];
            const ssize_t rn = ::pread(fd_, comp_pool[i].data(),
                                       j.comp_size, static_cast<off_t>(j.off));
            if (rn != static_cast<ssize_t>(j.comp_size))
                throw std::runtime_error("CompressedInputStream: short read on BGZF block");
        }

        // -------- PARALLEL DECOMPRESS: direct write into stash --------
        std::vector<size_t> produced(plan.size(), 0);
        std::vector<int>    status(plan.size(), 0); // 0 ok, 1 eof, -1 error

        tbb::parallel_for(size_t(0), plan.size(), [&](size_t i){
            const auto& j = plan[i];

            if (j.decomp_size == 0) { // EOF marker — no output
                produced[i] = 0;
                status[i]   = 1;
                return;
            }

        #ifdef USE_LIBDEFLATE
            size_t out_n = 0;
            // (optional) thread-local decompressor to avoid alloc/free overhead
            thread_local libdeflate_decompressor* tls_d = []{
                auto* p = libdeflate_alloc_decompressor();
                if (!p) std::abort();
                return p;
            }();
            auto st = libdeflate_gzip_decompress(
                tls_d,
                comp_pool[i].data(), j.comp_size,
                stash_base + j.dest_off, j.decomp_size,
                &out_n
            );
            if (st != LIBDEFLATE_SUCCESS || out_n != j.decomp_size) {
                status[i] = -1; return;
            }
        #else
            z_stream zs{};
            if (inflateInit2(&zs, 16 + MAX_WBITS) != Z_OK) { status[i] = -1; return; }
            zs.next_in   = const_cast<Bytef*>(reinterpret_cast<const Bytef*>(comp_pool[i].data()));
            zs.avail_in  = (uInt)j.comp_size;
            zs.next_out  = reinterpret_cast<Bytef*>(stash_base + j.dest_off);
            zs.avail_out = (uInt)j.decomp_size; // exact ISIZE bound
            const int ret = inflate(&zs, Z_FINISH);
            const size_t out_n = j.decomp_size - zs.avail_out;
            inflateEnd(&zs);
            if (ret != Z_STREAM_END || out_n != j.decomp_size) { status[i] = -1; return; }
        #endif

            produced[i] = j.decomp_size;
            status[i]   = 0;
        });

        // -------- SEQ VALIDATION + ADVANCE STREAM (no gaps) --------
        size_t appended = 0;
        for (size_t i = 0; i < plan.size(); ++i) {
            const auto& j = plan[i];

            if (j.decomp_size == 0) {
                // EOF marker: consume it in the compressed stream
                seen_eof_block_ = true;
                comp_off_ = j.off + j.comp_size;
                continue;
            }

            if (status[i] != 0 || produced[i] != j.decomp_size) {
                // Roll back unfilled tail beyond 'appended' and fail fast
                stash_.resize(old_size + appended);
                throw std::runtime_error("BGZF decompress failed or ISIZE mismatch (parallel direct-write)");
            }

            appended += j.decomp_size;
            // Block fully consumed in order; now advance compressed stream offset
            comp_off_ = j.off + j.comp_size;

            if (old_size + appended == stash_cap_) break; // stash is full
        }

        // After production: if still nothing available, break out
        if (stash_pos_ >= stash_.size()) break;
    }

    return static_cast<ssize_t>(total);
}

void CompressedInputStream::seekg(ssize_t pos) {
    // Use z_off_t for portability (zlib can be built with 64-bit z_off_t)
    z_off_t r = gzseek(fileDescriptor, static_cast<z_off_t>(pos), SEEK_SET);
    if (r < 0) {
        int err = Z_OK;
        const char* msg = gzerror(fileDescriptor, &err);
        throw std::runtime_error(std::string("gzseek failed: ") + (msg ? msg : "unknown"));
    }
}

const std::string& CompressedInputStream::getFilePath() { return path; }

// Return int so -1 is representable.
int CompressedInputStream::peek() {
    int c = gzgetc(fileDescriptor);
    if (c != -1) {
        // push back exactly what we read
        if (gzungetc(c, fileDescriptor) == -1) {
            int err = Z_OK;
            const char* msg = gzerror(fileDescriptor, &err);
            throw std::runtime_error(std::string("gzungetc failed: ") + (msg ? msg : "unknown"));
        }
    }
    return c; // -1 on EOF/error, else 0..255
}

bool CompressedInputStream::eof() {
    return gzeof(fileDescriptor) != 0;
}

bool CompressedInputStream::bad() {
    int err = Z_OK;
    (void)gzerror(fileDescriptor, &err);
    // Treat anything not OK or stream-end as bad.
    return err != Z_OK && err != Z_STREAM_END;
}

bool CompressedInputStream::good() {
    int err = Z_OK;
    (void)gzerror(fileDescriptor, &err);
    return (err == Z_OK) && !eof();
}

// Make close idempotent; destructor calls close().
void CompressedInputStream::close() {
    if (fileDescriptor) {
        gzclose(fileDescriptor);
        fileDescriptor = nullptr;
    }
}


// Probe BGZF block header/footer to get both compressed and decompressed sizes.
BGZFBlockInfo CompressedInputStream::probe_bgzf_block_size_(uint64_t off) {
    // --- Step 1: read minimal header ---
    uint8_t hdr[18];
    ssize_t n = ::pread(fd_, hdr, sizeof(hdr), static_cast<off_t>(off));
    if (n == 0) {
        // True EOF (no BGZF EOF marker)
        return {0, 0};
    }
    if (n < 12) throw std::runtime_error("Short GZIP header");
    if (hdr[0] != 0x1f || hdr[1] != 0x8b) throw std::runtime_error("Not a gzip stream");
    if ((hdr[3] & 0x04) == 0)             throw std::runtime_error("BGZF missing FEXTRA");

    const uint16_t xlen = le16_(hdr + 10);

    // --- Step 2: read FEXTRA to find BSIZE ---
    std::vector<uint8_t> extra(xlen);
    if (::pread(fd_, extra.data(), xlen, static_cast<off_t>(off + 12)) != static_cast<ssize_t>(xlen)) {
        throw std::runtime_error("Short GZIP extra");
    }

    std::size_t p = 0;
    uint16_t bsize = 0;
    while (p + 4 <= extra.size()) {
        const uint8_t  si1  = extra[p];
        const uint8_t  si2  = extra[p + 1];
        const uint16_t slen = le16_(extra.data() + p + 2);
        p += 4;
        if (p + slen > extra.size()) break;
        if (si1 == 'B' && si2 == 'C' && slen >= 2) {
            bsize = le16_(extra.data() + p); // BSIZE
            break;
        }
        p += slen;
    }
    if (bsize == 0) throw std::runtime_error("BGZF 'BC' subfield not found");

    // --- Step 3: compressed size = BSIZE + 1 ---
    const std::size_t comp_size = static_cast<std::size_t>(bsize) + 1u;

    // --- Step 4: read last 4 bytes for ISIZE (uncompressed size) ---
    uint8_t isize_bytes[4];
    off_t footer_off = static_cast<off_t>(off + comp_size - 4);
    if (::pread(fd_, isize_bytes, 4, footer_off) != 4) {
        throw std::runtime_error("Short GZIP footer");
    }
    const uint32_t isize = le32_(isize_bytes);

    // ISIZE is modulo 2^32, but BGZF ensures it's ≤ 64KB.
    if (isize > 65536) {
        throw std::runtime_error("Invalid BGZF ISIZE > 64KB");
    }

    return {comp_size, static_cast<std::size_t>(isize)};
}


bool CompressedInputStream::inflate_one_bgzf_block_() {
    auto [comp_size, decomp_size ] = probe_bgzf_block_size_(comp_off_);
    if (comp_size == 0) {
        // Reached physical EOF (no proper BGZF EOF marker)
        return false;
    }

    // Read the compressed block
    //tmp_comp_.resize(comp_size);
    const ssize_t rn = ::pread(fd_, tmp_comp_.data(), comp_size, static_cast<off_t>(comp_off_));
    if (rn != static_cast<ssize_t>(comp_size)) {
        throw std::runtime_error("Short read on BGZF block");
    }
    comp_off_ += comp_size;

    // Decompress to tmp_decomp_ (resize to max; shrink after)
    //tmp_decomp_.resize(65536);

#ifdef USE_LIBDEFLATE
        std::size_t out_n = 0;
        const auto st = libdeflate_gzip_decompress(
            decomp_,
            tmp_comp_.data(), tmp_comp_.size(),
            tmp_decomp_.data(), tmp_decomp_.size(),
            &out_n
        );
        if (st != LIBDEFLATE_SUCCESS) {
            throw std::runtime_error("libdeflate: decompress failed");
        }
#else
    // Reset zlib stream for a new gzip member
    if (inflateReset2(&zstrm_, 16 + MAX_WBITS) != Z_OK) {
        throw std::runtime_error("zlib inflateReset2 failed");
    }
    zstrm_.next_in   = const_cast<Bytef*>(reinterpret_cast<const Bytef*>(tmp_comp_.data()));
    zstrm_.avail_in  = static_cast<uInt>(tmp_comp_.size());
    zstrm_.next_out  = reinterpret_cast<Bytef*>(tmp_decomp_.data());
    zstrm_.avail_out = static_cast<uInt>(tmp_decomp_.size());

    int ret = inflate(&zstrm_, Z_FINISH);
    if (ret != Z_STREAM_END && ret != Z_OK) {
        throw std::runtime_error("zlib inflate failed");
    }
    const std::size_t out_n = tmp_decomp_.size() - zstrm_.avail_out;
#endif

    if (out_n == 0) {
        // BGZF EOF marker block (valid, zero-byte payload)
        //tmp_decomp_.clear();
        return false;
    }

    // Append to stash (compact first if needed)
    compact_stash_();
    const std::size_t old = stash_.size();
    stash_.resize(old + out_n);
    std::memcpy(stash_.data() + old, tmp_decomp_.data(), out_n);
    return true;
}

bool CompressedInputStream::inflate_bgzf_block_(BlockMeta& block) {
    // Read the compressed block
    //tmp_comp_.resize(comp_size);
    const ssize_t rn = ::pread(fd_, tmp_comp_.data() + block.comp_off, block.comp_size, static_cast<off_t>(comp_off_));
    if (rn != static_cast<ssize_t>(block.comp_size)) {
        throw std::runtime_error("Short read on BGZF block");
    }

#ifdef USE_LIBDEFLATE
        std::size_t out_n = 0;
        const auto st = libdeflate_gzip_decompress(
            decomp_,
            tmp_comp_.data() + block.comp_off, block.comp_size,
            tmp_decomp_.data() + block.decomp_off, block.decomp_size,
            &out_n
        );
        if (st != LIBDEFLATE_SUCCESS) {
            throw std::runtime_error("libdeflate: decompress failed");
        }
#else
    // Reset zlib stream for a new gzip member
    if (inflateReset2(&zstrm_, 16 + MAX_WBITS) != Z_OK) {
        throw std::runtime_error("zlib inflateReset2 failed");
    }
    zstrm_.next_in   = const_cast<Bytef*>(reinterpret_cast<const Bytef*>(tmp_comp_.data()));
    zstrm_.avail_in  = static_cast<uInt>(tmp_comp_.size());
    zstrm_.next_out  = reinterpret_cast<Bytef*>(tmp_decomp_.data());
    zstrm_.avail_out = static_cast<uInt>(tmp_decomp_.size());

    int ret = inflate(&zstrm_, Z_FINISH);
    if (ret != Z_STREAM_END && ret != Z_OK) {
        throw std::runtime_error("zlib inflate failed");
    }
    const std::size_t out_n = tmp_decomp_.size() - zstrm_.avail_out;
#endif

    if (out_n == 0) {
        // BGZF EOF marker block (valid, zero-byte payload)
        //tmp_decomp_.clear();
        return false;
    }

    // Append to stash (compact first if needed)
    //compact_stash_();
    const std::size_t old = stash_.size();
    stash_.resize(old + out_n);
    std::memcpy(stash_.data() + old, tmp_decomp_.data(), out_n);
    return true;
}

bool CompressedInputStream::get_bgzf_blocks_indexes(std::vector<BlockMeta> &blocks) {
    auto [comp_size, decomp_size] = probe_bgzf_block_size_(comp_off_);
    if (comp_size == 0) {
        return false;
    }

    blocks.emplace_back(BlockMeta{comp_off_, decomp_off_, comp_size, decomp_size});
    comp_off_ += comp_size;
    decomp_off_ += decomp_size;
    return true;
}

void CompressedInputStream::compact_stash_() {
    if (stash_pos_ == 0) return;
    const std::size_t tail = stash_.size() - stash_pos_;
    if (tail) {
        std::memmove(stash_.data(), stash_.data() + stash_pos_, tail);
        stash_.resize(tail);
    } else {
        stash_.clear();
    }
    stash_pos_ = 0;
}

CompressedInputStream::~CompressedInputStream() {
#ifdef USE_LIBDEFLATE
    if (decomp_) libdeflate_free_decompressor(decomp_);
#endif
    inflateEnd(&zstrm_);
    if (fd_ >= 0) ::close(fd_);
}
