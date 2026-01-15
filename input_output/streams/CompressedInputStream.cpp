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
    /*std::memset(&zstrm_, 0, sizeof(zstrm_));
    if (inflateInit2(&zstrm_, 16 + MAX_WBITS) != Z_OK) {
        ::close(fd_); fd_ = -1;
        throw std::runtime_error("zlib inflateInit2 failed");
    }*/

    decomp_ = libdeflate_alloc_decompressor();
    if (!decomp_) {
        //inflateEnd(&zstrm_);
        ::close(fd_); fd_ = -1;
        throw std::runtime_error("libdeflate_alloc_decompressor failed");
    }

    stash_.reserve(opts.stash_reserve_bytes);
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

//[[clang::optnone]]
// Put this helper in the .cpp (static) or as a private method.
static inline bool logical_end_reached_for_next_read(uint64_t comp_off,
                                                     uint64_t end_coff,
                                                     uint32_t end_uoff)
{
    // End is an EXCLUSIVE virtual offset [ .., (end_coff,end_uoff) )
    // If next block starts after end_coff => we're past end.
    if (comp_off > end_coff) return true;

    // If next block starts exactly at end_coff, then:
    // - if end_uoff == 0, end is at the start of this block => nothing to read.
    // - if end_uoff > 0, we must read this end block and clamp in it.
    if (comp_off == end_coff && end_uoff == 0) return true;

    return false;
}

//[[clang::optnone]]
ssize_t CompressedInputStream::read(std::string& buffer, ssize_t offset, ssize_t length)
{
    if (length < 0 || offset < 0)
        throw std::out_of_range("CompressedInputStream::read: negative offset/length");

    const size_t uoff = static_cast<size_t>(offset);
    const size_t ulen = static_cast<size_t>(length);
    if (uoff > buffer.size() || uoff + ulen > buffer.size())
        throw std::out_of_range("CompressedInputStream::read: buffer too small for requested span");

    if (comp_pool.empty())
        throw std::logic_error("CompressedInputStream: comp_pool not initialized (size==0)");

    const size_t MAX_BATCH_BLOCKS = comp_pool.size();
    const size_t MAX_BATCH_BYTES  = MAX_BATCH_BLOCKS * BGZF_BLOCK_SIZE;

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
        compact_stash_();

        // IMPORTANT: decide logical end BEFORE planning
        if (has_end_ && logical_end_reached_for_next_read(comp_off_, end_coff_, end_uoff_)) {
            seen_end_ = true;
            break;
        }

        if (seen_end_) break;
        if (seen_eof_block_) break;

        const size_t stash_room = (stash_cap_ > stash_.size()) ? (stash_cap_ - stash_.size()) : 0;
        if (stash_room == 0) break;

        const size_t need_bytes = std::min(
            std::max(want, static_cast<size_t>(BGZF_BLOCK_SIZE)),
            std::min(stash_room, MAX_BATCH_BYTES)
        );

        struct BlockPlan {
            uint64_t off;
            uint32_t comp_size;
            uint32_t decomp_size;
            size_t   dest_off;
        };

        std::vector<BlockPlan> plan;
        plan.reserve(std::min(MAX_BATCH_BLOCKS, (need_bytes + BGZF_BLOCK_SIZE - 1) / BGZF_BLOCK_SIZE));

        size_t   batch_bytes = 0;
        uint64_t off_plan    = comp_off_;
        bool     saw_phys_eof = false;
        bool     saw_logical_end = false;

        while (!saw_phys_eof && !saw_logical_end &&
               batch_bytes < need_bytes && plan.size() < MAX_BATCH_BLOCKS)
        {
            // Hard stop: don't plan blocks strictly after end_coff_
            if (has_end_ && off_plan > end_coff_) {
                saw_logical_end = true;
                break;
            }

            BGZFBlockInfo info = probe_bgzf_block_size_(off_plan);
            if (info.comp_size == 0) { saw_phys_eof = true; break; }

            const uint32_t comp_sz   = static_cast<uint32_t>(info.comp_size);
            const uint32_t decomp_sz = static_cast<uint32_t>(info.decomp_size);

            if (decomp_sz == 0) { // BGZF EOF marker
                plan.push_back(BlockPlan{ off_plan, comp_sz, 0u, batch_bytes });
                off_plan += comp_sz;
                saw_phys_eof = true;
                break;
            }

            if (batch_bytes + decomp_sz > stash_room)       break;
            if (batch_bytes + decomp_sz > MAX_BATCH_BYTES)  break;

            plan.push_back(BlockPlan{ off_plan, comp_sz, decomp_sz, batch_bytes });
            off_plan    += comp_sz;
            batch_bytes += decomp_sz;

            // If we planned the end block, stop planning further.
            if (has_end_ && plan.back().off == end_coff_) {
                saw_logical_end = true;
                break;
            }
        }

        // If nothing planned, set the correct termination flag and exit.
        if (plan.empty()) {
            if (saw_phys_eof)  seen_eof_block_ = true;
            if (saw_logical_end) seen_end_ = true;   // <<< THIS was missing / crucial
            break;
        }

        const size_t old_size = stash_.size();
        if (old_size + batch_bytes > stash_cap_)
            throw std::runtime_error("Planning error: batch would exceed stash capacity");

        stash_.resize(old_size + batch_bytes);
        uint8_t* const stash_base = stash_.data() + old_size;

        // -------- READ COMPRESSED (serial) --------
        for (size_t i = 0; i < plan.size(); ++i) {
            const auto& j = plan[i];
            const ssize_t rn = ::pread(fd_, comp_pool[i].data(), j.comp_size, static_cast<off_t>(j.off));
            if (rn < 0)
                throw std::runtime_error(std::string("pread failed: ") + std::strerror(errno));
            if (rn != static_cast<ssize_t>(j.comp_size))
                throw std::runtime_error("CompressedInputStream: short read on BGZF block");
        }

        // -------- PARALLEL DECOMPRESS --------
        std::vector<size_t> produced(plan.size(), 0);
        std::vector<int>    status(plan.size(), 0); // 0 ok, 1 eof, -1 error

        //tbb::parallel_for(size_t(0), plan.size(), [&](size_t i){
        for (size_t i = 0; i < plan.size(); ++i) {
            const auto& j = plan[i];
            if (j.decomp_size == 0) { status[i] = 1; produced[i] = 0; break; }

            size_t out_n = 0;
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
                status[i] = -1; produced[i] = out_n; break;
            }
            status[i] = 0;
            produced[i] = j.decomp_size;
        }//);

        // -------- VALIDATION + ADVANCE + CLAMP END --------
        size_t appended = 0;
        bool   hit_end_block = false;

        for (size_t i = 0; i < plan.size(); ++i) {
            const auto& j = plan[i];

            if (j.decomp_size == 0) {
                seen_eof_block_ = true;
                comp_off_ = j.off + j.comp_size;
                continue;
            }

            if (status[i] != 0 || produced[i] != j.decomp_size) {
                stash_.resize(old_size + appended);
                throw std::runtime_error("BGZF decompress failed or ISIZE mismatch");
            }

            appended += j.decomp_size;
            comp_off_ = j.off + j.comp_size;

            if (has_end_ && j.off == end_coff_) {
                // Clamp to EXCLUSIVE end_uoff_ inside this block
                const size_t keep_in_block = std::min<size_t>(static_cast<size_t>(end_uoff_), j.decomp_size);
                stash_.resize(old_size + j.dest_off + keep_in_block);
                seen_end_ = true;
                hit_end_block = true;
                break;
            }

            if (old_size + appended == stash_cap_) break;
        }

        if (!hit_end_block) {
            // Shrink to only what we validated (avoid keeping unvalidated tail)
            stash_.resize(old_size + appended);
        }

        // If after all this we still produced nothing, bail out safely.
        if (stash_pos_ >= stash_.size()) {
            if (seen_end_ || seen_eof_block_) break;
            // No progress would cause caller to spin: treat as EOF-ish.
            throw std::runtime_error("read(): no progress (unexpected)");
        }
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

// If you call seek_virtual(), you likely want to reset seen_end_ too.

static inline bool voff_ge(uint64_t coff, uint32_t uoff, uint64_t end_coff, uint32_t end_uoff) {
    return (coff > end_coff) || (coff == end_coff && uoff >= end_uoff);
}

// If you call seek_virtual(), you likely want to reset seen_end_ too.
void CompressedInputStream::seek_virtual(uint64_t voff) {
    const uint64_t coff = voff >> 16;
    const uint32_t uoff = static_cast<uint32_t>(voff & 0xFFFFu);

    comp_off_ = coff;
    seen_eof_block_ = false;

    stash_.clear();
    stash_pos_ = 0;

    // If end is set and the requested position is already at/after end, clamp immediately.
    if (has_end_ && voff_ge(coff, uoff, end_coff_, end_uoff_)) {
        seen_end_ = true;
        return;
    }
    seen_end_ = false;

    const bool produced = inflate_one_bgzf_block_();
    if (!produced) {
        throw std::runtime_error("seek_virtual: landed on BGZF EOF/empty block");
    }
    if (uoff > stash_.size()) {
        throw std::runtime_error("seek_virtual: uoff beyond decompressed block size");
    }
    stash_pos_ = static_cast<size_t>(uoff);

    // IMPORTANT: if end is inside THIS same block, clamp stash right away.
    if (has_end_ && coff == end_coff_) {
        const size_t keep = std::min<size_t>(static_cast<size_t>(end_uoff_), stash_.size());
        // keep is exclusive end inside this block; but we've already consumed uoff bytes logically
        // via stash_pos_. If stash_pos_ >= keep => empty stream.
        stash_.resize(keep);
        if (stash_pos_ >= stash_.size()) {
            stash_.clear();
            stash_pos_ = 0;
            seen_end_ = true;
        }
    }
}

void CompressedInputStream::set_virtual_end(uint64_t voff) {
    end_coff_ = (voff >> 16);
    end_uoff_ = static_cast<uint32_t>(voff & 0xFFFFu);
    has_end_  = true;
    seen_end_ = false;

    // If we are already beyond the end, clamp immediately.
    // (We can't compute exact current voff without tracking block-start, so we only handle the
    //  trivial case where stash is empty and comp_off_ already past end block.)
    if (comp_off_ > end_coff_) {
        seen_end_ = true;
        stash_.clear();
        stash_pos_ = 0;
    }
}

void CompressedInputStream::clear_virtual_end() {
    has_end_  = false;
    seen_end_ = false;
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
    tmp_comp_.resize(comp_size);
    const ssize_t rn = ::pread(fd_, tmp_comp_.data(), comp_size, static_cast<off_t>(comp_off_));
    if (rn != static_cast<ssize_t>(comp_size)) {
        throw std::runtime_error("Short read on BGZF block");
    }
    comp_off_ += comp_size;

    // Decompress to tmp_decomp_ (resize to max; shrink after)
    tmp_decomp_.resize(65536);

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
    if (out_n == 0) {
        // BGZF EOF marker block (valid, zero-byte payload)
        tmp_decomp_.clear();
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
    if (decomp_) libdeflate_free_decompressor(decomp_);
    if (fd_ >= 0) ::close(fd_);
}
