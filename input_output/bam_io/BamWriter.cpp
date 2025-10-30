#include "BamWriter.h"
#include "../formatter/BamHeaderFormatter.h"
#include "../formatter/BamAlignmentFormatter.h"
#include "../formatter/SamAlignmentFormatter.h"
#include "../formatter/SamHeaderFormatter.h"
#include "../../sam_data_structures/SamHeader.h"
#include <algorithm>
#include <stop_token>
#include <condition_variable>
#include <cstring>
#include <stdexcept>
#include <execution>



BamWriter::BamWriter(std::unique_ptr<OutputStream> sink) : out_(std::move(sink)) {}


// --- Write 16-bit little-endian integer
inline void put_u16_le(uint8_t*& p, uint16_t v) noexcept {
    *p++ = static_cast<uint8_t>(v & 0xFF);
    *p++ = static_cast<uint8_t>((v >> 8) & 0xFF);
}

// --- Write 32-bit little-endian unsigned integer
inline void put_u32_le(uint8_t*& p, uint32_t v) noexcept {
    *p++ = static_cast<uint8_t>(v & 0xFF);
    *p++ = static_cast<uint8_t>((v >> 8) & 0xFF);
    *p++ = static_cast<uint8_t>((v >> 16) & 0xFF);
    *p++ = static_cast<uint8_t>((v >> 24) & 0xFF);
}

// --- Write 32-bit little-endian *signed* integer
inline void put_i32_le(uint8_t*& p, int32_t v) noexcept {
    put_u32_le(p, static_cast<uint32_t>(v));
}

inline void BamWriter::push_u32(std::vector<uint8_t>& v, uint32_t x) {
    v.push_back(uint8_t(x));
    v.push_back(uint8_t(x >> 8));
    v.push_back(uint8_t(x >> 16));
    v.push_back(uint8_t(x >> 24));
}

inline void BamWriter::push_u16(std::vector<uint8_t>& v, uint16_t x) {
    v.push_back(uint8_t(x));
    v.push_back(uint8_t(x >> 8));
}

inline void BamWriter::push_u8 (std::vector<uint8_t>& v, uint8_t x) { v.push_back(x); }


uint32_t reg2bin(int32_t beg, int32_t end) {
    --end;
    if (beg>>14 == end>>14) return ((1<<15)-1)/7 + (beg>>14);
    if (beg>>17 == end>>17) return ((1<<12)-1)/7 + (beg>>17);
    if (beg>>20 == end>>20) return ((1<<9 )-1)/7 + (beg>>20);
    if (beg>>23 == end>>23) return ((1<<6 )-1)/7 + (beg>>23);
    if (beg>>26 == end>>26) return ((1<<3 )-1)/7 + (beg>>26);
    return 0;
}
uint32_t pack_cigar(uint32_t len, char op) {
    int code = 0;
    switch (op) {
        case 'M': code=0; break; case 'I': code=1; break; case 'D': code=2; break;
        case 'N': code=3; break; case 'S': code=4; break; case 'H': code=5; break;
        case 'P': code=6; break; case '=': code=7; break; case 'X': code=8; break;
        case 'B': code=9; break;
    }
    return (len << 4) | (code & 0xF);
}
inline uint8_t base_code(char b) {
    switch (b) {
        case '=': return 0;
        case 'A': case 'a': return 1;
        case 'C': case 'c': return 2;
        case 'G': case 'g': return 4;
        case 'T': case 't': return 8;
        default: return 15; // N/unknown
    }
}

inline void pack_seq(const char* seq, int len, uint8_t* dst)
{
    for (int i = 0; i < len; ++i) {
        const uint8_t code = base_code(seq[i]);
        if (i & 1)
            dst[i >> 1] = (dst[i >> 1] & 0xF0) | code;
        else
            dst[i >> 1] = (code << 4);
    }

    // If len is odd, last nibble (low 4 bits) remains undefined,
    // but this is fine — BAM readers ignore it.
}

//[[clang::optnone]]
void BamWriter::FormatHeader(Header& hdr) {
    std::call_once(header_once_, [&] {
        //if (header_written_) return;

        // 1) Produce SAM text from your formatter
        std::string sam_text;
        SamHeaderFormatter::FormatSamHeader(hdr, sam_text);  // ensure this returns SAM (text) header

        // 2) Build BAM header blob: magic + l_text + text + n_ref + refs
        std::vector<uint8_t> buf;
        buf.reserve(4096 + sam_text.size());

        // magic "BAM\1"
        buf.insert(buf.end(), {'B','A','M','\1'});

        // l_text
        push_u32(buf, (uint32_t)sam_text.size());
        // text
        buf.insert(buf.end(), sam_text.begin(), sam_text.end());

        // references
        const auto& seqs = hdr.getSequenceDictionary().getAll();
        push_u32(buf, (uint32_t)seqs.size());
        for (const auto& r : seqs) {
            const std::string& name = r.getSequenceName();
            push_u32(buf, (uint32_t)(name.size() + 1));
            buf.insert(buf.end(), name.begin(), name.end());
            buf.push_back('\0');
            push_u32(buf, (uint32_t)r.getSequenceLength());
        }

        out_->write(std::span<const uint8_t>(buf.data(), buf.size()));
        header_written_.store(true, std::memory_order_release);
    });
}

inline size_t record_wire_size(SamAlignments& r, size_t idx) {
    const size_t l_qname = r.GetQNAME(idx).size() + 1;  // NUL-terminated in BAM
    const size_t n_cigar = r.GetCIGAR(idx).size();
    const size_t l_seq   = r.GetSEQ(idx).size();

    size_t aux_len = r.TAGS.tagsBlobSize[idx];

    const size_t block_size =
        32
      + l_qname
      + 4 * n_cigar
      + ((l_seq + 1) / 2)
      + l_seq
      + aux_len;

    return 4 + block_size;
}

inline int32_t cigar_ref_length(std::span<CigarOperation> ops) noexcept {
    int32_t len = 0;
    for (auto& op : ops) {
        switch (op.Operation) {
            case 'M': case 'D': case 'N': case '=': case 'X': len += op.Length; break;
            default: break;
        }
    }
    return len;
}


inline void pack_seq_bytes(const char* s, int l, uint8_t* out) {
    // writes (l+1)/2 bytes into out
    for (int i = 0; i < l; ++i) {
        const uint8_t code = base_code(s[i]);
        if (i & 1) out[i >> 1] = (out[i >> 1] & 0xF0) | code;
        else       out[i >> 1] = (code << 4);
    }
}

inline int32_t cigar_ref_length(std::span<const CigarOperation> ops) {
    int32_t len = 0;
    for (auto& op : ops) {
        switch (op.Operation) {
            case 'M': case 'D': case 'N': case '=': case 'X': len += op.Length; break;
            default: break;
        }
    }
    return len > 0 ? len : 1; // at least 1 to avoid end==beg in unmapped oddities
}

// Correct BAM record writer
void write_record_into(uint8_t* dst, size_t bytes_expected,
                       SamAlignments& reads, size_t idx)
{
    uint8_t* p = dst;

    // reserve 4 bytes for block_size
    int32_t* block_size_ptr = reinterpret_cast<int32_t*>(p);
    p += 4;

    const bool unmapped = (reads.GetRNAME(idx) == SequenceRecord::RESERVED_SEQUENCE_NAMES['*']);

    const int32_t tid   = unmapped ? -1 : reads.GetRNAME(idx);
    const int32_t pos0  = unmapped ? -1 : std::max<int32_t>(0, reads.GetPOS(idx) - 1);
    const int32_t l_qn  = (int32_t)reads.GetQNAME(idx).size() + 1;
    const uint8_t mapq  = unmapped ? 0 : reads.GetMAPQ(idx);
    const uint16_t n_cig = unmapped ? 0
                                    : (uint16_t)reads.GetCIGAR(idx).size();
    const uint16_t flag  = reads.GetFLAG(idx);
    const int32_t  l_seq = (int32_t)reads.GetSEQ(idx).size();

    // *** Mate fields: do NOT zero just because 0x8 is set ***
    const bool mate_star = (reads.GetRNEXT(idx) == SequenceRecord::RESERVED_SEQUENCE_NAMES['*']);
    const int32_t mtid   = mate_star ? -1 : reads.GetRNEXT(idx);
    const int32_t mpos0  = mate_star ? -1 : std::max<int32_t>(0, reads.GetPNEXT(idx) - 1);

    // TLEN: keep whatever you parsed for paired reads; else 0
    const bool paired = (flag & 0x1) != 0;
    const int32_t tlen = paired ? reads.GetTLEN(idx) : 0;

    // bin only for mapped
    int32_t beg0 = pos0, end0 = pos0;
    if (!unmapped) {
        const auto cig = reads.GetCIGAR(idx);
        end0 = beg0 + cigar_ref_length(cig);
    }
    const uint16_t bin = unmapped ? 0 : (uint16_t)reg2bin(beg0, end0);

    // ---- 32-byte core in exact spec order ----
    put_i32_le(p, tid);           // refID
    put_i32_le(p, pos0);          // pos (0-based)
    *p++ = (uint8_t)l_qn;         // l_read_name
    *p++ = (uint8_t)mapq;         // mapq
    put_u16_le(p, bin);           // bin
    put_u16_le(p, n_cig);         // n_cigar_op
    put_u16_le(p, flag);          // flag
    put_i32_le(p, l_seq);         // l_seq  <-- you were missing this field
    put_i32_le(p, mtid);          // next_refID
    put_i32_le(p, mpos0);         // next_pos (0-based)
    put_i32_le(p, tlen);          // tlen

    // ---- read_name (NUL-terminated) ----
    const auto qn = reads.GetQNAME(idx);
    std::memcpy(p, qn.data(), qn.size());
    p[qn.size()] = '\0';
    p += qn.size() + 1;

    // ---- CIGAR ----
    if (!unmapped && n_cig) {
        auto* dst32 = reinterpret_cast<uint32_t*>(p);
        const auto cig = reads.GetCIGAR(idx);
        for (size_t j = 0; j < cig.size(); ++j) {
            dst32[j] = pack_cigar((uint32_t)cig[j].Length, cig[j].Operation);
        }
        p += n_cig * 4;
    }

    // ---- SEQ (packed) ----
    if (l_seq > 0) {
        pack_seq_bytes(reads.GetSEQ(idx).data(), l_seq, p);
        p += (l_seq + 1) >> 1;
    }

    // ---- QUAL (raw phred or 0xFF for missing) ----
    if (l_seq > 0) {
        const auto qv = reads.GetQUAL(idx);
        if (qv.empty()) {
            std::memset(p, 0xFF, l_seq);
        } else {
            std::memcpy(p, qv.data(), l_seq);
        }
        p += l_seq;
    }

    // ---- AUX (use offset + size!) ----
    size_t aux_off = 0, aux_sz = 0;
    // You must have these two arrays; if you only have size, add offset recording during parsing
    aux_off = reads.TAGS.tagsBlobOff[idx];   // <--- not tagsBlobSize
    aux_sz  = reads.TAGS.tagsBlobSize[idx];
    if (aux_sz) {
        const uint8_t* aux_ptr = reads.TAGS.tagArena.data() + aux_off;
        std::memcpy(p, aux_ptr, aux_sz);
        p += aux_sz;
    }

    // ---- block_size (payload after the int32 itself) ----
    const int32_t block_size = static_cast<int32_t>(p - (dst + 4));
    std::memcpy(block_size_ptr, &block_size, 4);

    if ((size_t)(p - dst) != bytes_expected) {
        throw std::runtime_error("BAM write_record_into(): size mismatch");
    }
}

static bool process_batch(const std::shared_ptr<OutputBatchBam>& bt,
                          SamAlignments& reads,
                          size_t start_idx,
                          size_t& n_emitted)
{
    n_emitted = 0;

    const size_t total = reads.GetAlnsNum();
    if (start_idx >= total) return false;

    const size_t N   = total - start_idx;
    const size_t cap = bt->capacity();
    if (cap == 0) return false;

    // Phase 1: size all candidate records
    std::vector<size_t> sizes(N);
    tbb::parallel_for(size_t(0), N, [&](size_t i) {
        sizes[i] = record_wire_size(reads, start_idx + i);
    });

    if (sizes[0] > cap) return false; // record too large for fixed batch

    // Prefix sums
    std::vector<size_t> offsets(N);
    std::exclusive_scan(std::execution::par, sizes.begin(), sizes.end(), offsets.begin(), size_t(0));

    // How many fit
    size_t lo = 1, hi = N + 1;
    while (lo < hi) {
        size_t mid = (lo + hi) >> 1;
        const size_t end_bytes = offsets[mid - 1] + sizes[mid - 1];
        if (end_bytes <= cap) lo = mid + 1;
        else                  hi = mid;
    }
    const size_t fit = lo - 1;
    if (fit == 0) return false;

    const size_t batch_bytes = offsets[fit - 1] + sizes[fit - 1];

    // Phase 2: write in parallel
    uint8_t* base = bt->data();
    tbb::parallel_for(size_t(0), fit, [&](size_t i) {
        write_record_into(base + offsets[i], sizes[i], reads, start_idx + i);
    });

    bt->seal(batch_bytes, fit);
    n_emitted = fit;
    return true;
}

//[[clang::optnone]]
void BamWriter::store_batch_(const std::shared_ptr<OutputBatchBam>& bt,
                             std::unique_ptr<OutputStream>& out) {
    const size_t n = bt->used();
    if (n) out->write(std::span<const uint8_t>(bt->data(), n));
    bt->reset(); // clean right after write
}

void BamWriter::BeginStreaming() {
    if (!header_written_) throw std::runtime_error("header must be written first");
    if (!out_)            throw std::runtime_error("no output stream");
    if (streaming_started_) return;

    writer_ex_ = nullptr;
    producer_ex_ = nullptr;
    have_work_ = wrote_ = stop_ = false;

    batch_process_  = std::make_shared<OutputBatchBam>(DEFAULT_BUFFER_SIZE);
    batch_transfer_ = std::make_shared<OutputBatchBam>(DEFAULT_BUFFER_SIZE);

    // spawn writer thread once
    writer_thread_ = std::jthread([this] (std::stop_token st) {
        try {
            for (;;) {
                std::unique_lock lk(m_);
                cv_.wait(lk, st, [&]{ return have_work_ || stop_; });
                if (st.stop_requested()) break;
                if (stop_ && !have_work_) break;
                have_work_ = false;
                auto bt = batch_transfer_;
                lk.unlock();

                store_batch_(bt, out_);

                lk.lock();
                wrote_ = true;
                lk.unlock();
                cv_.notify_one();
            }
        } catch (...) {
            writer_ex_ = std::current_exception();
            std::lock_guard<std::mutex> lk(m_);
            wrote_ = true;
            cv_.notify_one();
        }
    });

    streaming_started_ = true;
}

void BamWriter::WriteAlignmentBatch(SamAlignments& block) {
    FormatHeader(block.GetHdr()); // guarded by call_once
    if (!streaming_started_) BeginStreaming();

    try {
        // Keep filling the process buffer with this block’s alignments,
        // handing off to the writer whenever the buffer has data.
        size_t aln_idx = 0;
        const size_t n_alns = static_cast<size_t>(block.GetAlnsNum());

        // If process buffer is empty, ensure it’s reset
        if (batch_process_->used() == 0) {
            batch_process_->reset();
        }

        while (aln_idx < n_alns) {
            size_t n_emitted = 0;
            const bool any = process_batch(batch_process_, block, aln_idx, n_emitted);
            if (!any) break;               // nothing written for this block (e.g., no data)
            aln_idx += n_emitted;

            // If we produced *any* bytes, hand off immediately to overlap I/O
            if (batch_process_->used() > 0) {
                {
                    std::unique_lock lk(m_);
                    wrote_ = false;
                    std::swap(batch_process_, batch_transfer_);
                    have_work_ = true;
                    lk.unlock();
                    cv_.notify_one();
                }

                // We can start filling the other buffer right away
                batch_process_->reset();

                // Optional: if you want strict backpressure, wait here until write completes:
                {
                    std::unique_lock lk(m_);
                    cv_.wait(lk, [&]{ return wrote_ || writer_ex_ != nullptr; });
                }
                if (writer_ex_) std::rethrow_exception(writer_ex_);
            }
        }
    } catch (...) {
        producer_ex_ = std::current_exception();
        // propagate immediately
        if (producer_ex_) std::rethrow_exception(producer_ex_);
    }
}

void BamWriter::EndStreaming() {
    if (!streaming_started_) return;

    try {
        // Flush any remaining bytes in process buffer
        if (batch_process_ && batch_process_->used() > 0) {
            {
                std::unique_lock lk(m_);
                wrote_ = false;
                std::swap(batch_process_, batch_transfer_);
                have_work_ = true;
                lk.unlock();
                cv_.notify_one();
            }
            // Wait for the final write
            {
                std::unique_lock lk(m_);
                cv_.wait(lk, [&]{ return wrote_ || writer_ex_ != nullptr; });
            }
            if (writer_ex_) std::rethrow_exception(writer_ex_);
        }
    } catch (...) {
        producer_ex_ = std::current_exception();
    }

    // Shutdown writer thread
    {
        std::lock_guard lk(m_);
        stop_ = true;
    }
    writer_thread_.request_stop();
    cv_.notify_one();

    if (writer_thread_.joinable()) writer_thread_.join();

    streaming_started_ = false;

    if (producer_ex_) std::rethrow_exception(producer_ex_);
    if (writer_ex_)   std::rethrow_exception(writer_ex_);
}

void BamWriter::FormatAlignment(std::vector<SamAlignments>& blocks) {
    if (!header_written_) throw std::runtime_error("header must be written first");
    if (!out_)            throw std::runtime_error("no output stream");

    BeginStreaming();
    for (auto& blk : blocks) {
        WriteAlignmentBatch(blk);
    }
    EndStreaming();
}

void BamWriter::Close() {
    if (streaming_started_) {
        EndStreaming();
    }
    if (out_) out_->close();
}

BamWriter::~BamWriter() {
    // best-effort, but avoid double work:
    /*try {
        if (streaming_started_) EndStreaming();
    } catch (...) {}*/
}