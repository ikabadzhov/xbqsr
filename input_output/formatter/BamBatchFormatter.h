// BamBatchFormatter.h
#pragma once
#include <vector>
#include <span>
#include <cstdint>
#include <mutex>
#include <atomic>
#include <array>
#include <execution>
#include <cstring>
#include <stdexcept>
#include <algorithm>

#include <tbb/parallel_for.h>

#include "../../sam_data_structures/SamHeader.h"
#include "../../sam_data_structures/SamAlignments.h"
#include "../formatter/SamHeaderFormatter.h"
#include "../parser/IOBatch.h"

class BamBatchFormatter {
public:
    BamBatchFormatter() = default;

    // Produces BAM header bytes (magic + header text + ref dict).
    // Header is produced once per formatter instance (call_once).
    void FormatHeaderTo(Header& hdr, std::vector<uint8_t>& out);

    // Formats the whole SamAlignments batch into a contiguous byte blob:
    // [block_size(4) + record_payload] repeated for each record.
    void FormatAlignmentBatchTo(SamAlignments& reads, std::vector<uint8_t>& out);

    // Same, but writes into OutputBatchBam (grow-only). Seals used bytes + n_records.
    void FormatAlignmentBatchTo(SamAlignments& reads, OutputBatchBam& out);

    bool HeaderWritten() const { return header_written_.load(std::memory_order_acquire); }

private:
    // ---- little-endian put helpers (write into raw pointer) ----
    static inline void put_u16_le(uint8_t*& p, uint16_t v) noexcept {
        *p++ = static_cast<uint8_t>(v & 0xFF);
        *p++ = static_cast<uint8_t>((v >> 8) & 0xFF);
    }
    static inline void put_u32_le(uint8_t*& p, uint32_t v) noexcept {
        *p++ = static_cast<uint8_t>(v & 0xFF);
        *p++ = static_cast<uint8_t>((v >> 8) & 0xFF);
        *p++ = static_cast<uint8_t>((v >> 16) & 0xFF);
        *p++ = static_cast<uint8_t>((v >> 24) & 0xFF);
    }
    static inline void put_i32_le(uint8_t*& p, int32_t v) noexcept {
        put_u32_le(p, static_cast<uint32_t>(v));
    }

    // ---- vector append helpers ----
    static inline void push_u32(std::vector<uint8_t>& v, uint32_t x) {
        v.push_back(uint8_t(x));
        v.push_back(uint8_t(x >> 8));
        v.push_back(uint8_t(x >> 16));
        v.push_back(uint8_t(x >> 24));
    }
    static inline void push_u16(std::vector<uint8_t>& v, uint16_t x) {
        v.push_back(uint8_t(x));
        v.push_back(uint8_t(x >> 8));
    }
    static inline void push_u8(std::vector<uint8_t>& v, uint8_t x) { v.push_back(x); }

    static uint32_t reg2bin(int32_t beg, int32_t end) {
        --end;
        if (beg>>14 == end>>14) return ((1<<15)-1)/7 + (beg>>14);
        if (beg>>17 == end>>17) return ((1<<12)-1)/7 + (beg>>17);
        if (beg>>20 == end>>20) return ((1<<9 )-1)/7 + (beg>>20);
        if (beg>>23 == end>>23) return ((1<<6 )-1)/7 + (beg>>23);
        if (beg>>26 == end>>26) return ((1<<3 )-1)/7 + (beg>>26);
        return 0;
    }

    static uint32_t pack_cigar(uint32_t len, char op) {
        int code = 0;
        switch (op) {
            case 'M': code=0; break; case 'I': code=1; break; case 'D': code=2; break;
            case 'N': code=3; break; case 'S': code=4; break; case 'H': code=5; break;
            case 'P': code=6; break; case '=': code=7; break; case 'X': code=8; break;
            case 'B': code=9; break;
            default:  code=0; break;
        }
        return (len << 4) | (code & 0xF);
    }

    static inline uint8_t base_code(char b) {
        static constexpr std::array<uint8_t, 256> LUT = []{
            std::array<uint8_t,256> a{};
            a.fill(15);
            a[static_cast<unsigned char>('=')] = 0;
            a[static_cast<unsigned char>('A')] = 1; a[static_cast<unsigned char>('a')] = 1;
            a[static_cast<unsigned char>('C')] = 2; a[static_cast<unsigned char>('c')] = 2;
            a[static_cast<unsigned char>('G')] = 4; a[static_cast<unsigned char>('g')] = 4;
            a[static_cast<unsigned char>('T')] = 8; a[static_cast<unsigned char>('t')] = 8;
            return a;
        }();
        return LUT[static_cast<unsigned char>(b)];
    }

    // Packs DNA bases from s[0..l) into out[(l+1)/2] bytes.
    static inline void paq_seq(const char* s, int l, uint8_t* out) {
        int i = 0;
        for (; i + 1 < l; i += 2) {
            const uint8_t hi = base_code(s[i]);
            const uint8_t lo = base_code(s[i + 1]);
            out[i >> 1] = static_cast<uint8_t>((hi << 4) | (lo & 0x0F));
        }
        if (i < l) {
            const uint8_t hi = base_code(s[i]);
            out[i >> 1] = static_cast<uint8_t>(hi << 4);
        }
    }

    static inline int32_t cigar_ref_length(std::span<const CigarOperation> ops) {
        int32_t len = 0;
        for (auto& op : ops) {
            switch (op.Operation) {
                case 'M': case 'D': case 'N': case '=': case 'X': len += op.Length; break;
                default: break;
            }
        }
        return len > 0 ? len : 1;
    }

    static inline size_t record_wire_size(SamAlignments& r, size_t idx) {
        const size_t l_qname = r.GetQNAME(idx).size() + 1;  // NUL-terminated in BAM
        const size_t n_cigar = r.GetCIGAR(idx).size();
        const size_t l_seq   = r.GetSEQ(idx).size();
        const size_t aux_len = r.TAGS.tagsBlobSize[idx];

        const size_t block_size =
            32
          + l_qname
          + 4 * n_cigar
          + ((l_seq + 1) / 2)
          + l_seq
          + aux_len;

        return 4 + block_size; // includes block_size int32
    }

    static void write_record_into(uint8_t* dst, size_t bytes_expected,
                                  SamAlignments& reads, size_t idx)
    {
        uint8_t* p = dst;

        // reserve 4 bytes for block_size
        uint8_t* block_size_ptr = p;
        p += 4;

        const bool unmapped = (reads.GetRNAME(idx) == SequenceRecord::RESERVED_SEQUENCE_NAMES['*']);

        const int32_t tid   = unmapped ? -1 : reads.GetRNAME(idx);
        const int32_t pos0  = unmapped ? -1 : std::max<int32_t>(0, reads.GetPOS(idx) - 1);
        const int32_t l_qn  = (int32_t)reads.GetQNAME(idx).size() + 1;
        const uint8_t mapq  = unmapped ? 0 : reads.GetMAPQ(idx);
        const uint16_t n_cig = unmapped ? 0 : (uint16_t)reads.GetCIGAR(idx).size();
        const uint16_t flag  = reads.GetFLAG(idx);
        const int32_t  l_seq = (int32_t)reads.GetSEQ(idx).size();

        const bool mate_star = (reads.GetRNEXT(idx) == SequenceRecord::RESERVED_SEQUENCE_NAMES['*']);
        const int32_t mtid   = mate_star ? -1 : reads.GetRNEXT(idx);
        const int32_t mpos0  = mate_star ? -1 : std::max<int32_t>(0, reads.GetPNEXT(idx) - 1);

        const bool paired = (flag & 0x1) != 0;
        const int32_t tlen = paired ? reads.GetTLEN(idx) : 0;

        int32_t beg0 = pos0, end0 = pos0;
        if (!unmapped) {
            const auto cig = reads.GetCIGAR(idx);
            end0 = beg0 + cigar_ref_length(std::span<const CigarOperation>(cig.data(), cig.size()));
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
        put_i32_le(p, l_seq);         // l_seq
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
            paq_seq(reads.GetSEQ(idx).data(), l_seq, p);
            p += (l_seq + 1) >> 1;
        }

        // ---- QUAL ----
        if (l_seq > 0) {
            const auto qv = reads.GetQUAL(idx);
            if (qv.empty()) {
                std::memset(p, 0xFF, (size_t)l_seq);
            } else {
                std::memcpy(p, qv.data(), (size_t)l_seq);
            }
            p += (size_t)l_seq;
        }

        // ---- AUX (offset + size) ----
        const size_t aux_off = reads.TAGS.tagsBlobOff[idx];
        const size_t aux_sz  = reads.TAGS.tagsBlobSize[idx];
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

    std::once_flag header_once_;
    std::atomic<bool> header_written_{false};
};