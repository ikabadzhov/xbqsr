#include "BamAlignmentParser.h"
#include "../utils/ByteScanner.h"

// Parse exactly one BAM alignment from a complete record body (i.e., block_size bytes).
// Precondition: "bytes" is exactly the block body (NOT including the 4-byte block_size field).
void BamAlignmentParser::ParseBamAlignmentFromSlice(std::string_view bytes,
                                                    SamAlignments& reads,
                                                    size_t idx)
{

    // -------------------- PARSING --------------------
    ByteScanner bs;
    bs.reset(bytes);

    int32_t  refID        = bs.readInt32();
    int32_t  pos0         = bs.readInt32();
    uint8_t  l_read_name  = bs.readByte();
    uint8_t  mapq         = bs.readByte();
    uint16_t bin_mq_nl    = bs.readInt16(); (void)bin_mq_nl;
    uint16_t n_cigar_op   = bs.readUint16();
    uint16_t flag         = bs.readUint16();
    uint32_t l_seq        = bs.readUint32();
    int32_t  next_refID   = bs.readInt32();
    int32_t  next_pos0    = bs.readInt32();
    int32_t  tlen         = bs.readInt32();

    std::string qname = bs.readString(l_read_name);
    if (!qname.empty() && qname.back() == '\0') qname.pop_back();

    std::vector<CigarOperation> cigarVector(n_cigar_op);
    for (uint32_t i = 0; i < n_cigar_op; ++i) {
        uint32_t op = bs.readInt32();
        cigarVector[i] = CigarOperation{
            (op >> 4),
            cigarOps[static_cast<int>(op & 0x0F)]
        };
    }

    const size_t seq_bytes = (static_cast<size_t>(l_seq) + 1u) >> 1;
    std::string decoded_seq;
    decoded_seq.resize(l_seq);
    {
        size_t k = 0;
        for (size_t i = 0; i < seq_bytes; ++i) {
            uint8_t nibbles = bs.readByte();
            uint8_t first   = nibbles >> 4;
            uint8_t second  = nibbles & 0x0F;
            decoded_seq[k++] = seqEncoder[first];
            if (k < decoded_seq.size()) decoded_seq[k++] = seqEncoder[second];
        }
    }

    std::vector<char> qual(static_cast<size_t>(l_seq));
    if (l_seq) {
        const char* src = bytes.data() + bs.getIndex();
        std::memcpy(qual.data(), src, static_cast<size_t>(l_seq));
        bs.incrementIndex(static_cast<size_t>(l_seq));
    }

    // ---------------- SAM STRUCT UPDATES ----------------

    reads.RNAME[idx] = (refID < 0) ? SequenceRecord::RESERVED_SEQUENCE_NAMES['*'] : refID;
    if (next_refID < 0) next_refID = SequenceRecord::RESERVED_SEQUENCE_NAMES['*'];
    reads.RNEXT[idx] = next_refID;

    reads.POS[idx]   = pos0 + 1;
    reads.PNEXT[idx] = next_pos0 + 1;
    reads.MAPQ[idx]  = mapq;
    reads.FLAG[idx]  = flag;
    reads.TLEN[idx]  = tlen;

    if (!qname.empty()) {
        reads.TILE[idx] = computeTileInfo(qname);
    }

    reads.CIGAR.SetAt(cigarVector, idx);
    reads.SEQ.SetAt(decoded_seq, idx);
    reads.QUAL.SetAt(qual, idx);

    // ---- TAGS (parsing + writing counted as "updates") ----
    TagKeyAndType tagKeyAndType;
    while (bs.canRead(3)) {
        tagKeyAndType.setTagKey(TagKeyAndType::encode(bs.readString(2)));
        tagKeyAndType.setTagType(static_cast<char>(bs.readByte()));
        const char T = tagKeyAndType.getTagType();

        if (T == 'A') {
            reads.TAGS.InsertCharTagAt(idx, tagKeyAndType, bs.readByte());
        } else if (T == 'c') {
            reads.TAGS.InsertIntTagAt(idx, tagKeyAndType, static_cast<int8_t>(bs.readByte()));
        } else if (T == 'C') {
            reads.TAGS.InsertIntTagAt(idx, tagKeyAndType, static_cast<uint8_t>(bs.readByte()));
        } else if (T == 's') {
            reads.TAGS.InsertIntTagAt(idx, tagKeyAndType, bs.readInt16());
        } else if (T == 'S') {
            reads.TAGS.InsertIntTagAt(idx, tagKeyAndType, static_cast<uint16_t>(bs.readInt16()));
        } else if (T == 'i') {
            reads.TAGS.InsertIntTagAt(idx, tagKeyAndType, bs.readInt32());
        } else if (T == 'I') {
            reads.TAGS.InsertIntTagAt(idx, tagKeyAndType, static_cast<uint32_t>(bs.readInt32()));
        } else if (T == 'f') {
            reads.TAGS.InsertFloatTagAt(idx, tagKeyAndType, bs.readFloat());
        } else if (T == 'Z' || T == 'H') {
            reads.TAGS.InsertStringTagAt(idx, tagKeyAndType, bs.readUntil(0));
        } else if (T == 'B') {
            const char subtype = static_cast<char>(bs.readByte());
            tagKeyAndType.setSecondaryTagType(subtype);
            const int32_t count = bs.readInt32();
            size_t elem_size = (subtype=='c'||subtype=='C') ? 1 :
                               (subtype=='s'||subtype=='S') ? 2 : 4;
            size_t nbytes = (size_t)count * elem_size;
            const char* p = bs.data();

            if (subtype == 'c') {
                std::vector<int8_t> v(count);
                std::memcpy(v.data(), p, nbytes);
                reads.TAGS.InsertArrayInt8TagAt(idx, tagKeyAndType, v);
            } else if (subtype == 'C') {
                std::vector<uint8_t> v(count);
                std::memcpy(v.data(), p, nbytes);
                reads.TAGS.InsertArrayUint8TagAt(idx, tagKeyAndType, v);
            } else if (subtype == 's') {
                std::vector<int16_t> v(count);
                std::memcpy(v.data(), p, nbytes);
                reads.TAGS.InsertArrayInt16TagAt(idx, tagKeyAndType, v);
            } else if (subtype == 'S') {
                std::vector<uint16_t> v(count);
                std::memcpy(v.data(), p, nbytes);
                reads.TAGS.InsertArrayUint16TagAt(idx, tagKeyAndType, v);
            } else if (subtype == 'i') {
                std::vector<int32_t> v(count);
                std::memcpy(v.data(), p, nbytes);
                reads.TAGS.InsertArrayInt32TagAt(idx, tagKeyAndType, v);
            } else if (subtype == 'I') {
                std::vector<uint32_t> v(count);
                std::memcpy(v.data(), p, nbytes);
                reads.TAGS.InsertArrayUint32TagAt(idx, tagKeyAndType, v);
            } else { // 'f'
                std::vector<float> v(count);
                std::memcpy(v.data(), p, nbytes);
                reads.TAGS.InsertArrayFloatTagAt(idx, tagKeyAndType, v);
            }
            bs.incrementIndex(nbytes);
        } else {
            throw std::runtime_error("Unknown BAM tag type");
        }
    }
}

void BamAlignmentParser::ParseBamAlignmentFromSliceSimpleTypes(
    std::string_view bytes,
    SamAlignments &reads,
    size_t idx,
    VariableLenFieldsVar& fieldsOff,
    VariableLenFieldsVar& fieldsSizes
) {
    ByteScanner bs;
    bs.reset(bytes);

    // --- Fixed 32-byte BAM core block ---
    if (!bs.canRead(32))
        throw std::runtime_error("Truncated BAM core block");

    int32_t  refID        = bs.readInt32();
    int32_t  pos0         = bs.readInt32();
    uint8_t  l_read_name  = bs.readByte();
    uint8_t  mapq         = bs.readByte();
    uint16_t bin_mq_nl    = bs.readInt16(); (void)bin_mq_nl;
    uint16_t n_cigar_op   = bs.readUint16();
    uint16_t flag         = bs.readUint16();
    uint32_t l_seq        = bs.readUint32();
    int32_t  next_refID   = bs.readInt32();
    int32_t  next_pos0    = bs.readInt32();
    int32_t  tlen         = bs.readInt32();

    // --- QNAME ---
    if (!bs.canRead(l_read_name))
        throw std::runtime_error("Truncated QNAME");

    fieldsOff.qname   = bs.getIndex();
    fieldsSizes.qname = l_read_name;      // includes '\0'
    bs.incrementIndex(l_read_name);

    // --- CIGAR ---
    const size_t cigar_bytes = size_t(n_cigar_op) * sizeof(uint32_t);
    if (!bs.canRead(cigar_bytes))
        throw std::runtime_error("Truncated CIGAR");

    fieldsOff.cigar   = bs.getIndex();
    fieldsSizes.cigar = n_cigar_op;       // in ops
    bs.incrementIndex(cigar_bytes);

    // --- SEQ ---
    const size_t seq_bytes = (size_t(l_seq) + 1u) >> 1;
    if (!bs.canRead(seq_bytes))
        throw std::runtime_error("Truncated SEQ");

    fieldsOff.seq   = bs.getIndex();
    fieldsSizes.seq = l_seq;              // in bases
    bs.incrementIndex(seq_bytes);

    // --- QUAL ---
    if (!bs.canRead(size_t(l_seq)))
        throw std::runtime_error("Truncated QUAL");

    fieldsOff.qual   = bs.getIndex();
    fieldsSizes.qual = l_seq;
    bs.incrementIndex(size_t(l_seq));

    // --- AUX (whatever's left) ---
    const size_t aux_off  = bs.getIndex();
    const size_t aux_size = bytes.size() - aux_off;

    fieldsOff.aux   = aux_off;
    fieldsSizes.aux = aux_size;

    // --- Store into SoA ---

    reads.RNAME[idx] = (refID < 0)
        ? SequenceRecord::RESERVED_SEQUENCE_NAMES['*']
        : refID;

    if (next_refID < 0)
        next_refID = SequenceRecord::RESERVED_SEQUENCE_NAMES['*'];

    reads.RNEXT[idx] = next_refID;
    reads.POS[idx]   = pos0 + 1;
    reads.PNEXT[idx] = next_pos0 + 1;
    reads.MAPQ[idx]  = mapq;
    reads.FLAG[idx]  = flag;
    reads.TLEN[idx]  = tlen;
}

static const uint8_t* skip_payload(const uint8_t* q,
                                          const uint8_t* e,
                                          char T)
{
    auto adv = [&](size_t n) -> const uint8_t* {
        const uint8_t* r = q + n;
        return (r <= e) ? r : e;
    };

    switch (T) {
        // 1-byte scalar
        case 'A': case 'c': case 'C':
            return adv(1);

            // 2-byte scalar
        case 's': case 'S':
            return adv(2);

            // 4-byte scalar
        case 'i': case 'I': case 'f':
            return adv(4);

            // Zero-terminated string / hex
        case 'Z': case 'H': {
            const uint8_t* r = q;
            while (r < e && *r != 0) ++r;
            return (r < e) ? (r + 1) : e; // include NUL if present
        }
            // Array
        case 'B': {
            // Need at least subtype(1) + count(4)
            if (q + 5 > e) return e;

            const char subtype = char(*q++);
            uint32_t count = 0;
            std::memcpy(&count, q, 4); // BAM: little-endian
            q += 4;

            size_t elem_size;
            switch (subtype) {
                case 'c': case 'C': elem_size = 1; break;
                case 's': case 'S': elem_size = 2; break;
                case 'i': case 'I':
                case 'f':           elem_size = 4; break;
                default:            elem_size = 4; break; // safe fallback
            }

            // Overflow-safe: if count * elem_size would overflow, treat as truncated
            if (elem_size != 0 && count > SIZE_MAX / elem_size) {
                return e;
            }
            size_t bytes = size_t(count) * elem_size;

            const uint8_t* r = q + bytes;
            return (r <= e) ? r : e;
        }

        default:
            return e;
    }
}

void BamAlignmentParser::ParseBamAlignmentFromSliceComplexTypes(
    std::string_view bytes,
    SamAlignments &reads,
    size_t idx,
    VariableLenFieldsVar const& fieldsOff,
    VariableLenFieldsVar const& fieldsSizes,
    VariableLenFieldsVar const& outOff)
{
    ByteScanner bs; bs.reset(bytes);

    bs.setIndex(fieldsOff.qname);

    const uint32_t qname_bytes_including_nul = fieldsSizes.qname; // BAM l_read_name
    const size_t   qname_src_off = bs.getIndex();

    if (qname_bytes_including_nul) {
        const char *qsrc = bytes.data() + qname_src_off;

        // BAM requires a trailing NUL; be defensive if input is malformed.
        size_t qname_len = qname_bytes_including_nul;
        if (qname_len && qsrc[qname_len - 1] == '\0') {
            qname_len -= 1; // drop the NUL
        }

        // Copy exactly the visible characters of QNAME
        char* dst = reads.QNAME.Data() + outOff.qname;
        std::memcpy(dst, qsrc, qname_len);

        // Store span (offset + length without NUL)
        reads.QNAME.Ptrs()[idx] = StringField {
            outOff.qname,
            static_cast<uint16_t>(qname_len)
        };
    } else {
        // Empty name (shouldn't happen in valid BAM, but handle gracefully)
        reads.QNAME.Ptrs()[idx] = StringField{ outOff.qname, 0 };
    }

    reads.TILE[idx] = computeTileInfo(reads.GetQNAME(idx));

    // ---- CIGAR: decode directly into final buffer ----
    bs.setIndex(fieldsOff.cigar);
    CigarOperation* dst_cig = reads.CIGAR.Data() + outOff.cigar;
    for (uint32_t i = 0; i < fieldsSizes.cigar; ++i) {
        uint32_t op = bs.readInt32();
        dst_cig[i] = CigarOperation{
            (op >> 4),
            cigarOps[int(op & 0x0F)]
        };
    }
    reads.CIGAR.Ptrs()[idx] = SpanField{ outOff.cigar, static_cast<uint16_t>(fieldsSizes.cigar) };

    // ---- SEQ: decode packed nibbles directly into final buffer ----
    bs.setIndex(fieldsOff.seq);
    char* dst_seq = reads.SEQ.Data() + outOff.seq;
    const size_t l_seq    = fieldsSizes.seq;
    const size_t seqBytes = (l_seq + 1u) >> 1;
    size_t k = 0;
    for (size_t i = 0; i < seqBytes; ++i) {
        uint8_t nibbles = bs.readByte();
        dst_seq[k++] = seqEncoder[nibbles >> 4];
        if (k < l_seq) dst_seq[k++] = seqEncoder[nibbles & 0x0F];
    }
    reads.SEQ.Ptrs()[idx] = SpanField{ outOff.seq, static_cast<uint16_t>(l_seq) };

    // ---- QUAL: memcpy directly to final buffer ----
    bs.setIndex(fieldsOff.qual);
    std::memcpy(reads.QUAL.Data() + outOff.qual, bytes.data() + bs.getIndex(), l_seq);
    reads.QUAL.Ptrs()[idx] = SpanField{ outOff.qual, static_cast<uint16_t>(l_seq) };


    // ---- AUX: fixed-offset copy into arena (no insert/resize) ----
    OptionalFields& opt = reads.TAGS;
        if (fieldsSizes.aux) {
        std::memcpy(
            opt.tagArena.data() + outOff.aux,
            bytes.data() + fieldsOff.aux,
            fieldsSizes.aux
        );
        opt.tagsBlobOff[idx]  = outOff.aux;
        opt.tagsBlobSize[idx] = fieldsSizes.aux;

        // ---- RG extraction (Z string): map to int via header ----
        const uint8_t* p = reinterpret_cast<const uint8_t*>(opt.tagArena.data()) + outOff.aux;
        const uint8_t* e = p + fieldsSizes.aux;

        int rg_int = -1; // or your sentinel for "no RG"
        while (p + 3 <= e) {
            const char k1 = char(p[0]);
            const char k2 = char(p[1]);
            p += 2;
            const char T  = (p < e) ? char(*p++) : '\0';
            if (k1 == 'R' && k2 == 'G') {
                if (T == 'Z') {
                    const uint8_t* z = p;
                    while (z < e && *z != 0) ++z;
                    const size_t len = size_t(z - p);

                    std::string rg_id{ reinterpret_cast<const char*>(p), len };
                    rg_int = reads.GetHdr().getReadGroupRecord().getId(rg_id);

                    // advance past the NUL (if any)
                    p = (z < e) ? (z + 1) : e;
                    break; // assume at most one RG
                } else {
                    p = skip_payload(p, e, T);
                }
            } else {
                p = skip_payload(p, e, T);
            }
        }

        // store integer RG for this record
        opt.RG[idx] = rg_int;

    } else {
        opt.tagsBlobOff[idx]  = std::numeric_limits<size_t>::max();
        opt.tagsBlobSize[idx] = 0;
        opt.RG[idx] = -1;
    }
}