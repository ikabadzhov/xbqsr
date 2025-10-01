#include "BamAlignmentParser.h"
#include "../utils/ByteScanner.h"
#include "../../timer/ProfilerIO.h"

// Parse exactly one BAM alignment from a complete record body (i.e., block_size bytes).
// Precondition: 'bytes' is exactly the block body (NOT including the 4-byte block_size field).
void BamAlignmentParser::ParseBamAlignmentFromSlice(std::string_view bytes,
                                                    SamAlignments& reads,
                                                    size_t idx)
{
    //auto& profiler = ProfilerIO::instance();

    // -------------------- PARSING --------------------
    //profiler.start(IOStageId::Parsing);

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

    //profiler.stop(IOStageId::Parsing);

    // ---------------- SAM STRUCT UPDATES ----------------
    //profiler.start(IOStageId::SAMStructUpdates);

    reads.RNAME[idx] = (refID < 0) ? SequenceRecord::RESERVED_SEQUENCE_NAMES['*'] : refID;
    if (next_refID < 0) next_refID = SequenceRecord::RESERVED_SEQUENCE_NAMES['*'];
    reads.RNEXT[idx] = next_refID;

    reads.POS[idx]   = pos0 + 1;
    reads.PNEXT[idx] = next_pos0 + 1;
    reads.MAPQ[idx]  = mapq;
    reads.FLAG[idx]  = flag;
    reads.TLEN[idx]  = tlen;

    if (!qname.empty())
        reads.QNAME[idx] = reads.GetHdr().QnameMapper().addString(qname);
    reads.TILE[idx] = computeTileInfo(qname);

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
    //profiler.stop(IOStageId::SAMStructUpdates);
}

void BamAlignmentParser::ParseBamAlignmentFromSliceSimpleTypes(std::string_view bytes, SamAlignments &reads, size_t idx, VariableLenFieldsVar& fieldsOff, VariableLenFieldsVar& fieldsSizes)
{
//    auto& profiler = ProfilerIO::instance();

    // -------------------- PARSING --------------------
    ////profiler.start(IOStageId::Parsing);

    ByteScanner bs;
    bs.reset(bytes);

    // Core fixed fields
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

    // QNAME
    std::string qname = bs.readString(l_read_name);
    if (!qname.empty() && qname.back() == '\0') qname.pop_back();

    // CIGAR (n_cigar_op 32-bit words)
    const size_t cigar_off   = bs.getIndex();
    const size_t cigar_bytes = size_t(n_cigar_op) * sizeof(uint32_t);
    if (!bs.canRead(cigar_bytes)) throw std::runtime_error("Truncated CIGAR");
    fieldsOff.cigar   = cigar_off;      // byte offset within this record
    fieldsSizes.cigar = n_cigar_op;     // number of ops (not bytes)
    bs.incrementIndex(cigar_bytes);

    // SEQ (ceil(l_seq/2) bytes)
    const size_t seq_bytes = (size_t(l_seq) + 1u) >> 1;
    const size_t seq_off   = bs.getIndex();
    if (!bs.canRead(seq_bytes)) throw std::runtime_error("Truncated SEQ");
    fieldsOff.seq = seq_off;       // byte offset of packed sequence
    fieldsSizes.seq = l_seq;         // number of bases after decode
    bs.incrementIndex(seq_bytes);

    // QUAL (l_seq bytes)
    const size_t qual_off = bs.getIndex();
    if (!bs.canRead(size_t(l_seq))) throw std::runtime_error("Truncated QUAL");
    fieldsOff.qual   = qual_off;        // byte offset of qualities
    fieldsSizes.qual   = l_seq;           // number of qualities
    bs.incrementIndex(size_t(l_seq));

    // (optional) AUX start/len for raw-tag handling:

    const size_t aux_off  = static_cast<size_t>(bs.getIndex());
    const size_t aux_size = static_cast<size_t>(bytes.size()) - aux_off;

    fieldsOff.aux = aux_off;
    fieldsSizes.aux = aux_size; // safe since we validated earlier fields


    //profiler.stop(IOStageId::Parsing);

    // ---------------- SAM STRUCT UPDATES + TAG HANDLING ----------------
    //profiler.start(IOStageId::SAMStructUpdates);

    // Simple SoA fields
    reads.RNAME[idx] = (refID < 0) ? SequenceRecord::RESERVED_SEQUENCE_NAMES['*'] : refID;
    if (next_refID < 0) next_refID = SequenceRecord::RESERVED_SEQUENCE_NAMES['*'];
    reads.RNEXT[idx] = next_refID;

    reads.POS[idx]   = pos0 + 1;
    reads.PNEXT[idx] = next_pos0 + 1;
    reads.MAPQ[idx]  = mapq;
    reads.FLAG[idx]  = flag;
    reads.TLEN[idx]  = tlen;

    // If you want to skip QNAME interning here:
    if (!qname.empty()) {
        // reads.QNAME[idx] = reads.GetHdr().QnameMapper().addString(qname);
        reads.QNAME[idx] = 0; // as in your snippet
    }
    reads.TILE[idx] = computeTileInfo(qname);

    // ---- Optional fields (RAW TAG BLOB + RG extraction) ----
    // Adjust this to how you access your OptionalFields storage.
    OptionalFields& opt = reads.TAGS;  // <--- rename to your actual handle

    // ensure capacity (if you don't pre-size externally)
    if (opt.tagsBlobOff.size() <= idx) {
        opt.tagsBlobOff.resize(idx + 1, size_t{std::numeric_limits<size_t>::max()});
        opt.tagsBlobSize.resize(idx + 1, 0);
        // RG column should already be sized outside; otherwise:
        opt.RG.Resize(idx + 1);
    }

    if (aux_size == 0) {
        // No tags
        opt.tagsBlobOff[idx]  = std::numeric_limits<size_t>::max();
        opt.tagsBlobSize[idx] = 0;
    } else {
        // Save the raw tag blob
        const size_t arena_off = opt.tagArena.size();
        /*opt.tagArena.insert(opt.tagArena.end(),
                            reinterpret_cast<const std::byte*>(bytes.data() + aux_off),
                            reinterpret_cast<const std::byte*>(bytes.data() + aux_off + aux_size));*/
        opt.tagsBlobOff[idx]  = arena_off;
        opt.tagsBlobSize[idx] = aux_size;

        // Extract only RG (Z) from the just-copied blob
        const uint8_t* p = reinterpret_cast<const uint8_t*>(opt.tagArena.data()) + arena_off;
        const uint8_t* e = p + aux_size;

        auto skip_payload = [](const uint8_t* q, const uint8_t* e, char T) -> const uint8_t* {
            switch (T) {
                case 'A': case 'c': case 'C': return q + 1;
                case 's': case 'S':           return q + 2;
                case 'i': case 'I': case 'f': return q + 4;
                case 'Z': case 'H': {
                    // NUL-terminated
                    const uint8_t* r = q;
                    while (r < e && *r != 0) ++r;
                    return (r < e) ? (r + 1) : e; // include NUL if present
                }
                case 'B': {
                    if (q + 5 > e) return e;   // subtype + count
                    const char subtype = char(*q++);
                    (void)subtype;
                    int32_t count = 0;
                    std::memcpy(&count, q, 4); q += 4;
                    size_t elem = 0;
                    switch (subtype) {
                        case 'c': case 'C': elem = 1; break;
                        case 's': case 'S': elem = 2; break;
                        default:            elem = 4; break; // i,I,f
                    }
                    size_t bytes = size_t(count) * elem;
                    return (q + bytes <= e) ? (q + bytes) : e;
                }
                default: return e; // unknown → bail to end to avoid desync
            }
        };

        // scan for RG:Z
        while (p + 3 <= e) {
            const char k1 = char(p[0]);
            const char k2 = char(p[1]);
            p += 2;
            const char T  = char(*p++);
            if (k1 == 'R' && k2 == 'G') {
                if (T == 'Z') {
                    // Grab the Z string (without the trailing NUL)
                    const uint8_t* z = p;
                    while (z < e && *z != 0) ++z;
                    const size_t len = (z <= e) ? size_t(z - p) : size_t(0);

                    // Store RG as compressed/interned string (your CompressedStringCollection)
                    // I assume SetAt takes (std::string_view, idx). Adjust to your API:
                    //opt.RG.Set(idx, {reinterpret_cast<const char*>(p), len});

                    // Skip the Z payload (including NUL if it exists)
                    p = (z < e) ? (z + 1) : e;
                    // You can break here if you assume at most one RG per record
                    break;
                } else {
                    // RG exists but not a Z? Skip payload safely
                    p = skip_payload(p, e, T);
                }
            } else {
                // Not RG -> skip
                p = skip_payload(p, e, T);
            }
        }
    }

    //profiler.stop(IOStageId::SAMStructUpdates);
}

void BamAlignmentParser::ParseBamAlignmentFromSliceComplexTypes(std::string_view bytes, SamAlignments &reads, size_t idx, VariableLenFieldsVar& fieldsOff, VariableLenFieldsVar& fieldsSizes, VariableLenFieldsVar& outputOff)
{
    //auto& profiler = ProfilerIO::instance();

    // -------------------- PARSING --------------------
    //profiler.start(IOStageId::Parsing);

    ByteScanner bs;
    bs.reset(bytes);

    bs.setIndex(fieldsOff.name);

    size_t l_name = fieldsSizes.name;
    size_t l_seq = fieldsSizes.seq;
    size_t l_qual = fieldsSizes.qual;
    size_t l_cigar = fieldsSizes.cigar;

    std::string qname = bs.readString(l_name);
    if (!qname.empty() && qname.back() == '\0') qname.pop_back();

    bs.setIndex(fieldsOff.cigar);
    std::vector<CigarOperation> cigarVector(l_cigar);
    for (uint32_t i = 0; i < l_cigar; ++i) {
        uint32_t op = bs.readInt32();
        cigarVector[i] = CigarOperation{
            (op >> 4),
            cigarOps[static_cast<int>(op & 0x0F)]
        };
    }

    bs.setIndex(fieldsOff.seq);
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

    assert(l_seq == l_qual);
    bs.setIndex(fieldsOff.qual);
    std::vector<char> qual(static_cast<size_t>(l_seq));
    if (l_seq) {
        const char* src = bytes.data() + bs.getIndex();
        std::memcpy(qual.data(), src, static_cast<size_t>(l_seq));
        bs.incrementIndex(static_cast<size_t>(l_seq));
    }

    //profiler.stop(IOStageId::Parsing);

    // ---------------- SAM STRUCT UPDATES ----------------
    //profiler.start(IOStageId::SAMStructUpdates);

    if (!qname.empty())
        reads.QNAME[idx] = 0;//reads.GetHdr().QnameMapper().addString(qname);

    std::memcpy(reads.CIGAR.Data() + outputOff.cigar, cigarVector.data(), sizeof(CigarOperation)*l_cigar);
    std::memcpy(reads.SEQ.Data() + outputOff.seq, decoded_seq.data(), l_seq);
    std::memcpy(reads.QUAL.Data() + outputOff.qual, qual.data(), l_qual);
    reads.CIGAR.Ptrs()[idx] = SpanField{outputOff.cigar, static_cast<uint16_t>(l_cigar)};
    reads.SEQ.Ptrs()[idx] = SpanField{outputOff.seq, static_cast<uint16_t>(l_seq)};
    reads.QUAL.Ptrs()[idx] = SpanField{outputOff.qual, static_cast<uint16_t>(l_qual)};
    // reads.CIGAR.Data().SetAt(cigarVector, idx);
    // reads.SEQ.SetAt(decoded_seq, idx);
    // reads.QUAL.SetAt(qual, idx);

    //profiler.stop(IOStageId::SAMStructUpdates);
}