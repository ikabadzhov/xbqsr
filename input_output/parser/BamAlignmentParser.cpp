#include "BamAlignmentParser.h"
#include "../utils/ByteScanner.h"

// Parse exactly one BAM alignment from a complete record body (i.e., block_size bytes).
// Precondition: 'bytes' is exactly the block body (NOT including the 4-byte block_size field).
void BamAlignmentParser::ParseBamAlignmentFromSlice(std::string_view bytes,
                                                    SamAlignments& reads,
                                                    size_t idx)
{
    ByteScanner bs;
    bs.reset(bytes);

    // --- Core fixed fields (see BAM spec) ---
    int32_t  refID        = bs.readInt32();     // reference id
    int32_t  pos0         = bs.readInt32();     // 0-based
    uint8_t  l_read_name  = bs.readByte();
    uint8_t  mapq         = bs.readByte();
    uint16_t bin_mq_nl    = bs.readInt16();     // bin (hi 16 when split); we don't use it here
    (void)bin_mq_nl;
    uint16_t n_cigar_op   = bs.readUint16();
    uint16_t flag         = bs.readUint16();
    uint32_t l_seq        = bs.readUint32();
    int32_t  next_refID   = bs.readInt32();
    int32_t  next_pos0    = bs.readInt32();
    int32_t  tlen         = bs.readInt32();

    // Map ref IDs
    if (refID < 0)
        reads.RNAME[idx] = SequenceRecord::RESERVED_SEQUENCE_NAMES['*'];
    else
        reads.RNAME[idx] = refID;

    if (next_refID < 0)
        next_refID = SequenceRecord::RESERVED_SEQUENCE_NAMES['*'];
    reads.RNEXT[idx] = next_refID;

    // POS/PNEXT: BAM stores 0-based; if your internal is 1-based, +1 is correct.
    reads.POS[idx]  = pos0 + 1;
    reads.PNEXT[idx]= next_pos0 + 1;

    reads.MAPQ[idx] = mapq;
    reads.FLAG[idx] = flag;
    reads.TLEN[idx] = tlen;

    // --- QNAME (includes trailing NUL), then trim for storage ---
    std::string qname = bs.readString(l_read_name);
    if (!qname.empty() && qname.back() == '\0') qname.pop_back();
    if (!qname.empty())
        reads.QNAME[idx] = reads.GetHdr().QnameMapper().addString(qname);
    reads.TILE[idx] = computeTileInfo(qname);

    // --- CIGAR: n_cigar_op x 4 bytes ---
    std::vector<CigarOperation> cigarVector(n_cigar_op);
    for (uint32_t i = 0; i < n_cigar_op; ++i) {
        uint32_t op = bs.readInt32();
        cigarVector[i] = CigarOperation{
                static_cast<uint32_t>(op >> 4),            // length
                cigarOps[static_cast<int>(op & 0x0F)]      // op code
        };
    }
    reads.CIGAR.SetAt(cigarVector, idx);

    // --- SEQ: ceil(l_seq/2) bytes, high nibble first ---
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
            if (k < decoded_seq.size())
                decoded_seq[k++] = seqEncoder[second];
        }
    }
    reads.SEQ.SetAt(decoded_seq, idx);

    // --- QUAL: l_seq bytes (0xFF means missing). Copy directly. ---
    {
        std::vector<char> qual(static_cast<size_t>(l_seq));
        if (l_seq) {
            const char* src = bytes.data() + bs.getIndex();
            std::memcpy(qual.data(), src, static_cast<size_t>(l_seq));
            bs.incrementIndex(static_cast<size_t>(l_seq));
        }
        reads.QUAL.SetAt(qual, idx);
    }

    // --- Auxiliary tags (TLV stream) ---
    TagKeyAndType tagKeyAndType;
    while (bs.canRead(3)) { // key(2) + type(1)
        tagKeyAndType.setTagKey(TagKeyAndType::encode(bs.readString(2))); // 2-char key
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
            // NUL-terminated string (H is hex text)
            reads.TAGS.InsertStringTagAt(idx, tagKeyAndType, bs.readUntil(0));

        } else if (T == 'B') {
            // Array: subtype (1), count (int32 LE), then raw data
            const char subtype = static_cast<char>(bs.readByte());
            tagKeyAndType.setSecondaryTagType(subtype);

            const int32_t count = bs.readInt32(); // **int32**, not a byte
            if (count < 0) throw std::runtime_error("Negative B-array count");

            size_t elem_size = 0;
            switch (subtype) {
                case 'c': case 'C': elem_size = 1; break;
                case 's': case 'S': elem_size = 2; break;
                case 'i': case 'I': case 'f': elem_size = 4; break;
                default:
                    throw std::runtime_error("Unknown B-array subtype");
            }
            const size_t nbytes = static_cast<size_t>(count) * elem_size;
            if (!bs.canRead(nbytes)) throw std::runtime_error("Truncated B-array payload");

            const char* p = bs.data(); // pointer to payload start

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

            bs.incrementIndex(nbytes); // **advance after memcpy**
        } else {
            // Unknown tag type -> bail (avoid desync)
            throw std::runtime_error("Unknown BAM tag type");
        }
    }
}