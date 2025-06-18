#include "BamAlignmentParser.h"
#include "../utils/ByteScanner.h"


/*
 * This function is made in order to check if the buffer just read in BamUtilsIO::ParseSamAlignment
 * contains the EOF string specified in the SAM/BAM docs:
 * 1f 8b 08 04 00 00 00 00 00 ff 06 00 42 43 02 00 1b 00 03 00 00 00 00 00 00 00 00 00
 *
 * In some files, the end of the file can be found when all the sequence_size string are all '\0'.
 */
bool CheckEOFString(std::string_view buffer, size_t buffer_size, unsigned char* sequence, size_t sequence_size){
    const std::vector<uint8_t> uncompressedSEQ = {
            0x00, 0x00, 0x00, 0x00, 0x00, 0xFF, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00, 0x1B, 0x00,
            0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00
    };

    if (buffer.size() < uncompressedSEQ.size())
        return false; //Not enough data, the buffer wasn't able to fill up, so we finished reading the file

    for (size_t i = 0; i < uncompressedSEQ.size(); ++i) {
        if (static_cast<unsigned char>(uncompressedSEQ[i]) != static_cast<unsigned char>(buffer[i]))
            return false; // Mismatch
    }

    return true; // Matched

    /*for (size_t i = 0; i <= buffer_size - sequence_size; ++i) {
        bool found = true;
        for (size_t j = 0; j < sequence_size; ++j) {
            if (buffer[i + j] != sequence[j]) {
                found = false;
                break;
            }
        }
        if (found) {
            std::cout << "EOF sequence found in the buffer.\n";
            return true;
        }
    }

    return false;*/
}

/*
 * This function is designed to extract a single alignment from the BAM file. Upon reading the block, determined by the
 * uint32_t block_size field, it becomes crucial to refer to the documentation to comprehend the precise reading
 * methodology concerning the information's length and type. Within this process, strict adherence to the specifications
 * outlined in the documentation is crucial to avoid misinterpretation of data, such as encountering conflicting CIGAR
 * sequences. Moreover, the method respects the information indicating the lengths of subsequent fields to ensure
 * accurate data extraction. For instance, the POS (position) information, denoting the 0-based leftmost coordinate,
 * requires adjustment by subtracting one from the NEXT_POS value to align with the actual value read.
 **/
void BamAlignmentParser::ParseBamAlignmentFromSlice(std::string_view bytes, SamAlignments &reads, size_t idx) {
    ByteScanner bs;

    if (CheckEOFString(bytes, bytes.size(), const_cast<unsigned char *>(knownSequence), sizeof(knownSequence)))
        throw std::runtime_error("EOF condition encountered"); //return;
    bs.reset(bytes);

    // The RefID must be either STAR (unknown) or −1 ≤ refID < n ref
    int32_t temp_refID = bs.readInt32();
    if (temp_refID < 0) {
        reads.RNAME[idx] = SequenceRecord::RESERVED_SEQUENCE_NAMES['*'];
    }
    else {
        reads.RNAME[idx] = temp_refID;
    }

    //The RefID must be either STAR (unknown) or −1 ≤ refID < n ref
    /*if (temp_refID != std::to_string(42) && (RNAME_read_value < std::to_string(-1) &&
                                                   RNAME_read_value > std::to_string(reads.GetHdr().number_of_reference_sequences))) {
        std::cout << "refID not valid" << std::endl;
        return;
    }*/

    int32_t temp_pos = bs.readInt32();

    temp_pos += 1;
    reads.POS[idx] = temp_pos;

    uint8_t temp_lenght_read_name = bs.readByte();

    uint8_t temp_mapq = bs.readByte();

    reads.MAPQ[idx] = temp_mapq;

    //uint16_t temp_bin = bs.readInt16();

    // Note unmapped reads with POS 0 (which becomes −1 in BAM) therefore use reg2bin(−1, 0) which is computed as 4680
    //temp_bin -= 1; // Unused at the moment

    uint16_t temp_n_cigar_op = bs.readUint16();

    uint16_t temp_flag = bs.readUint16();
    reads.FLAG[idx] = temp_flag;

    uint32_t temp_length_seq = bs.readUint32();

    int32_t temp_next_refID = bs.readInt32();
    if (temp_next_refID < 0) {
        temp_next_refID = SequenceRecord::RESERVED_SEQUENCE_NAMES['*'];
    }
    /*
     * BAM and SAM parsing is different for this field. IF RNAME is "chrX" and the RNEXT is the same,
     * in SAM it is parsed with the sign "=" while in BAM it is the string is duplicated, so the final
     * value of the field will be "chrX". Gotta be careful when formatting.
    */
    reads.RNEXT[idx] = temp_next_refID;


    //The RefID must be either STAR (unknown) or −1 ≤ refID < n ref
    /*if (RNEXT_read_value != std::to_string(42) && (RNEXT_read_value < std::to_string(-1) &&
                                                   RNEXT_read_value > std::to_string(reads.GetHdr().number_of_reference_sequences))) {
        std::cout << "RNEXT not valid" << std::endl;
        return;
    }*/

    // 0-based leftmost coordinate, the real value is NEXT_POS-1
    int32_t temp_next_pos = bs.readInt32();
    temp_next_pos += 1;
    reads.PNEXT[idx] = temp_next_pos;

    int32_t temp_tlen = bs.readInt32();
    reads.TLEN[idx] = temp_tlen;

    std::string qname = bs.readString(temp_lenght_read_name);
    if (temp_lenght_read_name - 1 > 0) {
        reads.QNAME[idx] = reads.GetHdr().QnameMapper().addString(qname.substr(0, qname.size()-1));
    }
    reads.TILE[idx] = computeTileInfo(qname);

    //int index = readNameIndex + temp_lenght_read_name;
    std::vector<CigarOperation> cigarVector;
    cigarVector.resize(temp_n_cigar_op);

    for (uint32_t i = 0; i < temp_n_cigar_op; ++i) {
        uint32_t cigar = bs.readInt32();

        CigarOperation CigarIth = CigarOperation{
                static_cast<uint32_t>(cigar >> 4),
                cigarOps[static_cast<int>(0xF & cigar)]
        };
        cigarVector[i] = CigarIth;
    }
    //reads.CIGAR[idx] = cigarVector;
    reads.CIGAR.SetAt(cigarVector, idx);

    /* TRANSFORM VERSION OF NIBBLE SEQUENCE.
     * This version aims to read only once from the memory, reducing the (n) IO calls from the file.
     * The drawback found in VTune is that even though the complessity of the overall function drastically decreases,
     * the untracked (system and libz.so library ones) functions, such as func@0x8e44 requires way more CPU time
     * to complete their execution, adding up towards a greater elapsed time. This issue doesn't find a reference on the
     * QUAL vector parsing. Finally, the best combination is original version of nibble sequence parsing and transform
     * method applied to QUAL.
     * std::string decoded_sequence;
    */

    //ORIGINAL VERSION OF SEQUENCE PARSING
    std::string decoded_sequence;
    std::string nibble_sequence(temp_length_seq, 0);

    int32_t k = 0;
    for (size_t i = 0; i < temp_length_seq; i += 2) {
        uint8_t nibbles = bs.readByte();

        uint8_t first_nibble = nibbles >> 4;        //first 4 bits
        uint8_t second_nibble = nibbles & 0x0F;     //last 4 bits

        nibble_sequence[k] = seqEncoder[first_nibble]; // first nucleotide
        k++;
        if (i < temp_length_seq - 1) {
            nibble_sequence[k] = seqEncoder[second_nibble];// second nucleotide
            k++;
        }
    }
    reads.SEQ.SetAt(nibble_sequence, idx);

    std::vector<char> qualVector;
    qualVector.resize(temp_length_seq);

    //ORIGINAL VERSION OF QUALVECTOR PARSING
    std::transform(
            bytes.data() + bs.getIndex(),
            bytes.data() + bs.getIndex() + temp_length_seq,
            qualVector.begin(),
            [](char letter) {
                return letter;
            }
    );
    bs.incrementIndex(temp_length_seq);

    /*
     * Despite the memory efficiency of std::span, std::vector is chosen in this scenario.
     * This decision is due to the characteristics of std::span, which, although not owning the referenced data,
     * serves as a non-owning view over a continuous sequence of elements.
     * The usage of std::vector was necessitated because direct initialization with std::span requires an underlying
     * container or array-like structure.
     * Additional details are available in the documentation: https://en.cppreference.com/w/cpp/container/span.
    */
    //std::span<char> qualSpan(qualVector.data(), qualVector.size());

    reads.QUAL.SetAt(qualVector, idx);

    TagKeyAndType tagKeyAndType;
    //int64_t tagValueInt;
    //float tagValueFloat;
    std::vector<int32_t> tagValueArrayInt32;
    std::vector<uint32_t> tagValueArrayUInt32;
    std::vector<int16_t> tagValueArrayInt16;
    std::vector<uint16_t> tagValueArrayUInt16;
    std::vector<int8_t> tagValueArrayInt8;
    std::vector<uint8_t> tagValueArrayUInt8;
    std::vector<float> tagValueArrayFloat;

    while (bs.canRead(2)) {
        uint32_t elementCount = 0;
        tagKeyAndType.setTagKey(TagKeyAndType::encode(bs.readString(2)));
        tagKeyAndType.setTagType(static_cast<char>(bs.readByte()));

        if (tagKeyAndType.getTagType() == 'A') {
            char tagValue = bs.readByte();
            reads.TAGS.InsertCharTagAt(idx, tagKeyAndType, tagValue);
        } else if (tagKeyAndType.getTagType() == 'c' || tagKeyAndType.getTagType() == 'C') {
            reads.TAGS.InsertIntTagAt(idx, tagKeyAndType, static_cast<int8_t>(bs.readByte()));
        } else if (tagKeyAndType.getTagType() == 's' || tagKeyAndType.getTagType() == 'S') {
            reads.TAGS.InsertIntTagAt(idx, tagKeyAndType, bs.readInt16());
        } else if (tagKeyAndType.getTagType() == 'i' || tagKeyAndType.getTagType() == 'I') {
            reads.TAGS.InsertIntTagAt(idx, tagKeyAndType, bs.readInt32());
        } else if (tagKeyAndType.getTagType() == 'f') {
            reads.TAGS.InsertFloatTagAt(idx, tagKeyAndType, bs.readFloat());
        } else if (tagKeyAndType.getTagType() == 'Z' || tagKeyAndType.getTagType() == 'H') {
            reads.TAGS.InsertStringTagAt(idx, tagKeyAndType, bs.readUntil(0));
        } else if (tagKeyAndType.getTagType() == 'B') {
            tagKeyAndType.setSecondaryTagType(static_cast<char>(bs.readByte()));
            elementCount = bs.readByte();
            if (tagKeyAndType.getSecondaryTagType() == 'c') {
                tagValueArrayInt8.resize(elementCount);
                memcpy(tagValueArrayInt8.data(), bs.data(), elementCount * sizeof(int8_t));
                reads.TAGS.InsertArrayInt8TagAt(idx, tagKeyAndType, tagValueArrayInt8);
            } else if (tagKeyAndType.getSecondaryTagType() == 'C') {
                tagValueArrayUInt8.resize(elementCount);
                memcpy(tagValueArrayUInt8.data(), bs.data(), elementCount * sizeof(uint8_t));
                reads.TAGS.InsertArrayUint8TagAt(idx, tagKeyAndType, tagValueArrayUInt8);
            } else if (tagKeyAndType.getSecondaryTagType() == 's') {
                tagValueArrayInt16.resize(elementCount);
                memcpy(tagValueArrayInt16.data(), bs.data(), elementCount * sizeof(int16_t));
                reads.TAGS.InsertArrayInt16TagAt(idx, tagKeyAndType, tagValueArrayInt16);
            } else if (tagKeyAndType.getSecondaryTagType() == 'S') {
                tagValueArrayUInt16.resize(elementCount);
                memcpy(tagValueArrayUInt16.data(), bs.data(), elementCount * sizeof(uint16_t));
                reads.TAGS.InsertArrayUint16TagAt(idx, tagKeyAndType, tagValueArrayUInt16);
            } else if (tagKeyAndType.getSecondaryTagType() == 'i') {
                tagValueArrayInt32.resize(elementCount);
                memcpy(tagValueArrayInt32.data(), bs.data(), elementCount * sizeof(int32_t));
                reads.TAGS.InsertArrayInt32TagAt(idx, tagKeyAndType, tagValueArrayInt32);
            } else if (tagKeyAndType.getSecondaryTagType() == 'I') {
                tagValueArrayUInt32.resize(elementCount);
                memcpy(tagValueArrayUInt32.data(), bs.data(), elementCount * sizeof(uint32_t));
                reads.TAGS.InsertArrayUint32TagAt(idx, tagKeyAndType, tagValueArrayUInt32);
            } else if (tagKeyAndType.getSecondaryTagType() == 'f') {
                tagValueArrayFloat.resize(elementCount);
                memcpy(tagValueArrayFloat.data(), bs.data(), elementCount * sizeof(float));
                reads.TAGS.InsertArrayFloatTagAt(idx, tagKeyAndType, tagValueArrayFloat);
            }
        }
    }

}

