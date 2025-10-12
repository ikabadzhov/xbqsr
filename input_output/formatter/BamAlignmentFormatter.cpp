#include "BamAlignmentFormatter.h"
#include "../utils/ByteFormatter.h"


void BamAlignmentFormatter::FormatBamAlignment(SamAlignments &alns, ssize_t idx, std::string &formatted_line) {
    //Fixed size for the first 12 fields + variable sizes
    //int32_t blockSz = 32 + alns.GetHdr().QnameMapper().GetAt(alns.GetQNAME(idx)).size() + sizeof('\0');
    //blockSz += sizeof(int32_t) * alns.GetCIGAR(idx).size();
    //blockSz += sizeof(int8_t) * ((alns.GetSEQ(idx).size() + 1)/2);
    //blockSz += sizeof(char) * alns.GetQUAL(idx).size();

    ByteFormatter bf(formatted_line);
    //Field: Block size (uint32_t ) 1
    bf.appendInt32((uint32_t ) 1); //it will be replaced when the real value is defined
    //Field: redID
    bf.appendInt32(alns.GetRNAME(idx));
    //Field: pos
    bf.appendInt32(alns.GetPOS(idx)-1);
    bf.appendByte(alns.GetHdr().QnameMapper().GetAt(alns.GetQNAME(idx)).size() + 1); //?
    //Field: mapq
    bf.appendByte(alns.GetMAPQ(idx));
    //Field: bin, we don't use it
    bf.appendInt16(-1);
    //Field: n_cigar_op, to fix
    bf.appendInt16(alns.GetCIGAR(idx).size());
    //Field: flag
    bf.appendInt16(alns.GetFLAG(idx));
    //Field: l_seq
    bf.appendInt32(alns.GetSEQ(idx).size());
    //Field: next_refID
    bf.appendInt32(alns.GetRNEXT(idx));
    //Field: next_pos
    bf.appendInt32(alns.GetPNEXT(idx) - 1);
    //Field: tlen
    bf.appendInt32(alns.GetTLEN(idx));
    //Field: read_name
    //alns.GetQNAME(idx) returns the ID of the QNAME while this whole line gets the final name [read_name]
    bf.appendString(alns.GetHdr().QnameMapper().GetAt(alns.GetQNAME(idx)));
    bf.appendByte('\0');

    //Field: Cigar
    std::string cigarStr;
    //cigarToString(alns.GetCIGAR(idx), cigarStr);
    for(auto x:alns.GetCIGAR(idx)){
        bf.appendInt32(cigarStringToValue(std::to_string(x.Length) + x.Operation));
    }

    //Field: Seq
    std::string val((alns.GetSEQ(idx).size() + 1)/2, 0);

    for (size_t i = 0; i < alns.GetSEQ(idx).size(); i += 2) {
        bf.appendByte(encodeNucleotides(alns.GetSEQ(idx).data()[i], alns.GetSEQ(idx).data()[i+1]));
    }

    //Field: Qual
    val.resize(alns.GetQUAL(idx).size());
    for (size_t i = 0; i < alns.GetQUAL(idx).size(); i++) {
       val[i] = alns.GetQUAL(idx)[i];
    }
    bf.appendString(val);

    for (auto& tag : alns.TAGS.ALL_TAGS){

        //uint16_t tagKey = tag.second.getTagKey();
        std::string tagKeyAsString = tag.second.getTagKeyAsString();
        char tagType = tag.second.getTagType();
        char tagSecondaryType = tag.second.getSecondaryTagType();

        if (tagKeyAsString == "RG") { // TODO: Hide this inside the OptionalField class
            //bf.appendString(alns.TAGS.RG.toString(idx));
        } else {
            if (tagType == 'A') {
                if (alns.TAGS.GetCharTags(tagKeyAsString).isValid(idx)) {
                    bf.appendString(tagKeyAsString);
                    bf.appendByte(tagType);
                    bf.appendByte(static_cast<char>(alns.TAGS.GetCharTags(tagKeyAsString).GetAt(idx)));
                }
            } else if (tagType == 'c') {
                if (alns.TAGS.GetIntTags(tagKeyAsString).isValid(idx)) {
                    bf.appendString(tagKeyAsString);
                    bf.appendByte(tagType);
                    bf.appendByte(static_cast<int8_t>(alns.TAGS.GetIntTags(tagKeyAsString).GetAt(idx)));
                }
            } else if (tagType == 'C') {
                if (alns.TAGS.GetIntTags(tagKeyAsString).isValid(idx)) {
                    bf.appendString(tagKeyAsString);
                    bf.appendByte(tagType);
                    uint8_t newVal = static_cast<uint8_t>(alns.TAGS.GetIntTags(tagKeyAsString).GetAt(idx));
                    bf.appendByte(newVal);
                }
            } else if (tagType == 's') {
                if (alns.TAGS.GetIntTags(tagKeyAsString).isValid(idx)) {
                    bf.appendString(tagKeyAsString);
                    bf.appendByte(tagType);
                    bf.appendInt16(static_cast<int16_t>(alns.TAGS.GetIntTags(tagKeyAsString).GetAt(idx)));
                }
            } else if (tagType == 'S') {
                if (alns.TAGS.GetIntTags(tagKeyAsString).isValid(idx)) {
                    bf.appendString(tagKeyAsString);
                    bf.appendByte(tagType);
                    bf.appendInt16(static_cast<uint16_t>(alns.TAGS.GetIntTags(tagKeyAsString).GetAt(idx)));
                }
            } else if (tagType == 'i') {
                if (alns.TAGS.GetIntTags(tagKeyAsString).isValid(idx)) {
                    bf.appendString(tagKeyAsString);
                    bf.appendByte(tagType);
                    bf.appendInt32(static_cast<int32_t>(alns.TAGS.GetIntTags(tagKeyAsString).GetAt(idx)));
                }
            } else if (tagType == 'I') {
                if (alns.TAGS.GetIntTags(tagKeyAsString).isValid(idx)) {
                    bf.appendString(tagKeyAsString);
                    bf.appendByte(tagType);
                    bf.appendInt32(static_cast<uint32_t>(alns.TAGS.GetIntTags(tagKeyAsString).GetAt(idx)));
                }
            } else if (tagType == 'f') {
                if (alns.TAGS.GetFloatTags(tagKeyAsString).isValid(idx)) {
                    bf.appendString(tagKeyAsString);
                    bf.appendByte(tagType);
                    bf.appendFloat(static_cast<float>(alns.TAGS.GetFloatTags(tagKeyAsString).GetAt(idx)));
                }
            } else if (tagType == 'Z') {
                if (alns.TAGS.GetStringTags(tagKeyAsString).isValid(idx)) {
                    bf.appendString(tagKeyAsString);
                    bf.appendByte(tagType);
                    bf.appendString(alns.TAGS.GetStringTags(tagKeyAsString).GetAt(idx));
                    bf.appendByte('\0');
                }
            } else if (tagType == 'H') {
                if (alns.TAGS.GetStringTags(tagKeyAsString).isValid(idx)) {
                    bf.appendString(tagKeyAsString);
                    bf.appendByte(tagType);
                    bf.appendString(alns.TAGS.GetStringTags(tagKeyAsString).GetAt(idx));
                    bf.appendByte('\0');
                }
            } else if (tagType == 'B') {
                if (tagSecondaryType == 'c') {
                    if (alns.TAGS.GetArrayInt8Tags(tagKeyAsString).isValid(idx)) {
                        bf.appendString(tagKeyAsString);
                        bf.appendByte(tagType);
                        bf.appendByte(tagSecondaryType);
                        auto vec = alns.TAGS.GetArrayInt8Tags(tagKeyAsString).GetAt(idx);
                        bf.appendInt32(static_cast<int32_t>(vec.size()));
                        for (size_t i = 0; i < vec.size(); i++) {
                            bf.appendByte(static_cast<int8_t>(vec[i]));
                        }
                    }
                } else if (tagSecondaryType == 'C') {
                    if (alns.TAGS.GetArrayUint8Tags(tagKeyAsString).isValid(idx)) {
                        bf.appendString(tagKeyAsString);
                        bf.appendByte(tagType);
                        bf.appendByte(tagSecondaryType);
                        auto vec = alns.TAGS.GetArrayUint8Tags(tagKeyAsString).GetAt(idx);
                        bf.appendInt32(static_cast<int32_t>(vec.size()));
                        for (size_t i = 0; i < vec.size(); i++) {
                            bf.appendByte(static_cast<uint8_t>(vec[i]));
                        }
                    }
                } else if (tagSecondaryType == 's') {
                    if (alns.TAGS.GetArrayInt16Tags(tagKeyAsString).isValid(idx)) {
                        bf.appendString(tagKeyAsString);
                        bf.appendByte(tagType);
                        bf.appendByte(tagSecondaryType);
                        auto vec = alns.TAGS.GetArrayInt16Tags(tagKeyAsString).GetAt(idx);
                        bf.appendInt32(static_cast<int32_t>(vec.size()));
                        for (size_t i = 0; i < vec.size(); i++) {
                            bf.appendInt16(static_cast<int16_t>(vec[i]));
                        }
                    }
                } else if (tagSecondaryType == 'S') {
                    if (alns.TAGS.GetArrayUint16Tags(tagKeyAsString).isValid(idx)) {
                        bf.appendString(tagKeyAsString);
                        bf.appendByte(tagType);
                        bf.appendByte(tagSecondaryType);
                        auto vec = alns.TAGS.GetArrayUint16Tags(tagKeyAsString).GetAt(idx);
                        bf.appendInt32(static_cast<int32_t>(vec.size()));
                        for (size_t i = 0; i < vec.size(); i++) {
                            bf.appendInt16(static_cast<int16_t>(vec[i]));
                        }
                    }
                } else if (tagSecondaryType == 'i') {
                    if (alns.TAGS.GetArrayInt32Tags(tagKeyAsString).isValid(idx)) {
                        bf.appendString(tagKeyAsString);
                        bf.appendByte(tagType);
                        bf.appendByte(tagSecondaryType);
                        auto vec = alns.TAGS.GetArrayInt32Tags(tagKeyAsString).GetAt(idx);
                        bf.appendInt32(static_cast<int32_t>(vec.size()));
                        for (size_t i = 0; i < vec.size(); i++) {
                            bf.appendInt32(static_cast<int32_t>(vec[i]));
                        }
                    }
                } else if (tagSecondaryType == 'I') {
                    if (alns.TAGS.GetArrayUint32Tags(tagKeyAsString).isValid(idx)) {
                        bf.appendString(tagKeyAsString);
                        bf.appendByte(tagType);
                        bf.appendByte(tagSecondaryType);
                        auto vec = alns.TAGS.GetArrayUint32Tags(tagKeyAsString).GetAt(idx);
                        bf.appendInt32(static_cast<int32_t>(vec.size()));
                        for (size_t i = 0; i < vec.size(); i++) {
                            bf.appendInt32(static_cast<uint32_t>(vec[i]));
                        }
                    }
                } else if (tagSecondaryType == 'f') {
                    if (alns.TAGS.GetArrayFloatTags(tagKeyAsString).isValid(idx)) {
                        bf.appendString(tagKeyAsString);
                        bf.appendByte(tagType);
                        bf.appendByte(tagSecondaryType);
                        auto vec = alns.TAGS.GetArrayFloatTags(tagKeyAsString).GetAt(idx);
                        bf.appendInt32(static_cast<int32_t>(vec.size()));
                        for (size_t i = 0; i < vec.size(); i++) {
                            bf.appendFloat(static_cast<float>(vec[i]));
                        }
                    }
                }
            }
        }
    }
    //the size of the block is not calculated in the overall dimension
    int32_t new_block_size = formatted_line.size() - sizeof(int32_t);
    std::memcpy(&formatted_line[0], &new_block_size, sizeof(int32_t));
}

uint8_t BamAlignmentFormatter::encodeNucleotides(char nucleotide1, char nucleotide2) {
    uint8_t first_nibble = 0, second_nibble = 0;

    for (int i = 0; i < 16; ++i) {
        if (seqEncoder[i] == nucleotide1) {
            first_nibble = i;
            break;
        }
    }

    for (int i = 0; i < 16; ++i) {
        if (seqEncoder[i] == nucleotide2) {
            second_nibble = i;
            break;
        }
    }

    return (first_nibble << 4) | second_nibble;
}
