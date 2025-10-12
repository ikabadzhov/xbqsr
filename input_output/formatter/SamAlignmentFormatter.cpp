#include "SamAlignmentFormatter.h"

size_t SamAlignmentFormatter::GetSamAlignmentSize(SamAlignments& alns, ssize_t idx) {
    size_t tot_sz = alns.GetHdr().QnameMapper().GetAt(alns.GetQNAME(idx)).size();
    tot_sz +=1; // account for \t

    tot_sz += std::to_string(alns.GetFLAG(idx)).size();
    tot_sz +=1; // account for \t


    //line.append(alns.GetHdr().getSequenceDictionary().getSequenceName(alns.GetRNAME(idx)));
    tot_sz += alns.GetHdr().getSequenceDictionary().getSequenceName(alns.GetRNAME(idx)).size();
    tot_sz +=1; // account for \t

    //line.append(std::to_string(alns.GetPOS(idx)));
    tot_sz += std::to_string(alns.GetPOS(idx)).size();
    tot_sz +=1; // account for \t

    //line.append(std::to_string(alns.GetMAPQ(idx)));
    tot_sz += std::to_string(alns.GetMAPQ(idx)).size();
    tot_sz +=1; // account for \t

    std::string tmp_str;
    cigarToString(alns.GetCIGAR(idx), tmp_str);
    tot_sz += tmp_str.size();
    //line.append(tmp_str + "\t");
    tot_sz +=1; // account for \t

    if (alns.GetRNEXT(idx) == SequenceRecord::RESERVED_SEQUENCE_NAMES['=']) {
        //line.append("=\t");
        tot_sz +=2; // account for =\t
    } else if (alns.GetRNEXT(idx) == SequenceRecord::RESERVED_SEQUENCE_NAMES['*']) {
        //line.append("*\t");
        tot_sz +=2; // account for \t
    } else {
        if (alns.GetRNEXT(idx) == alns.GetRNAME(idx)) {
            //line.append("=\t");
            tot_sz +=2; // account for \t

        } else {
            tot_sz += alns.GetHdr().getSequenceDictionary().getSequenceName(alns.GetRNEXT(idx)).size();
            //line.append(1, '\t');
            tot_sz +=1; // account for \t
        }
    }
    tot_sz += std::to_string(alns.GetPNEXT(idx)).size();
    tot_sz += 1;

    //line.append(std::to_string(alns.GetTLEN(idx)));
    tot_sz += std::to_string(alns.GetTLEN(idx)).size();
    //line.append(1, '\t');
    tot_sz += 1;

    //size_t currentIndex = line.size();
    //line.resize(line.size() +  alns.GetSEQ(idx).size());
   // memcpy(line.data() + currentIndex, alns.GetSEQ(idx).data(), alns.GetSEQ(idx).size());
    tot_sz += alns.GetSEQ(idx).size();
    tot_sz += 1;

    //currentIndex = line.size();
    //line.resize(line.size() +  alns.GetQUAL(idx).size())
    tot_sz += alns.GetQUAL(idx).size();
    tot_sz += 1;

    std::string value;

    for (auto& tag : alns.TAGS.ALL_TAGS){

        std::string tagKeyAsString = tag.second.getTagKeyAsString();
        char tagType = tag.second.getTagType();

        if (tagType == 'c' || tagType == 'C' || tagType == 's' || tagType == 'S' || tagType == 'i' || tagType == 'I')
            tagType = 'i';

        char tagSecondaryType = tag.second.getSecondaryTagType();

        if (tagKeyAsString == "RG") {
            // it shouldn't need any check on the presence of RG tagValue
            tot_sz += alns.GetHdr().getReadGroupRecord().getString(alns.TAGS.RG[idx]).size();//alns.TAGS.RG.toString(idx).size();
            //line.append(1, '\t');
            tot_sz += 6;
        } else {
            if (tagType == 'i') {
                if (alns.TAGS.GetIntTags(tagKeyAsString).isValid(idx)) {
                    // The tag it's always 2 chars
                    //tot_sz += 2; //std::to_string(alns.TAGS.GetIntTags(tagKeyAsString).GetAt(idx)).size();
                    //line.append(1, '\t');
                    tot_sz += 6;
                    tot_sz += std::to_string(alns.TAGS.GetIntTags(tagKeyAsString).GetAt(idx)).size();
                }
            } else if (tagType == 'f') {
                if (alns.TAGS.GetFloatTags(tagKeyAsString).isValid(idx)) {
                    //value = std::to_string(alns.TAGS.GetFloatTags(tagKeyAsString).GetAt(idx));
                    //line.append(1, '\t');
                    tot_sz += 6;
                    //line.append(value);
                    tot_sz += std::to_string(alns.TAGS.GetFloatTags(tagKeyAsString).GetAt(idx)).size();
                }
            } else if (tagType == 'Z' || tagType == 'H') {
                if (alns.TAGS.GetStringTags(tagKeyAsString).isValid(idx)) {
                    //value = alns.TAGS.GetStringTags(tagKeyAsString).GetAt(idx);
                    //line.append(1, '\t');
                    //line.append(tagKeyAsString);
                    //line.append(1, ':');
                    //line.append(1, tagType);
                    //line.append(1, ':');
                    //tot_sz += (1+2+1+1+1);
                    tot_sz += 6;
                    tot_sz += alns.TAGS.GetStringTags(tagKeyAsString).GetAt(idx).size();
                    //line.append(value);
                }
            } else if (tagType == 'A') {
                if (alns.TAGS.GetCharTags(tagKeyAsString).isValid(idx)) {
                    value = alns.TAGS.GetCharTags(tagKeyAsString).GetAt(idx);
//                    line.append(1, '\t');
//                    line.append(tagKeyAsString);
//                    line.append(1, ':');
//                    line.append(1, tagType);
//                    line.append(1, ':')
                    tot_sz += (1+2+1+1+1);
                    tot_sz += alns.TAGS.GetCharTags(tagKeyAsString).GetAt(idx);
                    //line.append(value);
                }
            }
            else if (tagType == 'B') {
                if (tagSecondaryType == 'c') {
                    if (alns.TAGS.GetArrayInt8Tags(tagKeyAsString).isValid(idx)) {
                        tot_sz += (7+value.size());
//                        line.append(1, '\t');
//                        line.append(tagKeyAsString);
//                        line.append(1, ':');
//                        line.append(1, tagType);
//                        line.append(1, tagSecondaryType);
//                        line.append(1, ':');
//                        line.append(value);
                        auto vec = alns.TAGS.GetArrayInt8Tags(tagKeyAsString).GetAt(idx);
                        for (size_t vIdx = 0; vIdx < vec.size()-1; ++vIdx) {
//                            value.append(std::to_string(vec[vIdx]));
//                            value.append(",");
                            tot_sz += std::to_string(vec[vIdx]).size() + 1;
                        }
                        //value.append(std::to_string(vec.back()));
                        tot_sz += std::to_string(vec.back()).size();
                    }
                } else if (tagSecondaryType == 'C') {
                    if (alns.TAGS.GetArrayUint8Tags(tagKeyAsString).isValid(idx)) {
                        tot_sz += (7+value.size());
//                        line.append(1, '\t');
//                        line.append(tagKeyAsString);
//                        line.append(1, ':');
//                        line.append(1, tagType);
//                        line.append(1, tagSecondaryType);
//                        line.append(1, ':');
//                        line.append(value);
                        auto vec = alns.TAGS.GetArrayUint8Tags(tagKeyAsString).GetAt(idx);
                        for (size_t vIdx = 0; vIdx < vec.size()-1; ++vIdx) {
//                            value.append(std::to_string(vec[vIdx]));
//                            value.append(",");
                            tot_sz += std::to_string(vec[vIdx]).size() + 1;
                        }
                        //value.append(std::to_string(vec.back()));
                        tot_sz += std::to_string(vec.back()).size();
                    }
                } else if (tagSecondaryType == 's') {
                    if (alns.TAGS.GetArrayInt16Tags(tagKeyAsString).isValid(idx)) {
                        tot_sz += (7+value.size());
//                        line.append(1, '\t');
//                        line.append(tagKeyAsString);
//                        line.append(1, ':');
//                        line.append(1, tagType);
//                        line.append(1, tagSecondaryType);
//                        line.append(1, ':');
//                        line.append(value);
                        auto vec = alns.TAGS.GetArrayInt16Tags(tagKeyAsString).GetAt(idx);
                        for (size_t vIdx = 0; vIdx < vec.size()-1; ++vIdx) {
//                            value.append(std::to_string(vec[vIdx]));
//                            value.append(",");
                            tot_sz += std::to_string(vec[vIdx]).size() + 1;
                        }
                        //value.append(std::to_string(vec.back()));
                        tot_sz += std::to_string(vec.back()).size();
                    }
                } else if (tagSecondaryType == 'S') {
                    if (alns.TAGS.GetArrayUint16Tags(tagKeyAsString).isValid(idx)) {
                        tot_sz += (7+value.size());
//                        line.append(1, '\t');
//                        line.append(tagKeyAsString);
//                        line.append(1, ':');
//                        line.append(1, tagType);
//                        line.append(1, tagSecondaryType);
//                        line.append(1, ':');
//                        line.append(value);
                        auto vec = alns.TAGS.GetArrayUint16Tags(tagKeyAsString).GetAt(idx);
                        for (size_t vIdx = 0; vIdx < vec.size()-1; ++vIdx) {
//                            value.append(std::to_string(vec[vIdx]));
//                            value.append(",");
                            tot_sz += std::to_string(vec[vIdx]).size() + 1;
                        }
                        //value.append(std::to_string(vec.back()));
                        tot_sz += std::to_string(vec.back()).size();
                    }
                } else if (tagSecondaryType == 'i') {
                    if (alns.TAGS.GetArrayInt32Tags(tagKeyAsString).isValid(idx)) {
                        tot_sz += (7+value.size());
//                        line.append(1, '\t');
//                        line.append(tagKeyAsString);
//                        line.append(1, ':');
//                        line.append(1, tagType);
//                        line.append(1, tagSecondaryType);
//                        line.append(1, ':');
//                        line.append(value);
                        auto vec = alns.TAGS.GetArrayInt32Tags(tagKeyAsString).GetAt(idx);
                        for (size_t vIdx = 0; vIdx < vec.size()-1; ++vIdx) {
//                            value.append(std::to_string(vec[vIdx]));
//                            value.append(",");
                            tot_sz += std::to_string(vec[vIdx]).size() + 1;
                        }
                        //value.append(std::to_string(vec.back()));
                        tot_sz += std::to_string(vec.back()).size();
                    }
                } else if (tagSecondaryType == 'I') {
                    if (alns.TAGS.GetArrayUint32Tags(tagKeyAsString).isValid(idx)) {
                        tot_sz += (7+value.size());

//                        line.append(1, '\t');
//                        line.append(tagKeyAsString);
//                        line.append(1, ':');
//                        line.append(1, tagType);
//                        line.append(1, tagSecondaryType);
//                        line.append(1, ':');
//                        line.append(value);
                        auto vec = alns.TAGS.GetArrayUint32Tags(tagKeyAsString).GetAt(idx);
                        for (size_t vIdx = 0; vIdx < vec.size()-1; ++vIdx) {
//                            value.append(std::to_string(vec[vIdx]));
//                            value.append(",");
                            tot_sz += std::to_string(vec[vIdx]).size() + 1;
                        }
                        //value.append(std::to_string(vec.back()));
                        tot_sz += std::to_string(vec.back()).size();
                    }
                } else if (tagSecondaryType == 'f') {
                    if (alns.TAGS.GetArrayFloatTags(tagKeyAsString).isValid(idx)) {
                        tot_sz += (7+value.size());

//                        line.append(1, '\t');
//                        line.append(tagKeyAsString);
//                        line.append(1, ':');
//                        line.append(1, tagType);
//                        line.append(1, tagSecondaryType);
//                        line.append(1, ':');
//                        line.append(value);
                        auto vec = alns.TAGS.GetArrayFloatTags(tagKeyAsString).GetAt(idx);
                        for (size_t vIdx = 0; vIdx < vec.size()-1; ++vIdx) {
//                            value.append(std::to_string(vec[vIdx]));
//                            value.append(",");
                            tot_sz += std::to_string(vec[vIdx]).size() + 1;
                        }
                        //value.append(std::to_string(vec.back()));
                        tot_sz += std::to_string(vec.back()).size();
                    }
                }
            }
        }
    }
    return tot_sz;
   // line += "\n";
}

// Helper lambda to copy a std::string_view into out and update offset.
void copyStringView(char* out, std::string_view sv, size_t& offset) {
    std::memcpy(out + offset, sv.data(), sv.size());
    offset += sv.size();
};


//[[clang::optnone]]
void SamAlignmentFormatter::FormatSamAlignmentInBuffer(SamAlignments& alns, ssize_t idx, char* out) {
    size_t offset = 0;

    // 1. QNAME
    {
        std::string_view qname = alns.GetHdr().QnameMapper().GetAt(alns.GetQNAME(idx));
        copyStringView(out, qname, offset);
    }
    out[offset++] = '\t';

    // 2. FLAG
    {
        std::string flag = std::to_string(alns.GetFLAG(idx));
        copyStringView(out, flag, offset);
    }
    out[offset++] = '\t';

    // 3. RNAME
    {
        std::string_view rname = alns.GetHdr().getSequenceDictionary().getSequenceName(alns.GetRNAME(idx));
        copyStringView(out, rname, offset);
    }
    out[offset++] = '\t';

    // 4. POS
    {
        std::string pos = std::to_string(alns.GetPOS(idx));
        copyStringView(out, pos, offset);
    }
    out[offset++] = '\t';

    // 5. MAPQ
    {
        std::string mapq = std::to_string(alns.GetMAPQ(idx));
        copyStringView(out, mapq, offset);
    }
    out[offset++] = '\t';

    // 6. CIGAR
    {
        std::string tmp_str;
        cigarToString(alns.GetCIGAR(idx), tmp_str);
        // Append the CIGAR string followed by a tab.
        copyStringView(out, tmp_str, offset);
        out[offset++] = '\t';
    }

    // 7. RNEXT
    {
        if (alns.GetRNEXT(idx) == SequenceRecord::RESERVED_SEQUENCE_NAMES['=']) {
            copyStringView(out, "=\t", offset);
        } else if (alns.GetRNEXT(idx) == SequenceRecord::RESERVED_SEQUENCE_NAMES['*']) {
            copyStringView(out, "*\t", offset);
        } else {
            if (alns.GetRNEXT(idx) == alns.GetRNAME(idx)) {
                copyStringView(out, "=\t", offset);
            } else {
                copyStringView(out, alns.GetHdr().getSequenceDictionary().getSequenceName(alns.GetRNEXT(idx)), offset);
                //line.append(1, '\t');
                out[offset++] = '\t';
            }
        }
    }

    // 8. PNEXT
    {
        std::string pnext = std::to_string(alns.GetPNEXT(idx));
        copyStringView(out, pnext, offset);
    }
    out[offset++] = '\t';

    // 9. TLEN
    {
        std::string tlen = std::to_string(alns.GetTLEN(idx));
        copyStringView(out, tlen, offset);
    }
    out[offset++] = '\t';

    // 10. SEQ: copy the sequence as is.
    {
        size_t seqSize = alns.GetSEQ(idx).size();
        std::memcpy(out + offset, alns.GetSEQ(idx).data(), seqSize);
        offset += seqSize;
    }
    out[offset++] = '\t';

    // 11. QUAL: copy quality string and add 33 to each character.
    {
        size_t qualSize = alns.GetQUAL(idx).size();
        for (size_t i = 0; i < qualSize; i++) {
            // Adjust by 33 and store directly.
            out[offset + i] = alns.GetQUAL(idx)[i] + 33;
        }
        offset += qualSize;
    }

    // 12. Process TAG fields.
    // We'll use a temporary string 'value' for conversions.
    std::string value;
    for (auto& tag : alns.TAGS.ALL_TAGS) {
        uint16_t tagKey = tag.second.getTagKey();
        std::string tagKeyAsString = tag.second.getTagKeyAsString();
        char tagType = tag.second.getTagType();
        if (tagType == 'c' || tagType == 'C' || tagType == 's' ||
            tagType == 'S' || tagType == 'i' || tagType == 'I')
            tagType = 'i';
        char tagSecondaryType = tag.second.getSecondaryTagType();

        // Prepend a tab before each tag field.
        //out[offset++] = '\t';

        if (tagKeyAsString == "RG") {
            out[offset++] = '\t';
            // For the RG tag, get its value.
            value = alns.GetHdr().getReadGroupRecord().getString(alns.TAGS.RG[idx]);
            std::string decoded = TagKeyAndType::decode(tagKey);
            copyStringView(out, decoded, offset);
            out[offset++] = ':';
            out[offset++] = tagType;
            out[offset++] = ':';
            copyStringView(out, value, offset);
        } else {
            if (tagType == 'i') {
                if (alns.TAGS.GetIntTags(tagKeyAsString).isValid(idx)) {
                    out[offset++] = '\t';
                    value = std::to_string(alns.TAGS.GetIntTags(tagKeyAsString).GetAt(idx));
                    copyStringView(out, tagKeyAsString, offset);
                    out[offset++] = ':';
                    out[offset++] = tagType;
                    out[offset++] = ':';
                    copyStringView(out, value, offset);
                }
            } else if (tagType == 'f') {
                if (alns.TAGS.GetFloatTags(tagKeyAsString).isValid(idx)) {
                    out[offset++] = '\t';
                    value = std::to_string(alns.TAGS.GetFloatTags(tagKeyAsString).GetAt(idx));
                    copyStringView(out, tagKeyAsString, offset);
                    out[offset++] = ':';
                    out[offset++] = tagType;
                    out[offset++] = ':';
                    copyStringView(out, value, offset);
                }
            } else if (tagType == 'Z' || tagType == 'H') {
                if (alns.TAGS.GetStringTags(tagKeyAsString).isValid(idx)) {
                    out[offset++] = '\t';
                    value = alns.TAGS.GetStringTags(tagKeyAsString).GetAt(idx);
                    copyStringView(out, tagKeyAsString, offset);
                    out[offset++] = ':';
                    out[offset++] = tagType;
                    out[offset++] = ':';
                    copyStringView(out,value, offset);
                }
            } else if (tagType == 'A') {
                if (alns.TAGS.GetCharTags(tagKeyAsString).isValid(idx)) {
                    out[offset++] = '\t';
                    value = alns.TAGS.GetCharTags(tagKeyAsString).GetAt(idx);
                    copyStringView(out, tagKeyAsString, offset);
                    out[offset++] = ':';
                    out[offset++] = tagType;
                    out[offset++] = ':';
                    copyStringView(out, value, offset);
                }
            } else if (tagType == 'B') {
                // Array type: process based on secondary type.
                if (tagSecondaryType == 'c') {
                    if (alns.TAGS.GetArrayInt8Tags(tagKeyAsString).isValid(idx)) {
                        out[offset++] = '\t';
                        copyStringView(out, tagKeyAsString, offset);
                        out[offset++] = ':';
                        out[offset++] = tagType;
                        out[offset++] = tagSecondaryType;
                        out[offset++] = ':';
                        auto vec = alns.TAGS.GetArrayInt8Tags(tagKeyAsString).GetAt(idx);
                        for (size_t vIdx = 0; vIdx < vec.size() - 1; ++vIdx) {
                            std::string elem = std::to_string(vec[vIdx]);
                            copyStringView(out, elem, offset);
                            out[offset++] = ',';
                        }
                        copyStringView(out, std::to_string(vec.back()), offset);
                    }
                } else if (tagSecondaryType == 'C') {
                    if (alns.TAGS.GetArrayUint8Tags(tagKeyAsString).isValid(idx)) {
                        out[offset++] = '\t';
                        copyStringView(out, tagKeyAsString, offset);
                        out[offset++] = ':';
                        out[offset++] = tagType;
                        out[offset++] = tagSecondaryType;
                        out[offset++] = ':';
                        auto vec = alns.TAGS.GetArrayUint8Tags(tagKeyAsString).GetAt(idx);
                        for (size_t vIdx = 0; vIdx < vec.size(); ++vIdx) {
                            std::string elem = std::to_string(vec[vIdx]);
                            copyStringView(out, elem, offset);
                            if (vIdx < vec.size() - 1)
                                out[offset++] = ',';
                        }
                    }
                } else if (tagSecondaryType == 's') {
                    if (alns.TAGS.GetArrayInt16Tags(tagKeyAsString).isValid(idx)) {
                        out[offset++] = '\t';
                        copyStringView(out, tagKeyAsString, offset);
                        out[offset++] = ':';
                        out[offset++] = tagType;
                        out[offset++] = tagSecondaryType;
                        out[offset++] = ':';
                        auto vec = alns.TAGS.GetArrayInt16Tags(tagKeyAsString).GetAt(idx);
                        for (size_t vIdx = 0; vIdx < vec.size(); ++vIdx) {
                            std::string elem = std::to_string(vec[vIdx]);
                            copyStringView(out, elem, offset);
                            if (vIdx < vec.size() - 1)
                                out[offset++] = ',';
                        }
                    }
                } else if (tagSecondaryType == 'S') {
                    if (alns.TAGS.GetArrayUint16Tags(tagKeyAsString).isValid(idx)) {
                        out[offset++] = '\t';
                        copyStringView(out, tagKeyAsString, offset);
                        out[offset++] = ':';
                        out[offset++] = tagType;
                        out[offset++] = tagSecondaryType;
                        out[offset++] = ':';
                        auto vec = alns.TAGS.GetArrayUint16Tags(tagKeyAsString).GetAt(idx);
                        for (size_t vIdx = 0; vIdx < vec.size(); ++vIdx) {
                            std::string elem = std::to_string(vec[vIdx]);
                            copyStringView(out, elem, offset);
                            if (vIdx < vec.size() - 1)
                                out[offset++] = ',';
                        }
                    }
                } else if (tagSecondaryType == 'i') {
                    if (alns.TAGS.GetArrayInt32Tags(tagKeyAsString).isValid(idx)) {
                        out[offset++] = '\t';
                        copyStringView(out, tagKeyAsString, offset);
                        out[offset++] = ':';
                        out[offset++] = tagType;
                        out[offset++] = tagSecondaryType;
                        out[offset++] = ':';
                        auto vec = alns.TAGS.GetArrayInt32Tags(tagKeyAsString).GetAt(idx);
                        for (size_t vIdx = 0; vIdx < vec.size(); ++vIdx) {
                            std::string elem = std::to_string(vec[vIdx]);
                            copyStringView(out, elem, offset);
                            if (vIdx < vec.size() - 1)
                                out[offset++] = ',';
                        }
                    }
                } else if (tagSecondaryType == 'I') {
                    if (alns.TAGS.GetArrayUint32Tags(tagKeyAsString).isValid(idx)) {
                        out[offset++] = '\t';
                        copyStringView(out, tagKeyAsString, offset);
                        out[offset++] = ':';
                        out[offset++] = tagType;
                        out[offset++] = tagSecondaryType;
                        out[offset++] = ':';
                        auto vec = alns.TAGS.GetArrayUint32Tags(tagKeyAsString).GetAt(idx);
                        for (size_t vIdx = 0; vIdx < vec.size(); ++vIdx) {
                            std::string elem = std::to_string(vec[vIdx]);
                            copyStringView(out, elem, offset);
                            if (vIdx < vec.size() - 1)
                                out[offset++] = ',';
                        }
                    }
                } else if (tagSecondaryType == 'f') {
                    if (alns.TAGS.GetArrayFloatTags(tagKeyAsString).isValid(idx)) {
                        out[offset++] = '\t';
                        copyStringView(out, tagKeyAsString, offset);
                        out[offset++] = ':';
                        out[offset++] = tagType;
                        out[offset++] = tagSecondaryType;
                        out[offset++] = ':';
                        auto vec = alns.TAGS.GetArrayFloatTags(tagKeyAsString).GetAt(idx);
                        for (size_t vIdx = 0; vIdx < vec.size(); ++vIdx) {
                            std::string elem = std::to_string(vec[vIdx]);
                            copyStringView(out, elem, offset);
                            if (vIdx < vec.size() - 1)
                                out[offset++] = ',';
                        }
                    }
                }
            }
        }
    }
    // Finally, append newline.
    out[offset++] = '\n';
}

void SamAlignmentFormatter::FormatSamAlignment(SamAlignments& alns, ssize_t idx, std::string& line) {
    line.append(alns.GetHdr().QnameMapper().GetAt(alns.GetQNAME(idx)));
    line.append(1, '\t');

    line.append(std::to_string(alns.GetFLAG(idx)));
    line.append(1, '\t');

    line.append(alns.GetHdr().getSequenceDictionary().getSequenceName(alns.GetRNAME(idx)));
    line.append(1, '\t');

    line.append(std::to_string(alns.GetPOS(idx)));
    line.append(1, '\t');

    line.append(std::to_string(alns.GetMAPQ(idx)));
    line.append(1, '\t');

    std::string tmp_str;
    cigarToString(alns.GetCIGAR(idx), tmp_str);
    line.append(tmp_str + "\t");

    if (alns.GetRNEXT(idx) == SequenceRecord::RESERVED_SEQUENCE_NAMES['=']) {
        line.append("=\t");
    } else if (alns.GetRNEXT(idx) == SequenceRecord::RESERVED_SEQUENCE_NAMES['*']) {
        line.append("*\t");
    } else {
        if (alns.GetRNEXT(idx) == alns.GetRNAME(idx)) {
            line.append("=\t");
        } else {
            line.append(alns.GetHdr().getSequenceDictionary().getSequenceName(alns.GetRNEXT(idx)));
            line.append(1, '\t');
        }
    }
    line.append(std::to_string(alns.GetPNEXT(idx)));
    line.append(1, '\t');

    line.append(std::to_string(alns.GetTLEN(idx)));
    line.append(1, '\t');

    size_t currentIndex = line.size();
    line.resize(line.size() +  alns.GetSEQ(idx).size());
    memcpy(line.data() + currentIndex, alns.GetSEQ(idx).data(), alns.GetSEQ(idx).size());
    //for (size_t i = 0; i < l; i++) {
        //line.append(1, static_cast<char>(alns.GetSEQ(idx).Base(i)));
    //}

    line.append(1, '\t');

    currentIndex = line.size();
    line.resize(line.size() +  alns.GetQUAL(idx).size());
    memcpy(line.data() + currentIndex, alns.GetQUAL(idx).data(), alns.GetQUAL(idx).size());
    for (size_t i = 0; i < alns.GetQUAL(idx).size(); i++) {
        line[currentIndex+i] += 33;
    }

    std::string value;

    for (auto& tag : alns.TAGS.ALL_TAGS){

        uint16_t tagKey = tag.second.getTagKey();
        std::string tagKeyAsString = tag.second.getTagKeyAsString();
        char tagType = tag.second.getTagType();

        if (tagType == 'c' || tagType == 'C' || tagType == 's' || tagType == 'S' || tagType == 'i' || tagType == 'I')
            tagType = 'i';


        char tagSecondaryType = tag.second.getSecondaryTagType();

        if (tagKeyAsString == "RG") {
            // it shouldn't need any check on the presence of RG tagValue
            value = alns.GetHdr().getReadGroupRecord().getString(alns.TAGS.RG[idx]);
            line.append(1, '\t');
            line.append(TagKeyAndType::decode(tagKey));
            line.append(1, ':');
            line.append(1, tagType);
            line.append(1, ':');
            line.append(value);
        } else {
            if (tagType == 'i') {
                if (alns.TAGS.GetIntTags(tagKeyAsString).isValid(idx)) {
                    value = std::to_string(alns.TAGS.GetIntTags(tagKeyAsString).GetAt(idx));
                    line.append(1, '\t');
                    line.append(tagKeyAsString);
                    line.append(1, ':');
                    line.append(1, tagType);
                    line.append(1, ':');
                    line.append(value);
                }
            } else if (tagType == 'f') {
                if (alns.TAGS.GetFloatTags(tagKeyAsString).isValid(idx)) {
                    value = std::to_string(alns.TAGS.GetFloatTags(tagKeyAsString).GetAt(idx));
                    line.append(1, '\t');
                    line.append(tagKeyAsString);
                    line.append(1, ':');
                    line.append(1, tagType);
                    line.append(1, ':');
                    line.append(value);
                }
            } else if (tagType == 'Z' || tagType == 'H') {
                if (alns.TAGS.GetStringTags(tagKeyAsString).isValid(idx)) {
                    value = alns.TAGS.GetStringTags(tagKeyAsString).GetAt(idx);
                    line.append(1, '\t');
                    line.append(tagKeyAsString);
                    line.append(1, ':');
                    line.append(1, tagType);
                    line.append(1, ':');
                    line.append(value);
                }
            } else if (tagType == 'A') {
                if (alns.TAGS.GetCharTags(tagKeyAsString).isValid(idx)) {
                    value = alns.TAGS.GetCharTags(tagKeyAsString).GetAt(idx);
                    line.append(1, '\t');
                    line.append(tagKeyAsString);
                    line.append(1, ':');
                    line.append(1, tagType);
                    line.append(1, ':');
                    line.append(value);
                }
            }
            else if (tagType == 'B') {
                if (tagSecondaryType == 'c') {
                    if (alns.TAGS.GetArrayInt8Tags(tagKeyAsString).isValid(idx)) {
                        line.append(1, '\t');
                        line.append(tagKeyAsString);
                        line.append(1, ':');
                        line.append(1, tagType);
                        line.append(1, tagSecondaryType);
                        line.append(1, ':');
                        line.append(value);
                        auto vec = alns.TAGS.GetArrayInt8Tags(tagKeyAsString).GetAt(idx);
                        for (size_t vIdx = 0; vIdx < vec.size()-1; ++vIdx) {
                            value.append(std::to_string(vec[vIdx]));
                            value.append(",");
                        }
                        value.append(std::to_string(vec.back()));
                    }
                } else if (tagSecondaryType == 'C') {
                    if (alns.TAGS.GetArrayUint8Tags(tagKeyAsString).isValid(idx)) {
                        line.append(1, '\t');
                        line.append(tagKeyAsString);
                        line.append(1, ':');
                        line.append(1, tagType);
                        line.append(1, tagSecondaryType);
                        line.append(1, ':');
                        line.append(value);
                        auto vec = alns.TAGS.GetArrayUint8Tags(tagKeyAsString).GetAt(idx);
                        for (size_t vIdx = 0; vIdx < vec.size()-1; ++vIdx) {
                            value.append(std::to_string(vec[vIdx]));
                            value.append(",");
                        }
                        value.append(std::to_string(vec.back()));
                    }
                } else if (tagSecondaryType == 's') {
                    if (alns.TAGS.GetArrayInt16Tags(tagKeyAsString).isValid(idx)) {
                        line.append(1, '\t');
                        line.append(tagKeyAsString);
                        line.append(1, ':');
                        line.append(1, tagType);
                        line.append(1, tagSecondaryType);
                        line.append(1, ':');
                        line.append(value);
                        auto vec = alns.TAGS.GetArrayInt16Tags(tagKeyAsString).GetAt(idx);
                        for (size_t vIdx = 0; vIdx < vec.size()-1; ++vIdx) {
                            value.append(std::to_string(vec[vIdx]));
                            value.append(",");
                        }
                        value.append(std::to_string(vec.back()));
                    }
                } else if (tagSecondaryType == 'S') {
                    if (alns.TAGS.GetArrayUint16Tags(tagKeyAsString).isValid(idx)) {
                        line.append(1, '\t');
                        line.append(tagKeyAsString);
                        line.append(1, ':');
                        line.append(1, tagType);
                        line.append(1, tagSecondaryType);
                        line.append(1, ':');
                        line.append(value);
                        auto vec = alns.TAGS.GetArrayUint16Tags(tagKeyAsString).GetAt(idx);
                        for (size_t vIdx = 0; vIdx < vec.size()-1; ++vIdx) {
                            value.append(std::to_string(vec[vIdx]));
                            value.append(",");
                        }
                        value.append(std::to_string(vec.back()));
                    }
                } else if (tagSecondaryType == 'i') {
                    if (alns.TAGS.GetArrayInt32Tags(tagKeyAsString).isValid(idx)) {
                        line.append(1, '\t');
                        line.append(tagKeyAsString);
                        line.append(1, ':');
                        line.append(1, tagType);
                        line.append(1, tagSecondaryType);
                        line.append(1, ':');
                        line.append(value);
                        auto vec = alns.TAGS.GetArrayInt32Tags(tagKeyAsString).GetAt(idx);
                        for (size_t vIdx = 0; vIdx < vec.size()-1; ++vIdx) {
                            value.append(std::to_string(vec[vIdx]));
                            value.append(",");
                        }
                        value.append(std::to_string(vec.back()));
                    }
                } else if (tagSecondaryType == 'I') {
                    if (alns.TAGS.GetArrayUint32Tags(tagKeyAsString).isValid(idx)) {
                        line.append(1, '\t');
                        line.append(tagKeyAsString);
                        line.append(1, ':');
                        line.append(1, tagType);
                        line.append(1, tagSecondaryType);
                        line.append(1, ':');
                        line.append(value);
                        auto vec = alns.TAGS.GetArrayUint32Tags(tagKeyAsString).GetAt(idx);
                        for (size_t vIdx = 0; vIdx < vec.size()-1; ++vIdx) {
                            value.append(std::to_string(vec[vIdx]));
                            value.append(",");
                        }
                        value.append(std::to_string(vec.back()));
                    }
                } else if (tagSecondaryType == 'f') {
                    if (alns.TAGS.GetArrayFloatTags(tagKeyAsString).isValid(idx)) {
                        line.append(1, '\t');
                        line.append(tagKeyAsString);
                        line.append(1, ':');
                        line.append(1, tagType);
                        line.append(1, tagSecondaryType);
                        line.append(1, ':');
                        line.append(value);
                        auto vec = alns.TAGS.GetArrayFloatTags(tagKeyAsString).GetAt(idx);
                        for (size_t vIdx = 0; vIdx < vec.size()-1; ++vIdx) {
                            value.append(std::to_string(vec[vIdx]));
                            value.append(",");
                        }
                        value.append(std::to_string(vec.back()));
                    }
                }
            }
        }
    }
    line += "\n";
}
