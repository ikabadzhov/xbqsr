#include "SamAlignmentParser.h"

void SamAlignmentParser::ParseSamAlignmentFromSlice(std::string_view s, SamAlignments &reads, size_t idx) {
    StringScanner sc;
    sc.reset(s);
    SamAlignmentParser::ParseSamAlignment(sc, reads, idx);
}


// ParseSamAlignment parses a read alignment line in a SAM file and
// assign it to the alignment data structure. See
// http://samtools.github.io/hts-specs/SAMv1.pdf - Sections 1.4 and
// 1.5.
void SamAlignmentParser::ParseSamAlignment(StringScanner &sc, SamAlignments &reads, size_t idx) {
    std::string value;
    sc.doString(value);
    size_t qnameMappedIdx = reads.GetHdr().QnameMapper().addString(value);
    //reads.QNAME[idx] = qnameMappedIdx;
    //reads.TILE[idx] = computeTileInfo(value);

    reads.FLAG[idx] = static_cast<uint16_t>(sc.doUint());

    sc.doString(value);
    reads.RNAME[idx] = reads.GetHdr().getSequenceDictionary().getIdx(value);

    reads.POS[idx] = sc.doInt32();

    reads.MAPQ[idx] = static_cast<uint8_t>(sc.doUint());

    std::string cigar_string;
    sc.doString(cigar_string);
    std::vector<CigarOperation> cigar;
    ScanCigarString(cigar_string, cigar);
    reads.CIGAR.SetAt(cigar, idx);

    sc.doString(value);

    reads.RNEXT[idx] = reads.GetHdr().getSequenceDictionary().getIdx(value);

    reads.PNEXT[idx] = sc.doInt32();
    reads.TLEN[idx] = sc.doInt32();

    std::string seq;
    sc.doSeq(seq);
    reads.SEQ.SetAt(seq, idx);

    std::string qual;
    sc.readBytes(qual);
    std::vector<char> qualVec(qual.size());
    size_t i = 0;
    for (auto q: qual) {
        qualVec[i] = static_cast<char>(static_cast<int8_t>(q) - 33);
        i++;
    }
    reads.QUAL.SetAt(qualVec, idx);

    //Get all other field as optional fields
    //int64_t tagValueInt;
    //float tagValueFloat;
    std::vector<int32_t> tagValueArrayInt32;
    std::vector<uint32_t> tagValueArrayUInt32;
    std::vector<int16_t> tagValueArrayInt16;
    std::vector<uint16_t> tagValueArrayUInt16;
    std::vector<int8_t> tagValueArrayInt8;
    std::vector<uint8_t> tagValueArrayUInt8;
    std::vector<float> tagValueArrayFloat;

    std::string key(2, 0);

    while (sc.len() > 0) {
        tagValueArrayUInt32.clear();
        tagValueArrayInt16.clear();
        tagValueArrayUInt16.clear();
        tagValueArrayInt8.clear();
        tagValueArrayUInt8.clear();
        tagValueArrayFloat.clear();
        std::string types;
        std::string tagValue;

        TagKeyAndType tagKeyAndType;

        std::string fullTag;

        sc.readUntil('\t', fullTag);
        key[0] = fullTag[0];
        key[1] = fullTag[1];

        assert(fullTag[2] == ':');
        tagKeyAndType.setTagType(fullTag[3]);
        if (tagKeyAndType.getTagType() == 'B') {
            tagKeyAndType.setSecondaryTagType(fullTag[4]);
            assert(fullTag[5] == ':');
            tagValue = fullTag.substr(6);
        } else {
            tagValue = fullTag.substr(5);
            assert(fullTag[4] == ':');
        }

        tagKeyAndType.setTagKey(TagKeyAndType::encode(key));

        if (tagKeyAndType.getTagType() == 'A') {

            reads.TAGS.InsertCharTagAt(idx, tagKeyAndType, static_cast<char>(tagValue[0]));


        } else if (tagKeyAndType.getTagType() == 'c' || tagKeyAndType.getTagType() == 'C' ||
                   tagKeyAndType.getTagType() == 's' ||
                   tagKeyAndType.getTagType() == 'S' || tagKeyAndType.getTagType() == 'i' ||
                   tagKeyAndType.getTagType() == 'I') {

            reads.TAGS.InsertIntTagAt(idx, tagKeyAndType, static_cast<int64_t>(stoll(tagValue)));

        } else if (tagKeyAndType.getTagType() == 'f') {

            reads.TAGS.InsertFloatTagAt(idx, tagKeyAndType, static_cast<float>(stof(tagValue)));

        } else if (tagKeyAndType.getTagType() == 'Z') {

            reads.TAGS.InsertStringTagAt(idx, tagKeyAndType, tagValue);

        } else if (tagKeyAndType.getTagType() == 'H') {

            reads.TAGS.InsertArrayHexTagAt(idx, tagKeyAndType, tagValue);

        } else if (tagKeyAndType.getTagType() == 'B') {
            std::string arrayValue;
            StringScanner tagValueScanner;

            tagValueScanner.reset(tagValue);

            if (tagKeyAndType.getSecondaryTagType() == 'c') {
                tagValueArrayInt8.clear();
                while (tagValueScanner.len()) {
                    tagValueScanner.readUntil(',', arrayValue);
                    tagValueArrayInt8.emplace_back(static_cast<int8_t>(stol(arrayValue)));
                }
                reads.TAGS.InsertArrayInt8TagAt(idx, tagKeyAndType, tagValueArrayInt8);
            } else if (tagKeyAndType.getSecondaryTagType() == 'C') {
                tagValueArrayUInt8.clear();
                while (tagValueScanner.len()) {
                    tagValueScanner.readUntil(',', arrayValue);
                    tagValueArrayUInt8.emplace_back(static_cast<uint8_t>(stol(arrayValue)));
                }
                reads.TAGS.InsertArrayUint8TagAt(idx, tagKeyAndType, tagValueArrayUInt8);
            } else if (tagKeyAndType.getSecondaryTagType() == 's') {
                tagValueArrayInt16.clear();
                while (tagValueScanner.len()) {
                    tagValueScanner.readUntil(',', arrayValue);
                    tagValueArrayInt16.emplace_back(static_cast<int16_t>(stol(arrayValue)));
                }
                reads.TAGS.InsertArrayInt16TagAt(idx, tagKeyAndType, tagValueArrayInt16);
            } else if (tagKeyAndType.getSecondaryTagType() == 'S') {
                tagValueArrayUInt16.clear();
                while (tagValueScanner.len()) {
                    tagValueScanner.readUntil(',', arrayValue);
                    tagValueArrayUInt16.emplace_back(static_cast<uint16_t>(stol(arrayValue)));
                }
                reads.TAGS.InsertArrayUint16TagAt(idx, tagKeyAndType, tagValueArrayUInt16);
            } else if (tagKeyAndType.getSecondaryTagType() == 'i') {
                tagValueArrayInt32.clear();
                while (tagValueScanner.len()) {
                    tagValueScanner.readUntil(',', arrayValue);
                    tagValueArrayInt32.emplace_back(static_cast<int32_t>(stol(arrayValue)));
                }
                reads.TAGS.InsertArrayInt32TagAt(idx, tagKeyAndType, tagValueArrayInt32);
            } else if (tagKeyAndType.getSecondaryTagType() == 'I') {
                tagValueArrayUInt32.clear();
                while (tagValueScanner.len()) {
                    tagValueScanner.readUntil(',', arrayValue);
                    tagValueArrayUInt32.emplace_back(static_cast<uint32_t>(stol(arrayValue)));
                }
                reads.TAGS.InsertArrayUint32TagAt(idx, tagKeyAndType, tagValueArrayUInt32);
            } else if (tagKeyAndType.getSecondaryTagType() == 'f') {
                tagValueArrayFloat.clear();
                while (tagValueScanner.len()) {
                    tagValueScanner.readUntil(',', arrayValue);
                    tagValueArrayFloat.emplace_back(static_cast<float>(stof(arrayValue)));
                }
                reads.TAGS.InsertArrayFloatTagAt(idx, tagKeyAndType, tagValueArrayFloat);
            } else {
                std::cerr << "Panic: type not recognized" << std::endl;
                exit(-1);
            }
        } else {
            std::cerr << "Panic: type not recognized" << std::endl;
            exit(-1);
        }
    }
}