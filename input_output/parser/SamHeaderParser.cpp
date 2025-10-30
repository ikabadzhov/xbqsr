#include "SamHeaderParser.h"
#include "ParsedHeaderLine.h"

/*
 * ParseSamHeader parses a complete header in a SAM file. See
 * http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.3.
 */
void SamHeaderParser::ParseSamHeaderFromSlice(Header &hdr, const std::string& bytes, bool first) {
    ParsedHeaderLine parsedLine(bytes);
    if (parsedLine.isLineValid()) {
        switch (parsedLine.getHeaderRecordType()) {
            case HeaderRecordType::HD:
                ParseHDLine(hdr, parsedLine);
                break;
            case HeaderRecordType::SQ:
                ParseSQLine(hdr, parsedLine);
                break;
            case HeaderRecordType::RG:
                ParseRGLine(hdr, parsedLine);
                break;
            case HeaderRecordType::PG:
                ParsePGLine(hdr, parsedLine);
                break;
            case HeaderRecordType::CO:
                hdr.addComment(bytes);
                break;
            default:
                std::cerr << "Unrecognized header record type" << std::endl;
                exit(-1);
        }
    }
}



void SamHeaderParser::ParseHDLine(Header &hdr, ParsedHeaderLine& parsedLine) {
    parsedLine.requireTag("VN");
    std::string soString;
    std::string goString;

    if (parsedLine.containsKey("SO")) {
        soString = parsedLine.getValue("SO");
    }
    if (parsedLine.containsKey("GO")) {
        goString = parsedLine.getValue("GO");
    }
    AddAttributes(hdr.getHD(), parsedLine);
}

void SamHeaderParser::AddAttributes(SamHeaderRecord &record, ParsedHeaderLine &parsedHeaderLine) {
    parsedHeaderLine.keyValuePairs.for_each([&](const auto& key, const auto& value){
        record.putAttribute(key, value);
    });
}

void SamHeaderParser::ParsePGLine(Header &hdr, ParsedHeaderLine &parsedLine) {
    assert(parsedLine.getHeaderRecordType() == HeaderRecordType::PG);

    if (parsedLine.requireTag("ID")) {
        ProgramRecord programRecord(parsedLine.removeValue("ID"));
        AddAttributes(programRecord, parsedLine);
        hdr.addProgramRecord(programRecord);
    }
}

void SamHeaderParser::ParseRGLine(Header &hdr, ParsedHeaderLine &parsedLine) {
    assert(parsedLine.getHeaderRecordType() == HeaderRecordType::RG);

    if (parsedLine.requireTag("ID")) {
        parsedLine.requireTag("SM");
        ReadGroupRecord samReadGroupRecord(parsedLine.removeValue("ID"));
        AddAttributes(samReadGroupRecord, parsedLine);
        //auto predictedMedianInsertSize = ;
        if (samReadGroupRecord.containsAttribute("PI")) {
            try {
                const auto& s = samReadGroupRecord.getAttribute("PI");
                std::stoul(s, nullptr, 10);
            } catch (...) {
                throw std::runtime_error("PI is not numeric: " + samReadGroupRecord.getAttribute("PI"));
            }
        }

        /*std::string dateRunProduced = samReadGroupRecord.getAttribute("DT");
          if (!dateRunProduced.empty()) {
            try {
                date = this.mTagCodec.decodeDate(dateRunProduced);
            } catch (DateParser.InvalidDateException var7) {
                date = dateRunProduced;
                this.reportErrorParsingLine("DT tagValue value '" + dateRunProduced + "' is not parseable as a date", Type.INVALID_DATE_STRING, var7);
            }

            samReadGroupRecord.setAttribute("DT", date.toString());
        }*/

        hdr.addReadGroupRecord(samReadGroupRecord);
    }
}

void SamHeaderParser::ParseSQLine(Header &hdr, ParsedHeaderLine &parsedLine) {
    assert(parsedLine.getHeaderRecordType() == HeaderRecordType::SQ);

    if (parsedLine.requireTag("SN") && parsedLine.requireTag("LN")){
        std::string name = parsedLine.removeValue("SN");
        auto pos1 = name.find('\t');
        auto pos2 = name.find('\n');
        auto minPos = std::min(pos1,pos2);
        auto pos3 = name.find(' ');
        minPos = std::min(minPos, pos3);
        std::string lenStr = parsedLine.removeValue("LN");
        auto len = static_cast<int32_t>(strtol(lenStr.c_str(), nullptr, 10));
        SequenceRecord sequenceRecord(name.substr(0, minPos), len);
        AddAttributes(sequenceRecord, parsedLine);
        hdr.addSequenceRecord(sequenceRecord);
    }
}