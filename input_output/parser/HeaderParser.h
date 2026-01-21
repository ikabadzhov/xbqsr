#ifndef GENOMICPIPELINE_SAMHEADERPARSER_H
#define GENOMICPIPELINE_SAMHEADERPARSER_H

#include "ParsedHeaderLine.h"
#include "../../sam_data_structures/SamHeader.h"

namespace SamHeaderParser {
    void ParseSamHeaderFromSlice(Header &hdr, const std::string& bytes, bool first = false);
    void AddAttributes(SamHeaderRecord& record, ParsedHeaderLine& parsedHeaderLine);
    void ParseHDLine(Header &hdr, ParsedHeaderLine& parsedLine);
    void ParsePGLine(Header &hdr, ParsedHeaderLine& parsedLine);
    void ParseRGLine(Header &hdr, ParsedHeaderLine& parsedLine);
    void ParseSQLine(Header &hdr, ParsedHeaderLine& parsedLine);
};


#endif //GENOMICPIPELINE_SAMHEADERPARSER_H
