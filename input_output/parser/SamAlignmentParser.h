#ifndef GENOMICPIPELINE_SAMALIGNMENTPARSER_H
#define GENOMICPIPELINE_SAMALIGNMENTPARSER_H


#include "../../sam_data_structures/SamAlignments.h"
#include "../utils/StringScanner.h"

namespace SamAlignmentParser {
    void ParseSamAlignmentFromSlice(std::string_view s, SamAlignments &reads, size_t idx);
    void ParseSamAlignment(StringScanner &sc, SamAlignments &reads, size_t idx);
};


#endif //GENOMICPIPELINE_SAMALIGNMENTPARSER_H
