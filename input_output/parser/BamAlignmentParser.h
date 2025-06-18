#ifndef GENOMICPIPELINE_BAMALIGNMENTPARSER_H
#define GENOMICPIPELINE_BAMALIGNMENTPARSER_H

#include <string_view>
#include "../../sam_data_structures/SamAlignments.h"

class BamAlignmentParser {
public:
    static void ParseBamAlignmentFromSlice(std::string_view bytes, SamAlignments &reads, size_t idx);
};


#endif //GENOMICPIPELINE_BAMALIGNMENTPARSER_H
