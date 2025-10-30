#ifndef GENOMICPIPELINE_SAMALIGNMENTFORMATTER_H
#define GENOMICPIPELINE_SAMALIGNMENTFORMATTER_H


#include <string>
#include "../../sam_data_structures/SamAlignments.h"

class SamAlignmentFormatter {

public:
    static void FormatSamAlignment(SamAlignments& alns, ssize_t idx, std::string& line);

    static size_t GetSamAlignmentSize(SamAlignments &alns, ssize_t idx);

    static void FormatSamAlignmentInBuffer(SamAlignments &alns, ssize_t idx, char *out);
};


#endif //GENOMICPIPELINE_SAMALIGNMENTFORMATTER_H
