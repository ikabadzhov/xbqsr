#ifndef GENOMICPIPELINE_BAMALIGNMENTFORMATTER_H
#define GENOMICPIPELINE_BAMALIGNMENTFORMATTER_H


#include "../../sam_data_structures/SamAlignments.h"

class BamAlignmentFormatter {
public:
    static void FormatBamAlignment(SamAlignments& alns, ssize_t idx, std::string& formatted_line);
    static uint8_t encodeNucleotides(char nucleotide1, char nucleotide2);
};


#endif //GENOMICPIPELINE_BAMALIGNMENTFORMATTER_H
