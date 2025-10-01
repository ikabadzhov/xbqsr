#ifndef GENOMICPIPELINE_BAMALIGNMENTPARSER_H
#define GENOMICPIPELINE_BAMALIGNMENTPARSER_H

#include <string_view>
#include "../../sam_data_structures/SamAlignments.h"

class VariableLenFieldsVar {
public:
    size_t name;
    size_t cigar;
    size_t seq;
    size_t qual;
    size_t aux;
};



class BamAlignmentParser {
public:
    static void ParseBamAlignmentFromSlice(std::string_view bytes, SamAlignments &reads, size_t idx);

    static void ParseBamAlignmentFromSliceSimpleTypes(std::string_view bytes, SamAlignments &reads, size_t idx, VariableLenFieldsVar& fieldsOff, VariableLenFieldsVar& fieldsSizes);

    static void ParseBamAlignmentFromSliceComplexTypes(std::string_view bytes, SamAlignments &reads, size_t idx,
                                                VariableLenFieldsVar &fieldsOff, VariableLenFieldsVar& fieldsSizes, VariableLenFieldsVar& outputOff);
};


#endif //GENOMICPIPELINE_BAMALIGNMENTPARSER_H
