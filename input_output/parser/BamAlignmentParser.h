#ifndef GENOMICPIPELINE_BAMALIGNMENTPARSER_H
#define GENOMICPIPELINE_BAMALIGNMENTPARSER_H

#include <string_view>
#include "../../sam_data_structures/SamAlignments.h"

class VariableLenFieldsVar {
public:
    size_t qname{0};
    size_t cigar{0};
    size_t seq{0};
    size_t qual{0};
    size_t aux{0};

    VariableLenFieldsVar() = default;
    VariableLenFieldsVar(size_t q, size_t c, size_t s, size_t qu, size_t a)
        : qname(q), cigar(c), seq(s), qual(qu), aux(a) {}
};

class BamAlignmentParser {
public:
    static void ParseBamAlignmentFromSlice(std::string_view bytes, SamAlignments &reads, size_t idx);

    static void ParseBamAlignmentFromSliceSimpleTypes(std::string_view bytes, SamAlignments &reads, size_t idx, VariableLenFieldsVar& fieldsOff, VariableLenFieldsVar& fieldsSizes);

    static void ParseBamAlignmentFromSliceComplexTypes(std::string_view bytes, SamAlignments &reads, size_t idx, VariableLenFieldsVar const& fieldsOff, VariableLenFieldsVar const& fieldsSizes, VariableLenFieldsVar const& outOff);    // pass TAGS so we can write at fixed offsets
};


#endif //GENOMICPIPELINE_BAMALIGNMENTPARSER_H
