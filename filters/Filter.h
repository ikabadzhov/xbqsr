#ifndef GENOMICPIPELINE_FILTER_H
#define GENOMICPIPELINE_FILTER_H


#include <functional>
#include <memory>
#include <sstream>
#include "../sam_data_structures/SamHeader.h"
#include "../sam_data_structures/SamAlignments.h"

class Filter {

public:
    //static std::function<bool (Alignment& aln)> AddREFID(Header& header);
    static std::function<bool (SamAlignments&, size_t)> AddREFID(SamAlignments& reads);
    static std::function<bool (SamAlignments&, size_t)> AddNewTAG(TagKeyAndType tagKeyAndType, const std::string &tagValue);
    //static std::function<bool (Alignment& aln)> RemoveUnmappedReads(Header& header);
};


#endif //GENOMICPIPELINE_FILTER_H