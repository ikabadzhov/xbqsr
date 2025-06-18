#ifndef GENOMICPIPELINE_TMPFILEDS_H
#define GENOMICPIPELINE_TMPFILEDS_H

#include <string>
#include <vector>
#include <unordered_map>
#include "../../data_types/CompressedStringCollection.h"

class TmpFields {
public:
//    std::vector<std::string> libID;
    std::vector<int32_t> LIBID;
    std::vector<int32_t> REFID;
    std::vector<int32_t> nextREFID;
    std::vector<int32_t> adaptPos;
    std::vector<int32_t> adaptScore;
    CompressedStringCollection reads_group_covariate_;
    void Permute(std::vector<size_t>& indices);
};

#endif //GENOMICPIPELINE_TMPFILEDS_H
