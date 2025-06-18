#include <string>
#include <vector>
#include "TagsCollection.h"


void TagsCollection<std::string>::Resize(size_t sz) {
    ALN_SZ = sz;
    for (auto& [str, tagIdx] : tagToValueIndex) {
        tagValues[tagIdx].Resize(ALN_SZ);
    }
}

void TagsCollection<std::string>::AddTagIfNeeded(uint16_t TAG) {
    if (tagToValueIndex.count(TAG) <= 0) {
        //intToStr.emplace_back(TAG);
        tagValues.emplace_back();
        tagToValueIndex[TAG] = tagValues.size() - 1;
        tagValues[tagValues.size() - 1].Resize(ALN_SZ);
        //validity[tagValue.size() - 1].Resize(ALN_SZ);
    }
}

void TagsCollection<std::string>::AssignTagValue(size_t readIdx, uint16_t tagKey, const std::string& tagVal) {
    int32_t tagIdx = tagToValueIndex[tagKey];
//    if (readIdx >= tagValues.size()) {
//        throw std::out_of_range("String index not available");
//    }
    tagValues[tagIdx].setAt(readIdx, tagVal);
}

SparseVector<std::string>& TagsCollection<std::string>::GetTags(const std::string &tag) {
    auto tagKey = TagKeyAndType::encode(tag);
    int32_t tagIdx = tagToValueIndex[tagKey];
    return tagValues[tagIdx];
}

bool TagsCollection<std::string>::HasTAG(uint16_t TAG, size_t readIdx) {
    if(tagToValueIndex.contains(TAG)) {
        int32_t tagID = tagToValueIndex[TAG];
        return tagValues[tagID].isValid(readIdx);
    }
    return false;
}

void TagsCollection<std::string>::Permute(std::vector<size_t>& indices) {
    for (auto & tagValue : tagValues) {
        tagValue.Permute(indices);
    }
}

