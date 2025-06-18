#ifndef GENOMICPIPELINE_TAGSCOLLECTION_H
#define GENOMICPIPELINE_TAGSCOLLECTION_H

#include <vector>
#include <span>
#include <string>
#include <unordered_map>
#include "BitVector.h"
#include "../Utils.h"
#include "StringCollection.h"
#include "VectorCollection.h"
#include "SparseVector.h"
#include "../sam_data_structures/sam_fields/Tag.h"

/*
 * It's a collection storing all tags whose type is a simple type such as int, float, char.
 *
 * Sparse vector contains a validity bit and the vector containing the actual values
 * */
template<typename T>
class TagsCollection {
    size_t ALN_SZ = 0;

    std::vector<SparseVector<T>> tagValues;
    std::unordered_map<uint16_t, int32_t> tagToValueIndex;

public:
    void AddTagIfNeeded(uint16_t TAG);
    void AssignTagValue(size_t readIdx, uint16_t tagKey, T tagVal);

    SparseVector<T>& GetTags(const std::string& TAG);

    bool HasTAG(uint16_t TAG, size_t readIdx);
    void Resize(size_t sz);
    void Permute(std::vector<size_t>& indices);
};

template<typename T>
void TagsCollection<T>::Permute(std::vector<size_t> &indices) {
    for (size_t i = 0; i < tagValues.size(); i++) {
        tagValues[i].Permute(indices);
    }
}

template<typename T>
bool TagsCollection<T>::HasTAG(uint16_t TAG, size_t readIdx) {
    int32_t tagID = tagToValueIndex[TAG];
    return tagValues[tagID].is_valid(readIdx);
}


template<typename T>
void TagsCollection<T>::Resize(size_t sz) {
    ALN_SZ = sz;
    for (auto& [str, tagIdx] : tagToValueIndex) {
        tagValues[tagIdx].Resize(ALN_SZ);
    }
}

template<typename T>
void TagsCollection<T>::AssignTagValue(size_t readIdx, uint16_t tagKey, T val) {
    int32_t tagIdx = tagToValueIndex[tagKey];
    if (readIdx >= tagValues[tagIdx].Size()) {
        throw std::out_of_range("String index not available");
    }
    tagValues[tagIdx].setAt(readIdx, val);
}

template<typename T>
void TagsCollection<T>::AddTagIfNeeded(uint16_t TAG) {
    if (tagToValueIndex.count(TAG) <= 0) {
        tagValues.emplace_back();
        tagToValueIndex[TAG] = tagValues.size() - 1;
        tagValues.back().Resize(ALN_SZ);
    }
};


template<typename T>
SparseVector<T> &TagsCollection<T>::GetTags(const std::string &TAG) {
    auto tagKey = TagKeyAndType::encode(TAG);
    return tagValues[tagToValueIndex[tagKey]];
}

/*
 * Specialization for tags containing strings
 * */
template<>
class TagsCollection<std::string> {
    size_t ALN_SZ = 0;
    std::unordered_map<uint16_t, int32_t> tagToValueIndex;
    std::vector<SparseVector<std::string>> tagValues;

public:
    void Resize(size_t sz);
    void AddTagIfNeeded(uint16_t TAG);
    void AssignTagValue(size_t idx, uint16_t tagKey, const std::string& tagVal);
    SparseVector<std::string>& GetTags(const std::string& TAG);
    //Tag(std::string_view& TAG, int64_t val);
    bool HasTAG(uint16_t TAG, size_t readIdx);
    void Permute(std::vector<size_t>& indices);
};

/*
 * Class template for tags containing arrays
 * */
template<typename T>
class TagsCollectionArray {
    size_t ALN_SZ = 0;
    std::unordered_map<uint16_t, int32_t> tagToValueIndex;
    std::vector<SparseVector<std::vector<T>>> tagValues;

public:
    void Resize(size_t sz);
    void AddTagIfNeeded(uint16_t TAG);
    void AssignTagValue(size_t readIdx, uint16_t tagKey, const std::vector<T>& tagVal);
    SparseVector<std::vector<T>>& GetTags(const std::string& TAG);
    bool HasTAG(uint16_t TAG, size_t readIdx);
    void Permute(std::vector<size_t>& indices);
};

template<typename T>
void TagsCollectionArray<T>::AssignTagValue(size_t readIdx, uint16_t tagKey, const std::vector<T>& tagVal) {
    int32_t tagIdx = tagToValueIndex[tagKey];
    if (readIdx >= tagValues[tagIdx].Size()) {
        throw std::out_of_range("Read index not available");
    }
    tagValues[tagIdx].setAt(readIdx, tagVal);
}

template<typename T>
SparseVector<std::vector<T>>& TagsCollectionArray<T>::GetTags(const std::string& TAG) {
    auto tagKey = TagKeyAndType::encode(TAG);
    int32_t tagIdx = tagToValueIndex[tagKey];
    return tagValues[tagIdx];
}

template<typename T>
bool TagsCollectionArray<T>::HasTAG(uint16_t TAG, size_t readIdx) {
    int32_t tagID = tagToValueIndex[TAG];
    return tagValues[tagID].isValid(readIdx);
}

template<typename T>
void TagsCollectionArray<T>::Permute(std::vector<size_t> &indices) {
    for (size_t i = 0; i < tagValues.size(); i++) {
        tagValues[i].Permute(indices);
    }
}

template<typename T>
void TagsCollectionArray<T>::Resize(size_t sz) {
    ALN_SZ = sz;
    for (auto& [str, tagIdx] : tagToValueIndex) {
        tagValues[tagIdx].Resize(ALN_SZ);
    }
}

template<typename T>
void TagsCollectionArray<T>::AddTagIfNeeded(uint16_t TAG) {
    if (tagToValueIndex.count(TAG) <= 0) {
        tagValues.emplace_back();
        tagToValueIndex[TAG] = tagValues.size() - 1;
        tagValues.back().Resize(ALN_SZ);
    }
}

#endif //GENOMICPIPELINE_TAGSCOLLECTION_H
