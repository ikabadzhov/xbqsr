#ifndef GENOMICPIPELINE_VECTORCOLLECTION_H
#define GENOMICPIPELINE_VECTORCOLLECTION_H

#include <vector>
#include <span>
#include <cstdint>
#include <stdexcept>
#include <iostream>


struct SpanField {
public:
    size_t start;
    uint16_t len;
    explicit SpanField(size_t s = 0, uint16_t l = 0) : start(s), len(l) {}
};
struct SpanFieldLong {
public:
    size_t start;
    uint32_t len;
    explicit SpanFieldLong(size_t s = 0, uint32_t l = 0) : start(s), len(l) {}
};

template<typename T>
class VectorCollection {
    std::vector<T> data;
    std::vector<SpanField> ptrs;

public:
    void SetAt(std::vector<T> &str, size_t idx);
    void SetAt(std::span<T> str, size_t idx);
    void Resize(size_t sz);
    void Reserve(size_t sz);
    void ResizeData(size_t sz);
    void ReserveData(size_t sz);
    size_t Capacity();
    size_t CapacityData();
    std::span<T> operator[](size_t idx);
    void Permute(std::vector<size_t>& indexes);
    T* Data();
    std::vector<SpanField>& Ptrs();
    size_t size();
};

template<typename T>
std::vector<SpanField>& VectorCollection<T>::Ptrs() {
    return ptrs;
}

template<typename T>
size_t VectorCollection<T>::size() {
    return data.size();
}

template<typename T>
T* VectorCollection<T>::Data() {
    return data.data();
}

template<typename T>
void VectorCollection<T>::Resize(size_t sz) {
    ptrs.resize(sz);
}

template<typename T>
void VectorCollection<T>::Reserve(size_t sz) {
    ptrs.reserve(sz);
}

template<typename T>
void VectorCollection<T>::ResizeData(size_t sz) {
    data.resize(sz);
}

template<typename T>
void VectorCollection<T>::ReserveData(size_t sz) {
    data.reserve(sz);
}

template<typename T>
size_t VectorCollection<T>::Capacity() {
    return ptrs.capacity();
}

template<typename T>
size_t VectorCollection<T>::CapacityData() {
    return data.capacity();
}

template<typename T>
std::span<T> VectorCollection<T>::operator[](size_t idx) {
    if (idx >= ptrs.size()) {
        throw std::out_of_range("String index not available");
    }
    return {data.begin()+ptrs[idx].start, data.begin()+ptrs[idx].start+ptrs[idx].len};

}

template<typename T>
void VectorCollection<T>::SetAt(std::vector<T> &str, size_t idx) {
    size_t lastIdx = data.size();
    data.resize(data.size()+str.size());
    std::copy(str.begin(), str.end(), data.begin()+lastIdx);
    ptrs[idx].start = lastIdx;
    ptrs[idx].len = str.size();
}

template<typename T>
void VectorCollection<T>::SetAt(std::span<T> str, size_t idx) {
    size_t lastIdx = data.size();
    data.resize(data.size()+str.size());
    std::copy(str.begin(), str.end(), data.begin()+lastIdx);
    ptrs[idx].start = lastIdx;
    ptrs[idx].len = str.size();
}

template<typename T>
void VectorCollection<T>::Permute(std::vector<size_t>& indexes) {
    VectorCollection<T> sorted_collection;
    sorted_collection.Resize(indexes.size());

    size_t outIdx = 0;
    for (auto inIdx : indexes) {
        sorted_collection.SetAt(this->operator[](inIdx), outIdx);
        outIdx++;
    }

    // Swap collections
    data = std::move(sorted_collection.data);
    ptrs = std::move(sorted_collection.ptrs);
}


#endif //GENOMICPIPELINE_VECTORCOLLECTION_H
