#ifndef GENOMICPIPELINE_SPARSEVECTOR_H
#define GENOMICPIPELINE_SPARSEVECTOR_H

#include "BitVector.h"
#include "../Utils.h"

template<typename T>
class SparseVector {
    std::vector<T> values;
    BitVector validity;
public:
    void Resize(size_t sz);
    void Permute(std::vector<size_t>& indices);
    bool isValid(size_t idx);
    void setAt(size_t idx, T val);
    T GetAt(size_t idx);
    size_t Size();
};

template<typename T>
T SparseVector<T>::GetAt(size_t idx) {
    if (!validity.GetBit(idx)) {
        throw std::invalid_argument("Alignment " + std::to_string(idx) + " does not contain values for the give tag");
    }
    return values[idx];
}

template<typename T>
size_t SparseVector<T>::Size() {
    return values.size();
}

template<typename T>
bool SparseVector<T>::isValid(size_t idx) {
    return validity.GetBit(idx);
}

template<typename T>
void SparseVector<T>::setAt(size_t idx, T val) {
    values[idx] = val;
    validity.SetBit(idx, true);
}

template<typename T>
void SparseVector<T>::Permute(std::vector<size_t> &indices) {
    permute<>(indices, values);
    BitVector val(validity.SizeInBits(), true);
    for (size_t j = 0; j < indices.size(); j++) {
        val.SetBit(j, validity.GetBit(indices[j]));
    }
    validity = val;
}

template<typename T>
void SparseVector<T>::Resize(size_t sz) {
    values.resize(sz);
    validity.Resize(sz);
}


#endif //GENOMICPIPELINE_SPARSEVECTOR_H
