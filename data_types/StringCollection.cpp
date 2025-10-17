#include "StringCollection.h"
#include <stdexcept>
#include <iostream>


void StringCollection::SetAt(std::string_view str, size_t idx) {
    size_t lastIdx = data.size();
    data.resize(data.size()+str.size());
    std::copy(str.begin(), str.end(), data.begin()+lastIdx);
    ptrs[idx] = StringField(lastIdx, str.size());
}

std::string_view StringCollection::operator[](size_t idx) const {
    if (idx >= ptrs.size()) {
        throw std::out_of_range("String index not available");
    }
    return {data.begin()+ptrs[idx].start, data.begin()+ptrs[idx].start+ptrs[idx].len};
}

char * StringCollection::Data() {
    return data.data();
}

StringField * StringCollection::Ptrs() {
    return ptrs.data();
}

void StringCollection::Resize(size_t sz) {
    size = sz;
    ptrs.resize(sz);
}


void StringCollection::Permute(std::vector<size_t>& indexes) {
    StringCollection sorted_collection;
    //sorted_collection.Resize(indexes.size());
    sorted_collection.data.resize(this->data.size());
    sorted_collection.ptrs.resize(this->ptrs.size());
    sorted_collection.size = this->size;

    size_t j = 0;
    size_t lastIdx = 0;
    for (auto idx : indexes) {
        auto start = this->data.begin() + this->ptrs[idx].start;
        auto end = this->data.begin() + this->ptrs[idx].start + this->ptrs[idx].len;
        std::copy(start, end, sorted_collection.data.begin()+lastIdx);

        sorted_collection.ptrs[j].start = lastIdx;
        sorted_collection.ptrs[j].len = this->operator[](idx).size();

        lastIdx += this->operator[](idx).size();
        j++;
    }

    // Swap collections
    data = std::move(sorted_collection.data);
    ptrs = std::move(sorted_collection.ptrs);
}
