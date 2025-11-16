#ifndef GENOMICPIPELINE_STRINGCOLLECTION_H
#define GENOMICPIPELINE_STRINGCOLLECTION_H

#include <vector>
#include <string>
#include <string_view>
#include <cstdint>

struct StringField {
public:
    size_t start;
    uint32_t len;
    explicit StringField(size_t s = 0, uint32_t l = 0) : start(s), len(l) {}
};


class StringCollection {


public:
    std::string data;
    std::vector<StringField> ptrs;
    ssize_t size = 0;


    void SetAt(std::string_view str, size_t idx);
    void Resize(size_t sz);

    void Reserve(size_t sz);

    std::string_view operator[](size_t idx) const;
    char* Data();
    StringField* Ptrs();
    void Permute(std::vector<size_t>& indexes);


    void ResizeData(size_t sz) {
        data.resize(sz);
    }

    void ReserveData(size_t sz) {
        data.reserve(sz);
    }

    size_t Capacity() {
      return ptrs.capacity();
    }

    size_t CapacityData() {
        return data.capacity();
    }
};




#endif //GENOMICPIPELINE_STRINGCOLLECTION_H
