#include "CompressedStringCollection.h"
#include "../Utils.h"
#include <stdexcept>


CompressedStringCollection::CompressedStringCollection(std::string_view defaultString) {
    intToStr.emplace_back(defaultString);
    strToInt[std::string{defaultString}] = intToStr.size()-1;
}


void CompressedStringCollection::Set(size_t idx, const std::string &rg) {
    if (idx>=values.size()) {
        throw std::out_of_range ("CompressedStringCollection: Out Of Range");
    }
    if (rg.empty()) {
        values[idx] = EMPTY_STRING;
    } else {
        if (strToInt.count(rg) <= 0) {
            intToStr.emplace_back(rg);
            strToInt[rg] = static_cast<int32_t>(intToStr.size()) - 1;
        }
        values[idx] = strToInt[rg];
    }
}

int32_t CompressedStringCollection::Get(size_t idx) const {
    if (idx>=values.size()) {
        throw std::out_of_range ("CompressedStringCollection: Out Of Range");
    }
    return values[idx];
}

void CompressedStringCollection::Resize(size_t sz) {
    values.resize(sz);
}

std::string CompressedStringCollection::toString(size_t idx) const {
    if (idx >= values.size()) {
        throw std::out_of_range ("CompressedStringCollection: Out Of Range");
    }
    int32_t v = values[idx];
    return (v == EMPTY_STRING) ? "" : intToStr[v];
}

std::vector<int32_t> &CompressedStringCollection::toVec() {
    return values;
}

bool CompressedStringCollection::isEmpty(int32_t val) {
    return val == EMPTY_STRING;
}

int32_t CompressedStringCollection::StringToInt(std::string& str) {
    if (strToInt.count(str)<=0) {
        throw std::out_of_range ("CompressedStringCollection: Dictionary entry not there");
    }
    return strToInt[str];
}

std::string& CompressedStringCollection::IntToString(int32_t i) {
    if (i >= intToStr.size()) {
        throw std::out_of_range ("CompressedStringCollection: Dictionary entry not there");
    }
    return intToStr[i];
}

size_t CompressedStringCollection::NumEntriesTable() {
    return intToStr.size();
}

void CompressedStringCollection::Permute(std::vector<size_t> &indices) {
    permute(indices, values);
}
