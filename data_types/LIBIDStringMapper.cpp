#include "LIBIDStringMapper.h"
#include <stdexcept>

int32_t LIBIDStringMapper::addString(const std::string& str) {
    auto it = stringToId.find(str);
    if (it != stringToId.end()) {
        return it->second; // String already in map
    }

    // Insert new string
    idToString.emplace_back(str);
    stringToId[str] = nextId;
    return nextId++;
}

int32_t LIBIDStringMapper::getStringId(const std::string& str) const {
    if (str.empty()) return EMPTY_STRING;
    auto it = stringToId.find(str);
    if (it != stringToId.end()) {
        return it->second;
    }
    return STRING_NOT_FOUND;
}

std::string LIBIDStringMapper::getString(int32_t id) const {
    if (id == EMPTY_STRING) return "";
    if (id>= idToString.size()) {
        throw std::out_of_range ("LIBIDStringMapper: Out Of Range");
    }
    return idToString[id];
}


int32_t LIBIDStringMapper::size() const {
    return nextId;
}
