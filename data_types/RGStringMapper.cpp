#include <stdexcept>
#include "RGStringMapper.h"


int32_t RGStringMapper::addString(const std::string& str) {
    auto it = stringToId.find(str);
    if (it != stringToId.end()) {
        return it->second; // String already in map
    }

    // Insert new string
    idToString.emplace_back(str);
    stringToId[str] = nextId;
    return nextId++;
}

int32_t RGStringMapper::getStringId(const std::string& str) const {
    if (str.empty()) return EMPTY_STRING;

    auto it = stringToId.find(str);
    if (it != stringToId.end()) {
        return it->second;
    }
    return STRING_NOT_FOUND; // String not found
}

std::string RGStringMapper::getString(int32_t id) const {
    if (id == EMPTY_STRING) return "";
    if (id>= idToString.size()) {
        throw std::out_of_range ("RGStringMapper: Out Of Range");
    }
    return idToString[id];
}


int32_t RGStringMapper::size() const {
    return nextId;
}

bool RGStringMapper::isEmpty(int32_t val) {
    return val==EMPTY_STRING;
}

