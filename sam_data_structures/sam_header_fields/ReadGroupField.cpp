#include <stdexcept>
#include "ReadGroupField.h"

int32_t ReadGroupsDictionary::add(const ReadGroupRecord& readGroup) {
    auto it = stringToId.find(readGroup.getReadGroupId());
    if (it != stringToId.end()) {
        return it->second; // String already in map
    }
    readGroups.emplace_back(readGroup);
    idToString.emplace_back(readGroup.getReadGroupId());
    stringToId[readGroup.getReadGroupId()] = nextId;
    return nextId++;
}

int32_t ReadGroupsDictionary::getId(const std::string& str) const {
    if (str.empty()) return EMPTY_STRING;

    auto it = stringToId.find(str);
    if (it != stringToId.end()) {
        return it->second;
    }
    return STRING_NOT_FOUND; // String not found
}

std::string ReadGroupsDictionary::getString(int32_t id) const {
    //if (id == EMPTY_STRING) return "";
    if (id>= idToString.size()) {
        throw std::out_of_range ("ReadGroupsDictionary: Out Of Range");
    }
    return idToString[id];
}

ReadGroupRecord ReadGroupsDictionary::getReadGroupRecord(int32_t id) const {
    //if (id == EMPTY_STRING) return {};
    if (id >= idToString.size()) {
        throw std::out_of_range ("ReadGroupsDictionary: Out Of Range");
    }
    return readGroups[id];
}


int32_t ReadGroupsDictionary::size() const {
    return nextId;
}

bool ReadGroupsDictionary::isEmpty(int32_t id) {
    return id==EMPTY_STRING;
}

const std::vector<ReadGroupRecord> &ReadGroupsDictionary::getAll() const {
    return readGroups;
}

const std::string& ReadGroupRecord::getReadGroupId() const {
    return readGroupId;
}
