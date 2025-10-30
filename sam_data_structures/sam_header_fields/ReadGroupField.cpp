#include <stdexcept>
#include "ReadGroupField.h"

int32_t ReadGroupsDictionary::add(const ReadGroupRecord& rg) {
    const std::string& key = rg.getReadGroupId();
    if (auto it = stringToId.find(key); it != stringToId.end())
        return it->second;

    const int32_t id = nextId;

    readGroups.push_back(rg);
    stringToId.emplace(key, id);
    idToString.push_back(key);

    // Maintain RG -> PU map in lockstep with readGroups
    if (rg.containsAttribute("PU")) {
        const std::string& pu = rg.getAttribute("PU");     // <-- correct key
        const int32_t pu_id = pu_.add(pu);
        rgId2PuId_.push_back(pu_id);                       // single write here
    } else {
        rgId2PuId_.push_back(EMPTY_STRING);
    }

    ++nextId;
    return id;
}

int32_t ReadGroupsDictionary::getId(const std::string& s) const {
    if (s.empty()) return EMPTY_STRING;
    if (auto it = stringToId.find(s); it != stringToId.end())
        return it->second;
    return STRING_NOT_FOUND;
}

const std::string& ReadGroupsDictionary::getString(int32_t id) const {
    if (id < 0 || static_cast<size_t>(id) >= idToString.size())
        throw std::out_of_range("ReadGroupsDictionary::getString: id out of range");
    return idToString[static_cast<size_t>(id)];
}

const ReadGroupRecord& ReadGroupsDictionary::getReadGroupRecord(int32_t id) const {
    if (id < 0 || static_cast<size_t>(id) >= readGroups.size())
        throw std::out_of_range("ReadGroupsDictionary::getReadGroupRecord: id out of range");
    return readGroups[static_cast<size_t>(id)];
}

void ReadGroupsDictionary::setPUForRG(int32_t rg_id, int32_t pu_id) {
    if (rg_id < 0 || static_cast<size_t>(rg_id) >= rgId2PuId_.size())
        throw std::out_of_range("setPUForRG: rg id out of range");
    rgId2PuId_[static_cast<size_t>(rg_id)] = pu_id;
}

int32_t ReadGroupsDictionary::getPUForRG(int32_t rg_id) const {
    if (rg_id < 0 || static_cast<size_t>(rg_id) >= rgId2PuId_.size())
        return EMPTY_STRING;
    return rgId2PuId_[static_cast<size_t>(rg_id)];
}

const std::string& ReadGroupRecord::getReadGroupId() const {
    return readGroupId;
}

