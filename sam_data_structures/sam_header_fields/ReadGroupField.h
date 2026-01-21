#pragma once
#include <unordered_map>
#include <vector>
#include <string>

#include "SamHeaderRecord.h"
#include "../../data_types/StringDictionary.h"

class ReadGroupRecord : public SamHeaderRecord {
    std::string readGroupId;
public:
    explicit ReadGroupRecord(std::string readGroupId) : SamHeaderRecord("RG"), readGroupId(std::move(readGroupId)) {};
    const std::string& getReadGroupId() const;
};

class ReadGroupsDictionary {
    std::vector<ReadGroupRecord> readGroups;
    std::unordered_map<std::string, int> stringToId;
    std::vector<std::string> idToString;

    // Encoded PU dictionary and RG->PU mapping
    StringDictionary pu_;

    std::vector<int32_t> rgId2PuId_;

    int32_t nextId = 0;

public:
    static constexpr int EMPTY_STRING = -1;
    static constexpr int STRING_NOT_FOUND = -2;

    ReadGroupsDictionary() = default;

    int32_t add(const ReadGroupRecord& readGroup);
    int32_t getId(const std::string& str) const;

    const std::string&      getString(int32_t id) const;
    const ReadGroupRecord&  getReadGroupRecord(int32_t id) const;
    const std::vector<ReadGroupRecord>& getAll() const noexcept { return readGroups; }

    // RG -> PU mapping accessors
    void     setPUForRG(int32_t rg_id, int32_t pu_id);
    int32_t  getPUForRG(int32_t rg_id) const;

    int32_t size()  const noexcept { return static_cast<int32_t>(readGroups.size()); }
    bool    empty() const noexcept { return readGroups.empty(); }

    static constexpr bool isEmpty(int32_t id) noexcept { return id == EMPTY_STRING; }

    // Expose PU dictionary only if you need it elsewhere:
    const StringDictionary& platformUnits() const noexcept { return pu_; }
    StringDictionary&       platformUnits()       noexcept { return pu_; }
};
