#ifndef GENOMICPIPELINE_READGROUPFIELD_H
#define GENOMICPIPELINE_READGROUPFIELD_H

#include <unordered_map>
#include <vector>
#include <string>
#include "SamHeaderRecord.h"


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
    int32_t nextId = 0;

public:
    static const int EMPTY_STRING = -1;
    static const int STRING_NOT_FOUND = -2;

    int32_t add(const ReadGroupRecord& readGroup);
    int32_t getId(const std::string& str) const;
    std::string getString(int32_t id) const;
    ReadGroupRecord getReadGroupRecord(int32_t id) const;
    const std::vector<ReadGroupRecord>& getAll() const;
    int32_t size() const;

    static bool isEmpty(int32_t id);
};


#endif //GENOMICPIPELINE_READGROUPFIELD_H
