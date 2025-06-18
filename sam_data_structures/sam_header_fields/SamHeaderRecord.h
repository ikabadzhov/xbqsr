#ifndef GENOMICPIPELINE_SAMPROGRAMHEADER_H
#define GENOMICPIPELINE_SAMPROGRAMHEADER_H

#include <string>
#include "../../data_types/LinkedHashMap.h"


enum HeaderRecordType {
    HD, SQ, RG, PG, CO
};


class SamHeaderRecord {
protected:
    HeaderRecordType headerRecordTypes;
    LinkedHashMap<std::string, std::string> attributes;
public:
    explicit SamHeaderRecord(const std::string& recordType);

    void putAttribute(const std::string& key, const std::string& value);
    std::string getAttribute(const std::string& key) const;
    bool containsAttribute(const std::string& key) const;

    const LinkedHashMap<std::string, std::string>&  getAllAttributes() const;

    HeaderRecordType getRecordType();
    std::string getRecordTypeAsString();
};

#endif //GENOMICPIPELINE_SAMPROGRAMHEADER_H
