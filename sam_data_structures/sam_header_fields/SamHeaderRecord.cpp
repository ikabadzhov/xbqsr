#include <iostream>
#include "SamHeaderRecord.h"


const std::vector<std::string> recordTypesAsString = {
        "HD", "SQ", "RG", "PG", "CO"
};

const std::unordered_map<std::string, HeaderRecordType> recordTypesAsVal = {
        {"HD", HeaderRecordType::HD},
        {"SQ", HeaderRecordType::SQ},
        {"RG", HeaderRecordType::RG},
        {"PG", HeaderRecordType::PG},
        {"CO", HeaderRecordType::CO}
};


void SamHeaderRecord::putAttribute(const std::string &key, const std::string &value) {
    attributes.put(key, value);
}

std::string SamHeaderRecord::getAttribute(const std::string &key) const {
    return attributes.get(key);
}

const LinkedHashMap<std::string, std::string> &SamHeaderRecord::getAllAttributes() const {
    return attributes;
}

SamHeaderRecord::SamHeaderRecord(const std::string& recordType) {
    try {
        headerRecordTypes = recordTypesAsVal.at(recordType);
    } catch (std::exception& e) {
        std::cerr << "Invalid header type " << recordType << std::endl;
        exit(-1);
    }
}

HeaderRecordType SamHeaderRecord::getRecordType() {
    return headerRecordTypes;
}

std::string SamHeaderRecord::getRecordTypeAsString() {
    return recordTypesAsString[headerRecordTypes];
}

bool SamHeaderRecord::containsAttribute(const std::string &key) const {
    return attributes.containsKey(key);
}
