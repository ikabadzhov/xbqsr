#ifndef GENOMICPIPELINE_FILEMETADATA_H
#define GENOMICPIPELINE_FILEMETADATA_H

#include "SamHeaderRecord.h"

class FileMetadata : public SamHeaderRecord {
    // LinkedHashMap<std::string, std::string> attributes;
public:
    FileMetadata() : SamHeaderRecord("HD") {};
    void put(const std::string& key, const std::string& value);
    std::string get(const std::string &key);
};

#endif //GENOMICPIPELINE_FILEMETADATA_H
