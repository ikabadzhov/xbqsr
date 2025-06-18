#ifndef GENOMICPIPELINE_PARSEDHEADERLINE_H
#define GENOMICPIPELINE_PARSEDHEADERLINE_H

#include <string_view>
#include "../../data_types/LinkedHashMap.h"
#include "../../sam_data_structures/sam_header_fields/SamHeaderRecord.h"


class ParsedHeaderLine {
public:
    LinkedHashMap<std::string, std::string> keyValuePairs;
private:
    HeaderRecordType recordType;
    bool lineValid{};

    static std::vector<std::string> split(std::string_view line, char delim);
    static std::pair<std::string, std::string> splitKeyValue(std::string_view str, char delim);

public:
    explicit ParsedHeaderLine(std::string_view headerLine);

    HeaderRecordType getHeaderRecordType();
    std::string getValue(const std::string& key);
    bool containsKey(const std::string& key);
    std::string removeValue(const std::string& key);

    [[maybe_unused]] bool isLineValid() const;
    void validateSortOrderValue(std::vector<std::string>& value);

    bool requireTag(const std::string &tag);
};


#endif //GENOMICPIPELINE_PARSEDHEADERLINE_H
