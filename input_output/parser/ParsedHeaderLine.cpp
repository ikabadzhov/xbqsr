#include <iostream>
#include "ParsedHeaderLine.h"

ParsedHeaderLine::ParsedHeaderLine(std::string_view headerLine) {
    if (!headerLine.starts_with("@")) {
        throw std::invalid_argument("Line does not start with '@'");
    }

    auto fields = split(headerLine, '\t');
    try {
        SamHeaderRecord headerRecordType(fields[0].substr(1));
    } catch (const std::invalid_argument& e) {
        std::cerr << "Unrecognized header record type" << std::endl;
        return;
    }

    SamHeaderRecord headerRecord(fields[0].substr(1));
    recordType = headerRecord.getRecordType();
    if (recordType == HeaderRecordType::CO) {
        lineValid = true;
        return;
    }

    for (size_t i = 1; i < fields.size(); ++i) {
        try {
            auto keyAndValue = splitKeyValue(fields[i], ':');
            //std::string key = std::string(keyAndValue.first);
            //std::string value = std::string(keyAndValue.second);

            if (keyValuePairs.containsKey(keyAndValue.first) && keyValuePairs.get(keyAndValue.first) != keyAndValue.second) {
                std::cerr << "Key:value pair clashes with existing entry" << std::endl;
                continue;
            }

            if (keyAndValue.first.length() != 2) {
                std::cerr << "Key is not two characters" << std::endl;
                continue;
            }

            keyValuePairs.put(keyAndValue.first, keyAndValue.second);
        } catch (const std::invalid_argument& e) {
            std::cerr << e.what() << std::endl;
        }
    }
    lineValid = true;
}

std::vector<std::string> ParsedHeaderLine::split(std::string_view line, char delim) {
    std::vector<std::string> result;
    size_t start = 0;
    size_t end = line.find(delim);
    while (end != std::string_view::npos) {
        result.emplace_back(line.substr(start, end - start));
        start = end + 1;
        end = line.find(delim, start);
    }
    result.emplace_back(line.substr(start));
    return result;
}

std::pair<std::string, std::string> ParsedHeaderLine::splitKeyValue(std::string_view str, char delim) {
    size_t sep = str.find(delim);
    if (sep == std::string_view::npos) {
        throw std::invalid_argument("Invalid key:value pair");
    }
    return {std::string{str.substr(0, sep)}, std::string{str.substr(sep + 1)}};
}


bool ParsedHeaderLine::requireTag(const std::string& tag) {
    if (!keyValuePairs.containsKey(tag)) {
        std::cerr << "@" << recordType <<  " line missing " + tag + " tagValue" << std::endl;
        return false;
    } else {
        return true;
    }
}

[[maybe_unused]] bool ParsedHeaderLine::isLineValid() const {
    return lineValid;
}

void ParsedHeaderLine::validateSortOrderValue(std::vector<std::string> &value) {

}

bool ParsedHeaderLine::containsKey(const std::string &key) {
    return keyValuePairs.containsKey(key);
}

std::string ParsedHeaderLine::getValue(const std::string &key) {
    return keyValuePairs.get(key);
}

HeaderRecordType ParsedHeaderLine::getHeaderRecordType() {
    return recordType;
}

std::string ParsedHeaderLine::removeValue(const std::string &key) {
    std::string ret = keyValuePairs.get(key);
    keyValuePairs.soft_remove(key);
    return ret;
}

