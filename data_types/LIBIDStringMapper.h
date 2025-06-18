#ifndef GENOMICPIPELINE_LIBIDSTRINGMAPPER_H
#define GENOMICPIPELINE_LIBIDSTRINGMAPPER_H

#include <vector>
#include <string>
#include <algorithm>
#include <unordered_map>


class LIBIDStringMapper {
public:
    static const int EMPTY_STRING = -1;
    static const int DEFAULT_STRING = 0;
    static const int STRING_NOT_FOUND = -2;
    static const int INVALID_STRING = -3;

private:
    std::unordered_map<std::string, int32_t> stringToId = {{"undefinedLibrary", 0}};
    std::vector<std::string> idToString = {"undefinedLibrary"};
    int32_t nextId = 1;

    // Private constructor and deletion of copy/move constructors and assignment operators
    LIBIDStringMapper() = default;
    LIBIDStringMapper(const LIBIDStringMapper&) = delete;
    LIBIDStringMapper& operator=(const LIBIDStringMapper&) = delete;
    LIBIDStringMapper(LIBIDStringMapper&&) = delete;
    LIBIDStringMapper& operator=(LIBIDStringMapper&&) = delete;

public:
    static LIBIDStringMapper& getInstance() {
        static LIBIDStringMapper instance;
        return instance;
    }

    int32_t addString(const std::string& str);
    int32_t getStringId(const std::string& str) const;
    std::string getString(int32_t id) const;
    int32_t size() const;
};


#endif //GENOMICPIPELINE_LIBIDSTRINGMAPPER_H
