#ifndef GENOMICPIPELINE_RGSTRINGMAPPER_H
#define GENOMICPIPELINE_RGSTRINGMAPPER_H


#include <string>
#include <unordered_map>
#include <vector>

class RGStringMapper {
public:
    static const int EMPTY_STRING = -1;
    static const int STRING_NOT_FOUND = -2;

private:
    std::unordered_map<std::string, int> stringToId;
    std::vector<std::string> idToString;
    int32_t nextId = 0;

    // Private constructor and deletion of copy/move constructors and assignment operators
    RGStringMapper() = default;
    RGStringMapper(const RGStringMapper&) = delete;
    RGStringMapper& operator=(const RGStringMapper&) = delete;
    RGStringMapper(RGStringMapper&&) = delete;
    RGStringMapper& operator=(RGStringMapper&&) = delete;

public:
    static RGStringMapper& getInstance() {
        static RGStringMapper instance;
        return instance;
    }

    int32_t addString(const std::string& str);
    int32_t getStringId(const std::string& str) const;
    std::string getString(int32_t id) const;
    int32_t size() const;

    static bool isEmpty(int32_t i);
};



#endif //GENOMICPIPELINE_RGSTRINGMAPPER_H
