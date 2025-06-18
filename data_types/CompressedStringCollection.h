#ifndef GENOMICPIPELINE_COMPRESSEDSTRINGCOLLECTION_H
#define GENOMICPIPELINE_COMPRESSEDSTRINGCOLLECTION_H

#include <string>
#include <unordered_map>
#include <vector>

/**
 * Note:
 * RG must start from 0, as it is used as indexes in a multidimensinal array
 * On the other hand, LIBID must handle empty strings as undefinedLibrary.
 * */

const int32_t EMPTY_STRING = -1;
const int32_t DEFAULT_STRING = 0;
/**
 * Map strings into integer in the range [0/1,MAX_INT]
 * The value 0 may be reserved for the empty string, in
 * case an object can assume the value of an empty string.
 * Otherwise, the mapping start from 0, where 0 is the first
 * valid value.
 **/
class CompressedStringCollection {
    std::unordered_map<std::string, int32_t> strToInt;
    std::vector<std::string> intToStr;
    std::vector<int32_t> values;

public:
    explicit CompressedStringCollection(std::string_view defaultString);
    CompressedStringCollection() = default;

    void Set(size_t idx, const std::string& rg);
    int32_t Get(size_t idx) const;
    void Resize(size_t sz);
    size_t NumEntriesTable();
    std::string toString(size_t idx) const;
    int32_t StringToInt(std::string& str);
    std::string& IntToString(int32_t i);
    std::vector<int32_t>& toVec();
    static bool isEmpty(int32_t val);
    void Permute(std::vector<size_t>& indices);
};


#endif //GENOMICPIPELINE_COMPRESSEDSTRINGCOLLECTION_H
