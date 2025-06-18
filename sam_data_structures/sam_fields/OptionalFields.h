#ifndef GENOMICPIPELINE_OPTIONALFIELDS_H
#define GENOMICPIPELINE_OPTIONALFIELDS_H

#include <string>
#include <unordered_map>
#include <vector>
#include "../../data_types/CompressedStringCollection.h"
#include "../../data_types/TagsCollection.h"

class OptionalFields {
public:
    std::map<uint16_t, TagKeyAndType> ALL_TAGS;
    /*
    * Specific TAG
    * */
    CompressedStringCollection RG;

    /*
    * Generic TAG
    * */
    TagsCollection<char> CharTags;
    TagsCollection<int64_t> IntegerTags;
    TagsCollection<float> FloatTags;

    // For the current implementation StringTags contains String, ByteArray and NumericArray type
    TagsCollection<std::string> StringTags;
    TagsCollectionArray<uint8_t> ArrayTagsHex;

    TagsCollectionArray<int8_t> ArrayTags_c;
    TagsCollectionArray<uint8_t> ArrayTags_C;
    TagsCollectionArray<int16_t> ArrayTags_s;
    TagsCollectionArray<uint16_t> ArrayTags_S;
    TagsCollectionArray<int32_t> ArrayTags_i;
    TagsCollectionArray<uint32_t> ArrayTags_I;
    TagsCollectionArray<float> ArrayTags_f;

    void Resize(size_t sz);
    void Permute(std::vector<size_t>& indices);
    bool ContainTag(const std::string& tag);

    SparseVector<char>& GetCharTags(const std::string& tag);
    SparseVector<int64_t>& GetIntTags(const std::string& tag);
    SparseVector<float>& GetFloatTags(const std::string& tag);
    SparseVector<std::string>& GetStringTags(const std::string& tag);

    SparseVector<std::vector<uint8_t>>& GetArrayHexTags(const std::string& tag);

    SparseVector<std::vector<int8_t>>& GetArrayInt8Tags(const std::string& tag);
    SparseVector<std::vector<uint8_t>>& GetArrayUint8Tags(const std::string& tag);
    SparseVector<std::vector<int16_t>>& GetArrayInt16Tags(const std::string& tag);
    SparseVector<std::vector<uint16_t>>& GetArrayUint16Tags(const std::string& tag);
    SparseVector<std::vector<int32_t>>& GetArrayInt32Tags(const std::string& tag);
    SparseVector<std::vector<uint32_t>>& GetArrayUint32Tags(const std::string& tag);
    SparseVector<std::vector<float>>& GetArrayFloatTags(const std::string& tag);

    void InsertCharTagAt(size_t idx, TagKeyAndType& tag, char tagValue);
    void InsertIntTagAt(size_t idx, TagKeyAndType& tag, int64_t tagValue);
    void InsertFloatTagAt(size_t idx, TagKeyAndType& tag, float tagValue);
    void InsertStringTagAt(size_t idx, TagKeyAndType& tag, const std::string& tagValue);

    void InsertArrayHexTagAt(size_t idx, TagKeyAndType tag, const std::string& array);

    void InsertArrayInt8TagAt(size_t idx, TagKeyAndType tag, std::vector<int8_t>& array);
    void InsertArrayUint8TagAt(size_t idx, TagKeyAndType tag, std::vector<uint8_t>& array);
    void InsertArrayInt16TagAt(size_t idx, TagKeyAndType tag, std::vector<int16_t>& array);
    void InsertArrayUint16TagAt(size_t idx, TagKeyAndType tag, std::vector<uint16_t>& array);
    void InsertArrayInt32TagAt(size_t idx, TagKeyAndType tag, std::vector<int32_t>& array);
    void InsertArrayUint32TagAt(size_t idx, TagKeyAndType tag, std::vector<uint32_t>& array);
    void InsertArrayFloatTagAt(size_t idx, TagKeyAndType tag, std::vector<float>& array);

    static void addTag(OptionalFields& TAGS, size_t idx, TagKeyAndType tagKeyAndType, const std::string& tagValue);
};


#endif //GENOMICPIPELINE_OPTIONALFIELDS_H
