#include <iostream>
#include "OptionalFields.h"
#include "Tag.h"
#include "../../input_output/utils/StringScanner.h"

void OptionalFields::Resize(size_t sz) {
    RG.resize(sz);
    tagsBlobOff.resize(sz, 0);
    tagsBlobSize.resize(sz, 0);
    /*IntegerTags.Resize(sz);
    FloatTags.Resize(sz);
    StringTags.Resize(sz);
    CharTags.Resize(sz);

    ArrayTagsHex.Resize(sz);

    ArrayTags_c.Resize(sz);
    ArrayTags_C.Resize(sz);
    ArrayTags_s.Resize(sz);
    ArrayTags_S.Resize(sz);
    ArrayTags_i.Resize(sz);
    ArrayTags_I.Resize(sz);
    ArrayTags_f.Resize(sz);*/
}

void OptionalFields::Reserve(size_t sz) {
    RG.reserve(sz);
    tagsBlobOff.reserve(sz);
    tagsBlobSize.reserve(sz);
    /*IntegerTags.Resize(sz);
    FloatTags.Resize(sz);
    StringTags.Resize(sz);
    CharTags.Resize(sz);

    ArrayTagsHex.Resize(sz);

    ArrayTags_c.Resize(sz);
    ArrayTags_C.Resize(sz);
    ArrayTags_s.Resize(sz);
    ArrayTags_S.Resize(sz);
    ArrayTags_i.Resize(sz);
    ArrayTags_I.Resize(sz);
    ArrayTags_f.Resize(sz);*/
}

void OptionalFields::ReserveData(size_t sz) {
    tagArena.reserve(sz);
}

size_t OptionalFields::Capacity() {
    return tagsBlobSize.capacity();
}

size_t OptionalFields::CapacityData() {
    return tagArena.capacity();
}

void OptionalFields::Permute(std::vector<size_t> &indices) {
    permute(indices, this->RG);
    this->IntegerTags.Permute(indices);
    this->FloatTags.Permute(indices);
    this->StringTags.Permute(indices);
    this->CharTags.Permute(indices);
}

SparseVector<char> &OptionalFields::GetCharTags(const std::string& tag) {
    return CharTags.GetTags(tag);
}

SparseVector<int64_t> &OptionalFields::GetIntTags(const std::string& tag) {
    return IntegerTags.GetTags(tag);
}

SparseVector<float> &OptionalFields::GetFloatTags(const std::string& tag) {
    return FloatTags.GetTags(tag);
}

SparseVector<std::string>& OptionalFields::GetStringTags(const std::string& tag) {
    return StringTags.GetTags(tag);
}

SparseVector<std::vector<uint8_t>> &OptionalFields::GetArrayUint8Tags(const std::string& tag) {
    return ArrayTags_C.GetTags(tag);
}

SparseVector<std::vector<int8_t>> &OptionalFields::GetArrayInt8Tags(const std::string& tag) {
    return ArrayTags_c.GetTags(tag);
}

SparseVector<std::vector<float>> &OptionalFields::GetArrayFloatTags(const std::string& tag) {
    return ArrayTags_f.GetTags(tag);
}

SparseVector<std::vector<uint32_t>> &OptionalFields::GetArrayUint32Tags(const std::string& tag) {
    return ArrayTags_I.GetTags(tag);
}

SparseVector<std::vector<int32_t>> &OptionalFields::GetArrayInt32Tags(const std::string& tag) {
    return ArrayTags_i.GetTags(tag);
}

SparseVector<std::vector<int16_t>> &OptionalFields::GetArrayInt16Tags(const std::string& tag) {
    return ArrayTags_s.GetTags(tag);
}

SparseVector<std::vector<uint16_t>> &OptionalFields::GetArrayUint16Tags(const std::string& tag) {
    return ArrayTags_S.GetTags(tag);
}

SparseVector<std::vector<uint8_t>> &OptionalFields::GetArrayHexTags(const std::string& tag) {
    return ArrayTagsHex.GetTags(tag);
}

bool OptionalFields::ContainTag(const std::string &tag) {
    uint16_t tagKey = TagKeyAndType::encode(tag);
    return ALL_TAGS.contains(tagKey);
}

void OptionalFields::InsertCharTagAt(size_t idx, TagKeyAndType& tag, char tagValue) {
    if (ALL_TAGS.contains(tag.getTagKey()) <= 0) {
        ALL_TAGS[tag.getTagKey()] = tag;
    }
    CharTags.AddTagIfNeeded(tag.getTagKey());
    CharTags.AssignTagValue(idx, tag.getTagKey(), tagValue);
}

void OptionalFields::InsertIntTagAt(size_t idx, TagKeyAndType &tag, int64_t tagValue) {
    if (ALL_TAGS.contains(tag.getTagKey()) <= 0) {
        ALL_TAGS[tag.getTagKey()] = tag;
    }
    IntegerTags.AddTagIfNeeded(tag.getTagKey());
    IntegerTags.AssignTagValue(idx, tag.getTagKey(), tagValue);
}

void OptionalFields::InsertFloatTagAt(size_t idx, TagKeyAndType &tag, float tagValue) {
    if (ALL_TAGS.contains(tag.getTagKey()) <= 0) {
        ALL_TAGS[tag.getTagKey()] = tag;
    }
    FloatTags.AddTagIfNeeded(tag.getTagKey());
    FloatTags.AssignTagValue(idx, tag.getTagKey(), tagValue);
}

void OptionalFields::InsertStringTagAt(size_t idx, TagKeyAndType &tag, const std::string &tagValue) {
    /*if (ALL_TAGS.contains(tag.getTagKey()) <= 0) {
        ALL_TAGS[tag.getTagKey()] = tag;
    }
    if (tag.getTagKeyAsString() == "RG") {
        RG.Set(idx, tagValue);
    } else {
        StringTags.AddTagIfNeeded(tag.getTagKey());
        StringTags.AssignTagValue(idx, tag.getTagKey(), tagValue);
    }*/
}

void OptionalFields::InsertArrayHexTagAt(size_t idx, TagKeyAndType tag, const std::string &array) {
    if (ALL_TAGS.contains(tag.getTagKey()) <= 0) {
        ALL_TAGS[tag.getTagKey()] = tag;
    }
    StringTags.AddTagIfNeeded(tag.getTagKey());
    StringTags.AssignTagValue(idx, tag.getTagKey(), array);
}

void OptionalFields::InsertArrayInt8TagAt(size_t idx, TagKeyAndType tag, std::vector<int8_t> &array) {
    if (ALL_TAGS.contains(tag.getTagKey()) <= 0) {
        ALL_TAGS[tag.getTagKey()] = tag;
    }
    ArrayTags_c.AddTagIfNeeded(tag.getTagKey());
    ArrayTags_c.AssignTagValue(idx, tag.getTagKey(), array);
}

void OptionalFields::InsertArrayUint8TagAt(size_t idx, TagKeyAndType tag, std::vector<uint8_t> &array) {
    if (ALL_TAGS.contains(tag.getTagKey()) <= 0) {
        ALL_TAGS[tag.getTagKey()] = tag;
    }
    ArrayTags_C.AddTagIfNeeded(tag.getTagKey());
    ArrayTags_C.AssignTagValue(idx, tag.getTagKey(), array);
}

void OptionalFields::InsertArrayInt16TagAt(size_t idx, TagKeyAndType tag, std::vector<int16_t> &array) {
    if (ALL_TAGS.contains(tag.getTagKey()) <= 0) {
        ALL_TAGS[tag.getTagKey()] = tag;
    }
    ArrayTags_s.AddTagIfNeeded(tag.getTagKey());
    ArrayTags_s.AssignTagValue(idx, tag.getTagKey(), array);
}

void OptionalFields::InsertArrayUint16TagAt(size_t idx, TagKeyAndType tag, std::vector<uint16_t> &array) {
    if (ALL_TAGS.contains(tag.getTagKey()) <= 0) {
        ALL_TAGS[tag.getTagKey()] = tag;
    }
    ArrayTags_S.AddTagIfNeeded(tag.getTagKey());
    ArrayTags_S.AssignTagValue(idx, tag.getTagKey(), array);
}

void OptionalFields::InsertArrayInt32TagAt(size_t idx, TagKeyAndType tag, std::vector<int32_t> &array) {
    if (ALL_TAGS.contains(tag.getTagKey()) <= 0) {
        ALL_TAGS[tag.getTagKey()] = tag;
    }
    ArrayTags_i.AddTagIfNeeded(tag.getTagKey());
    ArrayTags_i.AssignTagValue(idx, tag.getTagKey(), array);
}

void OptionalFields::InsertArrayUint32TagAt(size_t idx, TagKeyAndType tag, std::vector<uint32_t> &array) {
    if (ALL_TAGS.contains(tag.getTagKey()) <= 0) {
        ALL_TAGS[tag.getTagKey()] = tag;
    }
    ArrayTags_I.AddTagIfNeeded(tag.getTagKey());
    ArrayTags_I.AssignTagValue(idx, tag.getTagKey(), array);
}

void OptionalFields::InsertArrayFloatTagAt(size_t idx, TagKeyAndType tag, std::vector<float> &array) {
    if (ALL_TAGS.contains(tag.getTagKey()) <= 0) {
        ALL_TAGS[tag.getTagKey()] = tag;
    }
    ArrayTags_f.AddTagIfNeeded(tag.getTagKey());
    ArrayTags_f.AssignTagValue(idx, tag.getTagKey(), array);
}


void OptionalFields::addTag(OptionalFields& TAGS, size_t idx, TagKeyAndType tagKeyAndType, const std::string &tagValue) {
    //int64_t tagValueInt;
    //float tagValueFloat;
    std::vector<int32_t> tagValueArrayInt32;
    std::vector<uint32_t> tagValueArrayUInt32;
    std::vector<int16_t> tagValueArrayInt16;
    std::vector<uint16_t> tagValueArrayUInt16;
    std::vector<int8_t> tagValueArrayInt8;
    std::vector<uint8_t> tagValueArrayUInt8;
    std::vector<float> tagValueArrayFloat;

    if (tagKeyAndType.getTagType() == 'A') {
        TAGS.InsertCharTagAt(idx, tagKeyAndType, static_cast<char>(tagValue[0]));

    } else if (tagKeyAndType.getTagType() == 'c' || tagKeyAndType.getTagType() == 'C' ||
               tagKeyAndType.getTagType() == 's' ||
               tagKeyAndType.getTagType() == 'S' || tagKeyAndType.getTagType() == 'i' ||
               tagKeyAndType.getTagType() == 'I') {

        TAGS.InsertIntTagAt(idx, tagKeyAndType, static_cast<int64_t>(stoll(tagValue)));

    } else if (tagKeyAndType.getTagType() == 'f') {

        TAGS.InsertFloatTagAt(idx, tagKeyAndType, static_cast<float>(stof(tagValue)));

    } else if (tagKeyAndType.getTagType() == 'Z') {

        TAGS.InsertStringTagAt(idx, tagKeyAndType, tagValue);

    } else if (tagKeyAndType.getTagType() == 'H') {

        TAGS.InsertArrayHexTagAt(idx, tagKeyAndType, tagValue);

    } else if (tagKeyAndType.getTagType() == 'B') {
        std::string arrayValue;
        StringScanner tagValueScanner;

        tagValueScanner.reset(tagValue);

        if (tagKeyAndType.getSecondaryTagType() == 'c') {
            tagValueArrayInt8.clear();
            while (tagValueScanner.len()) {
                tagValueScanner.readUntil(',', arrayValue);
                tagValueArrayInt8.emplace_back(static_cast<int8_t>(stol(arrayValue)));
            }
            TAGS.InsertArrayInt8TagAt(idx, tagKeyAndType, tagValueArrayInt8);
        } else if (tagKeyAndType.getSecondaryTagType() == 'C') {
            tagValueArrayUInt8.clear();
            while (tagValueScanner.len()) {
                tagValueScanner.readUntil(',', arrayValue);
                tagValueArrayUInt8.emplace_back(static_cast<uint8_t>(stol(arrayValue)));
            }
            TAGS.InsertArrayUint8TagAt(idx, tagKeyAndType, tagValueArrayUInt8);
        } else if (tagKeyAndType.getSecondaryTagType() == 's') {
            tagValueArrayInt16.clear();
            while (tagValueScanner.len()) {
                tagValueScanner.readUntil(',', arrayValue);
                tagValueArrayInt16.emplace_back(static_cast<int16_t>(stol(arrayValue)));
            }
            TAGS.InsertArrayInt16TagAt(idx, tagKeyAndType, tagValueArrayInt16);
        } else if (tagKeyAndType.getSecondaryTagType() == 'S') {
            tagValueArrayUInt16.clear();
            while (tagValueScanner.len()) {
                tagValueScanner.readUntil(',', arrayValue);
                tagValueArrayUInt16.emplace_back(static_cast<uint16_t>(stol(arrayValue)));
            }
            TAGS.InsertArrayUint16TagAt(idx, tagKeyAndType, tagValueArrayUInt16);
        } else if (tagKeyAndType.getSecondaryTagType() == 'i') {
            tagValueArrayInt32.clear();
            while (tagValueScanner.len()) {
                tagValueScanner.readUntil(',', arrayValue);
                tagValueArrayInt32.emplace_back(static_cast<int32_t>(stol(arrayValue)));
            }
            TAGS.InsertArrayInt32TagAt(idx, tagKeyAndType, tagValueArrayInt32);
        } else if (tagKeyAndType.getSecondaryTagType() == 'I') {
            tagValueArrayUInt32.clear();
            while (tagValueScanner.len()) {
                tagValueScanner.readUntil(',', arrayValue);
                tagValueArrayUInt32.emplace_back(static_cast<uint32_t>(stol(arrayValue)));
            }
            TAGS.InsertArrayUint32TagAt(idx, tagKeyAndType, tagValueArrayUInt32);
        } else if (tagKeyAndType.getSecondaryTagType() == 'f') {
            tagValueArrayFloat.clear();
            while (tagValueScanner.len()) {
                tagValueScanner.readUntil(',', arrayValue);
                tagValueArrayFloat.emplace_back(static_cast<float>(stof(arrayValue)));
            }
            TAGS.InsertArrayFloatTagAt(idx, tagKeyAndType, tagValueArrayFloat);
        } else {
            std::cerr << "Panic: type not recognized" << std::endl;
            exit(-1);
        }
    } else {
        std::cerr << "Panic: type not recognized" << std::endl;
        exit(-1);
    }
}

