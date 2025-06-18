#include <stdexcept>
#include "Tag.h"

uint16_t TagKeyAndType::getTagKey() const {
    return tagKey;
}

std::string TagKeyAndType::getTagKeyAsString() const {
    return decode(tagKey);
}

char TagKeyAndType::getTagType() const {
    return tagType;
}

char TagKeyAndType::getSecondaryTagType() const {
    return secondaryTagType;
}

uint16_t TagKeyAndType::encode(const std::string& tag) {
    return static_cast<uint16_t>(tag[1] << 8 | tag[0]);
}
std::string TagKeyAndType::decode(uint16_t tagKey) {
    return {char(tagKey & 255), char(tagKey >> 8 & 255)};
}

void TagKeyAndType::setTagKey(uint16_t key) {
    TagKeyAndType::tagKey = key;
}

void TagKeyAndType::setTagType(char type) {
    TagKeyAndType::tagType = type;
}

void TagKeyAndType::setSecondaryTagType(char type) {
    TagKeyAndType::secondaryTagType = type;
}
