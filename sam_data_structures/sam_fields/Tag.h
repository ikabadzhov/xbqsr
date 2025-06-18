#ifndef GENOMICPIPELINE_TAG_H
#define GENOMICPIPELINE_TAG_H

#include <string>
#include <vector>
#include <cstdint>

class TagKeyAndType {
    uint16_t tagKey{};
    char tagType{};
    char secondaryTagType{0};

public:
    [[nodiscard]] uint16_t getTagKey() const;
    [[nodiscard]] std::string getTagKeyAsString() const;
    [[nodiscard]] char getTagType() const;
    [[nodiscard]] char getSecondaryTagType() const;

    void setTagKey(uint16_t tagKey);
    void setTagType(char tagType);
    void setSecondaryTagType(char secondaryTagType);

    static uint16_t encode(const std::string& tag); // Encode string as key
    static std::string decode(uint16_t tagKey); // Encode string as key
};


#endif //GENOMICPIPELINE_TAG_H
