#ifndef GENOMICPIPELINE_BYTEFORMATTER_H
#define GENOMICPIPELINE_BYTEFORMATTER_H

#include <string>

class ByteFormatter {
    std::string& bytes;

public:
    explicit ByteFormatter(std::string& bytes) : bytes(bytes) {};
    void appendInt32(int32_t val);

    void appendString(const std::string &str);

    void appendByte(int8_t val);

    void appendInt16(int16_t val);

    void appendFloat(float val);
};


#endif //GENOMICPIPELINE_BYTEFORMATTER_H
