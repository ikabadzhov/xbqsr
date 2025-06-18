
#include <cstring>
#include "ByteFormatter.h"

void ByteFormatter::appendInt32(int32_t val) {
    //bytes.resize(sizeof(val), 0);
    bytes.resize(bytes.size() + sizeof(val));

    // Copy the int32_t value into the buffer
    std::memcpy(bytes.data() + bytes.size() - sizeof(val), &val, sizeof(val));
}

void ByteFormatter::appendString(const std::string &str) {
    bytes.append(str);
}

void ByteFormatter::appendByte(int8_t val) {
    //bytes.resize(sizeof(val), 0);
    bytes.resize(bytes.size() + sizeof(val));

    // Copy the byte value into the buffer
    std::memcpy(bytes.data() + bytes.size() - sizeof(val), &val, sizeof(val));
}

void ByteFormatter::appendInt16(int16_t val) {
    //bytes.resize(sizeof(val), 0);
    bytes.resize(bytes.size() + sizeof(val));

    // Copy the int16_t value into the buffer
    std::memcpy(bytes.data() + bytes.size() - sizeof(val), &val, sizeof(val));
}

void ByteFormatter::appendFloat(float val) {
    //bytes.resize(sizeof(val), 0);
    bytes.resize(bytes.size() + sizeof(val));

    // Copy the float value into the buffer
    std::memcpy(bytes.data() + bytes.size() - sizeof(val), &val, sizeof(val));
}
