#include "ByteScanner.h"

#include <stdexcept>
#include <cstring>

// Check if there aare enough data to read
void ByteScanner::checkAvailable(size_t size) const {
    if (!canRead(size)) {
        throw std::out_of_range("Not enough data to read");
    }
}

// Read a single byte and advance the index
uint8_t ByteScanner::readByte() {
    checkAvailable(sizeof(uint8_t));
    uint8_t value = static_cast<uint8_t>(bytes[index]);
    index++;
    return value;
}

// Read a 32-bit integer (4 bytes)
uint32_t ByteScanner::readUint32() {
    checkAvailable(sizeof(uint32_t));
    uint32_t value = static_cast<uint8_t>(bytes[index]) |
                     static_cast<uint8_t>(bytes[index + 1]) << 8 |
                     static_cast<uint8_t>(bytes[index + 2]) << 16 |
                     static_cast<uint8_t>(bytes[index + 3]) << 24;

    index += sizeof(uint32_t);
    return value;
}

// Read a 32-bit integer (4 bytes)
int32_t ByteScanner::readInt32() {
    checkAvailable(sizeof(uint32_t));
    int32_t value = static_cast<uint8_t>(bytes[index]) |
                     static_cast<uint8_t>(bytes[index + 1]) << 8 |
                     static_cast<uint8_t>(bytes[index + 2]) << 16 |
                     static_cast<uint8_t>(bytes[index + 3]) << 24;

    index += sizeof(int32_t);
    return value;
}

int16_t ByteScanner::readInt16() {
    checkAvailable(sizeof(int16_t));
    int16_t value = static_cast<uint8_t>(bytes[index]) |
                    static_cast<uint8_t>(bytes[index + 1]) << 8;

    index += sizeof(int16_t);
    return value;
}

uint16_t ByteScanner::readUint16() {
    checkAvailable(sizeof(uint16_t));
    uint16_t value = static_cast<uint8_t>(bytes[index]) |
                    static_cast<uint8_t>(bytes[index + 1]) << 8;

    index += sizeof(uint16_t);
    return value;
}

// Check if there are enough bytes left to read
bool ByteScanner::canRead(size_t size) const {
    return index + size <= bytes.size();
}

// Reset the scanner to the beginning
void ByteScanner::reset(std::string_view s) {
    bytes = s;
    index = 0;
}

std::string ByteScanner::readString(size_t len) {
    if (index + len > bytes.size()) {
        throw std::exception();
    }
    auto tmpIdx = index;
    index += len;
    return std::string{bytes.begin() + tmpIdx, len};
}

size_t ByteScanner::getIndex() {
    return index;
}

size_t ByteScanner::incrementIndex(size_t n) {
    return index += n;
}

size_t ByteScanner::setIndex(size_t pos) {
    return index = pos;
}

float ByteScanner::readFloat() {
    checkAvailable(sizeof(float));

    float value;
    std::memcpy(&value, bytes.data() + index, sizeof(float));

    index += sizeof(float);
    return value;
}

std::string ByteScanner::readUntil(char delim) {
    size_t i = index;
    size_t start = index;

    for ( ; i < this->bytes.size() && bytes[i] != delim; i++);
    size_t end = i;
    index = i+1;
    return {bytes.data() + start, bytes.data() + end};
}

const char *ByteScanner::data() {
    return bytes.data()+index;
}
