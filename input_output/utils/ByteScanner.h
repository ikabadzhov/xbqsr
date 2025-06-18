#ifndef GENOMICPIPELINE_BYTESCANNER_H
#define GENOMICPIPELINE_BYTESCANNER_H

#include <string_view>
#include <cstdint>

class ByteScanner {
    size_t index;
    std::string_view bytes;

    void checkAvailable(size_t size) const;

public:
    // TODO: BAM files always have data in little-endian format
    //  By default we are assuming the platform is little endian.
    //  To make it compatible with big endian platform, create a
    //  version of these functions that rearrange the byte order.
    //  The check can be done at compile time so that the exact
    //  functions are called

    void reset(std::string_view s);
    uint32_t readUint32();
    int32_t readInt32();
    int16_t readInt16();
    uint16_t readUint16();
    uint8_t readByte();
    float readFloat();
    std::string readString(size_t len);
    std::string readUntil(char delim);

    const char* data();

    [[nodiscard]] bool canRead(size_t size) const;

    size_t getIndex();
    size_t incrementIndex(size_t n);
};
#endif //GENOMICPIPELINE_BYTESCANNER_H
