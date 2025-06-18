#ifndef GENOMICPIPELINE_BUFFERWRITER_H
#define GENOMICPIPELINE_BUFFERWRITER_H

#include <string>
#include <memory>
#include <utility>
#include "OutputStream.h"
#include <cstring>


class BufferWriter {
    std::string buffer;
    size_t pos;   // Current read position in the buffer
    size_t count; // Number of valid characters in the buffer
    size_t size;

public:
    explicit BufferWriter(size_t size) : buffer(size, '\0'), pos(0), count(0) {
        this->size = size;
    }

    bool addSlice(std::string content);
    bool addNewLine();
    void flushContent();
    void getContent(size_t startPos, size_t length, std::string &value);
    size_t getLastCount();
    size_t getLastPos();
    size_t getSize();

    std::string_view getslice(char delim, bool& isDelimFound);
    std::string_view getslice(size_t numBytesToRead, size_t &numBytesRead);


    char peek();

    [[maybe_unused]] [[nodiscard]] bool has_data() const;
    [[maybe_unused]] [[nodiscard]] bool is_empty() const;

    size_t available_data();
    size_t remaining_space();
};


#endif //GENOMICPIPELINE_BUFFERREADER_H
