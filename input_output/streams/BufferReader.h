#pragma once
#include <string>
#include <memory>
#include "InputStream.h"

class BufferReader {
    std::string buffer;
    size_t pos;   // Current read position in the buffer
    size_t count_; // Number of valid characters in the buffer

public:
    explicit BufferReader(size_t size) : buffer(size, '\0'), pos(0), count_(0) {}
    void fill(std::unique_ptr<InputStream>& in);

    std::string_view getslice(char delim, bool& isDelimFound);
    std::string_view getslice(size_t numBytesToRead, size_t &numBytesRead);

    std::string_view getslicerange(size_t begin, size_t end);

    char peek();

    [[maybe_unused]] [[nodiscard]] bool has_data() const;
    [[maybe_unused]] [[nodiscard]] bool is_empty() const;

    size_t available_data();
    size_t count() {return count_;};
};
