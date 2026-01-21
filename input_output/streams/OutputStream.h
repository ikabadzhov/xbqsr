#pragma once

#include <cstdint>
#include <span>

class OutputStream {
public:

    virtual size_t write(std::span<const uint8_t> data) = 0;
    virtual void close() = 0;
    virtual void flush() = 0;
    virtual ~OutputStream() = default; // Virtual destructor
};