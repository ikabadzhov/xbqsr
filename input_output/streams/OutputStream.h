#ifndef GENOMICPIPELINE_OUTPUTSTREAM_H
#define GENOMICPIPELINE_OUTPUTSTREAM_H

#include <cstdint>
#include <span>
#include <string>

class OutputStream {
public:

    virtual size_t write(std::span<const uint8_t> data) = 0;
    virtual void close() = 0;
    virtual void flush() = 0;
    // Virtual destructor
    virtual ~OutputStream() = default;
};


#endif //GENOMICPIPELINE_OUTPUTSTREAM_H
