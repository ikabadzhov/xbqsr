#ifndef GENOMICPIPELINE_SAMSTREAM_H
#define GENOMICPIPELINE_SAMSTREAM_H

#include "PlainInputStream.h"
#include "CompressedInputStream.h"

class SamStream {
public:
    static ssize_t readBytes(PlainInputStream& stream, std::string& buffer, ssize_t offset, ssize_t length);
    static bool isBAMFile(PlainInputStream& fileStream);
};


#endif //GENOMICPIPELINE_SAMSTREAM_H
