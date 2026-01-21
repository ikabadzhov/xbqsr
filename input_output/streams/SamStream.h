#pragma once
#include "PlainInputStream.h"

class SamStream {
public:
    static ssize_t readBytes(PlainInputStream& stream, std::string& buffer, ssize_t offset, ssize_t length);
    static bool isBAMFile(PlainInputStream& fileStream);
};