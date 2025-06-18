#include <iostream>
#include "PlainInputStream.h"

void PlainInputStream::seekg(ssize_t pos) {
    in.seekg(pos);
}

ssize_t PlainInputStream::read(std::string& buffer, ssize_t offset, ssize_t length) {
    if (length == 0) {
        return 0;
    }
    in.read(buffer.data() + offset, length);
    return in.gcount();
}

const std::string& PlainInputStream::getFilePath() {
    return path;
}

char PlainInputStream::peek() {
    return static_cast<char>(in.peek());
}

bool PlainInputStream::eof() {
    return in.eof();
}

bool PlainInputStream::good() {
    return in.good();
}

bool PlainInputStream::bad() {
    return in.bad();
}

void PlainInputStream::close() {
    this->in.close();
}

PlainInputStream::~PlainInputStream() = default;
