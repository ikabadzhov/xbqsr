#include <string>
#include "PlainOutputStream.h"

bool PlainOutputStream::bad() {
    return false;
}

const std::string &PlainOutputStream::getFilePath() {
    return path;
}

bool PlainOutputStream::eof() {
    return false;
}

bool PlainOutputStream::good() {
    return false;
}

void PlainOutputStream::close(){
    this->out.close();
}

void PlainOutputStream::write(std::string &buffer, ssize_t n) {
    size_t sz = n>0?n:buffer.size();
    out.write(buffer.data(), sz);
    //out.flush();
    if (!out) {
        throw std::runtime_error("Error writing data to output stream");
    }
}


