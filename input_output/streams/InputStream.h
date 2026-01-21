#pragma once
#include <string>

class InputStream {
public:
    virtual const std::string& getFilePath() = 0;
    virtual ssize_t read(std::string& buffer, ssize_t offset, ssize_t length) = 0;
    virtual void seekg(ssize_t pos) = 0;
    virtual int peek() = 0;
    virtual bool good() = 0;
    virtual bool bad() = 0;
    virtual bool eof() = 0;
    virtual void close() = 0;
    virtual ~InputStream() = default;  // Virtual destructor
};