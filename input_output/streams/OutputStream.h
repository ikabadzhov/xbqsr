#ifndef GENOMICPIPELINE_OUTPUTSTREAM_H
#define GENOMICPIPELINE_OUTPUTSTREAM_H

#include <string>

class OutputStream {
public:
    virtual const std::string& getFilePath() = 0;

    virtual void write(std::string &buffer, ssize_t n) = 0;
    virtual bool good() = 0;
    virtual bool bad() = 0;
    virtual bool eof() = 0;
    virtual void close() = 0;
    // Virtual destructor
    virtual ~OutputStream() = default;
};


#endif //GENOMICPIPELINE_OUTPUTSTREAM_H
