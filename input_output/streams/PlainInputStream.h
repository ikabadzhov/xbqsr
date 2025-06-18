#ifndef GENOMICPIPELINE_PLAININPUTSTREAM_H
#define GENOMICPIPELINE_PLAININPUTSTREAM_H

#include <cstdio>
#include <string>
#include <iostream>
#include <fstream>
#include "InputStream.h"

class PlainInputStream : public InputStream {
    std::ifstream in;
    std::string path;
public:
    explicit PlainInputStream(const std::string& filename) : in(filename), path(filename) {
        if (!in.is_open()) {
            std::cerr << "Cannot open the file " << filename << std::endl;
            exit(-1);
        }
    };
    const std::string& getFilePath() override;
    ssize_t read(std::string& buffer, ssize_t offset, ssize_t length) override;
    void seekg(ssize_t pos) override;
    char peek() override;
    bool good() override;
    bool bad() override;
    bool eof() override;
    void close() override;
    ~PlainInputStream() override;
};


#endif //GENOMICPIPELINE_PLAININPUTSTREAM_H
