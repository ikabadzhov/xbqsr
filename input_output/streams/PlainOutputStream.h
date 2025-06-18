#ifndef GENOMICPIPELINE_PLAINOUTPUTSTREAM_H
#define GENOMICPIPELINE_PLAINOUTPUTSTREAM_H


#include <iostream>
#include <fstream>
#include "OutputStream.h"

class PlainOutputStream : public OutputStream {
    std::ofstream out;
    std::string path;
public:
    explicit PlainOutputStream(const std::string& filename) : out(filename), path(filename) {
        if (!out.is_open()) {
            std::cerr << "Cannot open the file " << filename << std::endl;
            exit(-1);
        }
    };
    bool bad() override;
    const std::string & getFilePath() override;
    bool eof() override;
    bool good() override;
    void write(std::string &buffer, ssize_t n) override;
    void close() override;
};


#endif //GENOMICPIPELINE_PLAINOUTPUTSTREAM_H
