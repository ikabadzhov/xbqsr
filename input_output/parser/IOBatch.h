#ifndef GENOMICPIPELINE_IOBATCH_H
#define GENOMICPIPELINE_IOBATCH_H


#include <string>
#include <utility>
#include "../streams/BufferReader.h"
#include "../streams/BufferWriter.h"

class IOBatch {

};

class InputBatchBam {
public:
    BufferReader buffer;
    std::string remaining;
    size_t n_processed = 0;
    bool header_to_process = true;

    explicit InputBatchBam(size_t size) : buffer(size) {};

};

class InputBatchSam {
public:
    BufferReader buffer;
    std::string remaining;
    size_t n_processed = 0;
    bool header_to_process = true;

    explicit InputBatchSam(size_t size) : buffer(size) {};

};

class OutputBatchSam {
public:
    std::string buffer;
    std::string remaining;
    size_t n_processed = 0;
    size_t n_loaded = 0;
    size_t n_written_in_buffer = 0;
    size_t from = 0;
    size_t to = 0;

    explicit OutputBatchSam(size_t size) : buffer(size, 0) {};

};

#endif //GENOMICPIPELINE_IOBATCH_H
