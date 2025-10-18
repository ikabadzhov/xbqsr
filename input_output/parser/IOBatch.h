#ifndef GENOMICPIPELINE_IOBATCH_H
#define GENOMICPIPELINE_IOBATCH_H


#include <string>
#include <utility>
#include <vector>

#include "../streams/BufferReader.h"

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

// Give this the same spirit as InputBatchBam, but for output.
class OutputBatchBam {
public:
    explicit OutputBatchBam(size_t cap) : buf_(cap), capacity_(cap) {}

    uint8_t* data()                  { return buf_.data(); }
    size_t   capacity()        const { return capacity_; }
    size_t   used()            const { return used_; }
    bool     empty()           const { return used_ == 0; }

    // Set the number of valid bytes and number of records for this batch
    void seal(size_t bytes, size_t nrec) { used_ = bytes; n_records_ = nrec; }

    // Reset to empty (caller will rewrite from offset 0)
    void reset() { used_ = 0; n_records_ = 0; }

private:
    std::vector<uint8_t> buf_;
    size_t capacity_{0};
    size_t used_{0};
    size_t n_records_{0};
};

#endif //GENOMICPIPELINE_IOBATCH_H
