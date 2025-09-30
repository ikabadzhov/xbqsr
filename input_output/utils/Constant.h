//
// Created by lorenzo on 3/11/24.
//

#ifndef GENOMICPIPELINE_CONSTANT_H
#define GENOMICPIPELINE_CONSTANT_H

#include <cstdio>
#include <cstdint>

//const ssize_t DEFAULT_BUFFER_SIZE = 2147483647; //2GB
const ssize_t DEFAULT_BUFFER_SIZE = (2UL << 30) - 1;//static_cast<size_t>(std::pow(2,30)); //1GB //10000000; //200MB
static const ssize_t BATCH_SIZE = (1UL << 20);

const std::string_view BAM_MAGIC = "BAM\u0001";

const char cigarOps[] = "MIDNSHP=X"; //CIGAR: op len<<4|op. ‘MIDNSHP=X’→‘012345678’
const char seqEncoder[] = "=ACMGRSVTWYHKDBN"; //4-bit encoded read: ‘=ACMGRSVTWYHKDBN’→ [0, 15].

const int32_t MAX_SEQ_LENGTH = 512;

const uint8_t knownSequence[] = {
        0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff,
        0x06, 0x00, 0x42, 0x43, 0x02, 0x00, 0x1b, 0x00, 0x03, 0x00,
        0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00
};

//const std::string uncompressedSEQ = "\000\000\000\000\000\377\006\000BC\002\000\033\000\003\000\000\000\000\000\000\000\000";

#endif //GENOMICPIPELINE_CONSTANT_H
