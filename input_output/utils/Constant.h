#ifndef GENOMICPIPELINE_CONSTANT_H
#define GENOMICPIPELINE_CONSTANT_H

#include <cstdio>
#include <cstdint>


const ssize_t DEFAULT_BUFFER_SIZE = 512 * (1UL << 20) - 1;
static const ssize_t DEFAULT_NUM_ALIGNMENTS = 10*(1UL << 20);

// Practical defaults for short-read workloads
// Note: they are used as defualt values in reserve method, not actually in resize.
const int32_t MAX_SEQ_LENGTH = 512;

constexpr int32_t MAX_QNAME_LENGTH   = 256;      // l_read_name is uint8_t
constexpr int32_t MAX_CIGAR_OPS      = 2 * MAX_SEQ_LENGTH;
constexpr int32_t MAX_QUAL_LENGTH      = MAX_SEQ_LENGTH;
constexpr int32_t MAX_TAGS_LENGTH      = 3 * MAX_SEQ_LENGTH;
constexpr int32_t MAX_SIMPLE_TYPES_SIZE = 19;

constexpr int32_t DEFAULT_CAPACITY_MARGIN = 1024;


const std::string_view BAM_MAGIC = "BAM\u0001";

const char cigarOps[] = "MIDNSHP=X"; //CIGAR: op len<<4|op. ‘MIDNSHP=X’→‘012345678’
const char seqEncoder[] = "=ACMGRSVTWYHKDBN"; //4-bit encoded read: ‘=ACMGRSVTWYHKDBN’→ [0, 15].

const uint8_t knownSequence[] = {
        0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff,
        0x06, 0x00, 0x42, 0x43, 0x02, 0x00, 0x1b, 0x00, 0x03, 0x00,
        0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00
};

//const std::string uncompressedSEQ = "\000\000\000\000\000\377\006\000BC\002\000\033\000\003\000\000\000\000\000\000\000\000";

#endif //GENOMICPIPELINE_CONSTANT_H
