#ifndef BQSR_BAMUTILSIO_H
#define BQSR_BAMUTILSIO_H
#include <iostream>
#include <cinttypes>
#include <memory>
#include "../streams/InputStream.h"

// --- Little-endian decoders (portable) ---
static inline uint32_t load_le_u32(const char* p) {
    return  (uint32_t)(uint8_t)p[0]
                               | ((uint32_t)(uint8_t)p[1] << 8)
    | ((uint32_t)(uint8_t)p[2] << 16)
    | ((uint32_t)(uint8_t)p[3] << 24);
}
static inline int32_t load_le_i32(const char* p) {
    return (int32_t)load_le_u32(p);
}

// --- Read exactly N bytes or throw (uses your InputStream::read API) ---
static inline void read_exact(std::unique_ptr<InputStream>& in, std::string& dst, size_t n) {
    dst.resize(n);
    size_t total = 0;
    while (total < n) {
        const ssize_t got = in->read(dst, (ssize_t)total, (ssize_t)(n - total));
        if (got <= 0) throw std::runtime_error("Unexpected EOF while reading BAM header");
        total += (size_t)got;
    }
}

static inline std::string read_exact_new(std::unique_ptr<InputStream>& in, size_t n) {
    std::string buf;
    read_exact(in, buf, n);
    return buf;
}

#endif //BQSR_BAMUTILSIO_H
