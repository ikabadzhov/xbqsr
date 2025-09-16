#include <cassert>
#include "CompressedInputStream.h"
#include <cstdint>
#include "../../zlib.h"
#include <stdexcept>
#include <string>

// ---- GZIP/BGZF header probes ------------------------------------------------

bool CompressedInputStream::isValidFile(PlainInputStream& stream) {
    // Read the first 18 bytes (enough for GZIP + BGZF extra header)
    std::string hdr(18, '\0');
    (void)stream.read(hdr, 0, 18);
    bool ok = isValidBlockHeader(hdr);
    stream.seekg(0); // reset to beginning
    return ok;
}

static inline unsigned U8(char c) { return static_cast<unsigned char>(c); }

bool CompressedInputStream::isValidBlockHeader(std::string& b) {
    if (b.size() < 18) return false;

    // GZIP magic + method
    if (U8(b[0]) != 0x1f || U8(b[1]) != 0x8b || U8(b[2]) != 0x08) return false;

    // FLG has FEXTRA bit set?
    const unsigned flg = U8(b[3]);
    if ((flg & 0x04) == 0) return false;

    // XLEN (LE) at bytes 10..11 should be 6 for BGZF
    const unsigned xlen = U8(b[10]) | (U8(b[11]) << 8);
    if (xlen < 6) return false;

    // BGZF subfield "BC" at 12..13, SLEN = 2 at 14..15
    if (U8(b[12]) != 'B' || U8(b[13]) != 'C') return false;
    const unsigned slen = U8(b[14]) | (U8(b[15]) << 8);
    if (slen != 2) return false;

    // If you want, you could read BSIZE at 16..17 here
    return true;
}

// ---- Core I/O ---------------------------------------------------------------

ssize_t CompressedInputStream::read(std::string& buffer, ssize_t offset, ssize_t length) {
    if (length < 0 || offset < 0 || static_cast<size_t>(offset) > buffer.size()) {
        throw std::out_of_range("Invalid offset/length");
    }
    if (static_cast<size_t>(offset + length) > buffer.size()) {
        throw std::out_of_range("Buffer too small");
    }

    ssize_t total = 0;
    while (total < length) {
        const unsigned want = static_cast<unsigned>(length - total);
        // gzread returns int: -1 on error, 0 on EOF, >0 bytes read
        int ret = gzread(fileDescriptor, &buffer[offset + total], want);
        if (ret < 0) {
            int err = Z_OK;
            const char* msg = gzerror(fileDescriptor, &err);
            throw std::runtime_error(std::string("gzread error: ") + (msg ? msg : "unknown"));
        }
        if (ret == 0) break; // EOF
        total += ret;
    }
    return total;
}

void CompressedInputStream::seekg(ssize_t pos) {
    // Use z_off_t for portability (zlib can be built with 64-bit z_off_t)
    z_off_t r = gzseek(fileDescriptor, static_cast<z_off_t>(pos), SEEK_SET);
    if (r < 0) {
        int err = Z_OK;
        const char* msg = gzerror(fileDescriptor, &err);
        throw std::runtime_error(std::string("gzseek failed: ") + (msg ? msg : "unknown"));
    }
}

const std::string& CompressedInputStream::getFilePath() { return path; }

// Return int so -1 is representable.
int CompressedInputStream::peek() {
    int c = gzgetc(fileDescriptor);
    if (c != -1) {
        // push back exactly what we read
        if (gzungetc(c, fileDescriptor) == -1) {
            int err = Z_OK;
            const char* msg = gzerror(fileDescriptor, &err);
            throw std::runtime_error(std::string("gzungetc failed: ") + (msg ? msg : "unknown"));
        }
    }
    return c; // -1 on EOF/error, else 0..255
}

bool CompressedInputStream::eof() {
    return gzeof(fileDescriptor) != 0;
}

bool CompressedInputStream::bad() {
    int err = Z_OK;
    (void)gzerror(fileDescriptor, &err);
    // Treat anything not OK or stream-end as bad.
    return err != Z_OK && err != Z_STREAM_END;
}

bool CompressedInputStream::good() {
    int err = Z_OK;
    (void)gzerror(fileDescriptor, &err);
    return (err == Z_OK) && !eof();
}

// Make close idempotent; destructor calls close().
void CompressedInputStream::close() {
    if (fileDescriptor) {
        gzclose(fileDescriptor);
        fileDescriptor = nullptr;
    }
}

CompressedInputStream::~CompressedInputStream() {
    // Do NOT call inflateEnd/deflateEnd when using gz* high-level API only.
    // Just close the gzFile once.
    try { close(); } catch (...) { /* no-throw destructor */ }
}