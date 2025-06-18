#include <cassert>
#include "CompressedInputStream.h"
#include <cstdint>
#include "../../zlib.h"

/*
 * It checks if the file has a valid GZIP header.
 * If the file is a GZIP then it can be either BAM or plain
 * text sam compressed in gzip format.
 * */
bool CompressedInputStream::isValidFile(PlainInputStream& stream) {
    std::string buffer(18, 0);
    stream.read(buffer, 0, 18);
    bool isValid = isValidBlockHeader(buffer);
    stream.seekg(std::ios::beg);
    return isValid;
}

bool CompressedInputStream::isValidBlockHeader(std::string& buffer) {
    return buffer[0] == 31 && (buffer[1] & 255) == 139 && (buffer[3] & 4) != 0 && buffer[10] == 6 && buffer[12] == 66 && buffer[13] == 67;
}


ssize_t CompressedInputStream::read(std::string& buffer, ssize_t offset, ssize_t length) {
    if (length < 0 || offset < 0 || offset > static_cast<ssize_t>(buffer.size())) {
        throw std::out_of_range("Invalid offset and length");
    }

    if (static_cast<size_t>(offset + length) > buffer.size()) {
        throw std::out_of_range("Invalid buffer size");
    }

    ssize_t totalBytesRead = 0;
    while (totalBytesRead < length) {
        //gzread support up to 2^32-1 before looping while reading [returns int instead of size_t which is ulong]
        uint32_t bytesRead = gzread(fileDescriptor, &buffer[offset + totalBytesRead], length - totalBytesRead);
        if (bytesRead <= 0) {
            // Handle end-of-file or error
            if (gzeof(fileDescriptor)) {
                break; // End of file reached
            } else {
                int err;
                const char* errMsg = gzerror(fileDescriptor, &err);
                if (err != Z_OK) {
                    throw std::runtime_error("Error reading gzip file: " + std::string(errMsg));
                }
            }
        }
        totalBytesRead += bytesRead;
    }
    return totalBytesRead;
}

void CompressedInputStream::seekg(ssize_t pos) {
    gzseek(fileDescriptor, pos, SEEK_SET);
}

const std::string &CompressedInputStream::getFilePath() {
    return path;
}

CompressedInputStream::~CompressedInputStream() {
    inflateEnd(&strm);
    gzclose(fileDescriptor);
}

char CompressedInputStream::peek() { // TODO: Test this
    int c = gzgetc(fileDescriptor);
    if (c != EOF) {
        gzungetc(c, fileDescriptor);
    }
    return static_cast<char>(c);
}

bool CompressedInputStream::eof() {
    return gzdirect(fileDescriptor) != 0 && gzeof(fileDescriptor);
}

bool CompressedInputStream::bad() {
    int err;
    gzerror(fileDescriptor, &err);
    return (err != Z_OK) && (err != Z_STREAM_END);
}

bool CompressedInputStream::good() {
    int err;
    gzerror(fileDescriptor, &err);
    return (err == Z_OK) && !gzeof(fileDescriptor);
}

void CompressedInputStream::close() {
    int retD = deflateEnd(&strm);
    assert(retD == Z_OK);
    int retGZ = gzclose(fileDescriptor);
    assert(retGZ == Z_OK);
}
