#ifndef GENOMICPIPELINE_COMPRESSEDINPUTSTREAM_H
#define GENOMICPIPELINE_COMPRESSEDINPUTSTREAM_H

#include <cstring>
#include <iostream>
#include "PlainInputStream.h"
#include "../../zlib.h"

class CompressedInputStream : public InputStream {
    gzFile fileDescriptor;
    std::string path;
    z_stream strm{};

    static bool isValidBlockHeader(std::string& buffer);

public:
    static bool isValidFile(PlainInputStream& stream);

    explicit CompressedInputStream(const std::string &input) : path(input) {
        fileDescriptor = gzopen(input.c_str(), "rb");
        memset(&strm, 0, sizeof(strm));

        // Use 32 for GZIP decoding, for more details go here: https://www.zlib.net/manual.html#inflateInit2
        if (inflateInit2(&strm, 32 + MAX_WBITS) != Z_OK) {
            std::cerr << "Error initializing zlib for decompression." << std::endl;
            gzclose(fileDescriptor);
            return;
        }
    }

    ~CompressedInputStream() override;

    ssize_t read(std::string& buffer, ssize_t offset, ssize_t length) override;
    void seekg(ssize_t pos) override;
    const std::string & getFilePath() override;
    int peek() override;
    bool good() override;
    bool bad() override;
    bool eof() override;
    void close() override;
};


#endif //GENOMICPIPELINE_COMPRESSEDINPUTSTREAM_H
