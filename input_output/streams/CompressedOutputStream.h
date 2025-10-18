#ifndef GENOMICPIPELINE_COMPRESSEDOUTPUTSTREAM_H
#define GENOMICPIPELINE_COMPRESSEDOUTPUTSTREAM_H

#include <string>
#include <cstring>
#include <iostream>
#include "../../zlib.h"
#include "OutputStream.h"

class CompressedOutputStream {
    gzFile fileDescriptor;
    std::string path;
    z_stream strm{};
    std::string filename;

public:
    explicit CompressedOutputStream(const std::string &filename) : path(filename) {
        this->filename = filename;
        fileDescriptor = gzopen(filename.c_str(), "wb");
        //fileDescriptor = gzopen("/media/pipeline-output/fileACaso.txt", "wb");
        memset(&strm, 0, sizeof(strm));

        if(fileDescriptor==nullptr)
            return;

        // Use 32 for GZIP decoding, for more details go here: https://www.zlib.net/manual.html#inflateInit2
        if(deflateInit(&strm, Z_DEFAULT_COMPRESSION) != Z_OK){
        //if (inflateInit2(&strm, 32 + MAX_WBITS) != Z_OK) {
            std::cerr << "Error initializing zlib for decompression." << std::endl;
            gzclose(fileDescriptor);
            return;
        }
    }

};


#endif //GENOMICPIPELINE_COMPRESSEDOUTPUTSTREAM_H
