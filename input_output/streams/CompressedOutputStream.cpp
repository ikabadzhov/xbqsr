#include <string>
#include <sstream>
#include "CompressedOutputStream.h"
#include "../bam_io/BamReader.h"

bool CompressedOutputStream::bad() {
    return false;
}

const std::string &CompressedOutputStream::getFilePath() {
    return path;
}

bool CompressedOutputStream::eof() {
    return false;
}

bool CompressedOutputStream::good() {
    return false;
}

/*gzclose will return Z_STREAM_ERROR if file is not valid, Z_ERRNO on a
   file operation error, Z_MEM_ERROR if out of memory, Z_BUF_ERROR if the
   last read ended in the middle of a gzip stream, or Z_OK on success.*/
void CompressedOutputStream::close() {
    //before closing it, we want to be sure all the information in the buffer has been flushed
    gzflush(fileDescriptor, Z_FULL_FLUSH);
    writeEOF();
    int retD = deflateEnd(&strm);
    assert(retD == Z_OK);
    int retGZ = gzclose(fileDescriptor);
    assert(retGZ == Z_OK);
}

void CompressedOutputStream::write(std::string &buffer, ssize_t n) {

    int sizeWritten = gzwrite(fileDescriptor, buffer.data(), buffer.size());

    /*
     * On error gzwrite() shall return a value less than or equal to 0.
     * Z_ERRNO: An underlying base library function has indicated an error. The global variable errno may be examined for further information.
     * Z_STREAM_ERROR: The stream is invalid, is not open for writing, or is in an invalid state.
     * Z_BUF_ERROR: no compression progress is possible (see deflate()).
     * Z_MEM_ERROR: Insufficient memory available to compress.
    */
    if (sizeWritten <=0) {
        const char *error = gzerror(fileDescriptor, nullptr);
        close();
        throw std::runtime_error("Error writing compressed data to file: " + std::string(error));
    }
}

void CompressedOutputStream::writeEOF(){
    gzwrite(fileDescriptor, knownSequence, sizeof(knownSequence));
}
