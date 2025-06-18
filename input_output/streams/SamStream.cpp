#include "SamStream.h"
#include "../bam_io/BamReader.h"

ssize_t SamStream::readBytes(PlainInputStream& stream, std::string& buffer, ssize_t offset, ssize_t length) {
    ssize_t bytesRead;
    bytesRead = stream.read(buffer, offset, length);
    return bytesRead;
}

bool SamStream::isBAMFile(PlainInputStream& fileStream){
    if (!CompressedInputStream::isValidFile(fileStream)) {
        return false;
    } else {
        CompressedInputStream cbi(fileStream.getFilePath());
        std::string bamMagicString(4, 0);
        cbi.read(bamMagicString, 0, 4);
        return BAM_MAGIC == bamMagicString;
    }
}