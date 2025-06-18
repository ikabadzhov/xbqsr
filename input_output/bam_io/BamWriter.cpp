#include "BamWriter.h"
#include "../formatter/BamHeaderFormatter.h"
#include "../formatter/BamAlignmentFormatter.h"
#include "../formatter/SamAlignmentFormatter.h"

void BamWriter::FormatHeader(Header &hdr) {
    std::string headerow;
    BamHeaderFormatter::FormatBamHeader(hdr, headerow);

    size_t startPos = 0;
    size_t bufferSize = buffer.getSize();

    while (startPos < headerow.size()) {
        size_t remainingSize = headerow.size() - startPos;
        size_t chunkSize = std::min(bufferSize, remainingSize);
        std::string partialContent = headerow.substr(startPos, chunkSize);
        startPos += chunkSize;

        buffer.addSlice(partialContent);
        this->out->write(partialContent, -1);
        buffer.flushContent();
    }

    if (buffer.getLastPos() > 0) {
        std::string remainedContent;
        buffer.getContent(0, buffer.getLastPos(), remainedContent);
        this->out->write(remainedContent, -1);
        buffer.flushContent();
    }
}

void BamWriter::FormatAlignment(SamAlignments &reads) {
    for(int idx=0; idx<reads.GetAlnsNum(); idx++) {
        std::string alignrow;
        BamAlignmentFormatter::FormatBamAlignment(reads, idx, alignrow);
        if(alignrow.size() < buffer.available_data()){
            buffer.addSlice(alignrow);
        }else {
            size_t startPos = buffer.getLastPos();
            size_t bufferSize = buffer.getSize();
            size_t remainingSize = bufferSize - startPos;

            buffer.addSlice(alignrow.substr(0, remainingSize));
            std::string partialContent;
            buffer.getContent(0, buffer.getLastPos(), partialContent);
            this->out->write(partialContent, -1);
            buffer.flushContent();

            if (alignrow.size() - remainingSize < buffer.available_data()){
                buffer.addSlice(alignrow.substr(remainingSize, alignrow.size()));
            }else{
                size_t chunkSize = alignrow.size()-remainingSize;
                std::string partialContent(&alignrow[remainingSize], chunkSize);
                this->out->write(partialContent, -1);
            }
        }
    }
    if(buffer.getLastPos()>0) {
        std::string remainedContent;
        buffer.getContent(0, buffer.getLastPos(), remainedContent);
        this->out->write(remainedContent, -1);
        buffer.flushContent();
    }
    Close();
}

void BamWriter::Close() {
    this->out->close();
}

BamWriter::~BamWriter() = default;

