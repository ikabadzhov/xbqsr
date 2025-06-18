#ifndef GENOMICPIPELINE_BAMREADER_H
#define GENOMICPIPELINE_BAMREADER_H


#include "../alignment_io/AlignmentIO.h"
#include "../../zlib.h"
#include "../streams/CompressedInputStream.h"
#include "../streams/BufferReader.h"
#include "../utils/Constant.h"


class BamReader : public AlignmentReader {

    BufferReader buffer;
    std::string_view getslice(size_t numBytesToRead, size_t &numBytesRead);

public:
    explicit BamReader(std::unique_ptr<InputStream> inputStream, ssize_t buffSz = DEFAULT_BUFFER_SIZE) : AlignmentReader(std::move(inputStream)),
                                                                                                         buffer(std::min(DEFAULT_BUFFER_SIZE, buffSz)) {
        in->seekg(static_cast<ssize_t>(BAM_MAGIC.size()));
    }

    void ParseHeader(Header& hdr) override;
    void Close() override;
    ssize_t ParseAlignment(SamAlignments& reads) override;
};


#endif //GENOMICPIPELINE_BAMREADER_H
