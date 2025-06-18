#ifndef GENOMICPIPELINE_BAMWRITER_H
#define GENOMICPIPELINE_BAMWRITER_H

#include "../alignment_io/AlignmentIO.h"
#include "../streams/OutputStream.h"
#include "../streams/BufferWriter.h"
#include "../utils/Constant.h"
#include <zlib.h>

class BamWriter : public AlignmentWriter {
    BufferWriter buffer;
    std::string_view getslice(char delim, bool& isDelimFound);

public:
    explicit BamWriter(std::unique_ptr<OutputStream> outputStream, ssize_t buffSz = DEFAULT_BUFFER_SIZE) :
                        AlignmentWriter(std::move(outputStream)),
                        buffer(std::min(DEFAULT_BUFFER_SIZE, buffSz)) {}

    void FormatHeader(Header& hdr) override;
    void FormatAlignment(SamAlignments& reads) override;
    void Close() override;

    ~BamWriter() override;
};


#endif //GENOMICPIPELINE_BAMWRITER_H
