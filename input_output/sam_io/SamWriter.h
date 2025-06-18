#ifndef GENOMICPIPELINE_SAMWRITER_H
#define GENOMICPIPELINE_SAMWRITER_H

#include <string>
#include "../alignment_io/AlignmentIO.h"
#include "../../sam_data_structures/SamAlignments.h"
#include "../../sam_data_structures/SamHeader.h"
#include "../streams/BufferWriter.h"
#include "../utils/Constant.h"

class SamWriter : public AlignmentWriter {
    BufferWriter buffer;
    std::string_view getslice(char delim, bool& isDelimFound);

public:
    explicit SamWriter(std::unique_ptr<OutputStream> outputStream, ssize_t buffSz = DEFAULT_BUFFER_SIZE) :
            AlignmentWriter(std::move(outputStream)),
            buffer(std::min(DEFAULT_BUFFER_SIZE, buffSz)) {}

    void FormatHeader(Header& hdr) override;
    void FormatAlignment(SamAlignments& reads) override;

    void Close() override;
    ~SamWriter() override;
};



#endif //GENOMICPIPELINE_SAMWRITER_H
