#ifndef GENOMICPIPELINE_SAMREADER_H
#define GENOMICPIPELINE_SAMREADER_H

#include <string>
#include "../alignment_io/AlignmentIO.h"
#include "../../sam_data_structures/SamHeader.h"
#include "../../sam_data_structures/SamAlignments.h"
#include "../utils/StringScanner.h"
#include "../streams/PlainInputStream.h"
#include "../streams/BufferReader.h"
#include "../utils/Constant.h"
#include "../parser/IOBatch.h"

class SamReader : public AlignmentReader {
    BufferReader buffer;

    std::string_view getslice(char delim, bool& isDelimFound);

public:
    explicit SamReader(std::unique_ptr<InputStream> inputStream, ssize_t buffSz = DEFAULT_BUFFER_SIZE) :
        AlignmentReader(std::move(inputStream)),
        buffer(std::min(DEFAULT_BUFFER_SIZE, buffSz)) {}

    void ParseHeader(Header& hdr) override;
    ssize_t ParseAlignment(std::vector<SamAlignments> &reads, std::shared_ptr<Header> &header) override;

    void Close() override;
};


#endif //GENOMICPIPELINE_SAMREADER_H
