#ifndef GENOMICPIPELINE_ALIGNMENTIO_H
#define GENOMICPIPELINE_ALIGNMENTIO_H

#include <cstdint>
#include <cstddef>
#include <vector>
#include <fstream>
#include <string>
#include <memory>
#include "../../sam_data_structures/SamHeader.h"
#include "../../sam_data_structures/SamAlignments.h"
#include "../bam_io/BaiReader.h"
#include "../formatter/BamBatchFormatter.h"
#include "../streams/InputStream.h"
#include "../streams/OutputStream.h"

class AlignmentReader {
protected:
    std::unique_ptr<InputStream> in; //TODO: Lets start with file only, then move to generic input source: file and/or stdin
    VOff end_voff_ = 0;
    bool has_end_  = false;

public:
    explicit AlignmentReader(std::unique_ptr<InputStream> inStream) : in(std::move(inStream)) {};
    virtual void ParseHeader(Header& hdr) = 0;

    virtual void ParseAlignmentStreamed(std::shared_ptr<Header> header, std::function<SamAlignments()> get_batch,
                                       std::function<void(SamAlignments&&)> on_batch) = 0;
    virtual void Close() = 0;
    virtual ~AlignmentReader() = default;
    bool isEOF();
    bool Good(); // Check the internal states of the stream
    bool Bad(); // It's not the opposite of good;
    virtual void SeekVirtual(VOff vstart, VOff vend) = 0;
};

class AlignmentWriter {
protected:

public:
    AlignmentWriter() = default;
    virtual void FormatHeader(Header& hdr) = 0;
    virtual void FormatAlignment(std::vector<SamAlignments> &reads) = 0;
    virtual void BeginStreaming() = 0;
    virtual void WriteAlignmentBatch(SamAlignments& block) = 0;
    virtual void WriteFormattedBatch(OutputBatchBam& block) = 0;
    virtual void WriteFormattedBatchDirect(OutputBatchBam& block) = 0;
    virtual void EndStreaming() = 0;
    virtual void Close() = 0;


    bool isEOF();
    bool Good(); // Check the internal states of the stream
    bool Bad();

    virtual ~AlignmentWriter() = 0;
};




#endif //GENOMICPIPELINE_ALIGNMENTIO_H
