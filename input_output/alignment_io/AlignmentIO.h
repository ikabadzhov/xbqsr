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
#include "../streams/InputStream.h"
#include "../streams/OutputStream.h"

class AlignmentReader {
protected:
    std::unique_ptr<InputStream> in; //TODO: Lets start with file only, then move to generic input source: file and/or stdin

public:
    explicit AlignmentReader(std::unique_ptr<InputStream> inStream) : in(std::move(inStream)) {};
    virtual void ParseHeader(Header& hdr) = 0;

    virtual ssize_t ParseAlignment(std::vector<SamAlignments> &reads, std::shared_ptr<Header> &header) = 0;
    virtual void Close() = 0;
    virtual ~AlignmentReader() = default;
    bool isEOF();
    bool Good(); // Check the internal states of the stream
    bool Bad(); // It's not the opposite of good;
};

class AlignmentWriter {
protected:
    std::unique_ptr<OutputStream> out; //TODO: Lets start with file only, then move to generic input source: file and/or stdout
public:
    explicit AlignmentWriter(std::unique_ptr<OutputStream> outStream) : out(std::move(outStream)) {};
    virtual void FormatHeader(Header& hdr) = 0;
    virtual void FormatAlignment(SamAlignments& reads) = 0;
    virtual void Close() = 0;

    bool isEOF();
    bool Good(); // Check the internal states of the stream
    bool Bad();

    virtual ~AlignmentWriter() = 0;
};




#endif //GENOMICPIPELINE_ALIGNMENTIO_H
