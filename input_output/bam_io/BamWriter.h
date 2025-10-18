#ifndef GENOMICPIPELINE_BAMWRITER_H
#define GENOMICPIPELINE_BAMWRITER_H

#include "../alignment_io/AlignmentIO.h"
#include "../streams/OutputStream.h"
#include "../../sam_data_structures/SamHeader.h"
#include "../utils/Constant.h"
#include <zlib.h>

class BamWriter : public AlignmentWriter {
public:
    // Construct with a concrete OutputStream (e.g., BGZFOutputStream)
    BamWriter(std::unique_ptr<OutputStream> sink);

    void FormatHeader(Header& hdr) override;
    void FormatAlignment(std::vector<SamAlignments> &reads) override;
    void Close() override;

private:

    static inline void push_u32(std::vector<uint8_t>& v, uint32_t x);
    static inline void push_u16(std::vector<uint8_t>& v, uint16_t x);
    static inline void push_u8 (std::vector<uint8_t>& v, uint8_t x);
private:
    std::unique_ptr<OutputStream> out_;
    bool header_written_ = false;

    // scratch buffer reused per record
    std::vector<uint8_t> recbuf_;
};

#endif //GENOMICPIPELINE_BAMWRITER_H
