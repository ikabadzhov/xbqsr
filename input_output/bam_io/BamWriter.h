#ifndef GENOMICPIPELINE_BAMWRITER_H
#define GENOMICPIPELINE_BAMWRITER_H

#include <condition_variable>

#include "../alignment_io/AlignmentIO.h"
#include "../streams/OutputStream.h"
#include "../../input_output/streams/OutputStream.h"
#include "../../sam_data_structures/SamHeader.h"
#include "../parser/IOBatch.h"
#include "../utils/Constant.h"
#include <zlib.h>


class BamWriter : public AlignmentWriter {
public:
    // Construct with a concrete OutputStream (e.g., BGZFOutputStream)
    BamWriter(std::unique_ptr<OutputStream> sink);
    // disallow copying if not already
    BamWriter(const BamWriter&) = delete;
    BamWriter& operator=(const BamWriter&) = delete;

    void FormatHeader(Header& hdr) override;
    void FormatAlignment(std::vector<SamAlignments> &reads) override;

    void BeginStreaming() override;
    void WriteAlignmentBatch(SamAlignments& block) override;
    void EndStreaming() override;

    void Close() override;
    ~BamWriter() override;

private:

    static inline void push_u32(std::vector<uint8_t>& v, uint32_t x);
    static inline void push_u16(std::vector<uint8_t>& v, uint16_t x);
    static inline void push_u8 (std::vector<uint8_t>& v, uint8_t x);

    // --- persistent streaming state ---
    std::shared_ptr<OutputBatchBam> batch_process_;
    std::shared_ptr<OutputBatchBam> batch_transfer_;

    bool have_work_ = false;
    bool wrote_     = false;
    bool stop_      = false;
    bool streaming_started_ = false;

    std::exception_ptr writer_ex_;
    std::exception_ptr producer_ex_;


    std::mutex m_;
    //std::condition_variable cv_;
    std::condition_variable_any cv_;

    std::jthread writer_thread_;

    // helper to write one filled buffer to out_
    static void store_batch_(const std::shared_ptr<OutputBatchBam>& bt,
                             std::unique_ptr<OutputStream>& out);

    std::unique_ptr<OutputStream> out_;
    //bool header_written_ = false;

    // scratch buffer reused per record
    std::vector<uint8_t> recbuf_;

    std::once_flag     header_once_;
    std::atomic<bool>  header_written_{false};
};

#endif //GENOMICPIPELINE_BAMWRITER_H
