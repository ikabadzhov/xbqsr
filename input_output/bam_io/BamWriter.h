// BamWriter.h
#pragma once
#include <memory>
#include <vector>
#include <mutex>
#include <condition_variable>
#include <stop_token>
#include <thread>
#include <exception>
#include <atomic>

#include "../parser/IOBatch.h"
#include "../formatter/BamBatchFormatter.h"
#include "../streams/OutputStream.h"
#include "../alignment_io/AlignmentIO.h"


class BamWriter : public AlignmentWriter {
public:
    explicit BamWriter(std::unique_ptr<OutputStream> sink);

    // Back-compatible API (same idea as your current implementation)
    void FormatHeader(Header& hdr) override;
    void BeginStreaming() override;
    void WriteAlignmentBatch(SamAlignments& block) override;

    void WriteFormattedBatchDirect(OutputBatchBam &batch) override;

    void WriteFormattedBatch(OutputBatchBam& batch) override;
    void EndStreaming() override;
    void FormatAlignment(std::vector<SamAlignments>& blocks) override;
    void Close() override;

    ~BamWriter();

private:
    // Writes bytes from bt to out, then resets bt.
    void store_batch_(const std::shared_ptr<OutputBatchBam>& bt, std::unique_ptr<OutputStream>& out);

private:
    std::unique_ptr<OutputStream> out_;

    // Formatting is now separated and reusable
    BamBatchFormatter formatter_;
    std::vector<uint8_t> header_buf_;

    // Old behavior preserved
    std::once_flag header_once_;
    std::atomic<bool> header_written_{false};
    bool streaming_started_{false};

    // Writer thread state (double buffer)
    std::shared_ptr<OutputBatchBam> batch_process_;
    std::shared_ptr<OutputBatchBam> batch_transfer_;

    std::jthread writer_thread_;
    std::mutex m_;
    std::condition_variable_any cv_;
    bool have_work_{false};
    bool wrote_{false};
    bool stop_{false};

    std::exception_ptr writer_ex_{nullptr};
    std::exception_ptr producer_ex_{nullptr};
};
