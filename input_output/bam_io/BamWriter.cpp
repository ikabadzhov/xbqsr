// BamWriter.cpp
#include "BamWriter.h"
#include "../../timer/TraceTimeline.h"
#include <stdexcept>

BamWriter::BamWriter(std::unique_ptr<OutputStream> sink) : out_(std::move(sink)) {}

void BamWriter::store_batch_(const std::shared_ptr<OutputBatchBam>& bt,
                            std::unique_ptr<OutputStream>& out) {
    const size_t n = bt->used();
    if (n) out->write(std::span<const uint8_t>(bt->data(), n));
    bt->reset();
}

void BamWriter::FormatHeader(Header& hdr) {
    std::call_once(header_once_, [&] {
        if (!out_) throw std::runtime_error("BamWriter::FormatHeader: no output stream");

        formatter_.FormatHeaderTo(hdr, header_buf_);
        if (!header_buf_.empty()) {
            out_->write(std::span<const uint8_t>(header_buf_.data(), header_buf_.size()));
        }
        header_written_.store(true, std::memory_order_release);
    });
}

void BamWriter::BeginStreaming() {
    if (!header_written_.load(std::memory_order_acquire))
        throw std::runtime_error("BamWriter::BeginStreaming: header must be written first");
    if (!out_) throw std::runtime_error("BamWriter::BeginStreaming: no output stream");
    if (streaming_started_) return;

    writer_ex_ = nullptr;
    producer_ex_ = nullptr;
    have_work_ = wrote_ = stop_ = false;

    batch_process_  = std::make_shared<OutputBatchBam>(/*cap*/ 1 << 20);
    batch_transfer_ = std::make_shared<OutputBatchBam>(/*cap*/ 1 << 20);

    writer_thread_ = std::jthread([this](std::stop_token st) {
        try {
            for (;;) {
                std::unique_lock lk(m_);
                cv_.wait(lk, st, [&]{ return have_work_ || stop_; });

                if (st.stop_requested()) break;
                if (stop_ && !have_work_) break;

                have_work_ = false;
                auto bt = batch_transfer_;
                lk.unlock();

                store_batch_(bt, out_);

                lk.lock();
                wrote_ = true;
                lk.unlock();
                cv_.notify_one();
            }
        } catch (...) {
            writer_ex_ = std::current_exception();
            std::lock_guard<std::mutex> lk(m_);
            wrote_ = true;
            cv_.notify_one();
        }
    });

    streaming_started_ = true;
}

void BamWriter::WriteAlignmentBatch(SamAlignments& block) {
    // Back-compatible: writer still takes SamAlignments and writes to the stream.
    FormatHeader(block.GetHdr());
    if (!streaming_started_) BeginStreaming();

    try {
        batch_process_->reset();

        // ---- formatting is now separated: formatter fills a flat byte batch ----
        formatter_.FormatAlignmentBatchTo(block, *batch_process_);

        if (batch_process_->used() > 0) {
            {
                std::unique_lock lk(m_);
                wrote_ = false;
                std::swap(batch_process_, batch_transfer_);
                have_work_ = true;
                lk.unlock();
                cv_.notify_one();
            }

            {
                std::unique_lock lk(m_);
                cv_.wait(lk, [&]{ return wrote_ || writer_ex_ != nullptr; });
            }
            if (writer_ex_) std::rethrow_exception(writer_ex_);
        }
    } catch (...) {
        producer_ex_ = std::current_exception();
        std::rethrow_exception(producer_ex_);
    }
}

void BamWriter::WriteFormattedBatchDirect(OutputBatchBam& batch) {
    if (!header_written_.load(std::memory_order_acquire))
        throw std::runtime_error("BamWriter::WriteFormattedBatchDirect: header must be written first");
    if (!out_) throw std::runtime_error("BamWriter::WriteFormattedBatchDirect: no output stream");

    const size_t n = batch.used();
    if (n == 0) return;

    out_->write(std::span<const uint8_t>(batch.data(), n));
}

void BamWriter::WriteFormattedBatch(OutputBatchBam &batch) {
    // Back-compatible: writer still takes SamAlignments and writes to the stream.
    //FormatHeader(block.GetHdr());
    if (!streaming_started_) BeginStreaming();

    try {
        const size_t n = batch.used();
        if (n == 0) return;

        batch_process_->reset();
        batch_process_->ensure_capacity(n);
        std::memcpy(batch_process_->data(), batch.data(), n);
        batch_process_->seal(n, batch.n_records()); // if you track it; else pass 0

        if (batch_process_->used() > 0) {
            {
                std::unique_lock lk(m_);
                wrote_ = false;
                std::swap(batch_process_, batch_transfer_);
                have_work_ = true;
                lk.unlock();
                cv_.notify_one();
            }

            {
                std::unique_lock lk(m_);
                cv_.wait(lk, [&]{ return wrote_ || writer_ex_ != nullptr; });
            }
            if (writer_ex_) std::rethrow_exception(writer_ex_);
        }
    } catch (...) {
        producer_ex_ = std::current_exception();
        std::rethrow_exception(producer_ex_);
    }
}

void BamWriter::EndStreaming() {
    if (!streaming_started_) return;

    try {
        // Flush any remaining bytes (should usually be empty because we hand off immediately)
        if (batch_process_ && batch_process_->used() > 0) {
            {
                std::unique_lock lk(m_);
                wrote_ = false;
                std::swap(batch_process_, batch_transfer_);
                have_work_ = true;
                lk.unlock();
                cv_.notify_one();
            }
            {
                std::unique_lock lk(m_);
                cv_.wait(lk, [&]{ return wrote_ || writer_ex_ != nullptr; });
            }
            if (writer_ex_) std::rethrow_exception(writer_ex_);
        }
    } catch (...) {
        producer_ex_ = std::current_exception();
    }

    {
        std::lock_guard lk(m_);
        stop_ = true;
    }
    writer_thread_.request_stop();
    cv_.notify_one();

    if (writer_thread_.joinable()) writer_thread_.join();
    streaming_started_ = false;

    if (producer_ex_) std::rethrow_exception(producer_ex_);
    if (writer_ex_)   std::rethrow_exception(writer_ex_);
}

void BamWriter::FormatAlignment(std::vector<SamAlignments>& blocks) {
    if (!out_) throw std::runtime_error("BamWriter::FormatAlignment: no output stream");

    // Keep old behavior: write header from first block (if any), then stream all.
    if (!blocks.empty()) FormatHeader(blocks.front().GetHdr());
    BeginStreaming();
    for (auto& blk : blocks) {
        WriteAlignmentBatch(blk);
    }
    EndStreaming();
}

void BamWriter::Close() {
    if (streaming_started_) {
        EndStreaming();
    }
    if (out_) out_->close();
}

BamWriter::~BamWriter() {
    // best-effort cleanup
    try {
        if (streaming_started_) EndStreaming();
    } catch (...) {}
}