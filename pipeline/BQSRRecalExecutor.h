// BQSRRecalExecutor.h
#pragma once
#include <memory>
#include <string>
#include <sycl/sycl.hpp>

#include "IExecutor.h"

// your project headers
#include "../sam_data_structures/SamAlignments.h"
#include "../sam_data_structures/SamHeader.h"
#include "../file_format/MappedFasta.h"
#include "../bqsr/BQSR.h"


struct RecalCtx {
    sycl::device dev_;
    sycl::queue  q_;
    std::shared_ptr<BaseRecalibratorTables> tables_;
    std::shared_ptr<BaseRecalibrator>       recal_;
};

class BQSRRecalExecutor final : public IExecutor {
public:
    BQSRRecalExecutor(const sycl::device& dev,
                      const BQSRTableConfig& table_cfg,
                      MappedFasta& ref,
                      KnownSites& known,
                      Header& hdr)
        : dev_(dev),
          q_(dev_)
    {
        // Each executor has its OWN tables + recalibrator (no sharing across threads/devices)
        tables_ = std::make_shared<BaseRecalibratorTables>(table_cfg, q_);
        recal_  = std::make_shared<BaseRecalibrator>(ref, known, q_, hdr, *tables_);
        recal_->initializeDeviceMemory();
    }

    void Process(SamAlignments& batch) override {
        recal_->Recalibrate(batch);
        // If Recalibrate enqueues async kernels and returns:
        // consider NOT waiting here (for overlap), but make sure you flush before merging.
        // q_.wait(); // (only if you need strict sync per-batch)
    }

    void Finalize() override {
        // Ensure all queued work is finished before export/merge
        q_.wait_and_throw();
        recal_->exportTables(*tables_);
        q_.wait_and_throw();
    }

    const char* Name() const { return "BQSRRecalExecutor"; }

    // Expose local tables for merging at the end
    BaseRecalibratorTables& Tables() { return *tables_; }
    const BaseRecalibratorTables& Tables() const { return *tables_; }

    BaseRecalibrator& Recal() { return *recal_; }
    const BaseRecalibrator& Recal() const { return *recal_; }

    const sycl::device& Device() const { return dev_; }
    const sycl::queue&  Queue()  const { return q_;  }

private:
    sycl::device dev_;
    sycl::queue  q_;

    std::shared_ptr<BaseRecalibratorTables> tables_;
    std::shared_ptr<BaseRecalibrator>       recal_;
};