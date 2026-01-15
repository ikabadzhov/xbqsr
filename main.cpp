#include <iostream>
#include <numeric>
#include <algorithm>
#include <chrono>
#include <optional>
#include <vector>
#include <memory>
#include <stdexcept>

#include <tbb/parallel_for.h>
#include <tbb/parallel_sort.h>

#include "Utils.h"
#include "input_output/alignment_io/AlignmentIO.h"
#include "input_output/bam_io/BamReader.h"
#include "sam_data_structures/SamHeader.h"
#include "filters/Filter.h"
#include "file_format/MappedFasta.h"
#include "bqsr/BQSR.h"
#include "sam_data_structures/SamAlignments.h"
#include "pipeline/BoundedChannel.h"
#include "pipeline/StreamingEngine.h"
#include "bam_data_structures/BamHeader.h"
#include "input_output/alignment_io/AlignmentReaderFactory.h"
#include "input_output/alignment_io/AlignmentWriterFactory.h"
#include "sam_data_structures/sam_fields/Tag.h"
#include "timer/Timer.h"

#include <sycl/queue.hpp>
#include <sycl/sycl.hpp>

#include "device/batch.h"
#include "input_output/bam_io/BaiReader.h"
#include "logger/Logger.h"

// NOTE: Executors headers no longer needed
// #include "pipeline/BQSRApplyExecutor.h"
// #include "pipeline/BQSRRecalExecutor.h"

#include "input_output/bam_io/BamWriter.h"
#include "input_output/formatter/BamBatchFormatter.h"
#include "sycl_parallel_library/intervalsmove/intervalmove.h"
#include "sycl_parallel_library/scan/scan.h"
#include "sycl_parallel_library/segmented_reduce/segreduce.h"
#include "sycl_parallel_library/segmented_scan/segscan.h"
#include "timer/ProfilerIO.h"

using namespace std;
using namespace chrono;

static inline int32_t le_i32(const uint8_t* p) {
    return int32_t(uint32_t(p[0]) |
                   (uint32_t(p[1]) << 8) |
                   (uint32_t(p[2]) << 16) |
                   (uint32_t(p[3]) << 24));
}

system_clock::time_point now() {
    return std::chrono::high_resolution_clock::now();
}

static std::mutex gpu_mtx;


size_t get_duration(std::chrono::system_clock::time_point start,
                    std::chrono::system_clock::time_point end) {
    return std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
}

// -------------------------
// Context structs (no IExecutor)
// -------------------------
struct RecalCtx {
    sycl::device dev_;
    sycl::queue  q_;
    std::shared_ptr<BaseRecalibratorTables> tables_;
    std::shared_ptr<BaseRecalibrator>       recal_;
};

struct ApplyCtx {
    sycl::device dev_;
    sycl::queue  q_;
    int quantize_level_{0};
    std::vector<uint8_t> sqqList_;
    std::shared_ptr<BaseRecalibratorTables> tables_;
};

//[[clang::optnone]]
int main(int argc, char **argv) {
    Arguments args = parse_arguments(argc, argv);

    // Enumerate all GPU devices available
    auto platforms = sycl::platform::get_platforms();
    std::vector<sycl::device> gpus;

    for (const auto& plat : platforms) {
        for (const auto& dev : plat.get_devices(sycl::info::device_type::gpu)) {
            gpus.push_back(dev);
        }
    }

    if (gpus.empty()) {
        std::cerr << "No GPU devices found!\n";
        return 1;
    }

    Logger::getLogger().log(LOG_INFO, "Available devices:");
    for (size_t i = 0; i < gpus.size(); ++i) {
        Logger::getLogger().log(LOG_INFO, "  [", i, "] ",
            gpus[i].get_info<sycl::info::device::name>(), " (",
            gpus[i].get_info<sycl::info::device::vendor>(), ")");
    }

    // Determine which devices to use
    std::vector<int> selected_devices;
    if (args.gpu_use_all) {
        for (size_t i = 2; i < gpus.size(); ++i) {
            selected_devices.push_back(static_cast<int>(i));
        }
    } else {
        selected_devices = args.gpu_devices;
    }

    // Validate device indices
    for (int dev_idx : selected_devices) {
        if (dev_idx < 0 || dev_idx >= static_cast<int>(gpus.size())) {
            Logger::getLogger().log(LOG_FATAL, "Invalid device index [", dev_idx,
                                    "]. Available devices: 0-", gpus.size() - 1);
            return 1;
        }
    }

    Logger::getLogger().log(LOG_INFO, "Selected device(s):");
    for (int dev_idx : selected_devices) {
        Logger::getLogger().log(LOG_INFO, "  [", dev_idx, "] ",
            gpus[dev_idx].get_info<sycl::info::device::name>(), " (",
            gpus[dev_idx].get_info<sycl::info::device::vendor>(), ")");
    }

    if (args.verbose) {
        Logger::getLogger().set_min_severity(LOG_DEBUG);
    } else {
        Logger::getLogger().set_min_severity(LOG_INFO);
    }

    std::vector<StreamingEngine::Stage> stages;
    std::vector<StreamingEngine::TransformFn> transforms;
    std::vector<StreamingEngine::FinalizeFn> finals;

    // ---------- Build slices / readers ----------
    const size_t n_producers = 3 * selected_devices.size();
    const std::string bai_path = find_bai_path(args.input_bam);

    auto bai = load_bai_linear(bai_path);
    std::vector<std::pair<VOff,VOff>> ranges =
        split_ranges_equal_by_compressed_span(bai.ioffsets, 512UL*1024*1024);

    //print_ranges_estimated_decompressed_over(args.input_bam, ranges);

    Logger::getLogger().log(LOG_INFO, "Partitioning slices");

    std::vector<std::vector<Slice>> slices(n_producers);
    for (size_t i = 0; i < ranges.size(); ++i) {
        slices[i % n_producers].emplace_back(ranges[i].first, ranges[i].second, i);
    }

    Logger::getLogger().log(LOG_INFO, "Opening " + std::to_string(n_producers) + " readers");

    std::vector<std::shared_ptr<AlignmentReader>> alignmentReader(n_producers);
    std::shared_ptr<AlignmentReader> alignmentForHeader = AlignmentReaderFactory::open(args.input_bam);

    tbb::parallel_for(0UL, n_producers, [&](size_t i) {
        alignmentReader[i] = AlignmentReaderFactory::open(args.input_bam);
        auto bam = dynamic_cast<BamReader*>(alignmentReader[i].get());
        if (!bam) throw std::runtime_error("AlignmentReader is not a BamReader");
        bam->SetSlices(std::move(slices[i]));
    });

    Logger::getLogger().log(LOG_INFO, "Parsing header");

    std::shared_ptr<Header> header = std::make_shared<Header>();
    alignmentForHeader->ParseHeader(*header);
    alignmentForHeader.reset();

    Logger::getLogger().log(LOG_INFO, "Building pipeline");
    header->buildReadGroupCovariateIndex();

    BQSRTableConfig bqsrTableConfig(
        header->getReadGroupCovariates().size(),
        bqsr::defaults::MAX_QUALITY_SCORE,
        bqsr::defaults::MAX_CYCLE,
        bqsr::defaults::CONTEXT_LEN
    );

    // ---------- Shared input-only objects ----------
    std::unique_ptr<MappedFasta> referenceGenome;
    std::unique_ptr<KnownSites>  knownSites;

    // ---------- Per-worker contexts ----------
    std::vector<RecalCtx> recal_ctxs;
    std::vector<ApplyCtx> apply_ctxs;

    std::shared_ptr<AlignmentWriter> alignmentWriter;

    // ---------- Recalibrate ----------
    if (args.is_recalibrate()) {
        referenceGenome = std::make_unique<MappedFasta>(args.reference, header);
        knownSites      = std::make_unique<KnownSites>(args.known_sites, header);

        recal_ctxs.reserve(selected_devices.size());

        for (int dev_idx : selected_devices) {
            sycl::device dev = gpus[dev_idx];

            RecalCtx ctx;
            ctx.dev_ = dev;

            // Optional: pass async handler to catch exceptions early
            ctx.q_ = sycl::queue(dev, [](sycl::exception_list el) {
                for (auto& e : el) {
                    try { std::rethrow_exception(e); }
                    catch (const std::exception& ex) {
                        std::fprintf(stderr, "[SYCL async] %s\n", ex.what());
                    }
                }
            });

            ctx.tables_ = std::make_shared<BaseRecalibratorTables>(bqsrTableConfig, ctx.q_);
            ctx.recal_  = std::make_shared<BaseRecalibrator>(*referenceGenome, *knownSites, ctx.q_, *header, *ctx.tables_);
            ctx.recal_->initializeDeviceMemory();

            recal_ctxs.push_back(std::move(ctx));
        }

        stages.emplace_back("BQSR-Recalibrate",
            [&contexts=recal_ctxs](SamAlignments& batch, size_t wid) {
                // device-local processing
                   // std::unique_lock<mutex> lk(gpu_mtx);
                contexts[wid].recal_->Recalibrate(batch);
                contexts[wid].q_.wait_and_throw();
            }
        );

        finals.emplace_back([&recal_ctxs, args, header] {
            for (auto& ctx : recal_ctxs) {
                ctx.q_.wait_and_throw();
                ctx.recal_->exportTables(*ctx.tables_);
                ctx.q_.wait_and_throw();
            }

            // merge into worker 0
            for (size_t i = 1; i < recal_ctxs.size(); ++i) {
                recal_ctxs[0].tables_->merge(*recal_ctxs[i].tables_);
            }

            recal_ctxs[0].tables_->FinalizeBQSRTables();
            recal_ctxs[0].tables_->PrintBQSRTables(args.recal_table, *header);
        });
    }

    // ---------- Apply ----------
    int quantize_level = 0;
    std::vector<uint8_t> sqqList;

    if (args.is_apply()) {
        alignmentWriter = AlignmentWriterFactory::open(args.output_bam);

        apply_ctxs.reserve(selected_devices.size());

        for (int dev_idx : selected_devices) {
            sycl::device dev = gpus[dev_idx];

            ApplyCtx ctx;
            ctx.dev_ = dev;

            ctx.q_ = sycl::queue(dev, [](sycl::exception_list el) {
                for (auto& e : el) {
                    try { std::rethrow_exception(e); }
                    catch (const std::exception& ex) {
                        std::fprintf(stderr, "[SYCL async] %s\n", ex.what());
                    }
                }
            });

            ctx.quantize_level_ = quantize_level;
            ctx.sqqList_        = sqqList;

            ctx.tables_ = std::make_shared<BaseRecalibratorTables>(bqsrTableConfig, ctx.q_);
            ctx.tables_->initializeDeviceMemory();
            ctx.tables_->LoadFromReport(args.recal_table, *header);

            apply_ctxs.push_back(std::move(ctx));
        }

        stages.emplace_back("BQSR-Apply",
            [&apply_ctxs](SamAlignments& batch, size_t wid) {
                auto& ctx = apply_ctxs[wid];

                // If ApplyBQSR expects a mutable copy, keep it local:
                auto sqqLocal = ctx.sqqList_;
                ctx.tables_->ApplyBQSR(batch, ctx.quantize_level_, sqqLocal);

                ctx.q_.wait_and_throw();
            }
        );

        alignmentWriter->FormatHeader(*header);
    }

    // ---------- Build and run engine ----------
    StreamingEngine::Config cfg;
    //cfg.channel_capacity = 4;

    // workers = number of devices used for the chosen mode
    if (args.is_recalibrate()) cfg.workers = recal_ctxs.size();
    if (args.is_apply())       cfg.workers = apply_ctxs.size();

    if (cfg.workers == 0) {
        throw std::runtime_error("No workers created (workers would be 0)");
    }

    Timer t;
    t.reset();

    StreamingEngine engine(cfg, std::move(stages), std::move(finals), alignmentWriter);

    // If your current StreamingEngine still takes `executors`, you must remove that parameter.
    engine.run_streaming(alignmentReader, header);

    Logger::getLogger().log(LOG_INFO, "Time to run all stages: " + std::to_string(t.elapsed()) + " sec");
    ProfilerIO::instance().report();
    return 0;
}