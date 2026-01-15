#include <iostream>
#include <numeric>
#include <algorithm>
#include <chrono>
#include <tbb/parallel_sort.h>
#include "Utils.h"
#include "input_output/alignment_io/AlignmentIO.h"
#include "input_output/bam_io/BamReader.h"
#include "sam_data_structures/SamHeader.h"
#include "filters/Filter.h"
#include "file_format/MappedFasta.h"
#include "bqsr/BQSR.h"
#include "sam_data_structures/SamAlignments.h"
#include "pipeline/StreamingChannel.h"
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
#include "pipeline/BQSRRecalExecutor.h"
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

size_t get_duration(std::chrono::system_clock::time_point start, std::chrono::system_clock::time_point end) {
    return std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
}

//[[clang::optnone]]
int main(int argc, char **argv) {
    // Parse all command-line arguments
    Arguments args = parse_arguments(argc, argv);

    /* Input parsing variable section */
    std::vector<SamAlignments> alignments;

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
        // Use all available devices
        for (size_t i = 2; i < gpus.size(); ++i) {
            selected_devices.push_back(static_cast<int>(i));
        }
    } else {
        // Use specified devices
        selected_devices = args.gpu_devices;
    }

    // Validate device indices
    for (int dev_idx : selected_devices) {
        if (dev_idx < 0 || dev_idx >= static_cast<int>(gpus.size())) {
            Logger::getLogger().log(LOG_FATAL, "Invalid device index [", dev_idx, "]. Available devices: 0-", gpus.size() - 1);
            return 1;
        }
    }

    Logger::getLogger().log(LOG_INFO, "Selected device(s):");
    for (int dev_idx : selected_devices) {
        Logger::getLogger().log(LOG_INFO, "  [", dev_idx, "] ",
            gpus[dev_idx].get_info<sycl::info::device::name>(), " (",
            gpus[dev_idx].get_info<sycl::info::device::vendor>(), ")");
    }

    // For now, use the first selected device (single-GPU support)
    sycl::device device = gpus[selected_devices[0]];

    sycl::queue queue(device);
    /**** Input parsing for reading Alignment Data *****/

    if (args.verbose) {
        Logger::getLogger().set_min_severity(LOG_DEBUG);
    } else {
        Logger::getLogger().set_min_severity(LOG_INFO);
    }

    std::vector<StreamingEngine::Stage> stages;
    std::vector<StreamingEngine::FinalizeFn> finals;

    size_t n_producers = selected_devices.size();
    const std::string bai_path = find_bai_path(args.input_bam);
    auto bai = load_bai_linear(bai_path);
    std::vector<std::pair<VOff,VOff>> ranges =
    split_ranges_equal_by_compressed_span(bai.ioffsets, 512UL*1024*1024);
    print_ranges_estimated_decompressed_over(args.input_bam, ranges);

    /* CompressedInputStream cis(args.input_bam);

    //const size_t n = (max_to_check == 0) ? slices.size() : std::min(slices.size(), max_to_check);

    for (size_t i = 0; i < ranges.size(); ++i) {
    const uint64_t beg = ranges[i].first;
    const uint64_t end = ranges[i].second;

    try {
        cis.set_virtual_end(end);
        cis.seek_virtual(beg);

        // Read the first 1 MiB of *decompressed* data from this slice
        const size_t probe_size = 1 << 20; // 1 MiB
        std::string buf(probe_size, '\0');

        const ssize_t got = cis.read(buf, 0, probe_size);
        if (got < 4) {
            std::ostringstream oss;
            oss << "Slice too small to contain a BAM record header (got=" << got << ")"
                << " at slice#" << i
                << " beg=(" << voff_coff(beg) << "," << voff_uoff(beg) << ")"
                << " end=(" << voff_coff(end) << "," << voff_uoff(end) << ")";
            throw std::runtime_error(oss.str());
        }

        const unsigned char* p =
            reinterpret_cast<const unsigned char*>(buf.data());

        const int32_t block_size =
            int32_t(uint32_t(p[0]) |
                    (uint32_t(p[1]) << 8) |
                    (uint32_t(p[2]) << 16) |
                    (uint32_t(p[3]) << 24));

        if (block_size <= 0 || block_size > (64 * 1024 * 1024)) {
            std::ostringstream oss;
            oss << "Invalid BAM block_size=" << block_size
                << " at slice#" << i
                << " beg=(" << voff_coff(beg) << "," << voff_uoff(beg) << ")"
                << " end=(" << voff_coff(end) << "," << voff_uoff(end) << ")";
            throw std::runtime_error(oss.str());
        }

        std::cout << "[SliceCheck #" << i << "] ok: block_size=" << block_size
                  << " beg=(" << voff_coff(beg) << "," << voff_uoff(beg) << ")"
                  << " end=(" << voff_coff(end) << "," << voff_uoff(end) << ")"
                  << "\n";
    } catch (...) {
        std::cerr << "[SliceCheck ERROR #" << i << "] "
                  << "beg=(" << voff_coff(beg) << "," << voff_uoff(beg) << ") "
                  << "end=(" << voff_coff(end) << "," << voff_uoff(end) << ")\n";
        throw;
    }
    }*/

    std::vector<std::vector<std::pair<VOff,VOff>>> slices(n_producers);
    for (size_t i = 0; i < ranges.size(); ++i) {
        slices[i % n_producers].emplace_back(ranges[i].first, ranges[i].second);
    }

    std::vector<std::shared_ptr<AlignmentReader>> alignmentReader(n_producers);
    std::shared_ptr<AlignmentReader> alignmentForHeader = AlignmentReaderFactory::open(args.input_bam);

    tbb::parallel_for (0UL, n_producers,[&] (size_t i) {
        alignmentReader[i] = AlignmentReaderFactory::open(args.input_bam);
        auto bam = dynamic_cast<BamReader*>(alignmentReader[i].get());
        bam->SetSlices(std::move(slices[i]));
    });

    std::shared_ptr<AlignmentWriter> alignmentWriter;

    std::shared_ptr<Header> header = std::make_shared<Header>();
    alignmentForHeader->ParseHeader(*header);
    alignmentForHeader.reset();

    BaiLinearIndex baiLinearIndex = load_bai_linear(bai_path);

    header->buildReadGroupCovariateIndex();

    BQSRTableConfig bqsrTableConfig(header->getReadGroupCovariates().size(), bqsr::defaults::MAX_QUALITY_SCORE, bqsr::defaults::MAX_CYCLE, bqsr::defaults::CONTEXT_LEN);
    std::shared_ptr<BaseRecalibratorTables> tables;

    std::unique_ptr<MappedFasta> referenceGenome;
    std::unique_ptr<KnownSites>  knownSites;
    std::shared_ptr<BaseRecalibrator> bqsrRecal;

    std::vector<std::unique_ptr<IExecutor>> executors;

    if (args.is_recalibrate()) {

        referenceGenome = std::make_unique<MappedFasta>(args.reference, header);
        knownSites = std::make_unique<KnownSites>(args.known_sites, header);

        executors.reserve(selected_devices.size());

        for (int dev_idx : selected_devices) {
            sycl::device dev = gpus[dev_idx];
            executors.emplace_back(
                std::make_unique<BQSRRecalExecutor>(
                    dev, bqsrTableConfig, *referenceGenome, *knownSites, *header
                )
            );
        }

        /*tables = std::make_shared<BaseRecalibratorTables>(bqsrTableConfig, queue);
        bqsrRecal = std::make_shared<BaseRecalibrator>(*referenceGenome, *knownSites, queue, *header, *tables);
        bqsrRecal->initializeDeviceMemory();
*/
        // Stage 1: Recalibrate (no output here)
        stages.emplace_back("BQSR-Recalibrate", [&executors](SamAlignments& batch, size_t wid){
            //recal->Recalibrate(batch);
            executors[wid]->Process(batch);
        });
        // Finalize: Print tables at the end of all batch processing
        finals.emplace_back([&] {
            for (int i = 0; i < executors.size(); i++) {
                executors[i]->Finalize();
            }
            auto* bqsrExecMain = dynamic_cast<BQSRRecalExecutor*>(executors[0].get());
            if (!bqsrExecMain) throw std::runtime_error("executor[0] is not BQSRRecalExecutor");

            for (int i = 1; i < executors.size(); i++) {
                auto* bqsrExec = dynamic_cast<BQSRRecalExecutor*>(executors[i].get());
                if (!bqsrExec) throw std::runtime_error("executor["+std::to_string(i)+"] is not BQSRRecalExecutor");
                bqsrExecMain->Tables().merge(bqsrExec->Tables());
            }
            //recal->exportTables(*table);

            bqsrExecMain->Tables().FinalizeBQSRTables();
            bqsrExecMain->Tables().PrintBQSRTables(args.recal_table, *header);
        });
    }

    int quantize_level = 0;
    std::vector<uint8_t> sqqList;

    /*if (args.is_apply()) {
        tables = std::make_unique<BaseRecalibratorTables>(bqsrTableConfig, queue);
        tables->initializeDeviceMemory();
        tables->LoadFromReport(args.recal_table, *header);
        alignmentWriter = AlignmentWriterFactory::open(args.output_bam);

        //bqsrRecal = std::make_unique<BaseRecalibrator>(*referenceGenome, *knownSites, device, *header, *tables);

        // Stage 1: Recalibrate (no output here)
        stages.emplace_back("BQSR-Apply", [table = tables, quantLevel=quantize_level, sqqCopy=sqqList](SamAlignments& batch){
            std::vector<uint8_t> sqqCopyLocal = sqqCopy;
            table->ApplyBQSR(batch, quantLevel, sqqCopyLocal);
        });

        // Stage 2: Write out batch on BAM file
        stages.emplace_back("Write as BAM", [writer = alignmentWriter](SamAlignments& batch){
            writer->WriteAlignmentBatch(batch);
        });

        finals.emplace_back([writer = alignmentWriter]{
            try {
                writer->Close();        // close file/FD
            } catch (...) {
                throw std::runtime_error("Error while closing output file");
            }
        });
    }*/

    // Build and run engine
    StreamingEngine::Config cfg;
    cfg.channel_capacity = 4;
    cfg.workers = executors.size();

    if (executors.empty()) {
        throw std::runtime_error("No executors created (workers would be 0)");
    }

    Timer t;
    t.reset();
    StreamingEngine engine(cfg, (std::move(stages)), (std::move(finals)));

    engine.run_streaming(alignmentReader, executors, header);
    //t.elapsed();
    Logger::getLogger().log(LOG_INFO, "Time to run all stages: " + std::to_string(t.elapsed()) + " sec");
    ProfilerIO::instance().report();
    return 0;
}

