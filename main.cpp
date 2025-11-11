#include <iostream>
#include <numeric>
#include <algorithm>
#include <chrono>
#include <tbb/parallel_sort.h>
#include "Utils.h"
#include "input_output/alignment_io/AlignmentIO.h"
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
#include "logger/Logger.h"
#include "sycl_parallel_library/intervalsmove/intervalmove.h"
#include "sycl_parallel_library/scan/scan.h"
#include "sycl_parallel_library/segmented_reduce/segreduce.h"
#include "sycl_parallel_library/segmented_scan/segscan.h"

using namespace std;
using namespace chrono;

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

    // Validate requested number of GPUs
    int n_gpus = args.gpus;
    if (n_gpus == -1) {
        // -1 means use all available GPUs
        n_gpus = gpus.size();
    } else if (n_gpus > static_cast<int>(gpus.size())) {
        Logger::getLogger().log(LOG_FATAL, "Requested ", n_gpus, " GPUs but only ", gpus.size(),
                  " GPU(s) available.");
        return 1;
    }

    //gpus.resize(n_gpus);

    Logger::getLogger().log(LOG_INFO, "Devices availables:");
    for (size_t i = 0; i < gpus.size(); ++i) {
        Logger::getLogger().log(LOG_INFO, "  [", i, "] ",
            gpus[i].get_info<sycl::info::device::name>(), " (",
            gpus[i].get_info<sycl::info::device::vendor>(), ")");
    }

    int selected_dev = 3;

    Logger::getLogger().log(LOG_INFO, "Selected device:");
    Logger::getLogger().log(LOG_INFO, "  [", selected_dev, "] ",
            gpus[selected_dev].get_info<sycl::info::device::name>(), " (",
            gpus[selected_dev].get_info<sycl::info::device::vendor>(), ")");

    gpus[0] = gpus[selected_dev];

    sycl::device device = gpus[0];

    sycl::queue queue(device);
    /**** Input parsing for reading Alignment Data *****/

    if (args.verbose) {
        Logger::getLogger().set_min_severity(LOG_DEBUG);
    } else {
        Logger::getLogger().set_min_severity(LOG_INFO);
    }

    std::vector<StreamingEngine::Stage> stages;
    std::vector<StreamingEngine::FinalizeFn> finals;

    std::shared_ptr<AlignmentReader> alignmentReader = AlignmentReaderFactory::open(args.input_bam);
    std::shared_ptr<AlignmentWriter> alignmentWriter;

    std::shared_ptr<Header> header = std::make_shared<Header>();
    alignmentReader->ParseHeader(*header);
    header->buildReadGroupCovariateIndex();

    BQSRTableConfig bqsrTableConfig(header->getReadGroupCovariates().size(), bqsr::defaults::MAX_QUALITY_SCORE, bqsr::defaults::MAX_CYCLE, bqsr::defaults::CONTEXT_LEN);
    std::shared_ptr<BaseRecalibratorTables> tables;

    std::unique_ptr<MappedFasta> referenceGenome;
    std::unique_ptr<KnownSites>  knownSites;
    std::shared_ptr<BaseRecalibrator> bqsrRecal;

    // Stage 0: AddREFID filter
    /*stages.emplace_back("AddREFID", [](SamAlignments& batch){
        auto addREFID = Filter::AddREFID(batch);
        for (size_t i = 0; i < batch.GetAlnsNum(); ++i) addREFID(batch, i);
    });*/

    if (args.is_recalibrate()) {
        referenceGenome = std::make_unique<MappedFasta>(args.reference, header);
        knownSites = std::make_unique<KnownSites>(args.known_sites, header);
        tables = std::make_shared<BaseRecalibratorTables>(bqsrTableConfig, queue);
        bqsrRecal = std::make_shared<BaseRecalibrator>(*referenceGenome, *knownSites, queue, *header, *tables);

        // Stage 1: Recalibrate (no output here)
        stages.emplace_back("BQSR-Recalibrate", [recal = bqsrRecal](SamAlignments& batch){
            recal->Recalibrate(batch);
        });
        // Finalize: Print tables at the end of all batch processing
        finals.emplace_back([&, recal = bqsrRecal, table = tables]{
            recal->exportTables(*table);
            table->FinalizeBQSRTables();
            table->PrintBQSRTables(args.recal_table, *header);
        });
    }

    int quantize_level = 0;
    std::vector<uint8_t> sqqList;

    if (args.is_apply()) {
        tables = std::make_unique<BaseRecalibratorTables>(bqsrTableConfig, queue);
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
    }

    // Build and run engine
    StreamingEngine::Config cfg;
    cfg.channel_capacity = 4;
    cfg.workers = 1;
    Timer t;
    t.reset();
    StreamingEngine engine(cfg, (std::move(stages)), (std::move(finals)));

    engine.run_streaming(*alignmentReader, header);
    //t.elapsed();
    Logger::getLogger().log(LOG_INFO, "Time to run all stages: " + std::to_string(t.elapsed()) + " sec");

    return 0;
}

