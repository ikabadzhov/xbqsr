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

    gpus.resize(n_gpus);

    Logger::getLogger().log(LOG_INFO, "Running on ", gpus.size(), " GPU devices:");
    for (size_t i = 0; i < gpus.size(); ++i) {
        Logger::getLogger().log(LOG_INFO, "  [", i, "] ",
            gpus[i].get_info<sycl::info::device::name>(), " (",
            gpus[i].get_info<sycl::info::device::vendor>(), ")");
    }

    // Select the first GPU (index 0)

    sycl::device device = gpus[0];
    /**** Input parsing for reading Alignment Data *****/

    Logger::getLogger().set_min_severity(LOG_DEBUG);

    std::shared_ptr<Header> header;
    std::unique_ptr<AlignmentReader> alignmentReader = AlignmentReaderFactory::open(args.input_bam);
    header = std::make_shared<Header>();
    alignmentReader->ParseHeader(*header);


    std::thread loader_sam([&]() {
        Timer t_in_sam;
        t_in_sam.reset();
        print_stage_header("Reading SAM File");
        alignmentReader->ParseAlignment(alignments, header);

        Logger::getLogger().log(LOG_INFO, "Time to import all alignments: " + std::to_string(t_in_sam.elapsed()) + " sec");
    });


    std::vector<std::unique_ptr<BaseRecalibratorTables>> localBaseRecalibratorTables(gpus.size());
    std::unique_ptr<BaseRecalibratorTables> globalBaseRecalibratorTables;

    BQSRTableConfig bqsrTableConfig(header->getReadGroupRecord().size(), bqsr::defaults::MAX_QUALITY_SCORE, bqsr::defaults::MAX_CYCLE, bqsr::defaults::CONTEXT_LEN);

    if (args.is_recalibrate()) {
        std::vector<std::unique_ptr<BaseRecalibrator>> baseRecals(gpus.size());

        std::thread loader_fasta([&]() {
            Timer t_in_ref;
            t_in_ref.reset();
            globalBaseRecalibratorTables = std::make_unique<BaseRecalibratorTables>(bqsrTableConfig);
            MappedFasta referenceGenome(args.reference, header);
            KnownSites known_sites(args.known_sites, header);

            tbb::parallel_for(0UL, gpus.size(), [&](size_t dev_idx) {
                localBaseRecalibratorTables[dev_idx] = std::make_unique<BaseRecalibratorTables>(bqsrTableConfig);
                baseRecals[dev_idx] = std::make_unique<BaseRecalibrator>(referenceGenome, known_sites, gpus[dev_idx], *header, *localBaseRecalibratorTables[dev_idx]);
            });
            Logger::getLogger().log(LOG_INFO, "Time to read reference genome: " + std::to_string(t_in_ref.elapsed()) + " sec");
        });
        loader_fasta.join();

        loader_sam.join();

        for (size_t i = 0; i < alignments.size(); i++) {
            auto addREFID = Filter::AddREFID(alignments[i]);
            for (size_t idx = 0; idx < alignments[i].GetAlnsNum(); idx++) {
                addREFID(alignments[i], idx);
            }
        }

        print_stage_header("BQSR - Recalibrate");

        Timer t_recal;
        t_recal.reset();

        size_t workload = (alignments.size() + n_gpus - 1) / n_gpus;

        tbb::parallel_for(0, n_gpus, [&](size_t dev_idx) {
            for (size_t i = workload * dev_idx; i < workload * dev_idx + workload; i++) {
                if (i >= alignments.size()) {
                    return;
                }
                baseRecals[dev_idx]->Recalibrate(gpus[dev_idx], alignments[i]);
            }
        });
        for (size_t dev_idx = 0; dev_idx < gpus.size(); dev_idx++) {
            baseRecals[dev_idx]->exportTables(*localBaseRecalibratorTables[dev_idx]);
            globalBaseRecalibratorTables->merge(*localBaseRecalibratorTables[dev_idx]);
        }

        globalBaseRecalibratorTables->FinalizeBQSRTables();

        Logger::getLogger().log(LOG_INFO, "Time to compute BQSR tables: " + std::to_string(t_recal.elapsed()) + " sec");

        globalBaseRecalibratorTables->PrintBQSRTables(args.recal_table, alignments[0].TAGS.RG, *header);
    }

    if (args.is_apply()) {
        Timer t_apply;
        t_apply.reset();

        tbb::parallel_for(0, n_gpus, [&](size_t dev_idx) {
            localBaseRecalibratorTables[dev_idx] = std::make_unique<BaseRecalibratorTables>(bqsrTableConfig);
            localBaseRecalibratorTables[dev_idx]->LoadFromReport(args.recal_table, *header);
        });

        Logger::getLogger().log(LOG_INFO, "Time to load BQSR table from file: " + std::to_string(t_apply.elapsed()) + " sec");

        loader_sam.join();

        print_stage_header("BQSR - Apply");

        for (size_t batchIdx = 0; batchIdx < alignments.size(); ++batchIdx) {
            auto addREFID = Filter::AddREFID(alignments[batchIdx]);
            for (size_t idx = 0; idx < alignments[batchIdx].GetAlnsNum(); idx++) {
                addREFID(alignments[batchIdx], idx);
            }
        }

        t_apply.reset();

        int quantize_level = 0;
        std::vector<uint8_t> sqqList;

        size_t workload = (alignments.size() + n_gpus - 1) / n_gpus;

        tbb::parallel_for(0, n_gpus, [&](size_t dev_idx) {
            for (size_t i = workload * dev_idx; i < workload * dev_idx + workload; i++) {
                if (i >= alignments.size()) {
                    return;
                }
                localBaseRecalibratorTables[dev_idx]->ApplyBQSR(gpus[dev_idx], alignments[i], quantize_level, sqqList);
            }
        });

        Logger::getLogger().log(LOG_INFO, "Time to apply BQSR: " + std::to_string(t_apply.elapsed()) + " sec");

        Timer t_out;
        t_out.reset();
        print_stage_header("Writing SAM File");

        std::unique_ptr<AlignmentWriter> wf = AlignmentWriterFactory::open(args.output_bam);

        wf->FormatHeader(*header);
        wf->FormatAlignment(alignments);

        Logger::getLogger().log(LOG_INFO, "Time to Write reads on disk: " + std::to_string(t_out.elapsed()) + " sec");
    }
    return 0;
}

