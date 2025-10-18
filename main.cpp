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
    Timer t;
    bool help{};
    namespace po = boost::program_options;
    po::options_description description("GenomicPipeline [options]" );
    init_options(description, help);

    po::command_line_parser parser{ argc, argv };
    parser.options( description );
    auto parsed_result = parser.run();
    po::variables_map vm;
    po::store( parsed_result, vm );
    po::notify( vm );

    if (help) {
        std::cerr << description << std::endl;
        return 0;
    }
    //bool filter_unmapped_reads = false;
    //bool filter_unmapped_reads_strict = false;
    //bool filter_mapping_quality = false;
    //bool filter_non_exact_mapping_reads = false;
    //bool filter_non_overlapping_reads = false;
    //bool replace_read_group = false;
    bool mark_duplicates = false;
    //bool mark_optical_duplicates = false;
    //bool optical_duplicates_pixel_distance = false;
    //bool remove_optional_fields = false;
    //bool keep_optional_fields = false;
    bool sorting = false;

    bool bqsr = false;
    bool applyBqsr = false;

    string input_filename;
    basic_string<char> output_filename;
    string metric_filename;
    string reference_filename;
    string recal_filename;
    string bed_filename;

    /* Input parsing variable section */
    std::vector<SamAlignments> alignments;

    std::vector<std::string> knownSites;

    string echo_command="gal ";

    //SortingOrder sorting_order = SortingOrder::Keep;

    if( vm.count("input") ){
        input_filename  = vm["input"].as<string>();
        echo_command+=("--input " + input_filename + " ");
    } else{
        std::cerr << "Input file not specified" << std::endl;
        std::cerr << description << std::endl;
        exit(-1);
    }

    if( vm.count("output") ){
        output_filename  = vm["output"].as<string>();
        echo_command+=("--output " + output_filename + " ");
    } else{
        std::cerr << "Output file not specified" << std::endl;
        std::cerr << description << std::endl;
        exit(-1);
    }

    if( vm.count("bqsr") ){
        if ( vm.count("reference") ) {
            recal_filename = vm["bqsr"].as<string>();
            echo_command += ("--bqsr " + recal_filename + " ");
            reference_filename = vm["reference"].as<string>();
            echo_command += ("--reference " + reference_filename + " ");
            bqsr = true;
        } else {
            std::cerr << "Reference file not specified for bqsr step." << std::endl;
            std::cerr << description << std::endl;
            exit(-1);
        }
    }

    if( vm.count("applyBqsr") ){
        if ( vm.count("reference") ) {
            recal_filename = vm["applyBqsr"].as<string>();
            echo_command += ("--bqsr " + recal_filename + " ");
            reference_filename = vm["reference"].as<string>();
            echo_command += ("--reference " + reference_filename + " ");
            applyBqsr = true;
        } else {
            std::cerr << "Reference file not specified for applyBqsr step." << std::endl;
            std::cerr << description << std::endl;
            exit(-1);
        }
    }

    if (vm.count("known-sites")) {
        std::string sites = vm["known-sites"].as<std::string>();
        echo_command += ("--known-sites ");
        std::string delimiter = ",";
        size_t pos = 0;
        std::string token;
        while ((pos = sites.find(delimiter)) != std::string::npos) {
            token = sites.substr(0, pos);
            knownSites.emplace_back(token);
            echo_command += (token + ",");
            sites.erase(0, pos + delimiter.length());
        }
        if (!sites.empty()) {
            echo_command += (sites);
            knownSites.emplace_back(sites);
        } else {
            echo_command[echo_command.size()-1] = ' ';
        }
    }
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

    std::cout << "Found " << gpus.size() << " GPU devices:\n";
    for (size_t i = 0; i < gpus.size(); ++i) {
        std::cout << "  [" << i << "] "
                  << gpus[i].get_info<sycl::info::device::name>() << " ("
                  << gpus[i].get_info<sycl::info::device::vendor>() << ")\n";
    }

    // Select the second GPU (index 1)
    size_t selected_idx = 0;  // zero-based index
    if (selected_idx >= gpus.size()) {
        std::cerr << "Requested GPU index " << selected_idx
                  << " but only " << gpus.size() << " GPUs available.\n";
        return 1;
    }

    sycl::device device = gpus[selected_idx];

    sycl::queue q(device);


    /**** Input parsing for reading Alignment Data *****/

    Logger::getLogger().set_min_severity(LOG_DEBUG);
    Logger::getLogger().log(LOG_INFO, "Running: " + echo_command);


    std::shared_ptr<Header> header;
    std::unique_ptr<AlignmentReader> alignmentReader = AlignmentReaderFactory::open(input_filename);
    header = std::make_shared<Header>();
    alignmentReader->ParseHeader(*header);


    std::thread loader_sam([&]() {
        t.reset();
        print_stage_header("Reading SAM File");
        size_t nAlns = alignmentReader->ParseAlignment(alignments, header);
        //assert(nAlns == alignments.GetAlnsNum());

        Logger::getLogger().log(LOG_INFO, "Time to import all alignments: " + std::to_string(t.elapsed()) + " sec");
    });

    std::unique_ptr<BaseRecalibrator> baseRecal;
    std::unique_ptr<BaseRecalibratorTables> baseRecalibratorTables;

    size_t maxCycle = 500;

    if (bqsr) {
        std::thread loader_fasta([&]() {
            baseRecalibratorTables = std::make_unique<BaseRecalibratorTables>(header->getReadGroupRecord().size(), 94,
                                                             maxCycle);
            MappedFasta referenceGenome;
            referenceGenome.Open(reference_filename);
            baseRecal = std::make_unique<BaseRecalibrator>(referenceGenome, knownSites, device, *header, *baseRecalibratorTables);
        });
        loader_fasta.join();
    }

    loader_sam.join();

    for (size_t batchIdx = 0; batchIdx < alignments.size(); ++batchIdx) {
        auto addREFID = Filter::AddREFID(alignments[batchIdx]);
        for (size_t idx = 0; idx < alignments[batchIdx].GetAlnsNum(); idx++) {
            addREFID(alignments[batchIdx], idx);
        }
    }

    if (bqsr) {
        print_stage_header("BQSR I - Recalibrate");
        t.reset();

        Logger::getLogger().log(LOG_INFO, "Time to read reference genome: " + std::to_string(t.elapsed()) + " sec");

        t.reset();

        for (size_t i = 0; i < alignments.size(); i++) {
            baseRecal->Recalibrate(device, alignments[i], maxCycle);
        }
        baseRecal->getTables(*baseRecalibratorTables);
        baseRecalibratorTables->FinalizeBQSRTables();

        Logger::getLogger().log(LOG_INFO, "Time to compute BQSR tables: " + std::to_string(t.elapsed()) + " sec");

        baseRecalibratorTables->PrintBQSRTables(recal_filename, alignments[0].TAGS.RG, alignments[0].GetHdr());
    }

    if (applyBqsr) {
        print_stage_header("BQSR II - Apply Recalibration");
        t.reset();

        int quantize_level = 0;
        std::vector<uint8_t> sqqList;

        BaseRecalibratorTables baseRecalibratorTables(alignments[0].GetHdr().getReadGroupRecord().size(), 94,
                                                      maxCycle);//TODO: Make the size as function of context size

        baseRecalibratorTables.ApplyBQSR(device, alignments[0], quantize_level, sqqList, maxCycle);


        Logger::getLogger().log(LOG_INFO, "Time to apply BQSR: " + std::to_string(t.elapsed()) + " sec");


        print_stage_header("Writing SAM File");

        std::unique_ptr<AlignmentWriter> wf = AlignmentWriterFactory::open(output_filename);

        Logger::getLogger().log(LOG_INFO, "Start writing output file");

        t.reset();

        wf->FormatHeader(alignments[0].GetHdr());

        wf->FormatAlignment(alignments);

        Logger::getLogger().log(LOG_INFO, "Time to Write reads on disk: " + std::to_string(t.elapsed()) + " sec");
    }
    return 0;
}

