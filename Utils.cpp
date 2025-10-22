#include "Utils.h"

#include <iomanip>

#include "logger/Logger.h"

using namespace std;
namespace po = boost::program_options;

void Fields(std::string &s, std::vector<std::string> &result) {
    size_t pos = 0;
    std::string token;
    while ((pos = s.find(' ')) != std::string::npos) {
        token = s.substr(0, pos);
        if(!token.empty()){
            //std::cout << token << std::endl;
            result.emplace_back(token);
        }
        s.erase(0, pos + 1);
    }
    result.emplace_back(s);
}


void init_options( po::options_description &description, bool &help ){

    description.add_options()( "help,h", po::bool_switch(&help), "Display help" )(
            "command", po::value<string>(), "Command: recalibrate or apply (positional)" )(
            "bam-in", po::value<string>(), "Input BAM file" )(
            "bam-out", po::value<string>(), "Output BAM file (apply command only)" )(
            "recal-out", po::value<string>(), "Recalibration table file" )(
            "recal", po::value<string>(), "Recalibration table file" )(
            "ref", po::value<string>(), "Reference FASTA file (recalibrate command only)" )(
            "known-sites", po::value<string>(), "Comma-delimited list of known variant sites (recalibrate command only)"
    );

}

std::pair<uint64_t, int> Uvarint(const std::vector<uint8_t>& buf) {
    uint64_t x = 0;
    int s = 0;
    for (int i = 0; i < buf.size(); ++i) {
        if (i == 9) {
// Catch byte reads past MaxVarintLen64.
// See issue https://golang.org/issues/41185
            return {0, -(i + 1)}; // overflow
        }
        if (buf[i] < 0x80) {
            if (i == 8 && buf[i] > 1) {
                return {0, -(i + 1)}; // overflow
            }
            return {x | uint64_t(buf[i]) << s, i + 1};
        }
        x |= uint64_t(buf[i] & 0x7f) << s;
        s += 7;
    }
    return {0, 0};
}

std::pair<int64_t, int> Varint(const std::vector<uint8_t>& buf) {
    auto [ux, n] = Uvarint(buf);
    int64_t x = int64_t(ux >> 1);
    if ((ux & 1) != 0) {
        x = ~x;
    }
    return {x, n};
}

void print_stage_header(const std::string& stage_name) {
    const int totalWidth = 74;  // Total width for the string inside the || (adjust as needed)

    assert(totalWidth > stage_name.size());

    auto stageLength = static_cast<int32_t>(stage_name.size());

    // Calculate padding for left and right
    int32_t leftPadding = (totalWidth - stageLength) / 2;
    int32_t rightPadding = totalWidth - stageLength - leftPadding;

    // Output the centered string using setw and manipulators
    std::stringstream sstr;
    sstr << "||"
              << std::setw(leftPadding + 1) << std::left << " "
              << stage_name
              << std::setw(rightPadding) << std::right << " "
              << "||";

    Logger::getLogger().log(LOG_INFO, "------------------------------------------------------------------------------");
    Logger::getLogger().log(LOG_INFO, "||                     SYCL Genomics Acceleration Library                    ||");
    Logger::getLogger().log(LOG_INFO, "||                              Version 0.0.0.1.                             ||");
    Logger::getLogger().log(LOG_INFO, sstr.str());
    Logger::getLogger().log(LOG_INFO, "-------------------------------------------------------------------------------");
}


void print_recalibrate_parameter() {
    std::cerr << "  --bam-in <file>       Input BAM file with aligned reads\n";
    std::cerr << "  --recal-out <file>    Output recalibration table file\n";
    std::cerr << "  --ref <file>          Reference genome FASTA file\n";
    std::cerr << "  --known-sites <file>  Known variants VCF file (comma-delimited for multiple)\n\n";
}

void print_apply_parameter() {
    std::cerr << "  --bam-in <file>       Input BAM file with original quality scores\n";
    std::cerr << "  --bam-out <file>      Output BAM file with recalibrated quality scores\n";
    std::cerr << "  --recal <file>    Recalibration table file (from recalibrate step)\n\n";
}


Arguments parse_arguments(int argc, char** argv) {
    bool help = false;
    po::options_description description("bqsr <command> [options]");
    po::positional_options_description positional;
    positional.add("command", 1);
    init_options(description, help);

    po::command_line_parser parser{argc, argv};
    parser.options(description).positional(positional);
    auto parsed_result = parser.run();
    po::variables_map vm;
    po::store(parsed_result, vm);
    po::notify(vm);

    // Parse command first
    string command;
    bool recalibrate = false;
    bool apply = false;

    if (vm.count("command")) {
        command = vm["command"].as<string>();
        if (command == "recalibrate") {
            recalibrate = true;
        } else if (command == "apply") {
            apply = true;
        } else if (command != "help" && command != "-h" && command != "--help") {
            std::cerr << "Unknown command: " << command << std::endl;
            std::cerr << "Must be 'recalibrate' or 'apply'" << std::endl;
            exit(-1);
        }
    } else {
        if (help) {
            std::cerr << "bqsr - Base Quality Score Recalibration\n\n";
            std::cerr << "Usage:\n";
            std::cerr << "  bqsr recalibrate [options]   Generate recalibration tables\n";
            std::cerr << "  bqsr apply [options]         Apply recalibration to reads\n";
            std::cerr << "  bqsr -h, --help              Display this help message\n\n";
            std::cerr << "================================================================================\n\n";
            std::cerr << "RECALIBRATE\n";
            std::cerr << "-----------\n";
            std::cerr << "bqsr recalibrate [options]\n\n";
            std::cerr << "Generate recalibration tables from aligned reads\n\n";
            std::cerr << "Required options:\n";
            print_recalibrate_parameter();
            std::cerr << "================================================================================\n\n";
            std::cerr << "APPLY\n";
            std::cerr << "-----\n";
            std::cerr << "bqsr apply [options]\n\n";
            std::cerr << "Apply recalibration tables to adjust base quality scores\n\n";
            std::cerr << "Required options:\n";
            print_apply_parameter();
            exit(0);
        }
        std::cerr << "Error: Command not specified\n\n";
        std::cerr << "bqsr - Base Quality Score Recalibration\n\n";
        std::cerr << "Usage:\n";
        std::cerr << "  bqsr recalibrate [options]   Generate recalibration tables\n";
        std::cerr << "  bqsr apply [options]         Apply recalibration to reads\n";
        std::cerr << "  bqsr -h, --help              Display this help message\n\n";
        std::cerr << "For command-specific help, use: bqsr <command> -h\n";
        exit(-1);
    }

    // Show command-specific help
    if (help) {
        std::cerr << "bqsr " << command << " [options]\n\n";
        if (recalibrate) {
            std::cerr << "Generate recalibration tables from aligned reads\n\n";
            std::cerr << "Required options:\n";
            print_recalibrate_parameter();
        } else if (apply) {
            std::cerr << "Apply recalibration tables to adjust base quality scores\n\n";
            std::cerr << "Required options:\n";
            print_apply_parameter();
        }
        exit(0);
    }

    Arguments args;
    args.command = command;

    // Parse common parameters
    if (vm.count("bam-in")) {
        args.input_bam = vm["bam-in"].as<string>();
    } else {
        std::cerr << "Error: Input BAM file not specified (--bam-in)\n" << std::endl;
        if (recalibrate) {
            std::cerr << "bqsr recalibrate [options]\n\n";
            std::cerr << "Generate recalibration tables from aligned reads\n\n";
            std::cerr << "Required options:\n";
            print_recalibrate_parameter();
        } else if (apply) {
            std::cerr << "bqsr apply [options]\n\n";
            std::cerr << "Apply recalibration tables to adjust base quality scores\n\n";
            std::cerr << "Required options:\n";
            print_apply_parameter();
        }
        exit(-1);
    }

    // Parse recalibrate-specific parameters
    if (recalibrate) {
        if (vm.count("recal-out")) {
            args.recal_table = vm["recal-out"].as<string>();
        } else {
            std::cerr << "Error: Output recalibration table not specified (--recal-out)\n" << std::endl;
            std::cerr << "bqsr recalibrate [options]\n\n";
            std::cerr << "Required options:\n";
            print_recalibrate_parameter();
            exit(-1);
        }

        if (vm.count("ref")) {
            args.reference = vm["ref"].as<string>();
        } else {
            std::cerr << "Error: Reference file not specified (--ref)\n" << std::endl;
            std::cerr << "bqsr recalibrate [options]\n\n";
            std::cerr << "Required options:\n";
            print_recalibrate_parameter();
            exit(-1);
        }
        // Parse known sites
        if (vm.count("known-sites")) {
            std::string sites = vm["known-sites"].as<std::string>();
            std::string delimiter = ",";
            size_t pos = 0;
            std::string token;
            while ((pos = sites.find(delimiter)) != std::string::npos) {
                token = sites.substr(0, pos);
                args.known_sites.emplace_back(token);
                sites.erase(0, pos + delimiter.length());
            }
            if (!sites.empty()) {
                args.known_sites.emplace_back(sites);
            }
        } else {
            std::cerr << "Error: VCF file not specified (--known-sites)\n" << std::endl;
            std::cerr << "bqsr recalibrate [options]\n\n";
            std::cerr << "Required options:\n";
            print_recalibrate_parameter();
            exit(-1);
        }
    }
    // Parse apply-specific parameters
    else if (apply) {
        if (vm.count("bam-out")) {
            args.output_bam = vm["bam-out"].as<string>();
        } else {
            std::cerr << "Error: Output BAM file not specified (--bam-out)\n" << std::endl;
            std::cerr << "bqsr apply [options]\n\n";
            std::cerr << "Required options:\n";
            print_apply_parameter();
            exit(-1);
        }

        if (vm.count("recal")) {
            args.recal_table = vm["recal"].as<string>();
        } else {
            std::cerr << "Error: Recalibration table not specified (--recal)\n" << std::endl;
            std::cerr << "bqsr apply [options]\n\n";
            std::cerr << "Required options:\n";
            print_apply_parameter();
            exit(-1);
        }
    }

    return args;
}

