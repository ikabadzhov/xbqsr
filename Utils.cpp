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
            "input", po::value<string>(), "arg")(
            "output", po::value<string>(), "arg")(
            "mark-duplicates", po::bool_switch()->default_value(false))(
            "metric-file", po::value<string>(), "file" )(
            "optical-duplicates-pixel-distance", po::value<int32_t>(), "nr" )(
            "remove-duplicates", po::bool_switch()->default_value(false), "" )(
            "sorting-order", po::value<string>(), "[ keep | unknown | unsorted | queryname | coordinate ]" )(
            "bqsr", po::value<string>(), " recal-file" )(
                "known-sites", po::value<string>(), " site-file" )(
            "reference", po::value<string>(), "fasta-file" )(
        	"contig-group-size", po::value<int32_t>(), "nr (sfm only)"
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

