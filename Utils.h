#ifndef GENOMICPIPELINE_UTILS_H
#define GENOMICPIPELINE_UTILS_H


#include <string>
#include <vector>
#include <boost/program_options.hpp>
#include <tbb/parallel_for.h>


// Command-line arguments structure
struct Arguments {
    std::string command;            // "recalibrate" or "apply"
    std::string input_bam;          // --bam-in
    std::string output_bam;         // --bam-out (for apply)
    std::string recal_table;        // --recal-out
    std::string reference;          // --ref (for recalibrate)
    std::vector<std::string> known_sites;  // --known-sites

    bool is_recalibrate() const { return command == "recalibrate"; }
    bool is_apply() const { return command == "apply"; }
};


template<typename T>
void permute(std::vector<size_t> &indices, std::vector<T> &values) {
    if (values.empty()) {
        return;
    }
    std::vector<T> sorted_vec(indices.size());
    tbb::parallel_for (size_t(0), indices.size(), [&](size_t i){
        sorted_vec[i] = values[indices[i]];
    });
    values = std::move(sorted_vec);
}

void Fields(std::string& s, std::vector<std::string>& result);

void init_options( boost::program_options::options_description &description, bool &help);

// Parse command-line arguments and return Arguments struct
Arguments parse_arguments(int argc, char** argv);

std::pair<uint64_t, int> Uvarint(const std::vector<uint8_t>& buf);

std::pair<int64_t, int> Varint(const std::vector<uint8_t>& buf);

void print_stage_header(const std::string& stage_name);

#endif //GENOMICPIPELINE_UTILS_H
