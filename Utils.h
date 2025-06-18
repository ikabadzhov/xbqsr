#ifndef GENOMICPIPELINE_UTILS_H
#define GENOMICPIPELINE_UTILS_H


#include <string>
#include <vector>
#include <boost/program_options.hpp>
#include <tbb/parallel_for.h>


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

std::pair<uint64_t, int> Uvarint(const std::vector<uint8_t>& buf);

std::pair<int64_t, int> Varint(const std::vector<uint8_t>& buf);

void print_stage_header(const std::string& stage_name);

#endif //GENOMICPIPELINE_UTILS_H
