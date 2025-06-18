#ifndef GENOMICPIPELINE_MAPPEDFASTA_H
#define GENOMICPIPELINE_MAPPEDFASTA_H

#include <string>
#include <vector>
#include <tbb/concurrent_unordered_map.h>

extern std::vector<uint8_t> ElfastaMagic;// = {0x31, 0xFA, 0x57, 0xA1};

class MappedFasta {
    tbb::concurrent_unordered_map<std::string, std::vector<uint8_t>> fasta;

public:
    std::vector<uint8_t>& Seq(std::string_view contig);
    tbb::concurrent_unordered_map<std::string, std::vector<uint8_t>>& getFullFasta() {
        return fasta;
    }
    void Open(const std::string &filename, bool toUpper = false, bool toN = false);
};


#endif //GENOMICPIPELINE_MAPPEDFASTA_H
