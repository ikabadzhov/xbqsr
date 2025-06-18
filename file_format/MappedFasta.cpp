#include "MappedFasta.h"
#include "../Utils.h"
#include <string_view>
#include <iostream>
#include <fstream>



constexpr size_t ASCII_SIZE = 256;

//const int MaxVarintLen64 = 10;

std::vector<uint8_t>& MappedFasta::Seq(std::string_view contig) {
    std::string s(contig);
    return this->fasta[s];
}

std::vector<char> iupacTable(ASCII_SIZE, 0);

// Initialize the IUPAC upper case normalization vector
std::vector<char> iupacUpperTable(ASCII_SIZE, 0);

void initializeIupacTables() {
    // Initialize the iupacTable for normal cases
    for (char c = 0; c < ASCII_SIZE; ++c) {
        if (c == 'A' || c == 'C' || c == 'G' || c == 'T' || c == 'a' || c == 'c' || c == 'g' || c == 't') {
            iupacTable[c] = c;
        } else if (std::isalpha(c)) {
            iupacTable[c] = 'N';
        }
    }

    // Initialize the iupacUpperTable for upper cases
    for (char c = 0; c < ASCII_SIZE; ++c) {
        if (c == 'A' || c == 'C' || c == 'G' || c == 'T') {
            iupacUpperTable[c] = c;
            iupacUpperTable[std::tolower(c)] = c;  // Set lower case to upper case
        } else if (std::isalpha(c)) {
            iupacUpperTable[c] = 'N';
        }
    }
}

char ToN(char base) {
    if (base < ASCII_SIZE) {
        return iupacTable[base];
    }
    return base;  // Return the base itself if it's out of range
}

char ToUpperAndN(char base) {
    if (base < ASCII_SIZE) {
        return iupacUpperTable[base];
    }
    return base;  // Return the base itself if it's out of range
}

void MappedFasta::Open(const std::string& filename, bool toUpper, bool toN) {
    std::ifstream file(filename, std::ios::binary);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open file " + filename);
    }

    std::string line;
    std::string contig;
    std::vector<unsigned char> seq;

    while (std::getline(file, line)) {
        if (line.empty())
            continue;

        if (line[0] == '>') {
            if (!contig.empty()) {
                fasta[contig] = seq;
                seq.clear();
            }
            std::stringstream ss(line);
            //auto pos = line.find_first_of('\t'); // find occurrences of multiple spaces of \t
            ss >> contig;//line.substr(1, pos);  // Extract contig name from header
            contig = contig.substr(1);
        } else {
            std::vector<unsigned char> lineData(line.begin(), line.end());
            if (toUpper) {
                std::transform(lineData.begin(), lineData.end(), lineData.begin(),
                               [](unsigned char c) { return std::toupper(c); });
            }
            if (toN) {
                // Example transformation, assuming iupacTable is a map
                for (int32_t j = 0; j < lineData[j]; j++) {
                    lineData[j] = ToN(lineData[j]);
                }
            }
            seq.insert(seq.end(), lineData.begin(), lineData.end());
        }
    }
    if (!contig.empty()) {
        fasta[contig] = seq;
    }

    if (!file.eof()) {
        throw std::runtime_error("Error reading the file " + filename);
    }
}