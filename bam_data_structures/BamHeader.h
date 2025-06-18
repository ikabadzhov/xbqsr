//
// Created by loren on 14/11/2023.
//

#ifndef GENOMICPIPELINE_BAMHEADER_H
#define GENOMICPIPELINE_BAMHEADER_H

#include <string>
#include <vector>
#include "../sam_data_structures/sam_fields/Cigar.h"

/* ENUM FOR INDIANESS */

enum class Endianness_old {
    LittleEndian,
    BigEndian,
    Unknown
};

/* STRUCTS */

struct ReferenceInformation_old {
    uint32_t length_name;
    std::string name;
    uint32_t length_sequence;
};

struct ReferenceSequence_old {
    std::string name;
    size_t length;
    ReferenceInformation_old referenceInformation;
};

struct HeaderRecord_old {
    uint32_t header_length;
    std::vector<ReferenceSequence_old> references;
    int n_references;
};


class BAMHeaderReader {
private:
    static void printField(const std::string& fieldName, const std::string& value);
    static std::vector<ReferenceSequence_old> parseHeader(const std::string& headerData);
public:
    static Endianness_old checkEndiannessOfHostMachine();
    static void readHeader(const std::string& bamFilePath);
    static void compressToBGZF(const std::string &inputFilename, const std::string &outputFilename);
    static void decompressBGZF(const std::string &inputFilename, const std::string &outputFilename);
};


#endif //GENOMICPIPELINE_BAMHEADER_H
