#ifndef GENOMICPIPELINE_BAMHEADERFORMATTER_H
#define GENOMICPIPELINE_BAMHEADERFORMATTER_H


#include <string>
#include "../../sam_data_structures/SamHeader.h"
#include "../streams/BufferWriter.h"

class BamHeaderFormatter {
public:
    static void FormatBamHeader(Header &hdr, std::string& line);
    static void FormatBamHeaderFromSlice(Header &hdr, BufferWriter &buffer);
};


#endif //GENOMICPIPELINE_BAMHEADERFORMATTER_H
