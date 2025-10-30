#ifndef GENOMICPIPELINE_SAMHEADERFORMATTER_H
#define GENOMICPIPELINE_SAMHEADERFORMATTER_H

#include <string>
#include "../../sam_data_structures/sam_header_fields/SequenceField.h"
#include "../../sam_data_structures/SamHeader.h"

class SamHeaderFormatter {
public:
    static void FormatSamHeader(Header &hdr, std::string& line);

private:
    static void writeHDLine(Header &hdr, std::string& line);
    static void getSQLine(const SequenceRecord& sequenceRecord, std::string& line);
    static void getPGLine(const ProgramRecord& program, std::string& line);
    static void getRGLine(const ReadGroupRecord &readGroupRecord, std::string& line);
    static void encodeTags(const SamHeaderRecord& sequenceRecord, std::vector<std::string>& fields, size_t offset);

    static void join(char sep, std::vector<std::string>& fields, std::string& line);

};


#endif //GENOMICPIPELINE_SAMHEADERFORMATTER_H
