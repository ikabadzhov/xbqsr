#ifndef GENOMICPIPELINE_STRINGSCANNER_H
#define GENOMICPIPELINE_STRINGSCANNER_H

#include <vector>
#include <string>
#include <boost/any.hpp>
#include "../../sam_data_structures/sam_fields/Sequence.h"


class StringScanner {
    size_t index;
    std::string_view data;

public:
    bool readUntil(char c, std::string &s);
    void readUntil2(char c1, char c2, std::string &s, char &b);
    bool readBytes(std::string &bytes);
    bool readByteUntil(char c, char &val);
    void reset(std::string_view s);
    size_t len();
    int32_t doInt32();
    uint64_t doUint();
    bool doSeq(std::string& seq);
    void doString(std::string &value);
};


#endif //GENOMICPIPELINE_STRINGSCANNER_H
