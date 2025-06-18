#ifndef GENOMICPIPELINE_SEQUENCE_H
#define GENOMICPIPELINE_SEQUENCE_H

#include <map>
#include <cstring>
#include <string>
#include <cstddef>
#include <numeric>
#include <vector>
#include "../../input_output/utils/Constant.h"


// In the following table, 'N' and all IUPAC codes are 0.
struct BaseToIntMap {
    std::vector<int32_t> baseToIntMap;
public:
    BaseToIntMap();
    int32_t operator[](char base);
};

struct SimpleBaseToBaseIndexTable {
    std::vector<int32_t> simpleBaseToBaseIndexTable;
public:
    SimpleBaseToBaseIndexTable();
    int32_t operator[](char base);
};

class BaseComplementTable {
    std::vector<char> baseComplementTable;
public:
    BaseComplementTable();
    char operator[](char base);
};

struct EncodeNucleotides {
    std::vector<int32_t> encodeNucleotides;
public:
    EncodeNucleotides(char nucleotide1, char nucleotide2);
    int32_t operator[](char base);
};

constexpr size_t BASE_MAP_SZ = 128;
// extern const int8_t baseToIntMap[BASE_MAP_SZ];
constexpr int8_t baseToIntMap[BASE_MAP_SZ] = {
        /* 0-31 */
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        /* 32-63*/
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1/* * */, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        /* 64-95 */
        0, 1/* A */, 0, 2/* C */, 0, 0, 0, 3/* G */, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 4/* T */, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        /* 96-111*/
        0, 1/* a */, 0, 2/* c */, 0, 0, 0, 3/* g */, 0, 0, 0, 0, 0, 0, 0, 0,
        /* 112-127*/
        0, 0, 0, 0, 4/* t*/, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
};

//extern const int8_t simpleBaseToBaseIndexTable[BASE_MAP_SZ];
constexpr int8_t simpleBaseToBaseIndexTable[BASE_MAP_SZ] = {
        /* 0-15 */
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
        /* 16-31 */
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
        /* 32-47 */
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0/* * */, -1, -1, -1, -1, -1,
        /* 48-63 */
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
        /* 64-79 */
        -1, 0/* A */, -1, 1/* C */, -1, -1, -1, 2/* G */, -1, -1, -1, -1, -1, -1, -1, -1,
        /* 80-95 */
        -1, -1, -1, -1, 3/* T */, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
        /* 96-111 */
        -1, 0/* a */, -1, 1/* c */, -1, -1, -1, 2/* g */, -1, -1, -1, -1, -1, -1, -1, -1,
        /* 112-127 */
        -1, -1, -1, -1, 3/* t */, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1
};

//extern const char baseComplementTable[BASE_MAP_SZ];

constexpr char baseComplementTable[BASE_MAP_SZ] = {
        /* 0-63 */ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27,
                   28, 29, 30, 31,
                   32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56,
                   57, 58, 59, 60, 61, 62, 63,
        /* 64-128 ('A' to '_') */
                   64, 'T' /* A -> T */, 66, 'G' /* C -> G */, 68, 69, 70, 'C' /* G-> C*/, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 'A' /* T -> A*/,
                   85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 'T' /* a -> T*/, 98, 'G' /* c -> G*/, 100, 101, 102, 'C' /* g -> C*/, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 'A' /* t -> A*/,
                   117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127};




class Sequence {
    std::string bytes;

public:
    Sequence() = default;

    explicit Sequence(size_t n);

    [[nodiscard]] size_t Len() const;

    [[nodiscard]] char Base(size_t index) const;
    void SetBase(size_t i, uint8_t base);
    void SetSequenceFromString(std::string_view sequence);
    const char* data() const;

    std::string_view String();

    std::string ParseString();
    void Append(char value);

    void AppendSlice(Sequence &m);

    void SliceInPlace(size_t low, size_t high);

    void Resize(size_t len);
    bool operator==(const Sequence &rhs) const;

    bool operator!=(const Sequence &rhs) const;

};


#endif //GENOMICPIPELINE_SEQUENCE_H