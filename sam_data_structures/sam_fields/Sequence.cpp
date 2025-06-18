#include <iostream>
#include <vector>
#include "Sequence.h"

/**
 * Base tables
 * */

//BaseToIntMap baseToIntMap;
//const int8_t baseToIntMap[BASE_MAP_SZ] = {
//        /* 0-39 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
//                   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, /* 40-47 */
//                   0, 1 /* '*' */, 0, 0, 0, 0, 0, /* 48-63 */
//        /* 64-95 ('A' to 'Z') */
//                   1 /* A */, 0, 2 /* C */, 0, 0, 0, 3 /* G */, 0, 0, 0, 0, 0, 0, 4 /* T */, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
//        /* 96-127 ('a' to 'z') */
//                   1 /* a */, 0, 2 /* c */, 0, 0, 0, 3 /* g */, 0, 0, 0, 0, 0, 0, 4 /* t */, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
//};

//SimpleBaseToBaseIndexTable simpleBaseToBaseIndexTable;
//const int8_t simpleBaseToBaseIndexTable[BASE_MAP_SZ] = {
//        /* 0-39 */ -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
//                   -1, -1, -1, -1, -1, -1, -1, 0 /* '*' */, -1, -1, -1, -1, -1, -1, /* 40-63 */
//        /* 64-95 ('A' to 'Z') */
//                   0 /* A */, -1, 1 /* C */, -1, -1, -1, 2 /* G */, -1, -1, -1, -1, -1, -1, 3 /* T */, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
//        /* 96-127 ('a' to 'z') */
//                   0 /* a */, -1, 1 /* c */, -1, -1, -1, 2 /* g */, -1, -1, -1, -1, -1, -1, 3 /* t */, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1
//};

//BaseComplementTable baseComplementTable;
//const char baseComplementTable[BASE_MAP_SZ] = {
//        /* 0-63 */ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27,
//                   28, 29, 30, 31,
//                   32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56,
//                   57, 58, 59, 60, 61, 62, 63,
//        /* 64-128 ('A' to '_') */
//                   64, 'T' /* A -> T */, 66, 'G' /* C -> G */, 68, 69, 70, 'C' /* G-> C*/, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 'A' /* T -> A*/,
//                   85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 'T' /* a -> T*/, 98, 'G' /* c -> G*/, 100, 101, 102, 'C' /* g -> C*/, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 'A' /* t -> A*/,
//                   117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127};
//


/*BaseToIntMap::BaseToIntMap() {
    baseToIntMap.resize(256, 0);
    baseToIntMap['a'] = 1;
    baseToIntMap['A'] = 1;
    baseToIntMap['*'] = 1;
    baseToIntMap['c'] = 2;
    baseToIntMap['C'] = 2;
    baseToIntMap['g'] = 3;
    baseToIntMap['G'] = 3;
    baseToIntMap['t'] = 4;
    baseToIntMap['T'] = 4;
}*/

/*
int32_t BaseToIntMap::operator[](char base) {
    return baseToIntMap[base];
}
*/

/*SimpleBaseToBaseIndexTable::SimpleBaseToBaseIndexTable() {
    simpleBaseToBaseIndexTable.resize(256, -1);
    simpleBaseToBaseIndexTable['A'] = 0;
    simpleBaseToBaseIndexTable['a'] = 0;
    simpleBaseToBaseIndexTable['*'] = 0;
    simpleBaseToBaseIndexTable['C'] = 1;
    simpleBaseToBaseIndexTable['c'] = 1;
    simpleBaseToBaseIndexTable['G'] = 2;
    simpleBaseToBaseIndexTable['g'] = 2;
    simpleBaseToBaseIndexTable['T'] = 3;
    simpleBaseToBaseIndexTable['t'] = 3;
}*/

/*int32_t SimpleBaseToBaseIndexTable::operator[](char base) {
    return simpleBaseToBaseIndexTable[base];
}

BaseComplementTable::BaseComplementTable() {
    baseComplementTable.resize(256);
    std::iota(baseComplementTable.begin(), baseComplementTable.end(), 0);
    baseComplementTable['A'] = 'T';
    baseComplementTable['a'] = 'T';
    baseComplementTable['C'] = 'G';
    baseComplementTable['c'] = 'G';
    baseComplementTable['G'] = 'C';
    baseComplementTable['g'] = 'C';
    baseComplementTable['T'] = 'A';
    baseComplementTable['t'] = 'A';
}*/
/*

char BaseComplementTable::operator[](char base) {
    return baseComplementTable[base];
}
*/


/**
 * Sequence class
 * */

// Make creates nibbles of the given length.
Sequence::Sequence(size_t n) {
    this->bytes.resize(n);
}

char Sequence::Base(size_t index) const {
    return bytes[index];
}

void Sequence::SetBase(size_t index, uint8_t base) {
    if(index >= bytes.size()) {
        std::cerr<<"index out of range"<<std::endl;
    }
    bytes[index] = base;
}


void Sequence::SetSequenceFromString(std::string_view sequence) {
    memcpy(bytes.data(), sequence.data(), bytes.size());
}


size_t Sequence::Len() const {
    return bytes.size();
}

void Sequence::Resize(size_t length) {
    bytes.resize(length);
}

bool Sequence::operator==(const Sequence &rhs) const {
    return bytes == rhs.bytes;
}

bool Sequence::operator!=(const Sequence &rhs) const {
    return !(rhs == *this);
}


std::string_view Sequence::String(){
    return bytes;
}

std::string Sequence::ParseString() {
    return bytes;
}

// SliceInPlace returns a subsequence of the given nibbles.
void Sequence::SliceInPlace(size_t low, size_t high) {
    if (high <= bytes.size() && low < high) {
        bytes.erase(bytes.begin() + high, bytes.end());
        bytes.erase(bytes.begin(), bytes.begin() + low);
        return;
    }
    exit(-1);
}

// AppendSlice appends the given nibbles.
void Sequence::AppendSlice(Sequence &m) {
    auto insertPosition = bytes.size();
    bytes.resize(bytes.size() + m.Len());
    std::copy(m.bytes.begin(), m.bytes.end(), bytes.begin() + insertPosition);
}

void Sequence::Append(char value) {
    bytes+=value;
}

const char *Sequence::data() const {
    return bytes.data();
}

