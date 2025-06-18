#ifndef GENOMICPIPELINE_CIGAR_H
#define GENOMICPIPELINE_CIGAR_H

#include <cstdint>
#include <string>
#include <iostream>
#include <vector>
#include <span>
#include <sycl/sycl.hpp>
#include "../../input_output/utils/Constant.h"
#include "../../data_types/VectorCollection.h"

const std::string cigarOperations = "MmIiDdNnSsHhPpXx=";

const size_t CIGAR_CONSUME_BASES_SZ = 128;


const int32_t CigarOperatorConsumesReadBases[CIGAR_CONSUME_BASES_SZ] = {
        /* 0-31 */
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        /* 32-63*/
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1/* = */, 0, 0,
        /* 64-79*/
        0, 0, 0, 0, 0, 0, 0, 0, 0, 1/* I */, 0, 0, 0, 1 /* M */, 0, 0,
        /* 80-95 */
        0, 0, 0, 1 /* S */, 0, 0, 0, 0, 1 /* X */, 0, 0, 0, 0, 0, 0, 0,
        /* 96-127*/
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
};
;/*={
        {'M',1},
        {'I',1},
        {'S',1},
        {'=',1},
        {'X',1}
        };*/

const int32_t CigarOperatorConsumesReferenceBases[CIGAR_CONSUME_BASES_SZ] = {
        /* 0-31 */
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        /* 32-63 */
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1/* = */, 0, 0,
        /* 64-79 */
        0, 0, 0, 0, 1/* D */, 0, 0, 0, 0, 0, 0, 0, 0, 1 /* M */, 1/* N */, 0,
        /* 80-95 */
        0, 0, 0, 0, 0, 0, 0, 0, 1 /* X */, 0, 0, 0, 0, 0, 0, 0,
        /* 96-127 */
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
};

;/*={
        {'M',1},
        {'D',1},
        {'N',1},
        {'=',1},
        {'X',1}
};*/

class CigarOperationsTable{
    std::vector<char> cigar_operations_table;
public:
    CigarOperationsTable();
    char &operator[](size_t i);
};

class CigarOperation {
public:
    uint32_t Length;
    char Operation;
    CigarOperation(uint32_t length = 0, char op = 0) : Length(length), Operation(op) {}

    SYCL_EXTERNAL bool operator==(const CigarOperation &rhs) const;

    SYCL_EXTERNAL bool operator!=(const CigarOperation &rhs) const;

};

CigarOperation new_cigar_operation( std::string& cigar, size_t i, size_t& j );
bool slow_scan_cigar_string(std::string& cigar, std::vector<CigarOperation>& out );
void  ScanCigarString( std::string& cigar_string, std::vector<CigarOperation>& out);

uint32_t ReadLengthFromCigar( const std::vector<CigarOperation>& cigar );
uint32_t ReadLengthFromCigar( std::span<CigarOperation> cigar );
uint32_t ReadLengthFromCigar(CigarOperation* cigar, SpanField cigarPtrs);

uint32_t ReferenceLengthFromCigar( const std::vector<CigarOperation>& cigar );
uint32_t ReferenceLengthFromCigar(CigarOperation* cigar, SpanField cigarPtrs);

void cigarToString(std::span<CigarOperation> cigar, std::string& res_string);
uint32_t cigarStringToValue(const std::string& cigarString);

extern std::unordered_map<std::string,std::vector<CigarOperation>> cigarCache;

#endif //GENOMICPIPELINE_CIGAR_H
