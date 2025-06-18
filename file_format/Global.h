#ifndef GENOMICPIPELINE_GLOBAL_H
#define GENOMICPIPELINE_GLOBAL_H

#include "../sam_data_structures/sam_fields/Cigar.h"
#include "../sam_data_structures/sam_fields/Sequence.h"

extern CigarOperationsTable cigarOperationsTable;

extern std::vector<char> nibbleToBase;



#endif //GENOMICPIPELINE_GLOBAL_H
