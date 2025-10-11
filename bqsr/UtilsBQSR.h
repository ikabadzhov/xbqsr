#ifndef UTILSBQSR_H
#define UTILSBQSR_H
#include <vector>
#include "../sam_data_structures/sam_fields/Cigar.h"
#include "../sam_data_structures/SamAlignments.h"


void deviceHardClipAdaptorSequence(CigarOperation *clipped_cigar, SpanField &cigar_ptr, int64_t& POS, char* sequence, SpanField& sequencePtrs, char* QUAL, SpanField& qualPtrs, size_t idx, int32_t RNAME, uint16_t FLAG, int32_t RNEXT, int32_t TLEN, int32_t PNEXT);
int64_t deviceComputeAdaptorBoundary(CigarOperation *cigar, SpanField &cigar_ptr, int64_t POS, int32_t RNAME, int32_t RNEXT, uint16_t FLAG, int32_t TLEN, int32_t PNEXT);
bool deviceIsInsideRead(CigarOperation *CIGAR, SpanField &cigar_ptr, int64_t POS, int64_t refCoord);
bool deviceHasWellDefinedFragmentSize(int32_t alnEnd, int64_t POS, int32_t RNAME, int32_t RNEXT, uint16_t FLAG, int32_t TLEN, int32_t PNEXT);
void deviceHardClip(CigarOperation *cigar, SpanField &cigar_ptr, SpanField& sequencePtrs, SpanField& qualPtrs, int64_t& POS, int32_t RNAME, uint16_t FLAG, int64_t start, int64_t stop);
void deviceHardClipByReferenceCoordinatesLeftTail(CigarOperation *clipped_cigar, SpanField &cigar_ptr, int64_t& POS, char* sequence, SpanField& sequencePtrs, char* QUAL, SpanField& qualPtrs, size_t idx, int32_t RNAME, uint16_t FLAG, int64_t refStop);
void deviceHardClipByReferenceCoordinatesRightTail(CigarOperation *clipped_cigar, SpanField &cigar_ptr, int64_t& POS, char* sequence, SpanField& sequencePtrs, char* QUAL, SpanField& qualPtrs, size_t idx, int32_t RNAME, uint16_t FLAG, int64_t refStart);
void deviceHardClipCigar(CigarOperation *cigar, SpanField &cigar_ptr, int64_t start, int64_t stop);
void deviceCleanHardClippedCigar(CigarOperation *cigar, SpanField &cigar_ptr);
std::pair<int64_t, CigarOperation> deviceGetReadCoordinateForReferenceCoordinate(CigarOperation* clipped_cigar, SpanField &clipped_cigar_ptr, int32_t softStart, int32_t refIndex);
std::pair<int32_t, CigarOperation> deviceComputeReadCoordinateForReferenceCoordinate(CigarOperation* cigar, SpanField &cigar_ptr, int32_t softStart, int32_t refCoord);
void deviceHardClipSoftClippedBases(CigarOperation* cigar, SpanField& cigar_ptr, int64_t& POS, char* sequence, SpanField& sequencePtrs, char* QUAL, SpanField& qualPtrs, int32_t RNAME, uint16_t FLAG);
void deviceHardClipCigar(CigarOperation* cigar, SpanField &cigar_ptr, int64_t start, int64_t stop);

#endif //UTILSBQSR_H
