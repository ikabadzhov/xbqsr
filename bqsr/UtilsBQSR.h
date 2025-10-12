#ifndef UTILSBQSR_H
#define UTILSBQSR_H
#include <vector>
#include "../sam_data_structures/sam_fields/Cigar.h"
#include "../sam_data_structures/SamAlignments.h"


SYCL_EXTERNAL void deviceHardClipAdaptorSequence(CigarOperation *clipped_cigar, SpanField &cigar_ptr, int32_t& POS, char* sequence, SpanField& sequencePtrs, char* QUAL, SpanField& qualPtrs, size_t idx, int32_t RNAME, uint16_t FLAG, int32_t RNEXT, int32_t TLEN, int32_t PNEXT);
SYCL_EXTERNAL int64_t deviceComputeAdaptorBoundary(CigarOperation *cigar, SpanField &cigar_ptr, int32_t POS, int32_t RNAME, int32_t RNEXT, uint16_t FLAG, int32_t TLEN, int32_t PNEXT);
SYCL_EXTERNAL bool deviceIsInsideRead(CigarOperation *CIGAR, SpanField &cigar_ptr, int32_t POS, int64_t refCoord);
SYCL_EXTERNAL bool deviceHasWellDefinedFragmentSize(int32_t alnEnd, int32_t POS, int32_t RNAME, int32_t RNEXT, uint16_t FLAG, int32_t TLEN, int32_t PNEXT);
SYCL_EXTERNAL void deviceHardClip(CigarOperation *cigar, SpanField &cigar_ptr, SpanField& sequencePtrs, SpanField& qualPtrs, int32_t& POS, int32_t RNAME, uint16_t FLAG, int64_t start, int64_t stop);
SYCL_EXTERNAL void deviceHardClipByReferenceCoordinatesLeftTail(CigarOperation *clipped_cigar, SpanField &cigar_ptr, int32_t& POS, char* sequence, SpanField& sequencePtrs, char* QUAL, SpanField& qualPtrs, size_t idx, int32_t RNAME, uint16_t FLAG, int64_t refStop);
SYCL_EXTERNAL void deviceHardClipByReferenceCoordinatesRightTail(CigarOperation *clipped_cigar, SpanField &cigar_ptr, int32_t& POS, char* sequence, SpanField& sequencePtrs, char* QUAL, SpanField& qualPtrs, size_t idx, int32_t RNAME, uint16_t FLAG, int64_t refStart);
SYCL_EXTERNAL void deviceHardClipCigar(CigarOperation *cigar, SpanField &cigar_ptr, int64_t start, int64_t stop);
SYCL_EXTERNAL void deviceCleanHardClippedCigar(CigarOperation *cigar, SpanField &cigar_ptr);
SYCL_EXTERNAL std::pair<int64_t, CigarOperation> deviceGetReadCoordinateForReferenceCoordinate(CigarOperation* clipped_cigar, SpanField &clipped_cigar_ptr, int32_t softStart, int32_t refIndex);
SYCL_EXTERNAL std::pair<int32_t, CigarOperation> deviceComputeReadCoordinateForReferenceCoordinate(CigarOperation* cigar, SpanField &cigar_ptr, int32_t softStart, int32_t refCoord);
SYCL_EXTERNAL void deviceHardClipSoftClippedBases(CigarOperation* cigar, SpanField& cigar_ptr, int32_t& POS, char* sequence, SpanField& sequencePtrs, char* QUAL, SpanField& qualPtrs, int32_t RNAME, uint16_t FLAG);
SYCL_EXTERNAL void deviceHardClipCigar(CigarOperation* cigar, SpanField &cigar_ptr, int64_t start, int64_t stop);

#endif //UTILSBQSR_H
