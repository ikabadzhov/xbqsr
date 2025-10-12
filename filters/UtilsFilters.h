#ifndef GENOMICPIPELINE_UTILSFILTERS_H
#define GENOMICPIPELINE_UTILSFILTERS_H

#include <cstdint>
#include <unordered_map>
#include <cstddef>
#include <span>
#include "../sam_data_structures/sam_fields/Cigar.h"
#include "../sam_data_structures/SamHeader.h"
#include "../sam_data_structures/SamAlignments.h"

using clippingTail = int32_t;

namespace ClippingTail {
    const clippingTail left = 0;
    const clippingTail right = 1;
}

enum ClippingRepresentation{
    SOFTCLIP_BASES, HARDCLIP_BASES
};


extern const int8_t operatorConsumesReadBases[128];
extern const int8_t operatorConsumesReferenceBases[128];

std::vector<CigarOperation> elementStradlessClippedRead(std::vector<CigarOperation>& newCigar, char op, int32_t relativeClippingPosition, std::int32_t clippedBases);

std::vector<CigarOperation> softClipEndOfRead(int32_t clipFrom, std::vector<CigarOperation>& cigars);

bool cigarContainsN(std::vector<CigarOperation>& cigarVec);
bool cigarContainsN(std::span<CigarOperation> cigarVec);
bool cigarContainsN(CigarOperation* cigar, SpanField& cigarPtrs);


bool alignmentAgreesWithHeader(Header& hdr, int32_t rname, int32_t pos);

SYCL_EXTERNAL bool isStrictUnmapped(bool isUnmapped, int32_t rname, int32_t pos);

SYCL_EXTERNAL bool isStrictNextUnmapped(bool isNextUnmapped, int32_t rnext, int32_t pnext);

bool hasWellDefinedFragmentSize(std::vector<CigarOperation>& CIGAR, int64_t POS, int32_t RNAME, int32_t RNEXT, uint16_t FLAG, int32_t TLEN, int32_t PNEXT);

int64_t computeAdaptorBoundary(std::vector<CigarOperation>& CIGAR, int64_t POS, int32_t RNAME, int32_t RNEXT, uint16_t FLAG, int32_t TLEN, int32_t PNEXT);

bool isInsideRead(std::vector<CigarOperation>& CIGAR, int64_t POS, int64_t refCoord);

bool isInsideDeletion(std::vector<CigarOperation>& cigar, int32_t offset);

int32_t nofAlignedBasesWithSoftClips(std::vector<CigarOperation>& cigar);

void hardClipAdaptorSequence(std::vector<CigarOperation>& clipped_cigar, int64_t& POS, char* sequence, SpanField& sequencePtrs, char* QUAL, SpanField& qualPtrs, size_t idx, int32_t RNAME, uint16_t FLAG, int32_t RNEXT, int32_t TLEN, int32_t PNEXT);

SYCL_EXTERNAL int64_t softStart(CigarOperation* cigar, int32_t cigar_len, int64_t pos);

SYCL_EXTERNAL int64_t softEnd(CigarOperation* cigar, int32_t cigar_len, int64_t pos, int64_t end);

void clipRead();

void clipByReferenceCoordinates(int64_t refStart, int64_t refStop, ClippingRepresentation clippingOp);

void hardClipByReferenceCoordinatesLeftTail(std::vector<CigarOperation>& clipped_cigar, int64_t& POS, char* sequence, SpanField& sequencePtrs, char* QUAL, SpanField& qualPtrs, size_t idx, int32_t RNAME, uint16_t FLAG, int64_t refStop);

void hardClipByReferenceCoordinatesRightTail(std::vector<CigarOperation>& clipped_cigar, int64_t& POS, char* sequence, SpanField& sequencePtrs, char* QUAL, SpanField& qualPtrs, size_t idx, int32_t RNAME, uint16_t FLAG, int64_t refStart);

std::pair<int32_t, CigarOperation> computeReadCoordinateForReferenceCoordinate(std::vector<CigarOperation>& clipped_cigar, int32_t softStart, int32_t refIndex);

std::pair<int32_t, CigarOperation> computeReadCoordinateForReferenceCoordinate(CigarOperation* cigar, SpanField cigar_ptr, int32_t softStart, int32_t refIndex);

SYCL_EXTERNAL std::pair<int64_t, CigarOperation> getReadCoordinateForReferenceCoordinate(CigarOperation* clipped_cigar, SpanField cigar_ptr, int32_t softStart, int32_t refCoord);

std::pair<int64_t, CigarOperation> getReadCoordinateForReferenceCoordinate(std::vector<CigarOperation>& clipped_cigar, int32_t softStart, int32_t refCoord);

int32_t calculateHardSoftOffset(std::vector<CigarOperation>& cigar);

int32_t calculateAlnStartShift(std::vector<CigarOperation>& cigar, std::vector<CigarOperation>& clippedCigar);

SYCL_EXTERNAL int32_t calculateHardClippingAlignmentShift(CigarOperation cigarOp, int cigarLength);

void hardClip(std::vector<CigarOperation>& cigar, char* SEQ, SpanField& sequencePtrs, char* QUAL, SpanField& qualPtrs, int64_t& POS, int32_t RNAME, uint16_t FLAG, int64_t start, int64_t stop);

void hardClipCigar(std::vector<CigarOperation>& cigar, std::vector<CigarOperation>& new_cigar, int64_t start, int64_t stop);

void cleanHardClippedCigar(std::vector<CigarOperation>& cigar);

void hardClipSoftClippedBases(std::vector<CigarOperation>& cigar, int64_t& POS, char* sequence, SpanField& sequencePtrs, char* QUAL, SpanField& qualPtrs, int32_t RNAME, uint16_t FLAG);




#endif //GENOMICPIPELINE_UTILSFILTERS_H
