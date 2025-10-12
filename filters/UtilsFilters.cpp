#include "UtilsFilters.h"

const int CONSUME_BASES_SZ = 128;

//TODO: Size can be reduced with proper offsets
const int8_t operatorConsumesReadBases[CONSUME_BASES_SZ] = {
        /* 0-63 */ 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0,
                   0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 1/* = */, 0, 0,
        /* 64-71 ('A' to 'Z') */
                   0, 0, 0, 0, 0, 0, 0, 0,
        /* 72-79 */
                   0, 1/* I */, 0, 0, 0, 1/* M */, 0, 0,
        /* 80-87*/
                   0, 0, 0, 1/* S */, 0, 0, 0, 0,
        /* 88-95 ('[' to '`') */
                   1/* X */, 0, 0, 0, 0, 0, 0, 0,
        /* 96-127 (lowercase and extended) */
                   0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0};


const int8_t operatorConsumesReferenceBases[CONSUME_BASES_SZ] = {
        /* 0-63 */ 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0,
                   0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 1/* = */, 0, 0,
        /* 64-71 ('A' to 'Z') */
                   0, 0, 0, 0, 1/* D */, 0, 0, 0,
        /* 72-79 */
                   0, 0, 0, 0, 0, 1/* M */, 1/* N */, 0,
        /* 80-87*/
                   0, 0, 0, 0, 0, 0, 0, 0,
        /* 88-95 ('[' to '`') */
                   1/* X */, 0, 0, 0, 0, 0, 0, 0,
        /* 96-127 (lowercase and extended) */
                   0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0};

//const std::unordered_map<char, bool> operatorConsumesReadBases = {{'M', true}, {'I', true}, {'S', true}, {'=', true}, {'X', true}};
//const std::unordered_map<char, bool> operatorConsumesReferenceBases = {{'M', true}, {'D', true}, {'N', true}, {'=', true}, {'X', true}};

std::vector<CigarOperation> elementStradlessClippedRead(std::vector<CigarOperation>& newCigar, char op, int32_t relativeClippingPosition, std::int32_t clippedBases) {
    if (operatorConsumesReadBases[static_cast<uint8_t>(op)]>0) {
        if (operatorConsumesReferenceBases[static_cast<uint8_t>(op)]>0) {
            if (relativeClippingPosition > 0) {
                newCigar.emplace_back(CigarOperation(relativeClippingPosition, op));
            }
        } else {
            clippedBases += relativeClippingPosition;
        }
    } else if (relativeClippingPosition != 0) {
        std::cerr << "Unexpected non-0 relative clipping position in CleanSam." << std::endl;
    }
    newCigar.emplace_back(CigarOperation(clippedBases, 'S'));
    return newCigar;
}

std::vector<CigarOperation> softClipEndOfRead(int32_t clipFrom, std::vector<CigarOperation>& cigars) {
    int32_t pos = 0;
    clipFrom--;
    std::vector<CigarOperation> newCigar;
    for (auto& op : cigars) {
        int32_t endPos = pos;
        if(CigarOperatorConsumesReadBases[static_cast<uint8_t>(op.Operation)]>0) {
            endPos += CigarOperatorConsumesReadBases[static_cast<uint8_t>(op.Operation)] * op.Length;
        }
        if (endPos < clipFrom) {
            newCigar.emplace_back(op);
        } else {
            auto clippedBases = ReadLengthFromCigar(cigars) + clipFrom;
            newCigar = elementStradlessClippedRead(newCigar, op.Operation, clipFrom-pos, clippedBases);
            break;
        }
        pos += endPos;
    }
    return newCigar;
}

bool cigarContainsN(std::vector<CigarOperation>& cigarVec) {
    for (auto& cigarOp : cigarVec){
        if (cigarOp.Operation == 'N') {
            return true;
        }
    }
    return false;
}

bool cigarContainsN(CigarOperation* cigar, SpanField& cigarPtrs) {
    for (size_t i = cigarPtrs.start; i < cigarPtrs.start + cigarPtrs.len; i++) {
        CigarOperation cigarOp = cigar[i];
        if (cigarOp.Operation == 'N') {
            return true;
        }
    }
    return false;
}

bool alignmentAgreesWithHeader(Header& hdr, int32_t rname, int32_t pos) {
    for (auto& sq : hdr.getSequenceDictionary().getAll()) {
        if (hdr.getSequenceDictionary().getIdx(sq.getSequenceName()) == rname) {
            return pos <= sq.getSequenceLength();
        }
    }
    return false;
}

SYCL_EXTERNAL bool isStrictUnmapped(bool isUnmapped, int32_t rname, int32_t pos) {
    return (isUnmapped || SequenceRecord::isEmpty(rname) || rname == SequenceRecord::RESERVED_SEQUENCE_NAMES['*'] || pos==0);
}

SYCL_EXTERNAL bool isStrictNextUnmapped(bool isNextUnmapped, int32_t rnext, int32_t pnext) {
    return (isNextUnmapped || SequenceRecord::UNAVAILABLE_SEQUENCE_INDEX == rnext || rnext == SequenceRecord::RESERVED_SEQUENCE_NAMES['*'] || pnext == 0);
}

bool hasWellDefinedFragmentSize(std::vector<CigarOperation>& cigar, int64_t POS, int32_t RNAME, int32_t RNEXT, uint16_t FLAG, int32_t TLEN, int32_t PNEXT) {
    if (TLEN == 0) {
        return false;
    }
    if (!SamAlignments::IsMultiple(FLAG)) {
        return false;
    }
    if (isStrictUnmapped(SamAlignments::IsUnmapped(FLAG), RNAME, POS) ||
        isStrictNextUnmapped(SamAlignments::IsNextUnmapped(FLAG), RNEXT, PNEXT)) {
        return false;
    }

    if (SamAlignments::IsReversed(FLAG) == SamAlignments::IsNextReversed(FLAG)) {
        return false;
    }

    if (SamAlignments::IsReversed(FLAG)) {
        int32_t alnEnd = SamAlignments::End(cigar, POS);
        return alnEnd > PNEXT;
    } else {
        return POS <= PNEXT + TLEN;
    }
}

int64_t computeAdaptorBoundary(std::vector<CigarOperation>& cigar, int64_t POS, int32_t RNAME, int32_t RNEXT, uint16_t FLAG, int32_t TLEN, int32_t PNEXT) {
    if (!hasWellDefinedFragmentSize(cigar, POS, RNAME, RNEXT, FLAG, TLEN, PNEXT)){
        return INT64_MIN;
    }
    if (SamAlignments::IsReversed(FLAG)) {
        return PNEXT - 1;
    } else {
        return POS + std::abs(TLEN);
    }
}


bool isInsideRead(std::vector<CigarOperation>& CIGAR, int64_t POS, int64_t refCoord) {
    return refCoord >= POS && refCoord <= SamAlignments::End(CIGAR, POS);
}

bool isInsideDeletion(std::vector<CigarOperation>& cigar, int32_t offset) {
    if (offset < 0) {
        return false;
    }

    int32_t pos = 0, prevPos = 0;

    for (auto& ce : cigar) {
        switch (ce.Operation){
            case 'I':
            case 'S':
            case 'D':
            case 'M':
            case '=':
            case 'X':
                prevPos = pos;
                pos += int32_t(ce.Length);
        }

        if (prevPos < offset && pos >= offset && ce.Operation == 'D') {
            return true;
        }
    }
    return false;
}

int32_t nofAlignedBasesWithSoftClips(std::vector<CigarOperation>& cigar) {
    int32_t n = 0;
    for (auto& ce : cigar) {
        switch (ce.Operation) {
            case 'M':
            case '=':
            case 'X':
            case 'S':
                n += static_cast<int32_t>(ce.Length);
        }
    }
    return n;
}

void hardClipAdaptorSequence(std::vector<CigarOperation>& clipped_cigar, int64_t& POS, char* sequence, SpanField& sequencePtrs, char* QUAL, SpanField& qualPtrs, size_t idx, int32_t RNAME, uint16_t FLAG, int32_t RNEXT, int32_t TLEN, int32_t PNEXT) {
    int64_t adaptorBoundary = computeAdaptorBoundary(clipped_cigar, POS, RNAME, RNEXT, FLAG, TLEN, PNEXT); // Doesn't modify anything
    if (adaptorBoundary == INT64_MIN || !isInsideRead(clipped_cigar, POS, adaptorBoundary)) {
        return;
    }
    if (SamAlignments::IsReversed(FLAG)) {
        hardClipByReferenceCoordinatesLeftTail(clipped_cigar, POS, sequence, sequencePtrs, QUAL, qualPtrs, idx, RNAME, FLAG, adaptorBoundary);
    } else {
        hardClipByReferenceCoordinatesRightTail(clipped_cigar, POS, sequence, sequencePtrs, QUAL, qualPtrs, idx, RNAME, FLAG, adaptorBoundary);
    }
}

int64_t softStart(CigarOperation* cigar, int32_t cigar_len, int64_t pos) {
    auto start = pos;
    for (size_t i = 0; i < cigar_len; i++) {
        CigarOperation op = cigar[i];
        if (op.Operation == 'S') {
            start -= op.Length;
        } else if (op.Operation != 'H') {
            break;
        }
    }
    return static_cast<int32_t>(start);
}

int64_t softEnd(CigarOperation* cigar, int32_t cigar_len, int64_t pos, int64_t end) {
    int64_t softEnd = end;
    for (int64_t i = static_cast<int64_t>(cigar_len) - 1; i >= 0; i--) {
        CigarOperation op = cigar[i];
        if (op.Operation == 'S') {
            softEnd += static_cast<int32_t>(op.Length);
        } else if (op.Operation != 'H') {
            return softEnd;
        }
    }
    return end;
}

void clipByReferenceCoordinates(int64_t refStart, int64_t refStop, ClippingRepresentation clippingOp = HARDCLIP_BASES) {
    return;
}

void hardClipByReferenceCoordinatesLeftTail(std::vector<CigarOperation>& clipped_cigar, int64_t& POS, char* sequence, SpanField& sequencePtrs, char* QUAL, SpanField& qualPtrs, size_t idx, int32_t RNAME, uint16_t FLAG, int64_t refStop) {

    std::pair<int32_t, CigarOperation> t = getReadCoordinateForReferenceCoordinate(clipped_cigar, softStart(clipped_cigar.data(), clipped_cigar.size(), POS), refStop);

    int64_t stop = t.first - (CigarOperatorConsumesReadBases[static_cast<uint8_t>(t.second.Operation)]>0 ? 0 : 1);
    if (stop < 0) {
        //std::cerr << "reference coordinate matches a non-existing base in read at index " << idx << std::endl;
        exit(-1);
        //throw sycl::exception(); // TODO: Add other non-exception mechanism for error handling
    }
    hardClip(clipped_cigar, sequence, sequencePtrs, QUAL, qualPtrs, POS, RNAME, FLAG, 0, stop);
}

void hardClipByReferenceCoordinatesRightTail(std::vector<CigarOperation>& clipped_cigar, int64_t& POS, char* sequence, SpanField& sequencePtrs, char* QUAL, SpanField& qualPtrs, size_t idx, int32_t RNAME, uint16_t FLAG, int64_t refStart) {
    std::pair<int64_t, CigarOperation> t = getReadCoordinateForReferenceCoordinate(clipped_cigar, softStart(clipped_cigar.data(), clipped_cigar.size(), POS), refStart);
    int64_t start = t.first;

    if (start < 0) {
        //std::cerr << "reference coordinate matches a non-existing base in read at index " << idx << std::endl;
        exit(100);
    }
    hardClip(clipped_cigar, sequence, sequencePtrs, QUAL, qualPtrs, POS, RNAME, FLAG, start, sequencePtrs.len - 1);
}

std::pair<int32_t, CigarOperation> computeReadCoordinateForReferenceCoordinate(std::vector<CigarOperation>& cigar, int32_t softStart, int32_t refCoord) {
    if (refCoord < softStart) {
        return {-1, CigarOperation{}};
    }
    int firstReadPosOfElement = 0;              //inclusive
    int firstRefPosOfElement = softStart;  //inclusive
    int lastReadPosOfElement = 0;               //exclusive
    int lastRefPosOfElement = softStart;   //exclusive
    for (size_t i = 0; i < cigar.size(); i++) {
        CigarOperation& op = cigar[i];
        firstReadPosOfElement = lastReadPosOfElement;
        firstRefPosOfElement = lastRefPosOfElement;
        lastReadPosOfElement += (operatorConsumesReadBases[static_cast<uint8_t>(op.Operation)] > 0 ? op.Length : 0);
        lastRefPosOfElement += operatorConsumesReferenceBases[static_cast<uint8_t>(op.Operation)] || op.Operation == 'S' ? op.Length : 0;

        if (firstRefPosOfElement <= refCoord && refCoord < lastRefPosOfElement) {   // refCoord falls within this cigar element
            int32_t readPosAtRefCoord = firstReadPosOfElement + (operatorConsumesReadBases[static_cast<uint8_t>(op.Operation)] > 0 ? (refCoord - firstRefPosOfElement) : 0);
            return {readPosAtRefCoord, op};
        }
    }
    return {-1, CigarOperation{}};
}

std::pair<int32_t, CigarOperation> computeReadCoordinateForReferenceCoordinate(CigarOperation* cigar, SpanField cigar_ptr, int32_t softStart, int32_t refCoord) {
    if (refCoord < softStart) {
        return {-1, CigarOperation{}};
    }
    int firstReadPosOfElement = 0;              //inclusive
    int firstRefPosOfElement = softStart;  //inclusive
    int lastReadPosOfElement = 0;               //exclusive
    int lastRefPosOfElement = softStart;   //exclusive
    for (size_t i = 0; i < cigar_ptr.len; i++) {
        CigarOperation& op = cigar[cigar_ptr.start + i];
        firstReadPosOfElement = lastReadPosOfElement;
        firstRefPosOfElement = lastRefPosOfElement;
        lastReadPosOfElement += (operatorConsumesReadBases[static_cast<uint8_t>(op.Operation)] > 0 ? op.Length : 0);
        lastRefPosOfElement += operatorConsumesReferenceBases[static_cast<uint8_t>(op.Operation)] || op.Operation == 'S' ? op.Length : 0;

        if (firstRefPosOfElement <= refCoord && refCoord < lastRefPosOfElement) {   // refCoord falls within this cigar element
            int32_t readPosAtRefCoord = firstReadPosOfElement + (operatorConsumesReadBases[static_cast<uint8_t>(op.Operation)] > 0 ? (refCoord - firstRefPosOfElement) : 0);
            return {readPosAtRefCoord, op};
        }
    }
    return {-1, CigarOperation{}};
}


std::pair<int32_t, bool> readStartsWithInsertion(std::span<CigarOperation> cigarVec) {
    for(auto& element : cigarVec){
        if (element.Operation=='I') {
            return {element.Length, true};
        } else if (element.Operation=='H' || element.Operation=='S') {
            continue;
        } else {
            return {-1, false};
        }
    }
    return {-1, false};
}

std::pair<int64_t, CigarOperation> getReadCoordinateForReferenceCoordinate(std::vector<CigarOperation>& clipped_cigar, int32_t softStart, int32_t refIndex) {
    std::pair<int64_t, CigarOperation> t = computeReadCoordinateForReferenceCoordinate(clipped_cigar, softStart, refIndex);
    return t;
}

std::pair<int64_t, CigarOperation> getReadCoordinateForReferenceCoordinate(CigarOperation* clipped_cigar, SpanField clipped_cigar_ptr, int32_t softStart, int32_t refIndex) {
    std::pair<int64_t, CigarOperation> t = computeReadCoordinateForReferenceCoordinate(clipped_cigar, clipped_cigar_ptr, softStart, refIndex);
    return t;
}

int32_t calculateHardSoftOffset(std::vector<CigarOperation>& cigar) {
    int32_t size = 0;
    int64_t i = 0;
    for ( ; i < cigar.size(); i++) {
        auto cigarOp = cigar[i];
        if (cigarOp.Operation == 'H') {
            size += static_cast<int32_t>(cigarOp.Length);
        } else {
            break;
        }
    }
    for ( ; i < cigar.size(); i++) {
        auto cigarOp = cigar[i];
        if (cigarOp.Operation == 'S') {
            size += static_cast<int32_t>(cigarOp.Length);
        } else {
            break;
        }
    }
    return size;
}

int32_t calculateAlnStartShift(std::vector<CigarOperation>& cigar, std::vector<CigarOperation>& clipped_cigar) {
    return calculateHardSoftOffset(clipped_cigar) - calculateHardSoftOffset(cigar);
}

int32_t calculateHardClippingAlignmentShift(CigarOperation cigarOp, int cigarLength) {
    switch (cigarOp.Operation) {
        case 'I':
            return -cigarLength;
        case 'D':
        case 'N':
            return static_cast<int32_t>(cigarOp.Length);
        default:
            return 0;
    }
}

/**
 * Modify:
 * SEQ
 * QUAL
 * CIGAR
 * */

void hardClip(std::vector<CigarOperation>& cigar, char* SEQ, SpanField& sequencePtrs, char* QUAL, SpanField& qualPtrs, int64_t& POS, int32_t RNAME, uint16_t FLAG, int64_t start, int64_t stop) {
    std::vector<CigarOperation> clippedCigar;
    hardClipCigar(cigar, clippedCigar, start, stop);
    auto readLength = sequencePtrs.len;
    auto newLength = readLength - (stop - start + 1);
    auto copyStart = sequencePtrs.start;
    if (start == 0) {
        copyStart = sequencePtrs.start + stop + 1;
    }

    sequencePtrs.start = copyStart;
    sequencePtrs.len = newLength;

    qualPtrs = SpanField(copyStart, newLength);
    auto tmp_cigar = cigar;
    cigar = clippedCigar;
    bool isUnmapped = SamAlignments::IsUnmapped(FLAG);
    if (start == 0 && !isStrictUnmapped(isUnmapped, RNAME, POS)) {
        POS += calculateAlnStartShift(tmp_cigar, clippedCigar);
    }
}


/*void hardClipCigar(std::vector<CigarOperation>& cigar, std::vector<CigarOperation>& new_cigar, int64_t start, int64_t stop) {
    if (cigar.empty()) return;

    // Write pointer for in-place rewrite
    size_t w = 0;

    auto consumes_read = [&](char op) -> int32_t {
        return CigarOperatorConsumesReadBases[static_cast<uint8_t>(op)] > 0 ? 1 : 0;
    };

    // FRONT CLIP: clip [0..stop] into leading 'H'
    if (start == 0) {
        int64_t read_index = 0;                 // read bases consumed so far
        int64_t totalHardClipCount = stop - start + 1;  // size of requested clip
        int64_t alignmentShift = 0;

        // Fold pre-existing leading H into totalHardClipCount
        size_t i = 0;
        while (i < cigar.size() && cigar[i].Operation == 'H') {
            totalHardClipCount += cigar[i].Length;
            ++i;
        }

        // Walk until we have clipped up to stop
        for (; i < cigar.size() && read_index <= stop; ++i) {
            const auto op = cigar[i];
            const int32_t op_len = op.Length;
            const int32_t r = consumes_read(op.Operation) ? op_len : 0;

            if (read_index + r == stop + 1) {
                // Boundary lands exactly at the end of this op
                alignmentShift += calculateHardClippingAlignmentShift(op, op_len);

                // Write the new leading H
                cigar[w++] = CigarOperation(static_cast<int32_t>(totalHardClipCount + alignmentShift), 'H');
                // Copy the remainder of the CIGAR as-is
                for (size_t k = i + 1; k < cigar.size(); ++k) cigar[w++] = cigar[k];
                cigar.resize(w);
                cleanHardClippedCigar(cigar);
                return;
            } else if (read_index + r > stop + 1) {
                // Boundary falls inside this op: split it
                const int32_t clip_here = static_cast<int32_t>(stop - read_index + 1);
                const int32_t remain   = op_len - clip_here;

                alignmentShift += calculateHardClippingAlignmentShift(op, clip_here);

                // Write the new leading H
                cigar[w++] = CigarOperation(static_cast<int32_t>(totalHardClipCount + alignmentShift), 'H');

                // Write the remainder of the boundary op (if any)
                if (remain > 0) cigar[w++] = CigarOperation(remain, op.Operation);

                // Copy the rest
                for (size_t k = i + 1; k < cigar.size(); ++k) cigar[w++] = cigar[k];
                cigar.resize(w);
                cleanHardClippedCigar(cigar);
                return;
            }

            // We entirely consumed this op in the clipped region
            read_index += r;
            alignmentShift += calculateHardClippingAlignmentShift(op, r);
        }

        // If we reach here, we clipped through (or up to) the end of CIGAR read span.
        cigar.clear();
        cigar.emplace_back(CigarOperation(static_cast<int32_t>(totalHardClipCount + alignmentShift), 'H'));
        cleanHardClippedCigar(cigar);
        return;
    }

    // BACK CLIP: keep prefix up to 'start', then convert suffix [start..end] to trailing 'H'
    {
        int64_t read_index = 0;
        int64_t totalHardClipCount = 0;
        int64_t alignmentShift = 0;

        size_t i = 0;

        // Copy ops strictly before 'start'; possibly truncate boundary op
        for (; i < cigar.size(); ++i) {
            const auto op = cigar[i];
            const int32_t op_len = op.Length;
            const int32_t r = consumes_read(op.Operation) ? op_len : 0;

            if (read_index + r < start) {
                // Entire op lies before 'start' → keep as-is
                cigar[w++] = op;
            } else if (read_index < start && read_index + r >= start) {
                // Boundary falls inside this op: keep left piece
                const int32_t keep_len = static_cast<int32_t>(start - read_index);

                if (op.Operation == 'H') {
                    // If it's H, we conceptually count the left part as kept;
                    // do not emit partial H before boundary (kept H stays as-is in previous ops).
                    // But usually H doesn't consume read bases; boundary won't land inside H
                    // unless start==read_index, which is handled by branches.
                } else if (keep_len > 0) {
                    cigar[w++] = CigarOperation(keep_len, op.Operation);
                }

                // The rest of this op (right piece) is part of the clipped suffix
                const int32_t clip_here = op_len - keep_len;
                alignmentShift += calculateHardClippingAlignmentShift(op, clip_here);

                ++i; // advance past boundary op
                break;
            } else {
                // read_index >= start: we're already in the clipped suffix
                break;
            }

            read_index += r;
        }

        // Now sum the alignment shift and collect trailing H across the suffix (including from i..end)
        for (; i < cigar.size(); ++i) {
            const auto op = cigar[i];
            if (op.Operation == 'H') {
                totalHardClipCount += op.Length;
            } else {
                const int32_t r = consumes_read(op.Operation) ? op.Length : 0;
                alignmentShift += calculateHardClippingAlignmentShift(op, r);
            }
        }

        // Append a single trailing H that represents the whole clipped suffix
        cigar[w++] = CigarOperation(static_cast<int32_t>(totalHardClipCount + alignmentShift), 'H');

        // Shrink and clean
        cigar.resize(w);
        cleanHardClippedCigar(cigar);
    }
}
*/

// ORIGINAL VERSION
void hardClipCigar(std::vector<CigarOperation>& cigar, std::vector<CigarOperation>& new_cigar, int64_t start, int64_t stop) {
    int64_t index = 0;
    auto totalHardClipCount = stop - start + 1;
    int64_t alignmentShift = 0;
    //size_t j = 0;
    new_cigar.clear();
    if (start == 0) {
        int64_t cigarOpIndex = 0;
        CigarOperation cigarOp(0, 0);
        for ( ; cigarOpIndex < cigar.size(); cigarOpIndex++) {
            cigarOp = cigar[cigarOpIndex];
            if (cigarOp.Operation != 'H') {
                break;
            }
            totalHardClipCount += static_cast<int32_t>(cigarOp.Length);
        }
        for (; index <= stop && cigarOpIndex < cigar.size(); cigarOpIndex++) {
            cigarOp = cigar[cigarOpIndex];
            int32_t cigarOpLength = cigarOp.Length;
            int32_t shift = 0;
            if (CigarOperatorConsumesReadBases[static_cast<uint8_t>(cigarOp.Operation)]>0) {
                shift = CigarOperatorConsumesReadBases[static_cast<uint8_t>(cigarOp.Operation)] * cigarOpLength;
            }
            if (index+shift == stop+1) {
                alignmentShift += calculateHardClippingAlignmentShift(cigarOp, cigarOpLength);
                new_cigar.emplace_back(CigarOperation(int32_t(totalHardClipCount + alignmentShift), 'H'));
                //j++;
            } else if (index+shift > stop+1) {
                auto lengthAfterChopping = cigarOpLength - (stop - index + 1);
                alignmentShift += calculateHardClippingAlignmentShift(cigarOp, stop-index+1);
                new_cigar.emplace_back(CigarOperation(int32_t(totalHardClipCount + alignmentShift), 'H'));
                //j++;
                new_cigar.emplace_back(CigarOperation(int32_t(lengthAfterChopping), cigarOp.Operation));
                //j++;
            }
            index += shift;
            alignmentShift += calculateHardClippingAlignmentShift(cigarOp, shift);
        }

        for (size_t i = cigarOpIndex; i < cigar.size(); i++){
            new_cigar.emplace_back(cigar[i]);
            //j++;
        }
    } else {
        int32_t cigarOpIndex = 0;
        for (; index < start && cigarOpIndex < cigar.size(); cigarOpIndex++) {
            auto cigarOp = cigar[cigarOpIndex];
            int32_t cigarOpLength = cigarOp.Length;

            int32_t shift = 0;
            if (CigarOperatorConsumesReadBases[static_cast<uint8_t>(cigarOp.Operation)]>0) {
                shift = CigarOperatorConsumesReadBases[static_cast<uint8_t>(cigarOp.Operation)] * cigarOpLength;
            }
            if (index+shift < start) {
                new_cigar.emplace_back(cigarOp);
                //j++;
            } else {
                auto lengthAfterChopping = start - index;
                alignmentShift += calculateHardClippingAlignmentShift(cigarOp, cigarOpLength-(start-index));
                if (cigarOp.Operation == 'H') {
                    totalHardClipCount += lengthAfterChopping;
                } else {
                    new_cigar.emplace_back(CigarOperation(int32_t(lengthAfterChopping), cigarOp.Operation));
                    //j++;
                }
            }
            index += shift;
        }
        for (; cigarOpIndex < cigar.size(); cigarOpIndex++) {
            auto cigarOp = cigar[cigarOpIndex];
            alignmentShift += calculateHardClippingAlignmentShift(cigarOp, int(cigarOp.Length));
            if (cigarOp.Operation == 'H') {
                totalHardClipCount += int32_t(cigarOp.Length);
            }
        }
        new_cigar.emplace_back(CigarOperation(int32_t(totalHardClipCount + alignmentShift), 'H'));
        //j++;
    }
//    newCigarPtrs.start = newCigarPtrs.start;
//    newCigarPtrs.len = j;
    cleanHardClippedCigar(new_cigar);
}

// VERSION TO BRING ON DEVICE
//[[clang::optnone]]
void hardClipCigar2(std::vector<CigarOperation>& cigar, std::vector<CigarOperation>& new_cigar, int64_t start, int64_t stop) {
    std::vector<CigarOperation> tmp_cigar(cigar.begin(), cigar.end());
    CigarOperation buff[4] = { CigarOperation(0,'\0'), CigarOperation(0,'\0'),
                               CigarOperation(0,'\0'), CigarOperation(0,'\0') };

    int j = 0;                       // number of new ops produced (0..2)
    int64_t index = 0;               // read-bases consumed so far
    auto totalHardClipCount = stop - start + 1;
    int64_t alignmentShift = 0;
    new_cigar.clear();

    if (start == 0) {
        // Work on a temporary copy; extend space once for safe shifts.
        const size_t cigar_sz = tmp_cigar.size();
        tmp_cigar.resize(cigar_sz + 16);

        int64_t cigarOpIndex = 0;
        CigarOperation cigarOp(0, 0);

        // Fold leading H into totalHardClipCount (identical to original)
        for (; cigarOpIndex < (int64_t)cigar_sz; ++cigarOpIndex) {
            cigarOp = tmp_cigar[(size_t)cigarOpIndex];
            if (cigarOp.Operation != 'H') break;
            totalHardClipCount += static_cast<int32_t>(cigarOp.Length);
        }

        // Walk until index > stop or we reach the end (identical stopping condition).
        for (; index <= stop && cigarOpIndex < (int64_t)cigar_sz; ++cigarOpIndex) {
            cigarOp = tmp_cigar[(size_t)cigarOpIndex];
            const int32_t cigarOpLength = cigarOp.Length;

            int32_t shift = 0;
            if (CigarOperatorConsumesReadBases[(uint8_t)cigarOp.Operation] > 0) {
                shift = CigarOperatorConsumesReadBases[(uint8_t)cigarOp.Operation] * cigarOpLength;
            }

            if (index + shift == stop + 1) {
                alignmentShift += calculateHardClippingAlignmentShift(cigarOp, cigarOpLength);
                buff[j] = CigarOperation(int32_t(totalHardClipCount + alignmentShift), 'H');
                ++j; // only one op added here
            } else if (index + shift > stop + 1) {
                const int32_t clipLen = int32_t(stop - index + 1);
                const int32_t lengthAfterChopping = cigarOpLength - clipLen;
                alignmentShift += calculateHardClippingAlignmentShift(cigarOp, clipLen);

                buff[j] = CigarOperation(int32_t(totalHardClipCount + alignmentShift), 'H');
                ++j;
                buff[j] = CigarOperation(int32_t(lengthAfterChopping), cigarOp.Operation);
                ++j; // two ops added here
            }

            index += shift;
            alignmentShift += calculateHardClippingAlignmentShift(cigarOp, shift);
            // NB: No manual break here — we mirror the original loop semantics exactly.
        }

        // Now splice buff[0..j-1] + tail(tmp_cigar[cigarOpIndex..cigar_sz-1]) into tmp_cigar[...]
        if (j < cigarOpIndex) {
            // Shift tail left: write buff at [0..j-1], then copy tail to follow.
            size_t writePos = 0;
            for (int k = 0; k < j; ++k) {
                tmp_cigar[writePos++] = buff[k];
            }
            for (size_t i = (size_t)cigarOpIndex; i < cigar_sz; ++i) {
                tmp_cigar[writePos++] = tmp_cigar[i];
            }
            tmp_cigar.resize(writePos);
        } else {
            // j >= cigarOpIndex → shift tail right by delta, then write buff at front.
            const size_t tailStart = (size_t)cigarOpIndex;
            const size_t tailLen   = cigar_sz - tailStart;
            const size_t delta     = (size_t)j - (size_t)cigarOpIndex; // how much more we insert than we removed

            // Move tail from end to begin to avoid clobbering.
            for (size_t n = tailLen; n-- > 0; ) {
                tmp_cigar[tailStart + delta + n] = tmp_cigar[tailStart + n];
            }
            // Write new ops at the very front
            for (int k = 0; k < j; ++k) {
                tmp_cigar[(size_t)k] = buff[k];
            }
            tmp_cigar.resize(j + tailLen);
        }

        // Return result exactly as original would do
        new_cigar.resize(tmp_cigar.size());
        std::copy(tmp_cigar.begin(), tmp_cigar.end(), new_cigar.begin());
    } else {
        // Suffix clip (start > 0): keep prefix up to 'start', then add a single trailing 'H'

        // Work on tmp_cigar, and give ourselves a small buffer once.
        size_t cigar_sz = tmp_cigar.size();
        tmp_cigar.resize(cigar_sz + 16);  // safe working space; we will shrink later

        int32_t cigarOpIndex = 0;
        size_t w = 0;                    // write head for kept/left-piece ops

        // 1) Copy ops strictly before the boundary; if boundary lands inside an op, keep its left piece
        for (; index < start && cigarOpIndex < (int32_t)cigar_sz; ++cigarOpIndex) {
            auto cigarOp = tmp_cigar[(size_t)cigarOpIndex];
            int32_t cigarOpLength = cigarOp.Length;

            int32_t shift = 0;
            if (CigarOperatorConsumesReadBases[(uint8_t)cigarOp.Operation] > 0) {
                shift = CigarOperatorConsumesReadBases[(uint8_t)cigarOp.Operation] * cigarOpLength;
            }

            if (index + shift < start) {
                // Entire op lies before the clip boundary → keep as-is
                tmp_cigar[w++] = cigarOp;
            } else {
                // Boundary falls inside this op
                auto lengthAfterChoppingLeft = (int32_t)(start - index); // portion to KEEP on the left

                // Accumulate alignment shift for the CLIPPED right part (same as original)
                alignmentShift += calculateHardClippingAlignmentShift(
                                      cigarOp,
                                      cigarOpLength - lengthAfterChoppingLeft);

                if (cigarOp.Operation == 'H') {
                    // Original behavior: if boundary lands inside H, add kept amount to totalHardClipCount
                    totalHardClipCount += lengthAfterChoppingLeft;
                } else if (lengthAfterChoppingLeft > 0) {
                    // Keep the left piece
                    tmp_cigar[w++] = CigarOperation(lengthAfterChoppingLeft, cigarOp.Operation);
                }
            }

            index += shift;
        }

        // 2) Everything from cigarOpIndex..end is fully clipped → accumulate shift and H lengths
        for (; cigarOpIndex < (int32_t)cigar_sz; ++cigarOpIndex) {
            auto cigarOp = tmp_cigar[(size_t)cigarOpIndex];
            alignmentShift += calculateHardClippingAlignmentShift(cigarOp, (int)cigarOp.Length);
            if (cigarOp.Operation == 'H') {
                totalHardClipCount += (int32_t)cigarOp.Length;
            }
        }

        // 3) Append the single trailing H summarizing the clipped suffix
        tmp_cigar[w++] = CigarOperation((int32_t)(totalHardClipCount + alignmentShift), 'H');

        // 4) Shrink to final size and export to new_cigar (to match original signature/behavior)
        tmp_cigar.resize(w);
        new_cigar.resize(tmp_cigar.size());
        std::copy(tmp_cigar.begin(), tmp_cigar.end(), new_cigar.begin());
    }
//    newCigarPtrs.start = newCigarPtrs.start;
//    newCigarPtrs.len = j;
    cleanHardClippedCigar(new_cigar);
}

void cleanHardClippedCigar(std::vector<CigarOperation>& cigar) {
    if (cigar.empty()) return;

    const size_t oldN = cigar.size();

    // 1) Scan leading H/D/N and accumulate their total length.
    int32_t lead_sum = 0;
    size_t lead_i = 0;
    for (; lead_i < oldN; ++lead_i) {
        char op = cigar[lead_i].Operation;
        if (op == 'H' || op == 'D' || op == 'N') {
            lead_sum += cigar[lead_i].Length;
        } else {
            break;
        }
    }

    // 2) Scan trailing H/D/N and accumulate their total length.
    int32_t trail_sum = 0;
    size_t trail_start = oldN; // index of first trailing element to drop/merge
    if (oldN > 0) {
        ssize_t back = static_cast<ssize_t>(oldN) - 1;
        for (; back >= 0; --back) {
            char op = cigar[(size_t)back].Operation;
            if (op == 'H' || op == 'D' || op == 'N') {
                trail_sum += cigar[(size_t)back].Length;
            } else {
                break;
            }
        }
        trail_start = (size_t)back + 1; // [trail_start .. oldN) are trailing H/D/N
    }

    // 3) If there is no middle (everything is H/D/N), output a single H(total)
    //    identical to the original forward+backward behavior.
    size_t writePos = 0;
    if (lead_i >= trail_start) {
        // All H/D/N (or middle empty): merge into one H with the total length.
        // Here, lead_sum == sum of all H/D/N because lead_i consumed everything.
        cigar[0] = CigarOperation(lead_sum, 'H');
        writePos = 1;
        cigar.resize(writePos);
        return;
    }

    // 4) Otherwise, compose: [optional leading H] + middle + [optional trailing H]

    // Optional leading H (only if there were leading H/D/N)
    if (lead_i > 0) {
        cigar[writePos++] = CigarOperation(lead_sum, 'H');
    }

    // Middle segment: copy [lead_i .. trail_start)
    // (forward copy is safe: dest index <= src index)
    for (size_t k = lead_i; k < trail_start; ++k) {
        cigar[writePos++] = cigar[k];
    }

    // Optional trailing H (only if there were trailing H/D/N)
    if (trail_start < oldN) {
        cigar[writePos++] = CigarOperation(trail_sum, 'H');
    }

    // Single final shrink
    cigar.resize(writePos);
}


void cleanHardClippedCigarOriginal(std::vector<CigarOperation>& cigar) {
    int32_t totalHardClip = 0;
    int32_t index = 0;
    //forward:
    for (; index < cigar.size(); index++) {
        auto element = cigar[index];
        if (element.Operation=='H' || element.Operation=='D' || element.Operation=='N') {
            totalHardClip += element.Length;
        } else {
            break;
        }
    }
    if (index > 0) {
        cigar[0] = CigarOperation(totalHardClip, 'H');
        auto first = cigar.cbegin() + 1;
        auto last = cigar.cbegin() + index;
        cigar.erase(first, last);
        //cigar = append(cigar[:1], cigar[index:]...);
    }
    totalHardClip = 0;
    index = cigar.size() - 1;

    //backward:
    for (; index >= 0; index--) {
        auto element = cigar[index];
        if (element.Operation=='H' || element.Operation=='D' || element.Operation=='N') {
            totalHardClip += static_cast<int32_t>(element.Length);
        } else {
            break;
        }
    }
    if (index < cigar.size()-1) {
        cigar.resize(index + 1, CigarOperation(0, ' '));
        cigar.emplace_back(totalHardClip, 'H');
    }
}

void hardClipSoftClippedBases(std::vector<CigarOperation>& cigar, int64_t& POS, char* sequence, SpanField& sequencePtrs, char* QUAL, SpanField& qualPtrs, int32_t RNAME, uint16_t FLAG) {
    int32_t readIndex = 0;
    int32_t cutLeft = -1;
    int32_t cutRight = -1;
    int32_t rightTail = false;
    for (size_t i = 0; i < cigar.size(); i++) {
        CigarOperation cigarOp = cigar[i];
        auto key = cigarOp.Operation;
        auto ln = cigarOp.Length;
        switch (key) {
            case 'S':
                if (rightTail) {
                    cutRight = readIndex;
                } else {
                    cutLeft = readIndex + ln - 1;
                }
                break;
            case 'H':
            default:
                rightTail = true;
                break;
        }
        if(CigarOperatorConsumesReadBases[static_cast<uint8_t>(key)]>0) {
            readIndex += CigarOperatorConsumesReadBases[static_cast<uint8_t>(key)] * ln;
        }
    }
    if (cutRight >= 0) {
        hardClip(cigar, sequence, sequencePtrs, QUAL, qualPtrs, POS, RNAME, FLAG, cutRight, sequencePtrs.len-1);
    }
    if (cutLeft >= 0) {
        hardClip(cigar, sequence, sequencePtrs, QUAL, qualPtrs, POS, RNAME, FLAG, 0, cutLeft);
    }
}
