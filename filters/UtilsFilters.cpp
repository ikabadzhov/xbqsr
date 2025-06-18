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

bool isStrictUnmapped(bool isUnmapped, int32_t rname, int32_t pos) {
    return (isUnmapped || SequenceRecord::isEmpty(rname) || rname == SequenceRecord::RESERVED_SEQUENCE_NAMES['*'] || pos==0);
}

bool isStrictNextUnmapped(bool isNextUnmapped, int32_t rnext, int32_t pnext) {
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

    //SEQ.SliceInPlace(copyStart, copyStart + newLength);
    sequencePtrs.start = copyStart;
    sequencePtrs.len = newLength;

    //QUAL = QUAL.substr(copyStart, newLength);
    qualPtrs = SpanField(copyStart, newLength);
    auto tmp_cigar = cigar;
    cigar = clippedCigar;
    bool isUnmapped = SamAlignments::IsUnmapped(FLAG);
    if (start == 0 && !isStrictUnmapped(isUnmapped, RNAME, POS)) {
        POS += calculateAlnStartShift(tmp_cigar, clippedCigar);
    }
}

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
    assert(new_cigar.size() < cigar.size() + 6);
}


void cleanHardClippedCigar(std::vector<CigarOperation>& cigar) {
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
