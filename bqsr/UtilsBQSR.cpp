#include "UtilsBQSR.h"
#include "../filters/UtilsFilters.h"


void deviceHardClipAdaptorSequence(CigarOperation *clipped_cigar, SpanField& cigar_ptr, int32_t &POS, char *sequence, SpanField &sequencePtrs, char *QUAL, SpanField &qualPtrs, size_t idx, int32_t RNAME, uint16_t FLAG, int32_t RNEXT, int32_t TLEN, int32_t PNEXT) {
    int64_t adaptorBoundary = deviceComputeAdaptorBoundary(clipped_cigar, cigar_ptr, POS, RNAME, RNEXT, FLAG, TLEN, PNEXT); // Doesn't modify anything
    if (adaptorBoundary == INT64_MIN || !deviceIsInsideRead(clipped_cigar, cigar_ptr, POS, adaptorBoundary)) {
        return;
    }
    if (SamAlignments::IsReversed(FLAG)) {
        deviceHardClipByReferenceCoordinatesLeftTail(clipped_cigar,cigar_ptr, POS, sequence, sequencePtrs, QUAL, qualPtrs, idx, RNAME, FLAG, adaptorBoundary);
    } else {
        deviceHardClipByReferenceCoordinatesRightTail(clipped_cigar, cigar_ptr, POS, sequence, sequencePtrs, QUAL, qualPtrs, idx, RNAME, FLAG, adaptorBoundary);
    }
}

int64_t deviceComputeAdaptorBoundary(CigarOperation *cigar, SpanField &cigar_ptr, int32_t POS,
                                     int32_t RNAME, int32_t RNEXT, uint16_t FLAG, int32_t TLEN, int32_t PNEXT) {

    int32_t alnEnd = SamAlignments::End(cigar, cigar_ptr, POS);

    if (!deviceHasWellDefinedFragmentSize(alnEnd, POS, RNAME, RNEXT, FLAG, TLEN, PNEXT)){
        return INT64_MIN;
    }
    if (SamAlignments::IsReversed(FLAG)) {
        return PNEXT - 1;
    } else {
        return POS + std::abs(TLEN);
    }
}

bool deviceIsInsideRead(CigarOperation *CIGAR, SpanField &cigar_ptr, int32_t POS, int64_t refCoord) {
    return refCoord >= POS && refCoord <= SamAlignments::End(CIGAR, cigar_ptr, POS);
}

bool deviceHasWellDefinedFragmentSize(int32_t alnEnd, int32_t POS, int32_t RNAME, int32_t RNEXT, uint16_t FLAG, int32_t TLEN, int32_t PNEXT) {
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
        //int32_t alnEnd = SamAlignments::End(alnEnd, POS);
        return alnEnd > PNEXT;
    } else {
        return POS <= PNEXT + TLEN;
    }
}

// Sum leading H then S lengths in [cigar_ptr.start .. start+len)
inline int32_t hardSoftOffset(const CigarOperation* cigar_buf, const SpanField& cigar_ptr) {
    const CigarOperation* a = cigar_buf + cigar_ptr.start;
    const int n = static_cast<int>(cigar_ptr.len);
    int32_t off = 0;
    int i = 0;

    // leading H
    while (i < n && a[i].Operation == 'H') {
        off += a[i].Length;
        ++i;
    }
    // then leading S (contiguous after the Hs)
    while (i < n && a[i].Operation == 'S') {
        off += a[i].Length;
        ++i;
    }
    return off;
}

int32_t deviceCalculateHardSoftOffset(CigarOperation* cigar, SpanField& cigar_ptr) {
    int32_t size = 0;
    int64_t i = 0;
    for ( ; i < cigar_ptr.len; i++) {
        auto cigarOp = cigar[cigar_ptr.start + i];
        if (cigarOp.Operation == 'H') {
            size += static_cast<int32_t>(cigarOp.Length);
        } else {
            break;
        }
    }
    for ( ; i < cigar_ptr.len; i++) {
        auto cigarOp = cigar[cigar_ptr.start + i];
        if (cigarOp.Operation == 'S') {
            size += static_cast<int32_t>(cigarOp.Length);
        } else {
            break;
        }
    }
    return size;
}

// Clips in-place (your existing in-place routine), then POS shift if clipping is left-side.
inline void deviceHardClip(
    CigarOperation* cigar_buf, SpanField& cigar_ptr,
    SpanField& sequencePtrs, SpanField& qualPtrs,
    int32_t &POS, int32_t RNAME, uint16_t FLAG,
    int64_t start, int64_t stop)
{
    auto readLength = sequencePtrs.len;
    auto newLength = readLength - (stop - start + 1);
    auto copyStart = sequencePtrs.start;

    if (start == 0) {
        copyStart = sequencePtrs.start + stop + 1;
    }

    sequencePtrs.start = copyStart;
    sequencePtrs.len = newLength;
    qualPtrs.start = copyStart;
    qualPtrs.len = newLength;

    auto off_before = deviceCalculateHardSoftOffset(cigar_buf, cigar_ptr);

    deviceHardClipCigar(cigar_buf, cigar_ptr, start, stop);

    auto off_after = deviceCalculateHardSoftOffset(cigar_buf, cigar_ptr);

    //auto tmp_cigar = cigar;
    //cigar = clippedCigar;
    bool isUnmapped = SamAlignments::IsUnmapped(FLAG);
    if (start == 0 && !isStrictUnmapped(isUnmapped, RNAME, POS)) {
        POS += off_after - off_before;
    }

}

void deviceHardClipByReferenceCoordinatesLeftTail(CigarOperation *clipped_cigar, SpanField &cigar_ptr,
                                                  int32_t &POS, char *sequence, SpanField &sequencePtrs, char *QUAL, SpanField &qualPtrs, size_t idx,
                                                  int32_t RNAME, uint16_t FLAG, int64_t refStop) {

    int32_t sftStart = softStart(clipped_cigar + cigar_ptr.start, cigar_ptr.len, POS);
    std::pair<int32_t, CigarOperation> t = deviceGetReadCoordinateForReferenceCoordinate(clipped_cigar, cigar_ptr, sftStart, refStop);

    int64_t stop = t.first - (CigarOperatorConsumesReadBases[static_cast<uint8_t>(t.second.Operation)]>0 ? 0 : 1);
    if (stop < 0) {
        //exit(-1);
        //abort();
        assert(false && "Negative stop coordinate in hard clip");
    }
    deviceHardClip(clipped_cigar, cigar_ptr, sequencePtrs, qualPtrs, POS, RNAME, FLAG, 0, stop);
}

void deviceHardClipByReferenceCoordinatesRightTail(CigarOperation *clipped_cigar, SpanField &cigar_ptr,
                                                   int32_t &POS, char *sequence, SpanField &sequencePtrs, char *QUAL, SpanField &qualPtrs, size_t idx,
                                                   int32_t RNAME, uint16_t FLAG, int64_t refStart) {
    auto sftStart = softStart(clipped_cigar + cigar_ptr.start, cigar_ptr.len, POS);
    std::pair<int64_t, CigarOperation> t = deviceGetReadCoordinateForReferenceCoordinate(clipped_cigar, cigar_ptr, sftStart, refStart);
    int64_t start = t.first;

    if (start < 0) {
        //std::cerr << "reference coordinate matches a non-existing base in read at index " << idx << std::endl;
        //exit(100);
        assert(false && "Negative stop coordinate in hard clip");
        //abort();
    }
    deviceHardClip(clipped_cigar, cigar_ptr, sequencePtrs, qualPtrs, POS, RNAME, FLAG, start, sequencePtrs.len - 1);
}


std::pair<int64_t, CigarOperation> deviceGetReadCoordinateForReferenceCoordinate(CigarOperation *clipped_cigar,
    SpanField &clipped_cigar_ptr, int32_t softStart, int32_t refIndex) {
    std::pair<int64_t, CigarOperation> t = deviceComputeReadCoordinateForReferenceCoordinate(clipped_cigar, clipped_cigar_ptr, softStart, refIndex);
    return t;
}

void deviceCleanHardClippedCigar(CigarOperation* cigar, SpanField& cigar_ptr) {
    if (cigar_ptr.len == 0) return;
    auto* acc = cigar + cigar_ptr.start;
    const size_t oldN = cigar_ptr.len;

    // 1) Scan leading H/D/N and accumulate their total length.
    int32_t lead_sum = 0;
    size_t lead_i = 0;
    for (; lead_i < oldN; ++lead_i) {
        char op = acc[lead_i].Operation;
        if (op == 'H' || op == 'D' || op == 'N') {
            lead_sum += acc[lead_i].Length;
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
            char op = acc[(size_t)back].Operation;
            if (op == 'H' || op == 'D' || op == 'N') {
                trail_sum += acc[(size_t)back].Length;
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
        acc[0] = CigarOperation(lead_sum, 'H');
        writePos = 1;
        cigar_ptr.len = writePos;
        //cigar.resize(writePos);
        return;
    }

    // 4) Otherwise, compose: [optional leading H] + middle + [optional trailing H]

    // Optional leading H (only if there were leading H/D/N)
    if (lead_i > 0) {
        acc[writePos++] = CigarOperation(lead_sum, 'H');
    }

    // Middle segment: copy [lead_i .. trail_start)
    // (forward copy is safe: dest index <= src index)
    for (size_t k = lead_i; k < trail_start; ++k) {
        acc[writePos++] = acc[k];
    }

    // Optional trailing H (only if there were trailing H/D/N)
    if (trail_start < oldN) {
        acc[writePos++] = CigarOperation(trail_sum, 'H');
    }

    // Single final shrink
    cigar_ptr.len = writePos;
    //cigar.resize(writePos);
}

void deviceHardClipCigar(CigarOperation* cigar, SpanField &cigar_ptr, int64_t start, int64_t stop) {
    //std::vector<CigarOperation> tmp_cigar(cigar + cigar_ptr.start, cigar + cigar_ptr.start + cigar_ptr.len);
    auto* acc_cigar = cigar + cigar_ptr.start;
    CigarOperation buff[4] = { CigarOperation(0,'\0'), CigarOperation(0,'\0'),
                               CigarOperation(0,'\0'), CigarOperation(0,'\0') };

    int j = 0;                       // number of new ops produced (0..2)
    int64_t index = 0;               // read-bases consumed so far
    auto totalHardClipCount = stop - start + 1;
    int64_t alignmentShift = 0;
    //new_cigar.clear();

    if (start == 0) {
        // Work on a temporary copy; extend space once for safe shifts.
        const size_t cigar_sz = cigar_ptr.len;//.size();
        //tmp_cigar.resize(cigar_sz + 16);

        int64_t cigarOpIndex = 0;
        CigarOperation cigarOp(0, 0);

        // Fold leading H into totalHardClipCount (identical to original)
        for (; cigarOpIndex < (int64_t)cigar_sz; ++cigarOpIndex) {
            cigarOp = acc_cigar[(size_t)cigarOpIndex];
            if (cigarOp.Operation != 'H') break;
            totalHardClipCount += static_cast<int32_t>(cigarOp.Length);
        }

        // Walk until index > stop or we reach the end (identical stopping condition).
        for (; index <= stop && cigarOpIndex < (int64_t)cigar_sz; ++cigarOpIndex) {
            cigarOp = acc_cigar[(size_t)cigarOpIndex];
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
                acc_cigar[writePos++] = buff[k];
            }
            for (size_t i = (size_t)cigarOpIndex; i < cigar_sz; ++i) {
                acc_cigar[writePos++] = acc_cigar[i];
            }
            //acc_cigar.resize(writePos);
            cigar_ptr.len = writePos;
        } else {
            // j >= cigarOpIndex → shift tail right by delta, then write buff at front.
            const size_t tailStart = (size_t)cigarOpIndex;
            const size_t tailLen   = cigar_sz - tailStart;
            const size_t delta     = (size_t)j - (size_t)cigarOpIndex; // how much more we insert than we removed

            // Move tail from end to begin to avoid clobbering.
            for (size_t n = tailLen; n-- > 0; ) {
                acc_cigar[tailStart + delta + n] = acc_cigar[tailStart + n];
            }
            // Write new ops at the very front
            for (int k = 0; k < j; ++k) {
                acc_cigar[(size_t)k] = buff[k];
            }
            //tmp_cigar.resize(j + tailLen);
            cigar_ptr.len = (j + tailLen);
        }

        // Return result exactly as original would do
        //new_cigar.resize(tmp_cigar.size());

        //std::copy(tmp_cigar.begin(), tmp_cigar.end(), new_cigar.begin());
    } else {
        // Suffix clip (start > 0): keep prefix up to 'start', then add a single trailing 'H'

        // Work on tmp_cigar, and give ourselves a small buffer once.
        size_t cigar_sz = cigar_ptr.len;
        //tmp_cigar.resize(cigar_sz + 16);  // safe working space; we will shrink later

        int32_t cigarOpIndex = 0;
        size_t w = 0;                    // write head for kept/left-piece ops

        // 1) Copy ops strictly before the boundary; if boundary lands inside an op, keep its left piece
        for (; index < start && cigarOpIndex < (int32_t)cigar_sz; ++cigarOpIndex) {
            auto cigarOp = acc_cigar[(size_t)cigarOpIndex];
            int32_t cigarOpLength = cigarOp.Length;

            int32_t shift = 0;
            if (CigarOperatorConsumesReadBases[(uint8_t)cigarOp.Operation] > 0) {
                shift = CigarOperatorConsumesReadBases[(uint8_t)cigarOp.Operation] * cigarOpLength;
            }

            if (index + shift < start) {
                // Entire op lies before the clip boundary → keep as-is
                acc_cigar[w++] = cigarOp;
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
                    acc_cigar[w++] = CigarOperation(lengthAfterChoppingLeft, cigarOp.Operation);
                }
            }

            index += shift;
        }

        // 2) Everything from cigarOpIndex..end is fully clipped → accumulate shift and H lengths
        for (; cigarOpIndex < (int32_t)cigar_sz; ++cigarOpIndex) {
            auto cigarOp = acc_cigar[(size_t)cigarOpIndex];
            alignmentShift += calculateHardClippingAlignmentShift(cigarOp, (int)cigarOp.Length);
            if (cigarOp.Operation == 'H') {
                totalHardClipCount += (int32_t)cigarOp.Length;
            }
        }

        // 3) Append the single trailing H summarizing the clipped suffix
        acc_cigar[w++] = CigarOperation((int32_t)(totalHardClipCount + alignmentShift), 'H');

        // 4) Shrink to final size and export to new_cigar (to match original signature/behavior)
        //tmp_cigar.resize(w);
        cigar_ptr.len = w;
        //new_cigar.resize(tmp_cigar.size());
       // std::copy(tmp_cigar.begin(), tmp_cigar.end(), new_cigar.begin());
    }
    //std::vector<CigarOperation> tmp_cigar(cigar + cigar_ptr.start, cigar + cigar_ptr.start + cigar_ptr.len);
    deviceCleanHardClippedCigar(cigar, cigar_ptr);
    //std::memcpy(cigar + cigar_ptr.start, tmp_cigar.data(), tmp_cigar.size() * sizeof(CigarOperation));
    //cigar_ptr.len = tmp_cigar.size();
}


void deviceHardClipSoftClippedBases(CigarOperation *cigar, SpanField &cigar_ptr, int32_t &POS, char *sequence,
    SpanField &sequencePtrs, char *QUAL, SpanField &qualPtrs, int32_t RNAME, uint16_t FLAG) {
    int32_t readIndex = 0;
    int32_t cutLeft = -1;
    int32_t cutRight = -1;
    int32_t rightTail = false;
    for (size_t i = cigar_ptr.start; i < cigar_ptr.start + cigar_ptr.len; i++) {
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
        deviceHardClip(cigar, cigar_ptr, sequencePtrs, qualPtrs, POS, RNAME, FLAG, cutRight, sequencePtrs.len-1);
    }
    if (cutLeft >= 0) {
        deviceHardClip(cigar, cigar_ptr, sequencePtrs, qualPtrs, POS, RNAME, FLAG, 0, cutLeft);
    }
}

std::pair<int32_t, CigarOperation> deviceComputeReadCoordinateForReferenceCoordinate(CigarOperation *cigar,
                                                                                     SpanField &cigar_ptr, int32_t softStart, int32_t refCoord) {
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
