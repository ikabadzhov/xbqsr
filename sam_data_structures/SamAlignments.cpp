#include "SamAlignments.h"

#include <future>

void SamAlignments::Resize(size_t n) {
    N_ALNS = n;
    PNEXT.resize(n);
    TLEN.resize(n);
    POS.resize(n);
    FLAG.resize(n);
    MAPQ.resize(n);
    QNAME.resize(n);
    TILE.resize(n);
    RNAME.resize(n);
    RNEXT.resize(n);
    SEQ.Resize(n);
    QUAL.Resize(n);
    CIGAR.Resize(n);
    TAGS.Resize(n);
    //UnusedTags.resize(n);
}

void SamAlignments::Permute(std::vector<size_t>& indices) {

    //auto pnext_future = std::async(std::launch::async, permute, std::ref(indices), PNEXT);
    // permute(indices, PNEXT);
    // permute(indices, TLEN);
    // permute(indices, POS);
    // permute(indices, FLAG);
    // permute(indices, MAPQ);
    // permute(indices, QNAME);
    // permute(indices, TILE);
    // permute(indices, RNAME);
    // permute(indices, RNEXT);
    // SEQ.Permute(indices);
    // QUAL.Permute(indices);
    // CIGAR.Permute(indices);
    // TAGS.Permute(indices);
    // TMP.Permute(indices);
    std::vector<std::future<void>> futures;

    futures.emplace_back(std::async(std::launch::async, permute<int32_t>, std::ref(indices), std::ref(PNEXT)));
    futures.emplace_back(std::async(std::launch::async, permute<int32_t>, std::ref(indices), std::ref(TLEN)));
    futures.emplace_back(std::async(std::launch::async, permute<int32_t>, std::ref(indices), std::ref(POS)));
    futures.emplace_back(std::async(std::launch::async, permute<uint16_t>, std::ref(indices), std::ref(FLAG)));
    futures.emplace_back(std::async(std::launch::async, permute<uint8_t>, std::ref(indices), std::ref(MAPQ)));
    futures.emplace_back(std::async(std::launch::async, permute<int32_t>, std::ref(indices), std::ref(QNAME)));
    futures.emplace_back(std::async(std::launch::async, permute<tileInfo>, std::ref(indices), std::ref(TILE)));
    futures.emplace_back(std::async(std::launch::async, permute<int32_t>, std::ref(indices), std::ref(RNAME)));
    futures.emplace_back(std::async(std::launch::async, permute<int32_t>, std::ref(indices), std::ref(RNEXT)));

    futures.emplace_back(std::async(std::launch::async, [&]() {
        SEQ.Permute(indices);
    }));

    futures.emplace_back(std::async(std::launch::async, [&]() {
        QUAL.Permute(indices);
    }));

    futures.emplace_back(std::async(std::launch::async, [&]() {
        CIGAR.Permute(indices);
    }));

    futures.emplace_back(std::async(std::launch::async, [&]() {
        TAGS.Permute(indices);
    }));

    futures.emplace_back(std::async(std::launch::async, [&]() {
        TMP.Permute(indices);
    }));

    // Wait for all asynchronous tasks to complete.
    for (auto& fut : futures) {
        fut.get();
    }
}

void SamAlignments::ResizeAdaptPos(size_t sz) {
    this->TMP.adaptPos.resize(sz);
}

void SamAlignments::ResizeAdaptScore(size_t sz) {
    this->TMP.adaptScore.resize(sz);
}

Header& SamAlignments::GetHdr() {
    return *hdr;
}

size_t SamAlignments::GetAlnsNum() const {
    return N_ALNS;
}

void SamAlignments::SetAlnsNum(size_t n) {
    N_ALNS = n;
}

ssize_t SamAlignments::GetQNAME(size_t idx) const {
    return QNAME[idx];
}

ssize_t SamAlignments::GetRNAME(size_t idx) const {
    return RNAME[idx];
}

const int32_t SamAlignments::GetTLEN(size_t idx) const {
    return TLEN[idx];
}

int32_t SamAlignments::GetPNEXT(size_t idx) const {
    return PNEXT[idx];
}

ssize_t SamAlignments::GetRNEXT(size_t idx) const {
    return RNEXT[idx];
}

int64_t SamAlignments::GetPOS(size_t idx) const {
    return POS[idx];
}

uint16_t SamAlignments::GetFLAG(size_t idx) const {
    return FLAG[idx];
}
std::span<char> SamAlignments::GetSEQ(size_t idx) {
    return SEQ[idx];
}

uint8_t SamAlignments::GetMAPQ(size_t idx) {
    return MAPQ[idx];
}

std::span<char> SamAlignments::GetQUAL(size_t idx) {
    return QUAL[idx];
}

std::span<CigarOperation> SamAlignments::GetCIGAR(size_t idx) {
    return CIGAR[idx];
}


/*int32_t SamAlignments::GetRG(size_t idx) const {
    return TAGS.RG.Get(idx);
}

std::string SamAlignments::GetStringRG(size_t idx) const {
    return TAGS.RG.toString(idx);
}*/




// IsMultiple checks for template having multiple segments in
// sequencing. See http://samtools.github.io/hts-specs/SAMv1.pdf -
// Section 1.4.2.
bool SamAlignments::IsMultiple(uint16_t readFlag) { return (readFlag & Multiple) != 0; }

// IsProper checks for each segment being properly aligned according
// to the aligner. See http://samtools.github.io/hts-specs/SAMv1.pdf -
// Section 1.4.2.
bool SamAlignments::IsProper(uint16_t readFlag) { return (readFlag & Proper) != 0; }

// IsUnmapped checks for segment unmapped. See
// http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.4.2.
bool SamAlignments::IsUnmapped(uint16_t readFlag) { return (readFlag & Unmapped) != 0; }

// IsNextUnmapped checks for next segment in the template
// unmapped. See http://samtools.github.io/hts-specs/SAMv1.pdf -
// Section 1.4.2.
bool SamAlignments::IsNextUnmapped(uint16_t readFlag) { return (readFlag & NextUnmapped) != 0; }

// IsReversed checks for SEQ being reversed complemented. See
// http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.4.2.
SYCL_EXTERNAL bool SamAlignments::IsReversed(uint16_t readFlag) { return (readFlag & Reversed) != 0; }

// IsNextReversed check for SEQ of the next segment in the template
// being reverse complemented. See
// http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.4.2.
bool SamAlignments::IsNextReversed(uint16_t readFlag) { return (readFlag & NextReversed) != 0; }

// IsFirst checks for being the first segment in the template. See
// http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.4.2.
bool SamAlignments::IsFirst(uint16_t readFlag) { return (readFlag & First) != 0; }

// IsLast checks for being the last segment in the template. See
// http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.4.2.
bool SamAlignments::IsLast(uint16_t readFlag) { return (readFlag & Last) != 0; }

// IsSecondary checks for secondary alignment. See
// http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.4.2.
bool SamAlignments::IsSecondary(uint16_t readFlag) { return (readFlag & Secondary) != 0; }

// IsQCFailed checks for not passing filters, such as platform/vendor
// quality controls. See http://samtools.github.io/hts-specs/SAMv1.pdf
// - Section 1.4.2.
bool SamAlignments::IsQCFailed(uint16_t readFlag) { return (readFlag & QCFailed) != 0; }

// IsDuplicate checks for PCR or optical duplicate. See
// http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.4.2.
bool SamAlignments::IsDuplicate(uint16_t readFlag) { return (readFlag & Duplicate) != 0; }

// IsSupplementary checks for supplementary alignment. See
// http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.4.2.
bool SamAlignments::IsSupplementary(uint16_t readFlag) { return (readFlag & Supplementary) != 0; }

// FlagEvery checks for every bit set in the given flag being also set
// in this->FLAG.
bool SamAlignments::FlagEvery(uint16_t readFlag, uint16_t flag) { return (readFlag & flag) == flag; }

// FlagSome checks for some bits set in the given flag being also set
// in this->FLAG.
bool SamAlignments::FlagSome(uint16_t readFlag, uint16_t flag) { return (readFlag & flag) != 0; }

// FlagNotEvery checks for not every bit set in the given flag being
// also set in this->FLAG.
bool SamAlignments::FlagNotEvery(uint16_t readFlag, uint16_t flag) { return (readFlag & flag) != flag; }

// FlagNotAny checks for not any bit set in the given flag being also
// set in this->FLAG.
bool SamAlignments::FlagNotAny(uint16_t readFlag, uint16_t flag) {
    return (readFlag & flag) == 0;
}

bool SamAlignments::IsTrueFragment(uint16_t flag) {
    return (flag & (Multiple | NextUnmapped)) != Multiple;
}

// Is this alignment definitely part of pair?
bool SamAlignments::IsTruePair(uint16_t flag) {
    return (flag & (Multiple | NextUnmapped)) == Multiple;
}

int32_t SamAlignments::End(CigarOperation* cigar, SpanField cigarPtr, int32_t pos) {
    int32_t length=0;
    for (size_t i = cigarPtr.start; i < (cigarPtr.start + cigarPtr.len); i++) {
        CigarOperation op = cigar[i];
        if (CigarOperatorConsumesReferenceBases[static_cast<uint8_t>(op.Operation)] > 0) {
            length += static_cast<int32_t>(CigarOperatorConsumesReferenceBases[static_cast<uint8_t>(op.Operation)] * op.Length);
        }
    }
    return pos + length - 1;
}

int32_t SamAlignments::End(const std::vector<CigarOperation>& cigar, int32_t pos) {
    int32_t length=0;
    for (size_t i = 0; i < (cigar.size()); i++) {
        CigarOperation op = cigar[i];
        if (CigarOperatorConsumesReferenceBases[static_cast<uint8_t>(op.Operation)] > 0) {
            length += static_cast<int32_t>(CigarOperatorConsumesReferenceBases[static_cast<uint8_t>(op.Operation)] * op.Length);
        }
    }
    return pos + length - 1;
}


void SamAlignments::SetFLAG(size_t idx, uint16_t flag) {
    FLAG[idx] = flag;
}

void SamAlignments::SetLIBID(size_t idx, int32_t libid) {
  TMP.LIBID[idx]=libid;
}

void SamAlignments::SetAdaptPOS(size_t idx, int32_t pos) {
    TMP.adaptPos[idx] = pos;
}

void SamAlignments::SetAdaptSCORE(size_t idx, int32_t score) {
    TMP.adaptScore[idx] = score;
}

uint16_t modFlag(uint16_t flag) {

    if ( (flag&Multiple) == 0 ) {
        flag = flag & ~NextUnmapped;
        flag = flag & ~NextReversed;
    }
    if( (flag&Unmapped) != 0 ) {
        flag = flag & ~Reversed;
    }
    if( (flag&NextUnmapped) != 0 ) {
        flag = flag & ~NextReversed;
    }
    return flag;
}

// CoordinateLess compares two alignments according to their
// coordinate. See http://samtools.github.io/hts-specs/SAMv1.pdf -
// Section 1.3, Tag @HD, SO.
bool CoordinateLess(size_t idx1, size_t idx2, SamAlignments& reads) {
    int32_t refid1 = reads.TMP.REFID[idx1];
    int32_t refid2 = reads.TMP.REFID[idx2];

    if( refid1 < refid2 ) {
        return refid1 >= 0;
    }
    if( refid2 < refid1 ) {
        return refid2 < 0;
    }
    if ( reads.POS[idx1] < reads.POS[idx2] ) {
        return true;
    }
    if( reads.POS[idx1] > reads.POS[idx2] ) {
        return false;
    }
    if (SamAlignments::IsReversed(reads.FLAG[idx1]) != SamAlignments::IsReversed(reads.FLAG[idx2]) ) {
        return !SamAlignments::IsReversed(reads.FLAG[idx1]);
    }
    if (QnameStringMapper::EMPTY_STRING != reads.QNAME[idx1] && QnameStringMapper::EMPTY_STRING != reads.QNAME[idx2]) {
        if( reads.QNAME[idx1] < reads.QNAME[idx2] ) {
            return true;
        }
        if ( reads.QNAME[idx1] > reads.QNAME[idx2] ) {
            return false;
        }
    }

    uint16_t flag1 = modFlag(reads.FLAG[idx1]);
    uint16_t flag2 = modFlag(reads.FLAG[idx2]);

    if( flag1 < flag2 ) {
        return true;
    }
    if( flag1 > flag2 ){
        return false;
    }
    if( reads.MAPQ[idx1] < reads.MAPQ[idx2] ) {
        return true;
    }
    if( reads.MAPQ[idx1] > reads.MAPQ[idx2] ) {
        return false;
    }
    if(SamAlignments::IsMultiple(reads.FLAG[idx1]) && SamAlignments::IsMultiple(reads.FLAG[idx2]) ) {
        int32_t nextRefid1 = reads.TMP.nextREFID[idx1];
        int32_t nextRefid2 = reads.TMP.nextREFID[idx2];

        if( nextRefid1 < nextRefid2 ) {
            return true; // no special treatment of negative values!
        }
        if( nextRefid1 > nextRefid2 ) {
            return false; // no special treatment of negative values!
        }
        if( reads.PNEXT[idx1] < reads.PNEXT[idx2] ) {
            return true;
        }
        if( reads.PNEXT[idx1] > reads.PNEXT[idx2] ) {
            return false;
        }
    }
    return reads.TLEN[idx1] < reads.TLEN[idx2];
}