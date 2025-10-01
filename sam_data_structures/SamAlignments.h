#ifndef GENOMICPIPELINE_SAMALIGNMENTS_H
#define GENOMICPIPELINE_SAMALIGNMENTS_H

#include <span>
#include "sam_fields/Sequence.h"
#include "sam_fields/Cigar.h"
#include "sam_fields/OptionalFields.h"
#include "sam_fields/TmpFields.h"
#include "../data_types/VectorCollection.h"
#include "SamHeader.h"
#include "../data_types/StringCollection.h"
#include "../tileinfo/tileInfo.h"

// Bit values for the FLAG field in the Alignment struct. See
// http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.4.2.

class tileInfo;
// Template having multiple segments in sequencing.
const uint16_t Multiple = 0x1;

// Each segment properly aligned according to the aligner.
const uint16_t Proper = 0x2;

// Segment unmapped.
const uint16_t Unmapped = 0x4;

// Next segment in the template unmapped.
const uint16_t NextUnmapped = 0x8;

// SEQ being reversed complemented.
const uint16_t Reversed = 0x10;

// SEQ of the next segment in the template being reverse
// complemented.
const uint16_t NextReversed = 0x20;

// The first segment in the template.
const uint16_t First = 0x40;

// The last segment in the template.
const uint16_t Last = 0x80;

// Secondary alignment.
const uint16_t Secondary = 0x100;

// Not passing filters, such as platform/vendor quality controls.
const uint16_t QCFailed = 0x200;

// PCR or optical duplicate.
const uint16_t Duplicate = 0x400;

// Supplementary alignment.
const uint16_t Supplementary = 0x800;

class SamAlignments {
private:
    /**
     * Header fields
     **/
    std::shared_ptr<Header> hdr;
    size_t N_ALNS{};

    /**
     * Mandatory fields
     **/
public:
    explicit SamAlignments(std::shared_ptr<Header> h) : hdr(std::move(h)) {}

    std::vector<int32_t> QNAME; // TODO: For the moment òet's keep it int32_t because the empty string is represented with -1, but later change everything in uint32_t
    std::vector<tileInfo> TILE;
    std::vector<int32_t> RNAME;
    std::vector<int32_t> RNEXT;

    VectorCollection<char> QUAL;
    VectorCollection<CigarOperation> CIGAR;

    std::vector<int32_t> TLEN;
    std::vector<int32_t> PNEXT;
    std::vector<int32_t> POS; // TODO: Check if it's enough int32
    std::vector<uint16_t> FLAG;
    std::vector<uint8_t> MAPQ;
    VectorCollection<char> SEQ;

    OptionalFields TAGS; // Keep ony the fields that will be actively use during the following stages
    //std::vector<UnusedOptionalFields> UnusedTags; // Store all others TAGS for compatibility purpose only

    TmpFields TMP;

public:

    void Resize(size_t n);

    [[nodiscard]] Header& GetHdr();
    size_t GetAlnsNum() const;
    void SetAlnsNum(size_t n);

    void Permute(std::vector<size_t>& indices);

    void ResizeAdaptPos(size_t sz);
    void ResizeAdaptScore(size_t sz);

    /**
     * Getters for a single value
     * */
    [[nodiscard]] ssize_t GetQNAME(size_t idx) const;
    [[nodiscard]] ssize_t GetRNAME(size_t idx) const;
    [[nodiscard]] const int32_t GetTLEN(size_t idx) const;
    [[nodiscard]] int32_t GetPNEXT(size_t idx) const;
    [[nodiscard]] ssize_t GetRNEXT(size_t idx) const;
    [[nodiscard]] int64_t GetPOS(size_t idx) const;
    [[nodiscard]] uint16_t GetFLAG(size_t idx) const;
    [[nodiscard]] uint8_t GetMAPQ(size_t idx);
    [[nodiscard]] std::span<char> GetSEQ(size_t idx);
    [[nodiscard]] std::span<char> GetQUAL(size_t idx);
    [[nodiscard]] std::span<CigarOperation> GetCIGAR(size_t idx);

    [[nodiscard]] std::string_view GetSR(size_t idx) const;

    void SetFLAG(size_t idx, uint16_t flag);
    void SetAdaptPOS(size_t idx, int32_t pos);
    void SetAdaptSCORE(size_t idx, int32_t score);
    void SetLIBID(size_t idx, int32_t libid);

    SYCL_EXTERNAL static bool IsMultiple(uint16_t readFlag);
    SYCL_EXTERNAL static bool IsUnmapped(uint16_t readFlag);
    SYCL_EXTERNAL static bool IsNextUnmapped(uint16_t readFlag);
    SYCL_EXTERNAL static bool IsNextReversed(uint16_t readFlag);
    SYCL_EXTERNAL static bool IsReversed(uint16_t readFlag);
    static bool IsLast(uint16_t readFlag);
    static bool IsSecondary(uint16_t readFlag);
    static bool IsQCFailed(uint16_t readFlag);
    static bool IsSupplementary(uint16_t readFlag);
    static bool IsDuplicate(uint16_t readFlag);
    static bool FlagEvery(uint16_t readFlag, uint16_t flag);
    static bool FlagSome(uint16_t readFlag, uint16_t flag);
    SYCL_EXTERNAL static bool FlagNotAny(uint16_t readFlag, uint16_t flag);
    static bool FlagNotEvery(uint16_t readFlag, uint16_t flag);
    static bool IsProper(uint16_t readFlag);
    SYCL_EXTERNAL static bool IsFirst(uint16_t readFlag);
    SYCL_EXTERNAL static int32_t End(CigarOperation* cigar, SpanField cigarPtr, int32_t pos);
    SYCL_EXTERNAL static int32_t End(const std::vector<CigarOperation>& cigar, int32_t pos);
    static bool IsTrueFragment(uint16_t flag);
    static bool IsTruePair(uint16_t flag);
};

bool CoordinateLess(size_t idx1, size_t idx2, SamAlignments& reads);

#endif //GENOMICPIPELINE_SAMALIGNMENTS_H
