#ifndef GENOMICPIPELINE_MAPPEDFASTA_H
#define GENOMICPIPELINE_MAPPEDFASTA_H

#include <string>
#include <vector>
#include <memory>
#include <unordered_map>
#include <tbb/concurrent_unordered_map.h>

#include "../intervals/Interval.h"
#include "../sam_data_structures/SamHeader.h"

// forward-declare helper if you keep it somewhere else
// struct SpanFieldLong { uint64_t start; uint64_t len; };

extern std::vector<uint8_t> ElfastaMagic;

class MappedFasta {
    // per-contig reference sequence before flattening
    tbb::concurrent_unordered_map<std::string, std::vector<uint8_t>> fasta_;

    // header to map contig name -> index
    std::shared_ptr<Header> hdr_;

    // flattened view:
    //  - flat_fasta_ptrs_[chrIdx] = {start,len} into flat_fasta_
    //  - flat_fasta_[start ... start+len)
    std::vector<SpanFieldLong> flat_fasta_ptrs_;
    std::vector<uint8_t>       flat_fasta_;

public:
    MappedFasta(const std::string &filename,
                std::shared_ptr<Header> hdr,
                bool toUpper = false,
                bool toN     = false);

    // non-copyable (optional, but usually good for big buffers)
    MappedFasta(const MappedFasta&) = delete;
    MappedFasta& operator=(const MappedFasta&) = delete;

    // movable (optional)
    MappedFasta(MappedFasta&&) = default;
    MappedFasta& operator=(MappedFasta&&) = default;

    // return per-contig vector (host map view)
    std::vector<uint8_t>& Seq(std::string_view contig);

    // flattened data for GPU upload
    std::vector<uint8_t>&       getFullFasta()      { return flat_fasta_; }
    std::vector<SpanFieldLong>& getFullFastaPtrs()  { return flat_fasta_ptrs_; }

    const std::vector<uint8_t>&       getFullFasta()     const { return flat_fasta_; }
    const std::vector<SpanFieldLong>& getFullFastaPtrs() const { return flat_fasta_ptrs_; }
};


class KnownSites {
    // original, per-contig intervals before flattening
    std::unordered_map<std::string, std::vector<Interval>> knownIntervals_;

    // flattened view, by contig index
    std::vector<Interval>      flat_known_intervals_;
    std::vector<SpanFieldLong> flat_known_intervals_ptrs_;

    std::shared_ptr<Header> hdr_;

public:
    KnownSites(std::vector<std::string> &knownSites,
               std::shared_ptr<Header> hdr);

    // non-copyable (optional)
    KnownSites(const KnownSites&) = delete;
    KnownSites& operator=(const KnownSites&) = delete;

    // movable
    KnownSites(KnownSites&&) = default;
    KnownSites& operator=(KnownSites&&) = default;

    // per-contig pre-flatten view (only if still populated)
    std::vector<Interval>& getKnownIntervals(std::string_view contig);
    // pointer metadata per contig
    SpanFieldLong          getKnownIntervalsPtr(std::string_view contig);

    // flattened global views (for GPU upload)
    std::vector<Interval>&      getIntervals()        { return flat_known_intervals_; }
    std::vector<SpanFieldLong>& getIntervalsPtrs()    { return flat_known_intervals_ptrs_; }

    const std::vector<Interval>&      getIntervals()     const { return flat_known_intervals_; }
    const std::vector<SpanFieldLong>& getIntervalsPtrs() const { return flat_known_intervals_ptrs_; }
};

#endif // GENOMICPIPELINE_MAPPEDFASTA_H