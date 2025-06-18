#ifndef GENOMICPIPELINE_INTERVAL_H
#define GENOMICPIPELINE_INTERVAL_H

#include <vector>
#include <span>
#include <cstdint>
#include <unordered_map>
#include <string_view>
#include "../data_types/VectorCollection.h"
#include <sycl/sycl.hpp>

// Interval is a generic struct with a start and an end position.
class Interval {
public:
    uint32_t Start;
    uint32_t End;

    bool operator<(const Interval &rhs) const;
    bool operator>(const Interval &rhs) const;
    bool operator<=(const Interval &rhs) const;
    bool operator>=(const Interval &rhs) const;

    explicit Interval(uint32_t start = 0, uint32_t end = 0) : Start(start), End(end) {};
    bool Extend( Interval& interval2 );
    static void ImportIntervals(const std::string& filename, std::unordered_map<std::string, std::vector<Interval>>& intervals);
    static void ParallelSortByStart(std::vector<Interval>& intervals );
    static std::vector<Interval> ParallelFlatten(std::vector<Interval>& intervals);
    static std::vector<Interval> Flatten(std::vector<Interval>& intervals);

    static void ImportIntervalsFromBed(const std::string& filename, std::unordered_map<std::string, std::vector<Interval>> &intervals);

    static void ImportIntervalsFromCsv(const std::string& filename, std::unordered_map<std::string, std::vector<Interval>> &intervals);
};

SYCL_EXTERNAL std::pair<size_t, size_t> Intersect(Interval* intervals, SpanFieldLong& intervalsPtrs, uint32_t start, uint32_t end);



#endif //GENOMICPIPELINE_INTERVAL_H
