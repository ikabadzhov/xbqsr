#include <tbb/parallel_sort.h>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <charconv>
#include <zlib.h>
#include "Interval.h"


// SortByStart sorts a slice of Interval by Start position.
void SortByStart(std::vector<Interval>& intervals) {
    sort(intervals.begin(), intervals.end(), [](Interval& e1, Interval& e2) {return e1.Start < e2.Start;});
}

// ParallelSortByStart sorts a slice of Interval by Start position using
// a parallel stable sort.
void Interval::ParallelSortByStart( std::vector<Interval>& intervals ) {
    tbb::parallel_sort(intervals.begin(), intervals.end());
}

// Extend makes interval1 larger if it overlaps with interval2,
// by storing max(interval1.End, interval2.End) in interval1.End;
// otherwise, interval1 remains unchanged.
// Returns true if the two intervals overlap, false otherwise.
// interval2.Start >= interval1.Start must be true before
// calling Extend.
bool Interval::Extend( Interval& interval2 ) {
    if( interval2.Start > this->End ){
        return false;
    }
    if( interval2.End > this->End ){
        this->End = interval2.End;
    }
    return true;
}

// Flatten merges overlapping intervals into larger intervals.
// intervals must be sorted by Start before calling Flatten.
// The resulting slice is sorted by Start, and no two
// intervals in the result overlap with each other.
// The result shares memory with the intervals argument.
std::vector<Interval> Flatten(std::vector<Interval>& intervals) {
    for( size_t i=0, n = intervals.size()-1; i < n; i++ ){
        if( intervals[i].Extend(intervals[i+1]) ){
            n++;
            for( size_t j = i + 1; j < n; j++ ){
                if( !intervals[i].Extend(intervals[j]) ){
                    i++;
                    intervals[i] = intervals[j];
                }
            }
            return { intervals.begin(),intervals.begin()+(i+1) };//[:i+1]};
        }
    }
    return intervals;
}

const int32_t parallelFlattenGrainSize = 0x1000;

// ParallelFlatten merges overlapping intervals into larger intervals,
// using a parallel algorithm.
// intervals must be sorted by Start before calling Flatten.
// The resulting slice is sorted by Start, and no two
// intervals in the result overlap with each other.
// The result shares memory with the intervals argument.
std::vector<Interval> ParallelFlatten(std::vector<Interval>& intervals) {
    if( intervals.size() < parallelFlattenGrainSize ){
        return Flatten(intervals);
    }
    size_t half = intervals.size() >> 1;

    std::vector<Interval> left(intervals.begin(), intervals.begin()+half);//[:half]);
    std::vector<Interval> right(intervals.begin()+half, intervals.end());//[:half]);
    ParallelFlatten(left);
    ParallelFlatten(right);
    while( left[left.size()-1].Extend(right[0]) ){
        right.erase(right.cbegin());
    }
    left.insert(left.end(), right.begin(), right.end());
    return left;
}

// Overlap determines whether the given start/end range overlaps
// with any of the given intervals.
// intervals must be Flattened and sorted by Start.
bool Overlap(std::vector<Interval>& intervals, int32_t start, int32_t end) {
    for( size_t left=0, right = intervals.size()-1; left <= right; ) {
        size_t mid = (left + right) / 2;
        size_t intervalStart = intervals[mid].Start;
        size_t intervalEnd = intervals[mid].End;
        if( intervalStart > end-1 ){
            right = mid - 1;
        } else if( intervalEnd <= start-1 ){
            left = mid + 1;
        } else {
            return true;
        }
    }
    return false;
}

// Intersect returns a slice of all intervals that overlap with the
// given start/end range.
// intervals must be Flattened and sorted by Start.
// The result shares memory with the intervals argument.
std::pair<size_t, size_t>  Intersect(Interval* intervals, SpanFieldLong& intervalsPtrs, uint32_t start, uint32_t end) {

    auto start_it = std::lower_bound(intervals + intervalsPtrs.start, intervals + intervalsPtrs.start + intervalsPtrs.len, start, [](const auto& interval, uint32_t value) {
        return interval.End < value;
    });
    auto end_it = std::upper_bound(intervals + intervalsPtrs.start, intervals + intervalsPtrs.start + intervalsPtrs.len, end, [](uint32_t value, const auto& interval) {
        return value < interval.Start;
    });
    size_t idx_start = std::distance(intervals + intervalsPtrs.start, start_it);
    size_t idx_end = std::distance(intervals + intervalsPtrs.start, end_it);
    return {intervalsPtrs.start + idx_start, intervalsPtrs.start + idx_end};
}

// Helper to read arbitrary-length line from gzFile (supports plain and gzipped files)
static bool gz_getline(gzFile gz, std::string& line, std::vector<char>& buf) {
    line.clear();
    while (true) {
        if (!gzgets(gz, buf.data(), static_cast<int>(buf.size()))) {
            if (!line.empty()) {
                if (line.back() == '\r') {
                    line.pop_back();
                }
                return true;
            }
            return false;
        }
        line.append(buf.data());
        if (!line.empty() && line.back() == '\n') {
            line.pop_back();
            if (!line.empty() && line.back() == '\r') {
                line.pop_back();
            }
            return true;
        }
    }
}

// ImportIntervals loads intervals.
void Interval::ImportIntervals(const std::string& filename, std::unordered_map<std::string, std::vector<Interval>>& intervals) {
    auto ends_with = [](const std::string& str, const std::string& suffix) {
        if (str.size() < suffix.size()) return false;
        return std::equal(suffix.rbegin(), suffix.rend(), str.rbegin(),
                          [](char a, char b) {
                              return std::tolower(static_cast<unsigned char>(a)) == std::tolower(static_cast<unsigned char>(b));
                          });
    };

    if (ends_with(filename, ".bed") || ends_with(filename, ".bed.gz")) {
        ImportIntervalsFromBed(filename, intervals);
    } else if (ends_with(filename, ".vcf") || ends_with(filename, ".vcf.gz")) {
        ImportIntervalsFromCsv(filename, intervals);
    } else {
        std::cerr << "Warning: Unrecognized interval file format: " << filename << std::endl;
    }
}

bool Interval::operator<(const Interval &rhs) const {
    return Start < rhs.Start;
}

bool Interval::operator>(const Interval &rhs) const {
    return rhs < *this;
}

bool Interval::operator<=(const Interval &rhs) const {
    return !(rhs < *this);
}

bool Interval::operator>=(const Interval &rhs) const {
    return !(*this < rhs);
}

std::vector<Interval> Interval::ParallelFlatten(std::vector<Interval>& intervals) {
    //Not implemented yet
    return Interval::Flatten(intervals);
}

std::vector<Interval> Interval::Flatten(std::vector<Interval>& intervals) {
    for (size_t i = 0, n = intervals.size()-1; i < n; i++) {
        if (intervals[i].Extend(intervals[i+1])) {
            n++;
            for (size_t j = i + 1; j < n; j++) {
                if (!intervals[i].Extend(intervals[j])) {
                    i++;
                    intervals[i] = intervals[j];
                }
            }
            return {intervals.begin(), intervals.begin()+i+1};
        }
    }
    return intervals;
}

void Interval::ImportIntervalsFromBed(const std::string& filename, std::unordered_map<std::string, std::vector<Interval>> &intervals) {
    gzFile file = gzopen(filename.c_str(), "rb");
    if (!file) {
        std::cerr << "Error in reading " << filename << " file" << std::endl;
        exit(-1);
    }

    std::vector<char> buf(65536);
    std::string line;
    std::vector<std::string> fields;
    fields.resize(3);

    while (gz_getline(file, line, buf)) {
        if (line.empty() || line[0] == '#' || line.rfind("browser", 0) == 0 || line.rfind("track", 0) == 0) {
            continue;
        }

        size_t startPos = 0;
        auto endPos = std::string::npos;
        bool parse_ok = true;

        for (size_t i = 0; i < 3; i++) {
            if (i < 2) {
                endPos = line.find('\t', startPos);
                if (endPos == std::string::npos) {
                    parse_ok = false;
                    break;
                }
                fields[i] = line.substr(startPos, endPos - startPos);
                startPos = endPos + 1;
            } else {
                endPos = line.find('\t', startPos);
                if (endPos != std::string::npos) {
                    fields[i] = line.substr(startPos, endPos - startPos);
                } else {
                    fields[i] = line.substr(startPos);
                }
            }
        }

        if (!parse_ok) {
            std::cerr << "Error in parsing bed line: " << line << std::endl;
            gzclose(file);
            exit(-1);
        }

        uint32_t pos = std::stoul(fields[1]);
        uint32_t end = std::stoul(fields[2]);
        intervals[fields[0]].emplace_back(pos, end);
    }

    gzclose(file);
}

void Interval::ImportIntervalsFromCsv(const std::string& filename, std::unordered_map<std::string, std::vector<Interval>> &intervals) {
    gzFile file = gzopen(filename.c_str(), "rb");
    if (!file) {
        std::cerr << "Error in reading " << filename << " file" << std::endl;
        exit(-1);
    }

    std::vector<char> buf(65536);
    std::string line;
    std::vector<std::string> fields(4);

    bool found = false;
    size_t fieldCount = 0;

    while (gz_getline(file, line, buf)) {
        if (line.empty()) continue;
        if (line.size() >= 2 && line[0] == '#' && line[1] == '#') {
            continue;
        } else if (line[0] == '#') {
            size_t startPos = 0;
            auto endPos = std::string::npos;
            while ((endPos = line.find('\t', startPos)) != std::string::npos) {
                std::string fieldKey = line.substr(startPos, endPos - startPos);
                if (fieldKey == "END") {
                    found = true;
                    break;
                }
                fieldCount++;
                startPos = endPos + 1;
            }
            if (!found) {
                std::string fieldKey = line.substr(startPos);
                if (fieldKey == "END") {
                    found = true;
                }
            }
            continue;
        } else {
            size_t startPos = 0;
            auto endPos = std::string::npos;
            size_t nFieldsToScan = fields.size();
            if (found) {
                nFieldsToScan = fieldCount + 1;
                fields.resize(nFieldsToScan);
            }
            for (size_t i = 0; i < nFieldsToScan; i++) {
                if ((endPos = line.find('\t', startPos)) != std::string::npos) {
                    fields[i] = line.substr(startPos, endPos - startPos);
                    startPos = endPos + 1;
                } else {
                    if (i == nFieldsToScan - 1) {
                        fields[i] = line.substr(startPos);
                    } else {
                        std::cerr << "Error in parsing vcf file: " << filename << std::endl;
                        gzclose(file);
                        exit(-1);
                    }
                }
            }
        }
        uint32_t pos = std::stoul(fields[1]);
        uint32_t refLen = static_cast<uint32_t>(fields[3].size());
        uint32_t end = 0;
        if (found) {
            end = std::stoul(fields[fieldCount]);
        } else {
            end = pos - 1 + refLen;
        }
        intervals[fields[0]].emplace_back(pos, end);
    }

    gzclose(file);
}
