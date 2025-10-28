#include "MappedFasta.h"
#include "../Utils.h"

#include <cctype>
#include <cassert>
#include <fstream>
#include <sstream>
#include <algorithm>

// ----- IUPAC helpers -----

static constexpr size_t ASCII_SIZE = 256;

// maps char -> {A,C,G,T,a,c,g,t stay same, other letters -> 'N'}
static char iupacToN[ASCII_SIZE];

// maps char -> uppercased A/C/G/T or 'N'
static char iupacToUpperOrN[ASCII_SIZE];

static void initIupacTables() {
    for (size_t c = 0; c < ASCII_SIZE; ++c) {
        char ch = static_cast<char>(c);

        // default:
        iupacToN[c]           = ch;
        iupacToUpperOrN[c]    = ch;

        const bool isAlpha = std::isalpha(static_cast<unsigned char>(ch)) != 0;
        if (isAlpha) {
            // normalize to uppercase for second table
            char up = static_cast<char>(std::toupper(static_cast<unsigned char>(ch)));

            // allow only A/C/G/T to persist; everything else -> N
            if (up == 'A' || up == 'C' || up == 'G' || up == 'T') {
                iupacToN[c]           = ch;
                iupacToUpperOrN[c]    = up;
            } else {
                iupacToN[c]           = 'N';
                iupacToUpperOrN[c]    = 'N';
            }
        }
    }
}

// returns 'N' for ambiguity but preserves case if it's canonical
static inline char ToN(char base) {
    return iupacToN[static_cast<unsigned char>(base)];
}

// returns uppercased A/C/G/T or 'N' for others
static inline char ToUpperAndN(char base) {
    return iupacToUpperOrN[static_cast<unsigned char>(base)];
}

// ----- MappedFasta -----

std::vector<uint8_t>& MappedFasta::Seq(std::string_view contig) {
    std::string key(contig);
    return fasta_[key];
}

MappedFasta::MappedFasta(const std::string &filename,
                         std::shared_ptr<Header> hdr,
                         bool toUpper,
                         bool toN)
    : hdr_(std::move(hdr))
{
    static bool tablesInit = false;
    if (!tablesInit) {
        initIupacTables();
        tablesInit = true;
    }

    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open fasta file " + filename);
    }

    std::string line;
    std::string contigName;
    std::vector<uint8_t> seq; // accumulate bases for current contig

    while (std::getline(file, line)) {
        if (line.empty()) {
            continue;
        }

        if (line[0] == '>') {
            // flush previous contig
            if (!contigName.empty()) {
                fasta_[contigName] = std::move(seq);
                seq.clear();
            }

            // parse new header line, e.g. ">chr1 blah blah"
            std::stringstream ss(line);
            std::string raw;
            ss >> raw; // raw might be ">chr1"
            if (raw.empty() || raw[0] != '>') {
                throw std::runtime_error("Malformed FASTA header line: " + line);
            }
            contigName = raw.substr(1); // drop '>'

        } else {
            // sequence line
            std::vector<uint8_t> lineData(line.begin(), line.end());

            if (toUpper && toN) {
                // collapse to canonical uppercase A/C/G/T or 'N'
                for (auto &ch : lineData) {
                    ch = static_cast<uint8_t>(ToUpperAndN(static_cast<char>(ch)));
                }
            } else if (toUpper) {
                for (auto &ch : lineData) {
                    ch = static_cast<uint8_t>(std::toupper(static_cast<unsigned char>(ch)));
                }
            } else if (toN) {
                for (auto &ch : lineData) {
                    ch = static_cast<uint8_t>(ToN(static_cast<char>(ch)));
                }
            }
            seq.insert(seq.end(), lineData.begin(), lineData.end());
        }
    }

    // flush last contig
    if (!contigName.empty()) {
        fasta_[contigName] = std::move(seq);
    }

    if (!file.eof()) {
        throw std::runtime_error("Error while reading FASTA file " + filename);
    }

    // ---- Flatten ----
    auto &seqDict = hdr_->getSequenceDictionary();
    const size_t nContigs = seqDict.getAll().size();

    flat_fasta_ptrs_.resize(nContigs);

    size_t totalBases = 0;
    // first pass: figure out lengths
    for (size_t chrIdx = 0; chrIdx < nContigs; ++chrIdx) {
        const std::string &chrName = seqDict.getSequenceName(chrIdx);

        auto it = fasta_.find(chrName);
        if (it != fasta_.end()) {
            const size_t L = it->second.size();
            flat_fasta_ptrs_[chrIdx] = SpanFieldLong{0, static_cast<unsigned int>(L)};
            totalBases += L;
        } else {
            flat_fasta_ptrs_[chrIdx] = SpanFieldLong{0, 0};
        }
    }

    flat_fasta_.resize(totalBases);

    // prefix-scan assign start offsets
    {
        size_t running = 0;
        for (size_t chrIdx = 0; chrIdx < nContigs; ++chrIdx) {
            auto &span = flat_fasta_ptrs_[chrIdx];
            span.start = running;
            running += span.len;
        }
    }

    // copy into flat buffer
    for (size_t chrIdx = 0; chrIdx < nContigs; ++chrIdx) {
        const std::string &chrName = seqDict.getSequenceName(chrIdx);
        auto it = fasta_.find(chrName);
        if (it == fasta_.end()) {
            continue;
        }

        const auto &bases = it->second;
        const auto &span  = flat_fasta_ptrs_[chrIdx];
        assert(span.len == bases.size());

        std::copy(bases.begin(), bases.end(),
                  flat_fasta_.begin() + span.start);
    }

    // if you don't need the map form anymore:
    // fasta_.clear();
    // fasta_.rehash(0);
}

// ----- KnownSites -----

KnownSites::KnownSites(std::vector<std::string> &knownSites,
                       std::shared_ptr<Header> hdr)
    : hdr_(std::move(hdr))
{
    // Load all intervals into knownIntervals_[contig]
    for (const std::string &fname : knownSites) {
        std::unordered_map<std::string, std::vector<Interval>> tmp;
        Interval::ImportIntervals(fname, tmp);

        for (auto &kv : tmp) {
            const std::string &chrom = kv.first;
            auto &ivals = kv.second;

            auto &dstVec = knownIntervals_[chrom];
            const size_t oldSz = dstVec.size();
            dstVec.reserve(oldSz + ivals.size());
            std::move(ivals.begin(), ivals.end(), std::back_inserter(dstVec));
        }
    }

    // sort + flatten each contig's intervals, then replace
    for (auto &kv : knownIntervals_) {
        auto &ivals = kv.second;
        Interval::ParallelSortByStart(ivals);
        ivals = Interval::ParallelFlatten(ivals);
    }

    // Flatten across contigs in header order
    auto &seqDict = hdr_->getSequenceDictionary();
    const size_t nContigs = seqDict.getAll().size();

    flat_known_intervals_ptrs_.resize(nContigs);

    size_t totalIntervals = 0;
    // first pass: lengths
    for (size_t chrIdx = 0; chrIdx < nContigs; ++chrIdx) {
        const std::string &chrName = seqDict.getSequenceName(chrIdx);

        auto it = knownIntervals_.find(chrName);
        if (it != knownIntervals_.end()) {
            const size_t L = it->second.size();
            flat_known_intervals_ptrs_[chrIdx] =
                SpanFieldLong{0, static_cast<unsigned int>(L)};
            totalIntervals += L;
        } else {
            flat_known_intervals_ptrs_[chrIdx] = SpanFieldLong{0, 0};
        }
    }

    flat_known_intervals_.resize(totalIntervals);

    // prefix-scan for starts
    {
        size_t running = 0;
        for (size_t chrIdx = 0; chrIdx < nContigs; ++chrIdx) {
            auto &span = flat_known_intervals_ptrs_[chrIdx];
            span.start = running;
            running += span.len;
        }
    }

    // copy intervals contig by contig into flat_known_intervals_
    for (size_t chrIdx = 0; chrIdx < nContigs; ++chrIdx) {
        const std::string &chrName = seqDict.getSequenceName(chrIdx);

        auto it = knownIntervals_.find(chrName);
        if (it == knownIntervals_.end()) {
            continue;
        }

        auto &ivals = it->second;
        const auto &span = flat_known_intervals_ptrs_[chrIdx];
        assert(span.len == ivals.size());

        std::copy(ivals.begin(), ivals.end(),
                  flat_known_intervals_.begin() + span.start);
    }

    // you can now drop the map if you only want flat data:
    // knownIntervals_.clear();
    // knownIntervals_.rehash(0);
}

// return pre-flatten intervals for a contig (host helper)
std::vector<Interval>& KnownSites::getKnownIntervals(std::string_view contig) {
    std::string key(contig);
    return knownIntervals_[key];
}

// return SpanFieldLong for a specific contig
SpanFieldLong KnownSites::getKnownIntervalsPtr(std::string_view contig) {
    auto &seqDict = hdr_->getSequenceDictionary();
    if (!seqDict.containsSequence(std::string(contig))) {
        return SpanFieldLong{0,0};
    }
    const size_t idx = seqDict.getIdx(std::string(contig));
    return flat_known_intervals_ptrs_[idx];
}