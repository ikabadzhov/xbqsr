#ifndef BAIREADER_H
#define BAIREADER_H

#include <cstdint>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <vector>
#include <algorithm>
#include <filesystem>
#include <tbb/parallel_sort.h>

using VOff = uint64_t;

struct BaiLinearIndex {
    // all ioffsets for all references in one flat vector
    std::vector<VOff> ioffsets;
};

static inline bool voff_valid(VOff v) {
    if (v == 0) return false;
    const uint32_t uoff = uint32_t(v & 0xFFFFu);
    // uoff must fit within a BGZF block
    if (uoff > 65535u) return false;
    return true;
}

static uint32_t read_u32(FILE* f) {
    uint32_t x;
    if (fread(&x, 4, 1, f) != 1) throw std::runtime_error("BAI read_u32 failed");
    return x;
}
static uint64_t read_u64(FILE* f) {
    uint64_t x;
    if (fread(&x, 8, 1, f) != 1) throw std::runtime_error("BAI read_u64 failed");
    return x;
}

using VOff = uint64_t;

static inline uint64_t voff_coff(VOff v) { return v >> 16; }
static inline uint32_t voff_uoff(VOff v) { return uint32_t(v & 0xFFFFu); }

static std::vector<std::pair<VOff, VOff>>
split_ranges_equal_by_compressed_span(const std::vector<VOff>& ioffsets, int n_producers)
{
    if (n_producers <= 0) throw std::runtime_error("n_producers must be > 0");
    if ((int)ioffsets.size() < n_producers + 1)
        throw std::runtime_error("Not enough BAI ioffsets to split into ranges");

    // ioffsets must be sorted+unique already (your loader does that)
    const uint64_t coff_min = voff_coff(ioffsets.front());
    const uint64_t coff_max = voff_coff(ioffsets.back());
    const uint64_t span     = (coff_max > coff_min) ? (coff_max - coff_min) : 1;

    std::vector<std::pair<VOff, VOff>> ranges;
    ranges.reserve(n_producers);

    size_t last_idx = 0;
    for (int p = 1; p < n_producers; ++p) {
        const uint64_t target = coff_min + (span * uint64_t(p)) / uint64_t(n_producers);

        // Find first ioffset whose coff >= target, starting from last_idx
        auto it = std::lower_bound(ioffsets.begin() + last_idx, ioffsets.end(), target,
            [](VOff v, uint64_t coff_t) { return voff_coff(v) < coff_t; });

        size_t idx = (it == ioffsets.end()) ? (ioffsets.size() - 1) : size_t(it - ioffsets.begin());
        if (idx <= last_idx) idx = std::min(last_idx + 1, ioffsets.size() - 1);

        ranges.emplace_back(ioffsets[last_idx], ioffsets[idx]);
        last_idx = idx;
    }

    ranges.emplace_back(ioffsets[last_idx], ioffsets.back());

    // Remove empty/degenerate ranges
    ranges.erase(std::remove_if(ranges.begin(), ranges.end(),
        [](auto& r){ return r.first >= r.second; }), ranges.end());

    if ((int)ranges.size() == 0)
        throw std::runtime_error("Range split produced no ranges");

    return ranges;
}

static BaiLinearIndex load_bai_linear(const std::string& bai_path) {
    FILE* f = std::fopen(bai_path.c_str(), "rb");
    if (!f) throw std::runtime_error("Cannot open BAI: " + bai_path);

    char magic[4];
    if (fread(magic, 1, 4, f) != 4) throw std::runtime_error("BAI read magic failed");
    if (!(magic[0]=='B' && magic[1]=='A' && magic[2]=='I' && magic[3]==1))
        throw std::runtime_error("Not a BAI file (magic mismatch)");

    uint32_t n_ref = read_u32(f);

    BaiLinearIndex out;
    out.ioffsets.reserve(200000); // human-ish guess

    for (uint32_t r = 0; r < n_ref; ++r) {
        // bins
        uint32_t n_bin = read_u32(f);
        for (uint32_t b = 0; b < n_bin; ++b) {
            (void)read_u32(f);          // bin id
            uint32_t n_chunk = read_u32(f);
            // skip chunks: beg+end each
            for (uint32_t c = 0; c < n_chunk; ++c) {
                (void)read_u64(f);      // chunk beg
                (void)read_u64(f);      // chunk end
            }
        }

        // linear index
        uint32_t n_intv = read_u32(f);
        for (uint32_t i = 0; i < n_intv; ++i) {
            VOff voff = read_u64(f);
            if (voff_valid(voff)) out.ioffsets.push_back(voff);
        }

        // optional metadata pseudo-bin + trailing n_no_coor exists, but we can ignore:
        // if present, it's stored as a bin=37450 in n_bin above, already skipped.
    }

    std::fclose(f);

    // sort + unique; keep only record-aligned split candidates
    tbb::parallel_sort(out.ioffsets.begin(), out.ioffsets.end());
    out.ioffsets.erase(std::unique(out.ioffsets.begin(), out.ioffsets.end()), out.ioffsets.end());

    // remove pathological tiny offsets if needed
    // (e.g. offsets inside header). Usually first offsets are okay; but you may want to drop < first_alignment_voff.
    return out;
}

static std::string find_bai_path(const std::string& bam_path) {
    namespace fs = std::filesystem;

    fs::path p(bam_path);

    fs::path bai1 = p;
    bai1 += ".bai";          // file.bam.bai

    fs::path bai2 = p;
    bai2.replace_extension(".bai"); // file.bai

    if (fs::exists(bai1)) return bai1.string();
    if (fs::exists(bai2)) return bai2.string();

    throw std::runtime_error("BAI index not found for BAM: " + bam_path);
}

#endif //BAIREADER_H
