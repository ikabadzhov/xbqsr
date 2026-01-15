#ifndef BAIREADER_H
#define BAIREADER_H

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <tbb/parallel_sort.h>

using VOff = uint64_t;

// ------------------------------------------------------------
// VOff helpers
// ------------------------------------------------------------
static inline uint64_t voff_coff(VOff v) { return v >> 16; }
static inline uint32_t voff_uoff(VOff v) { return uint32_t(v & 0xFFFFu); }
static inline VOff make_voff(uint64_t coff, uint32_t uoff) {
  return (VOff(coff) << 16) | VOff(uoff);
}

class Slice {
public:
    VOff start;
    VOff end;
    ssize_t id;
};

static inline bool voff_valid(VOff v) {
  if (v == 0)
    return false;
  const uint32_t uoff = voff_uoff(v);
  if (uoff > 65535u)
    return false;
  return true;
}

// ------------------------------------------------------------
// Little-endian helpers
// ------------------------------------------------------------
static inline uint16_t rd_le16(const uint8_t *p) {
  return uint16_t(p[0]) | (uint16_t(p[1]) << 8);
}
static inline uint32_t rd_le32(const uint8_t *p) {
  return uint32_t(p[0]) | (uint32_t(p[1]) << 8) | (uint32_t(p[2]) << 16) |
         (uint32_t(p[3]) << 24);
}

// ------------------------------------------------------------
// BAI IO helpers
// ------------------------------------------------------------
static uint32_t read_u32(FILE *f) {
  uint32_t x;
  if (fread(&x, 4, 1, f) != 1)
    throw std::runtime_error("BAI read_u32 failed");
  return x;
}
static uint64_t read_u64(FILE *f) {
  uint64_t x;
  if (fread(&x, 8, 1, f) != 1)
    throw std::runtime_error("BAI read_u64 failed");
  return x;
}

// ------------------------------------------------------------
// Linear index
// ------------------------------------------------------------
struct BaiLinearIndex {
  std::vector<VOff> ioffsets; // flattened across references
};

// Tracks the last BAI path we loaded, so the splitter can infer BAM path
// without changing signatures (your main calls load_bai_linear() first).
static inline std::string &_last_loaded_bai_path() {
  static std::string s;
  return s;
}

static BaiLinearIndex load_bai_linear(const std::string &bai_path) {
  _last_loaded_bai_path() = bai_path;

  FILE *f = std::fopen(bai_path.c_str(), "rb");
  if (!f)
    throw std::runtime_error("Cannot open BAI: " + bai_path);

  char magic[4];
  if (fread(magic, 1, 4, f) != 4)
    throw std::runtime_error("BAI read magic failed");
  if (!(magic[0] == 'B' && magic[1] == 'A' && magic[2] == 'I' && magic[3] == 1))
    throw std::runtime_error("Not a BAI file (magic mismatch)");

  uint32_t n_ref = read_u32(f);

  BaiLinearIndex out;
  out.ioffsets.reserve(200000);

  for (uint32_t r = 0; r < n_ref; ++r) {
    // bins (skip)
    uint32_t n_bin = read_u32(f);
    for (uint32_t b = 0; b < n_bin; ++b) {
      (void)read_u32(f); // bin id
      uint32_t n_chunk = read_u32(f);
      for (uint32_t c = 0; c < n_chunk; ++c) {
        (void)read_u64(f); // chunk beg
        (void)read_u64(f); // chunk end
      }
    }

    // linear index
    uint32_t n_intv = read_u32(f);
    for (uint32_t i = 0; i < n_intv; ++i) {
      VOff voff = read_u64(f);
      if (voff_valid(voff))
        out.ioffsets.push_back(voff);
    }

    // metadata pseudo-bin (37450) is inside bins above; we are skipping it here
    // on purpose. Unmapped/no-coor will be handled by appending a tail slice to
    // EOF (see splitter).
  }

  std::fclose(f);

  tbb::parallel_sort(out.ioffsets.begin(), out.ioffsets.end());
  out.ioffsets.erase(std::unique(out.ioffsets.begin(), out.ioffsets.end()),
                     out.ioffsets.end());
  return out;
}

// ------------------------------------------------------------
// Find BAI for BAM (unchanged)
// ------------------------------------------------------------
static std::string find_bai_path(const std::string &bam_path) {
  namespace fs = std::filesystem;

  fs::path p(bam_path);

  fs::path bai1 = p;
  bai1 += ".bai"; // file.bam.bai

  fs::path bai2 = p;
  bai2.replace_extension(".bai"); // file.bai

  if (fs::exists(bai1))
    return bai1.string();
  if (fs::exists(bai2))
    return bai2.string();

  throw std::runtime_error("BAI index not found for BAM: " + bam_path);
}

// Infer BAM path from the last loaded BAI path (so we don't change signatures).
static std::string _infer_bam_path_from_last_bai() {
  namespace fs = std::filesystem;
  fs::path bai(_last_loaded_bai_path());
  if (bai.empty()) {
    throw std::runtime_error(
        "split_ranges_equal_by_compressed_span: no BAI path recorded. "
        "Did you call load_bai_linear() first?");
  }

  std::string s = bai.string();
  // Case 1: file.bam.bai  -> file.bam
  if (s.size() >= 4 && s.substr(s.size() - 4) == ".bai") {
    // remove ".bai"
    s.resize(s.size() - 4);
    // Now it is either ".../file.bam" (good) or ".../file" (if original was
    // file.bai)
  }

  fs::path p(s);
  if (p.extension() == ".bam")
    return p.string();

  // Case 2: file.bai came from bam_path replace_extension(".bai")
  // We try adding ".bam"
  fs::path bam = p;
  bam.replace_extension(".bam");
  if (fs::exists(bam))
    return bam.string();

  // Fallback: maybe the bam had the same basename and lives next to bai
  // (If this fails, user can still pass bam explicitly to print function;
  // slicing won't work)
  throw std::runtime_error("Could not infer BAM path from BAI path: " +
                           bai.string());
}

// ------------------------------------------------------------
// BGZF parsing: read one block {cbytes, ubytes} and advance to next block.
// Robustly handles optional gzip fields to avoid desync.
// ------------------------------------------------------------
static std::pair<uint32_t, uint32_t> read_bgzf_block_cbytes_ubytes(FILE *f) {
  uint8_t hdr[10];
  if (fread(hdr, 1, 10, f) != 10)
    throw std::runtime_error("EOF reading gzip header");
  if (hdr[0] != 31 || hdr[1] != 139 || hdr[2] != 8)
    throw std::runtime_error("Not gzip at expected BGZF block start");

  const uint8_t flg = hdr[3];
  if ((flg & 0x04) == 0)
    throw std::runtime_error("Not BGZF (FEXTRA flag not set)");

  // XLEN
  uint8_t xlen_buf[2];
  if (fread(xlen_buf, 1, 2, f) != 2)
    throw std::runtime_error("EOF reading XLEN");
  uint16_t xlen = rd_le16(xlen_buf);

  // Extra
  std::vector<uint8_t> extra(xlen);
  if (xlen && fread(extra.data(), 1, xlen, f) != xlen)
    throw std::runtime_error("EOF reading extra");

  // Find BC/BSIZE
  bool found_bc = false;
  uint16_t bsize = 0;
  size_t i = 0;
  while (i + 4 <= extra.size()) {
    uint8_t si1 = extra[i + 0];
    uint8_t si2 = extra[i + 1];
    uint16_t slen = rd_le16(&extra[i + 2]);
    i += 4;
    if (i + slen > extra.size())
      break;

    if (si1 == 'B' && si2 == 'C' && slen == 2) {
      bsize = rd_le16(&extra[i]); // BSIZE = total_block_size - 1
      found_bc = true;
      break;
    }
    i += slen;
  }
  if (!found_bc)
    throw std::runtime_error("Not BGZF (missing BC field)");

  uint32_t block_cbytes = uint32_t(bsize) + 1u;

  // bytes consumed so far in this member
  uint32_t consumed = 10u + 2u + uint32_t(xlen);

  // Optional gzip fields (rare in BGZF, but handle to avoid desync)
  if (flg & 0x08) { // FNAME
    int c;
    do {
      c = fgetc(f);
      if (c == EOF)
        throw std::runtime_error("EOF in FNAME");
      consumed++;
    } while (c != 0);
  }
  if (flg & 0x10) { // FCOMMENT
    int c;
    do {
      c = fgetc(f);
      if (c == EOF)
        throw std::runtime_error("EOF in FCOMMENT");
      consumed++;
    } while (c != 0);
  }
  if (flg & 0x02) { // FHCRC
    uint8_t hcrc[2];
    if (fread(hcrc, 1, 2, f) != 2)
      throw std::runtime_error("EOF in FHCRC");
    consumed += 2;
  }

  if (block_cbytes < consumed + 8u)
    throw std::runtime_error("Corrupt BGZF block (BSIZE too small)");

  // Skip payload to footer (CRC32+ISIZE)
  uint32_t to_skip = block_cbytes - consumed - 8u;
  if (to_skip) {
    if (fseek(f, long(to_skip), SEEK_CUR) != 0)
      throw std::runtime_error("fseek failed skipping payload");
  }

  uint8_t foot[8];
  if (fread(foot, 1, 8, f) != 8)
    throw std::runtime_error("EOF reading gzip footer");
  uint32_t isize = rd_le32(&foot[4]);
  if (isize > 65536u)
    throw std::runtime_error("Suspicious ISIZE > 65536");

  return {block_cbytes, isize};
}

// ------------------------------------------------------------
// Build a light BGZF index: block coff + ubytes + prefix_ubytes.
// We scan from the first relevant coff to EOF once.
// ------------------------------------------------------------
struct _BgzfLightIndex {
  std::vector<uint64_t> coff;          // block starts
  std::vector<uint32_t> ubytes;        // ISIZE per block
  std::vector<uint64_t> prefix_ubytes; // prefix sum at block start
};

static _BgzfLightIndex _build_bgzf_light_index(const std::string &bam_path,
                                               uint64_t coff_start) {
  FILE *f = std::fopen(bam_path.c_str(), "rb");
  if (!f)
    throw std::runtime_error("Cannot open BAM: " + bam_path);

#if defined(_WIN32)
  if (_fseeki64(f, (long long)coff_start, SEEK_SET) != 0) {
    std::fclose(f);
    throw std::runtime_error("_fseeki64 failed seeking BAM");
  }
#else
  if (fseeko(f, (off_t)coff_start, SEEK_SET) != 0) {
    std::fclose(f);
    throw std::runtime_error("fseeko failed seeking BAM");
  }
#endif

  _BgzfLightIndex idx;
  idx.coff.reserve(1 << 20);
  idx.ubytes.reserve(1 << 20);
  idx.prefix_ubytes.reserve(1 << 20);

  uint64_t prefix = 0;

  while (true) {
    // current coff
#if defined(_WIN32)
    long long pos_ll = _ftelli64(f);
    if (pos_ll < 0) {
      std::fclose(f);
      throw std::runtime_error("_ftelli64 failed");
    }
    uint64_t coff = uint64_t(pos_ll);
#else
    off_t pos = ftello(f);
    if (pos < 0) {
      std::fclose(f);
      throw std::runtime_error("ftello failed");
    }
    uint64_t coff = uint64_t(pos);
#endif

    // EOF check
    int c = std::fgetc(f);
    if (c == EOF)
      break;
    std::ungetc(c, f);

    auto [cbytes, ubytes] = read_bgzf_block_cbytes_ubytes(f);

    idx.coff.push_back(coff);
    idx.ubytes.push_back(ubytes);
    idx.prefix_ubytes.push_back(prefix);
    prefix += uint64_t(ubytes);
    (void)cbytes; // file position already advanced
  }

  std::fclose(f);
  return idx;
}

static VOff _compute_eof_voff_from_light_index(const _BgzfLightIndex &bg,
                                               uint64_t coff_start_global) {
  // coff_start_global is unused in this form; kept for clarity.
  if (bg.coff.empty())
    throw std::runtime_error("BGZF index empty (scan produced no blocks)");
  size_t last = bg.coff.size() - 1;
  return make_voff(bg.coff[last],
                   bg.ubytes[last]); // end of last block (uoff = ISIZE)
}

// Snap helpers for over-estimating “full blocks touched” between two voffs.
static size_t _snap_beg_block_idx(const _BgzfLightIndex &bg,
                                  uint64_t coff_beg) {
  // want last block start <= coff_beg
  auto it = std::upper_bound(bg.coff.begin(), bg.coff.end(), coff_beg);
  if (it == bg.coff.begin())
    return 0;
  return size_t((it - bg.coff.begin()) - 1);
}
static size_t _snap_end_block_idx(const _BgzfLightIndex &bg,
                                  uint64_t coff_end) {
  // want first block start >= coff_end
  auto it = std::lower_bound(bg.coff.begin(), bg.coff.end(), coff_end);
  if (it == bg.coff.end())
    return bg.coff.size() - 1;
  return size_t(it - bg.coff.begin());
}

static uint64_t _overestimated_span_bytes_full_blocks(const _BgzfLightIndex &bg,
                                                      VOff beg, VOff end) {
  uint64_t cb = voff_coff(beg);
  uint64_t ce = voff_coff(end);
  if (ce < cb)
    throw std::runtime_error("range end before start");

  size_t ib = _snap_beg_block_idx(bg, cb);
  size_t ie = _snap_end_block_idx(bg, ce);

  // bytes from start of ib to end of ie (inclusive)
  uint64_t start_prefix = bg.prefix_ubytes[ib];
  uint64_t end_prefix = bg.prefix_ubytes[ie] + uint64_t(bg.ubytes[ie]);
  if (end_prefix < start_prefix)
    return 0;
  return end_prefix - start_prefix;
}

// ------------------------------------------------------------
// (1)+(2) Improved splitter, same signature.
// It returns many ranges so each estimated decompressed size <= limit.
// Adds an extra tail slice to EOF to include unmapped/no-coor reads.
// ------------------------------------------------------------
static std::vector<std::pair<VOff, VOff>>
split_ranges_equal_by_compressed_span(const std::vector<VOff> &ioffsets,
                                      size_t max_size) {
  if (max_size <= 0)
    throw std::runtime_error("max_size must be > 0");
  if (ioffsets.size() < 2)
    throw std::runtime_error("Not enough BAI ioffsets to create ranges");

  // Interpret as "max slice size" (MiB) to satisfy the new requirement without
  // signature changes.
  const uint64_t limit_bytes = max_size;

  // Need BAM path to estimate decompressed bytes (block ISIZE). Infer from last
  // loaded BAI path.
  const std::string bam_path = _infer_bam_path_from_last_bai();

  // Build a light BGZF index starting at the first ioffset coff.
  const uint64_t coff_min = voff_coff(ioffsets.front());
  _BgzfLightIndex bg = _build_bgzf_light_index(bam_path, coff_min);

  // EOF virtual offset (for unmapped/no-coor tail)
  const VOff eof_voff = _compute_eof_voff_from_light_index(bg, coff_min);

  // We will create ranges using ioffsets as candidate cut points.
  // Overestimate sizes by counting full BGZF blocks touched.
  std::vector<std::pair<VOff, VOff>> ranges;
  ranges.reserve(1024);

  size_t start_idx = 0;
  while (start_idx + 1 < ioffsets.size()) {
    // If even the smallest next step exceeds limit, we still must make
    // progress. We'll choose the largest end_idx such that span <= limit (or
    // start_idx+1 if none).
    size_t lo = start_idx + 1;
    size_t hi = ioffsets.size() - 1;
    size_t best = lo;

    // If the limit is tiny, we may only be able to take one step.
    // Binary search for the farthest end within limit (monotonic in end_idx).
    while (lo <= hi) {
      size_t mid = lo + (hi - lo) / 2;

      uint64_t span = _overestimated_span_bytes_full_blocks(
          bg, ioffsets[start_idx], ioffsets[mid]);
      if (span <= limit_bytes) {
        best = mid;
        lo = mid + 1;
      } else {
        if (mid == 0)
          break;
        hi = mid - 1;
      }
    }

    // Ensure progress
    if (best <= start_idx)
      best = start_idx + 1;

    ranges.emplace_back(ioffsets[start_idx], ioffsets[best]);
    start_idx = best;
  }

  // If the last produced end does not reach EOF, append tail to capture
  // unmapped/no-coor.
  if (!ranges.empty()) {
    VOff last_end = ranges.back().second;
    if (voff_coff(last_end) < voff_coff(eof_voff) ||
        (voff_coff(last_end) == voff_coff(eof_voff) &&
         voff_uoff(last_end) < voff_uoff(eof_voff))) {
      ranges.emplace_back(last_end, eof_voff);
    }
  } else {
    // Degenerate: if ranges empty, make a single range to EOF
    ranges.emplace_back(ioffsets.front(), eof_voff);
  }

  // Remove empty/degenerate ranges
  ranges.erase(std::remove_if(ranges.begin(), ranges.end(),
                              [](auto &r) { return r.first >= r.second; }),
               ranges.end());

  if (ranges.empty())
    throw std::runtime_error("Range split produced no ranges");

  return ranges;
}

// ------------------------------------------------------------
// Existing estimator/printing (unchanged interface)
// Uses per-range scanning; fine for debug printing.
// ------------------------------------------------------------

// Over-estimated decompressed bytes in [beg_voff, end_voff):
// counts full BGZF blocks touched between coff(beg) and coff(end), inclusive of
// start/end blocks.
static uint64_t estimate_range_decompressed_bytes_over(FILE *f, VOff beg_voff,
                                                       VOff end_voff) {
  const uint64_t coff_b = voff_coff(beg_voff);
  const uint64_t coff_e = voff_coff(end_voff);

  if (coff_e < coff_b)
    throw std::runtime_error("range end before start");

#if defined(_WIN32)
  if (_fseeki64(f, (long long)coff_b, SEEK_SET) != 0)
    throw std::runtime_error("_fseeki64 failed");
#else
  if (fseeko(f, (off_t)coff_b, SEEK_SET) != 0)
    throw std::runtime_error("fseeko failed");
#endif

  uint64_t total = 0;

  while (true) {
#if defined(_WIN32)
    long long pos_ll = _ftelli64(f);
    if (pos_ll < 0)
      throw std::runtime_error("_ftelli64 failed");
    uint64_t coff_cur = uint64_t(pos_ll);
#else
    off_t pos = ftello(f);
    if (pos < 0)
      throw std::runtime_error("ftello failed");
    uint64_t coff_cur = uint64_t(pos);
#endif

    if (coff_cur > coff_e)
      break;

    auto [cbytes, ubytes] = read_bgzf_block_cbytes_ubytes(f);
    (void)cbytes;
    total += uint64_t(ubytes);

    if (coff_cur == coff_e)
      break;
  }

  return total;
}

static void print_ranges_estimated_decompressed_over(
    const std::string &bam_path,
    const std::vector<std::pair<VOff, VOff>> &ranges) {
  if (ranges.empty()) {
    std::cout << "[Ranges] none\n";
    return;
  }

  FILE *f = std::fopen(bam_path.c_str(), "rb");
  if (!f)
    throw std::runtime_error("Cannot open BAM: " + bam_path);

  auto mib = [](uint64_t b) { return double(b) / (1024.0 * 1024.0); };

  uint64_t total = 0;
  uint64_t min_u = UINT64_MAX, max_u = 0;

  std::cout << std::fixed << std::setprecision(2);
  std::cout << "[Ranges] count=" << ranges.size()
            << " (decompressed estimate = sum(ISIZE) of BGZF blocks touched; "
               "start/end blocks counted fully)\n";

  for (size_t i = 0; i < ranges.size(); ++i) {
    VOff beg = ranges[i].first;
    VOff end = ranges[i].second;

    uint64_t u = estimate_range_decompressed_bytes_over(f, beg, end);

    total += u;
    min_u = std::min(min_u, u);
    max_u = std::max(max_u, u);

    std::cout << "  #" << i << " est_u=" << mib(u) << " MiB"
              << " beg=(" << voff_coff(beg) << "," << voff_uoff(beg) << ")"
              << " end=(" << voff_coff(end) << "," << voff_uoff(end) << ")"
              << "\n";
  }

  double avg = double(total) / double(ranges.size());
  std::cout << "[Ranges] total=" << mib(total) << " MiB"
            << " avg=" << mib((uint64_t)avg) << " MiB"
            << " min=" << mib(min_u) << " MiB"
            << " max=" << mib(max_u) << " MiB" << std::endl;

  std::fclose(f);
}

#endif // BAIREADER_H