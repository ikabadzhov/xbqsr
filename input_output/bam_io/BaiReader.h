#pragma once
#include <algorithm>
#include <atomic>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fcntl.h>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#include <mutex>
#include <condition_variable>
#include <thread>

#include <errno.h>
#include <unistd.h>

#include <tbb/parallel_sort.h>

#include "../../logger/Logger.h"
#include "../../timer/Timer.h"

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
    VOff start{0};
    VOff end{0};
    ssize_t id{0};

    Slice() = default;
    Slice(VOff s, VOff e, ssize_t i) : start(s), end(e), id(i) {}
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
inline uint16_t rd_le16(const uint8_t *p) {
  return uint16_t(p[0]) | (uint16_t(p[1]) << 8);
}
inline uint32_t rd_le32(const uint8_t *p) {
  return uint32_t(p[0]) | (uint32_t(p[1]) << 8) | (uint32_t(p[2]) << 16) |
         (uint32_t(p[3]) << 24);
}

inline uint32_t read_u32(FILE *f) {
  uint32_t x;
  if (fread(&x, 4, 1, f) != 1)
    throw std::runtime_error("BAI read_u32 failed");
  return x;
}

inline uint64_t read_u64(FILE *f) {
  uint64_t x;
  if (fread(&x, 8, 1, f) != 1)
    throw std::runtime_error("BAI read_u64 failed");
  return x;
}

struct BaiLinearIndex {
  std::vector<VOff> ioffsets; // flattened across references
};

// Tracks the last BAI path we loaded, so the splitter can infer BAM path
// without changing signatures (your main calls load_bai_linear() first).
inline std::string &_last_loaded_bai_path() {
  static std::string s;
  return s;
}

inline BaiLinearIndex load_bai_linear(const std::string &bai_path) {
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

inline std::string find_bai_path(const std::string &bam_path) {
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
inline std::string _infer_bam_path_from_last_bai() {
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
  fs::path bam = p;
  bam.replace_extension(".bam");
  if (fs::exists(bam))
    return bam.string();

  throw std::runtime_error("Could not infer BAM path from BAI path: " +
                           bai.string());
}

// ------------------------------------------------------------
// BGZF parsing: read one block {cbytes, ubytes} and advance to next block.
// Robustly handles optional gzip fields to avoid desync.
// ------------------------------------------------------------
inline std::pair<uint32_t, uint32_t> read_bgzf_block_cbytes_ubytes(FILE *f) {
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

// ------------------------------------------------------------
// Minimal patch: robust tri-state header parser from memory
// ------------------------------------------------------------
enum class _ParseStatus : uint8_t { OK, NEED_MORE, CORRUPT };

struct _BlockInfo {
  uint32_t c_size = 0;
  uint32_t u_size = 0;
  _ParseStatus st = _ParseStatus::CORRUPT;
};

static inline _BlockInfo _parse_bgzf_header_mem_robust(const uint8_t* p, size_t remaining) {
  _BlockInfo out;

  // Need at least gzip(10) + XLEN(2)
  if (remaining < 12) {
    out.st = _ParseStatus::NEED_MORE;
    return out;
  }

  if (p[0] != 0x1f || p[1] != 0x8b || p[2] != 0x08) {
    out.st = _ParseStatus::CORRUPT;
    return out;
  }

  const uint8_t flg = p[3];
  if ((flg & 0x04) == 0) { // FEXTRA required for BGZF
    out.st = _ParseStatus::CORRUPT;
    return out;
  }

  // XLEN
  const uint16_t xlen = uint16_t(p[10]) | (uint16_t(p[11]) << 8);
  if (remaining < 12ull + xlen) {
    out.st = _ParseStatus::NEED_MORE;
    return out;
  }

  // Scan extra for "BC" + SLEN==2 to read BSIZE
  size_t off = 12;
  const size_t extra_end = 12ull + xlen;
  bool found_bc = false;
  uint16_t bsize = 0;

  while (off + 4 <= extra_end) {
    const uint8_t si1 = p[off + 0];
    const uint8_t si2 = p[off + 1];
    const uint16_t slen = uint16_t(p[off + 2]) | (uint16_t(p[off + 3]) << 8);
    off += 4;

    if (off + slen > extra_end) { // malformed extra
      out.st = _ParseStatus::CORRUPT;
      return out;
    }

    if (si1 == 'B' && si2 == 'C' && slen == 2) {
      bsize = uint16_t(p[off + 0]) | (uint16_t(p[off + 1]) << 8);
      found_bc = true;
      break;
    }
    off += slen;
  }

  if (!found_bc) {
    out.st = _ParseStatus::CORRUPT;
    return out;
  }

  const uint32_t c_size = uint32_t(bsize) + 1u;

  // Practical BGZF bounds: blocks are <= 64KiB.
  if (c_size < 26u || c_size > 65536u) {
    out.st = _ParseStatus::CORRUPT;
    return out;
  }

  if (remaining < c_size) {
    out.st = _ParseStatus::NEED_MORE;
    return out;
  }

  const uint32_t isize =
      uint32_t(p[c_size - 4]) |
      (uint32_t(p[c_size - 3]) << 8) |
      (uint32_t(p[c_size - 2]) << 16) |
      (uint32_t(p[c_size - 1]) << 24);

  if (isize > 65536u) {
    out.st = _ParseStatus::CORRUPT;
    return out;
  }

  out.c_size = c_size;
  out.u_size = isize;
  out.st = _ParseStatus::OK;
  return out;
}

inline _BgzfLightIndex _build_bgzf_light_index(const std::string& bam_path,
                                               uint64_t coff_start) {
  // ---- O_DIRECT constraints (typical) ----
  constexpr size_t kAlign = 4096;                 // common logical block size
  constexpr size_t kBufSize = (size_t)1 << 30;    // 1 GiB, divisible by 4096

  int fd = ::open(bam_path.c_str(), O_RDONLY | O_DIRECT);
  if (fd < 0) {
    throw std::runtime_error("Cannot open BAM with O_DIRECT: " + bam_path +
                             " (errno=" + std::to_string(errno) + ")");
  }

  // Align read offset down to kAlign for O_DIRECT.
  const uint64_t aligned_start = (coff_start / kAlign) * kAlign;
  const uint64_t skip_first = coff_start - aligned_start;

  // Allocate two aligned buffers
  void* buf0 = nullptr;
  void* buf1 = nullptr;
  if (::posix_memalign(&buf0, kAlign, kBufSize) != 0 || buf0 == nullptr) {
    ::close(fd);
    throw std::runtime_error("posix_memalign failed for buf0");
  }
  if (::posix_memalign(&buf1, kAlign, kBufSize) != 0 || buf1 == nullptr) {
    ::free(buf0);
    ::close(fd);
    throw std::runtime_error("posix_memalign failed for buf1");
  }

  struct Slot {
    uint8_t* data = nullptr;
    size_t   nbytes = 0;         // bytes actually read
    uint64_t file_off = 0;       // file offset corresponding to data[0]
    bool     ready = false;
    bool     eof = false;
  };

  Slot slots[2];
  slots[0].data = static_cast<uint8_t*>(buf0);
  slots[1].data = static_cast<uint8_t*>(buf1);

  std::mutex m;
  std::condition_variable cv_ready;
  std::condition_variable cv_free;
  bool slot_free[2] = {true, true};
  bool done = false;

  std::atomic<bool> loader_error{false};
  std::string loader_errmsg;

  uint64_t next_off = aligned_start;

  // Async loader thread: fills slots in alternating fashion
  std::thread loader([&] {
    int cur = 0;
    while (true) {
      // Wait for a free slot
      {
        std::unique_lock<std::mutex> lk(m);
        cv_free.wait(lk, [&]{ return done || slot_free[cur]; });
        if (done) return;
        slot_free[cur] = false;
      }

      slots[cur].file_off = next_off;

      // pread with O_DIRECT: buffer + size + offset must be aligned
      ssize_t r = ::pread(fd, slots[cur].data, kBufSize, (off_t)next_off);
      if (r < 0) {
        loader_error = true;
        loader_errmsg = "pread failed (errno=" + std::to_string(errno) + ")";
        {
          std::lock_guard<std::mutex> lk(m);
          slots[cur].nbytes = 0;
          slots[cur].eof = true;
          slots[cur].ready = true;
        }
        cv_ready.notify_one();
        return;
      }

      slots[cur].nbytes = (size_t)r;
      slots[cur].eof = (r == 0);
      slots[cur].ready = true;

      next_off += (uint64_t)r;

      cv_ready.notify_one();

      if (r == 0) return;

      cur ^= 1;
    }
  });

  _BgzfLightIndex idx;
  idx.coff.reserve(1 << 20);
  idx.ubytes.reserve(1 << 20);
  idx.prefix_ubytes.reserve(1 << 20);

  uint64_t prefix = 0;

  // ---- Minimal patch: robust carry handling across buffers ----
  std::vector<uint8_t> carry;
  carry.reserve(65536);
  uint64_t carry_off = 0; // file offset corresponding to carry[0]

  int consume = 0;
  bool first_chunk = true;

  while (true) {
    {
      std::unique_lock<std::mutex> lk(m);
      cv_ready.wait(lk, [&]{ return slots[consume].ready; });
    }

    if (loader_error) {
      {
        std::lock_guard<std::mutex> lk(m);
        done = true;
      }
      cv_free.notify_all();
      if (loader.joinable()) loader.join();
      ::free(buf0);
      ::free(buf1);
      ::close(fd);
      throw std::runtime_error("Async loader error: " + loader_errmsg);
    }

    const size_t n = slots[consume].nbytes;
    const uint64_t file_off = slots[consume].file_off;
    const bool eof = slots[consume].eof;

    const uint8_t* base = slots[consume].data;
    const uint8_t* ptr  = base;
    const uint8_t* end  = base + n;

    // Apply initial skip once (due to aligned_start)
    if (first_chunk && skip_first > 0) {
      if (skip_first >= n) {
        ptr = end;
      } else {
        ptr += skip_first;
      }
      first_chunk = false;
    } else {
      first_chunk = false;
    }

    // 1) If we have a carry (partial block), complete it robustly.
    // 1) If we have a carry (partial block), complete exactly ONE block.
    // IMPORTANT: never slurp arbitrary bytes into carry, or you'll drop data and desync.
    if (!carry.empty()) {
      while (true) {
        // Step A: ensure we can read XLEN (need 12 bytes)
        if (carry.size() < 12) {
          const size_t need = 12 - carry.size();
          const size_t take = std::min<size_t>(need, (size_t)(end - ptr));
          if (take == 0) break; // wait next chunk
          carry.insert(carry.end(), ptr, ptr + take);
          ptr += take;
          continue;
        }

        // Basic sanity (like B): if carry doesn't look like gzip+FEXTRA -> corruption/desync
        if (!(carry[0] == 0x1f && carry[1] == 0x8b && carry[2] == 0x08) ||
            ((carry[3] & 0x04) == 0)) {
          throw std::runtime_error("Corrupt BGZF stream at buffer boundary (bad magic/flags)");
        }

        const uint16_t xlen = uint16_t(carry[10]) | (uint16_t(carry[11]) << 8);
        const size_t need_extra = 12ull + xlen;
        if (need_extra > 65536ull) {
          throw std::runtime_error("Corrupt BGZF: XLEN implies >64KiB header");
        }

        // Step B: ensure extra is present
        if (carry.size() < need_extra) {
          const size_t need = need_extra - carry.size();
          const size_t take = std::min<size_t>(need, (size_t)(end - ptr));
          if (take == 0) break; // wait next chunk
          carry.insert(carry.end(), ptr, ptr + take);
          ptr += take;
          continue;
        }

        // Step C: now we can parse BC/BSIZE (might still need more for whole block)
        _BlockInfo h = _parse_bgzf_header_mem_robust(carry.data(), carry.size());
        if (h.st == _ParseStatus::CORRUPT) {
          throw std::runtime_error("Corrupt BGZF stream at buffer boundary (bad extra/BC/BSIZE)");
        }

        // If parser still says NEED_MORE here, it can only mean: need full block bytes (c_size)
        // But in that case h.c_size will still be 0 (by our struct init). So we must compute c_size
        // by parsing again after ensuring enough bytes — easiest: attempt to grow to 64KiB cap safely.
        if (h.st == _ParseStatus::NEED_MORE) {
          // We don't know c_size yet (parser didn't produce it), so append a small amount and retry.
          // This path is extremely rare; choose a small step to avoid slurping multiple blocks.
          const size_t step = 256; // tiny step; overhead negligible
          const size_t take = std::min<size_t>(step, (size_t)(end - ptr));
          if (take == 0) break;
          carry.insert(carry.end(), ptr, ptr + take);
          ptr += take;
          if (carry.size() > 65536u) {
            throw std::runtime_error("Carry exceeded BGZF max block size (desync/corrupt BAM?)");
          }
          continue;
        }

        // Step D: ensure we have the whole block (exactly)
        if (carry.size() < h.c_size) {
          const size_t need = (size_t)h.c_size - carry.size();
          const size_t take = std::min<size_t>(need, (size_t)(end - ptr));
          if (take == 0) break; // wait next chunk
          carry.insert(carry.end(), ptr, ptr + take);
          ptr += take;
          continue;
        }

        // Now carry contains >= full block; with the steps above it should be == full block.
        // If it's larger (shouldn't happen), keep leftover bytes (do NOT drop them).
        idx.coff.push_back(carry_off);
        idx.ubytes.push_back(h.u_size);
        idx.prefix_ubytes.push_back(prefix);
        prefix += uint64_t(h.u_size);

        if (carry.size() == h.c_size) {
          carry.clear();
        } else {
          // Keep leftover (rare, but safe)
          std::vector<uint8_t> leftover;
          leftover.assign(carry.begin() + h.c_size, carry.end());
          carry.swap(leftover);
          carry_off += uint64_t(h.c_size);
          // And loop again: the leftover might contain more blocks or partial bytes
          continue;
        }

        break; // carry cleared; proceed to parse current chunk from ptr
      }
    }

    // 2) Parse full blocks in current chunk.
    while (ptr < end) {
      const size_t rem = (size_t)(end - ptr);

      _BlockInfo h = _parse_bgzf_header_mem_robust(ptr, rem);
      if (h.st == _ParseStatus::OK) {
        const uint64_t coff = file_off + uint64_t(ptr - base);
        idx.coff.push_back(coff);
        idx.ubytes.push_back(h.u_size);
        idx.prefix_ubytes.push_back(prefix);
        prefix += uint64_t(h.u_size);

        ptr += h.c_size;
        continue;
      }

      if (h.st == _ParseStatus::NEED_MORE) {
        // Save tail as carry (partial header or partial block)
        carry_off = file_off + uint64_t(ptr - base);
        carry.assign(ptr, end);
        ptr = end;
        break;
      }

      // CORRUPT at this position means we lost sync (or BAM is corrupt).
      throw std::runtime_error("Corrupt BGZF stream while scanning (bad magic/extra/BC)");
    }

    // Mark slot free and move on
    {
      std::lock_guard<std::mutex> lk(m);
      slots[consume].ready = false;
      slots[consume].nbytes = 0;
      slots[consume].eof = false;
      slot_free[consume] = true;
    }
    cv_free.notify_one();

    if (eof) break;

    consume ^= 1;
  }

  // EOF: carry must be empty or a complete final block (rare). Otherwise truncated.
  if (!carry.empty()) {
    _BlockInfo h = _parse_bgzf_header_mem_robust(carry.data(), carry.size());
    if (h.st == _ParseStatus::OK && carry.size() == h.c_size) {
      idx.coff.push_back(carry_off);
      idx.ubytes.push_back(h.u_size);
      idx.prefix_ubytes.push_back(prefix);
      prefix += uint64_t(h.u_size);
      carry.clear();
    } else {
      ::free(buf0);
      ::free(buf1);
      {
        std::lock_guard<std::mutex> lk(m);
        done = true;
      }
      cv_free.notify_all();
      if (loader.joinable()) loader.join();
      ::close(fd);
      throw std::runtime_error("BAM appears truncated: incomplete BGZF block at EOF");
    }
  }

  // Stop loader if still running (normally it ends on EOF)
  {
    std::lock_guard<std::mutex> lk(m);
    done = true;
  }
  cv_free.notify_all();
  if (loader.joinable()) loader.join();

  ::free(buf0);
  ::free(buf1);
  ::close(fd);

  return idx;
}

inline VOff _compute_eof_voff_from_light_index(const _BgzfLightIndex &bg,
                                               uint64_t coff_start_global) {
  (void)coff_start_global;
  if (bg.coff.empty())
    throw std::runtime_error("BGZF index empty (scan produced no blocks)");
  size_t last = bg.coff.size() - 1;
  return make_voff(bg.coff[last],
                   bg.ubytes[last]); // end of last block (uoff = ISIZE)
}

// Snap helpers for over-estimating “full blocks touched” between two voffs.
inline size_t _snap_beg_block_idx(const _BgzfLightIndex &bg,
                                  uint64_t coff_beg) {
  auto it = std::upper_bound(bg.coff.begin(), bg.coff.end(), coff_beg);
  if (it == bg.coff.begin())
    return 0;
  return size_t((it - bg.coff.begin()) - 1);
}
inline size_t _snap_end_block_idx(const _BgzfLightIndex &bg,
                                  uint64_t coff_end) {
  auto it = std::lower_bound(bg.coff.begin(), bg.coff.end(), coff_end);
  if (it == bg.coff.end())
    return bg.coff.size() - 1;
  return size_t(it - bg.coff.begin());
}

inline uint64_t _overestimated_span_bytes_full_blocks(const _BgzfLightIndex &bg,
                                                      VOff beg, VOff end) {
  uint64_t cb = voff_coff(beg);
  uint64_t ce = voff_coff(end);
  if (ce < cb)
    throw std::runtime_error("range end before start");

  size_t ib = _snap_beg_block_idx(bg, cb);
  size_t ie = _snap_end_block_idx(bg, ce);

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
inline std::vector<std::pair<VOff, VOff>>
split_ranges_equal_by_compressed_span(const std::string bam_path, const std::vector<VOff> &ioffsets,
                                      size_t max_size) {
  if (max_size <= 0)
    throw std::runtime_error("max_size must be > 0");
  if (ioffsets.size() < 2)
    throw std::runtime_error("Not enough BAI ioffsets to create ranges");

  Timer t;
  t.reset();

  const uint64_t limit_bytes = max_size;

  const uint64_t coff_min = voff_coff(ioffsets.front());
  _BgzfLightIndex bg = _build_bgzf_light_index(bam_path, coff_min);

  Logger::getLogger().log(LOG_DEBUG,"_build_bgzf_light_index: " + std::to_string(t.elapsed()));
  t.reset();

  const VOff eof_voff = _compute_eof_voff_from_light_index(bg, coff_min);
  Logger::getLogger().log(LOG_DEBUG,"_compute_eof_voff_from_light_index: " + std::to_string(t.elapsed()));
  t.reset();

  std::vector<std::pair<VOff, VOff>> ranges;
  ranges.reserve(1024);

  size_t start_idx = 0;
  while (start_idx + 1 < ioffsets.size()) {
    size_t lo = start_idx + 1;
    size_t hi = ioffsets.size() - 1;
    size_t best = lo;

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

    if (best <= start_idx)
      best = start_idx + 1;

    ranges.emplace_back(ioffsets[start_idx], ioffsets[best]);
    start_idx = best;
  }

  if (!ranges.empty()) {
    VOff last_end = ranges.back().second;
    if (voff_coff(last_end) < voff_coff(eof_voff) ||
        (voff_coff(last_end) == voff_coff(eof_voff) &&
         voff_uoff(last_end) < voff_uoff(eof_voff))) {
      ranges.emplace_back(last_end, eof_voff);
    }
  } else {
    ranges.emplace_back(ioffsets.front(), eof_voff);
  }

  Logger::getLogger().log(LOG_DEBUG,"while loop: " + std::to_string(t.elapsed()));
  t.reset();

  ranges.erase(std::remove_if(ranges.begin(), ranges.end(),
                              [](auto &r) { return r.first >= r.second; }),
               ranges.end());

  if (ranges.empty())
    throw std::runtime_error("Range split produced no ranges");

  Logger::getLogger().log(LOG_DEBUG,"erase empty/degenerate ranges: " + std::to_string(t.elapsed()));
  t.reset();
  return ranges;
}

inline std::vector<std::pair<VOff, VOff>>
split_ranges_from_ioffsets_adjacent(const std::vector<VOff>& ioffsets) {
  if (ioffsets.size() < 2)
    throw std::runtime_error("Not enough BAI ioffsets to create ranges");

  std::vector<std::pair<VOff,VOff>> ranges;
  ranges.reserve(ioffsets.size() - 1);

  for (size_t i = 0; i + 1 < ioffsets.size(); ++i) {
    VOff a = ioffsets[i];
    VOff b = ioffsets[i + 1];
    if (a < b) ranges.emplace_back(a, b);
  }

  if (ranges.empty())
    throw std::runtime_error("Range split produced no ranges");
  return ranges;
}

// ------------------------------------------------------------
// Existing estimator/printing (unchanged interface)
// Uses per-range scanning; fine for debug printing.
// ------------------------------------------------------------
inline uint64_t estimate_range_decompressed_bytes_over(FILE *f, VOff beg_voff,
                                                       VOff end_voff) {
  const uint64_t coff_b = voff_coff(beg_voff);
  const uint64_t coff_e = voff_coff(end_voff);

  if (coff_e < coff_b)
    throw std::runtime_error("range end before start");

  if (fseeko(f, (off_t)coff_b, SEEK_SET) != 0)
    throw std::runtime_error("fseeko failed");

  uint64_t total = 0;

  while (true) {
    off_t pos = ftello(f);
    if (pos < 0)
      throw std::runtime_error("ftello failed");
    uint64_t coff_cur = uint64_t(pos);

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

inline void print_ranges_estimated_decompressed_over(
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