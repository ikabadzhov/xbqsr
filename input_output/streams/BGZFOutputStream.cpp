#include "BGZFOutputStream.h"

#include <algorithm>
#include <cstring>
#include <fcntl.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>

// Canonical 28-byte BGZF EOF block
static const uint8_t BGZF_EOF[28] = {
    31,139,8,4,0,0,0,0,0,255,  6,0,  'B','C', 2,0,  27,0,  3,0,  0,0,0,0,  0,0,0,0
};
//                      XLEN     SI    SLEN     BSIZE   CRC    ISIZE
// BSIZE = total_block_size - 1 = 28 - 1 = 27

int BGZFOutputStream::open_write_trunc_(const std::string& path) {
    int fd = ::open(path.c_str(), O_CREAT | O_TRUNC | O_WRONLY | O_CLOEXEC, 0644);
    if (fd < 0) throw std::system_error(errno, std::generic_category(), "open(" + path + ")");
    return fd;
}

BGZFOutputStream::BGZFOutputStream(const std::string& path,
                                   int compression_level,
                                   unsigned threads)
    : path_(path),
      level_(std::clamp(compression_level, 1, 12))
{
    if (threads == 0) threads = 1;
    fd_ = open_write_trunc_(path_);
    start_workers_(std::min<unsigned>(threads, 16));

    cur_ = std::make_unique<BgzfTask>();
    cur_->seq = next_seq_++;
    cur_->un.resize(UN_TARGET_);
    cur_->un_len = 0;
}

BGZFOutputStream::~BGZFOutputStream() {
    try { close(); } catch (...) {}
}

inline void BGZFOutputStream::put_u16le(std::vector<uint8_t>& v, uint16_t x) {
    v.push_back(uint8_t(x & 0xFF));
    v.push_back(uint8_t(x >> 8));
}
inline void BGZFOutputStream::put_u32le(std::vector<uint8_t>& v, uint32_t x) {
    v.push_back(uint8_t(x & 0xFF));
    v.push_back(uint8_t((x >> 8) & 0xFF));
    v.push_back(uint8_t((x >> 16) & 0xFF));
    v.push_back(uint8_t((x >> 24) & 0xFF));
}

void BGZFOutputStream::build_bgzf_block(const BgzfTask& t, std::vector<uint8_t>& out) {
    const size_t total = BGZF_HDR_SIZE + t.cz_len + BGZF_TRL_SIZE;
    if (total > BGZF_MAX_TOTAL) {
        throw std::runtime_error("BGZF block > 64KiB (unexpected). Reduce UN_TARGET_.");
    }
    out.clear();
    out.reserve(total);

    // gzip header (ID1,ID2,CM=8, FLG=FEXTRA, MTIME=0, XFL=0, OS=255)
    out.push_back(0x1f); out.push_back(0x8b); out.push_back(0x08);
    out.push_back(0x04);
    out.push_back(0x00); out.push_back(0x00); out.push_back(0x00); out.push_back(0x00);
    out.push_back(0x00);
    out.push_back(0xff);

    // XLEN=6
    put_u16le(out, 6);
    // BGZF extra: 'BC', SLEN=2, BSIZE=total-1
    out.push_back('B'); out.push_back('C');
    put_u16le(out, 2);
    uint16_t bsize = static_cast<uint16_t>(total - 1);
    put_u16le(out, bsize);

    // payload (compressed)
    out.insert(out.end(), t.cz.begin(), t.cz.begin() + t.cz_len);

    // trailer: CRC32 (uncompressed), ISIZE (uncompressed size mod 2^32)
    put_u32le(out, t.crc);
    put_u32le(out, static_cast<uint32_t>(t.un_len & 0xFFFFFFFFu));
}

void BGZFOutputStream::start_workers_(unsigned n) {
    workers_.reserve(n);
    for (unsigned i = 0; i < n; ++i) {
        workers_.emplace_back(&BGZFOutputStream::worker_loop_, this);
    }
}

void BGZFOutputStream::stop_workers_() {
    // push poison tasks (one per worker)
    {
        std::lock_guard lk(m_);
        for (size_t i = 0; i < workers_.size(); ++i) {
            auto t = std::make_unique<BgzfTask>();
            t->poison = true;
            t->seq = UINT64_MAX;
            q_in_.push_back(std::move(t));
        }
    }
    cv_.notify_all();
    for (auto& th : workers_) th.join();
    workers_.clear();
}

void BGZFOutputStream::worker_loop_() {
    std::unique_ptr<libdeflate_compressor, CompDeleter> comp(
            libdeflate_alloc_compressor(level_));
    if (!comp) throw std::runtime_error("libdeflate_alloc_compressor failed");

    std::unique_ptr<BgzfTask> task;
    while (true) {
        {
            std::unique_lock lk(m_);
            cv_.wait(lk, [&]{ return !q_in_.empty(); });
            task = std::move(q_in_.front());
            q_in_.pop_front();
        }
        if (task->poison) break;

        // compress
        const size_t bound = libdeflate_deflate_compress_bound(comp.get(), task->un_len);
        task->cz.resize(bound);
        task->crc = libdeflate_crc32(0, task->un.data(), task->un_len);
        size_t got = libdeflate_deflate_compress(comp.get(),
                                                 task->un.data(), task->un_len,
                                                 task->cz.data(), task->cz.size());
        if (got == 0) throw std::runtime_error("libdeflate_deflate_compress failed");
        task->cz_len = got;

        // enqueue into out_map_ keyed by seq
        {
            std::lock_guard lk(m_);
            out_map_.emplace(task->seq, std::move(task));
        }
        cv_.notify_all();    // notify drainers
        drain_and_write_ready_();
    }
}

void BGZFOutputStream::drain_and_write_ready_() {
    static thread_local std::vector<uint8_t> block_bytes;

    for (;;) {
        std::unique_ptr<BgzfTask> next;
        {
            std::lock_guard lk(m_);
            auto it = out_map_.find(next_to_write_);
            if (it == out_map_.end()) break; // next seq not ready yet
            next = std::move(it->second);
            out_map_.erase(it);
        }
        build_bgzf_block(*next, block_bytes);
        fd_write_all_(std::span<const uint8_t>(block_bytes.data(), block_bytes.size()));
        ++next_to_write_;
    }
}

void BGZFOutputStream::enqueue_full_task_() {
    // Submit current block to workers
    {
        std::lock_guard lk(m_);
        ++in_flight_;
        q_in_.push_back(std::move(cur_));
    }
    cv_.notify_one();

    // Allocate a fresh current task
    cur_ = std::make_unique<BgzfTask>();
    cur_->seq = next_seq_++;
    cur_->un.resize(UN_TARGET_);
    cur_->un_len = 0;
}

size_t BGZFOutputStream::write(std::span<const uint8_t> data) {
    if (closed_) throw std::runtime_error("BGZFOutputStream::write on closed stream");
    size_t written = 0;
    const uint8_t* p = data.data();
    size_t n = data.size();

    while (n > 0) {
        size_t space = UN_TARGET_ - cur_->un_len;
        if (space == 0) {
            enqueue_full_task_();
            space = UN_TARGET_;
        }
        size_t take = std::min(space, n);
        std::memcpy(cur_->un.data() + cur_->un_len, p, take);
        cur_->un_len += take;
        p += take;
        n -= take;
        written += take;
    }
    return written;
}

void BGZFOutputStream::flush() {
    if (closed_) return;
    // Submit partial block if it has data
    if (cur_ && cur_->un_len > 0) {
        enqueue_full_task_();
    }

    // Wait until all submitted to workers
    {
        std::unique_lock lk(m_);
        cv_.wait(lk, [&]{ return q_in_.empty(); });
    }

    // Strong barrier: wait until all in-flight tasks are written
    for (;;) {
        {
            std::lock_guard lk(m_);
            if (in_flight_ == 0 && out_map_.empty()) break;
        }
        drain_and_write_ready_();
        std::this_thread::yield();
    }
}

void BGZFOutputStream::close() {
    if (closed_) return;

    // Submit any remaining uncompressed data as a last task
    if (cur_ && cur_->un_len > 0) enqueue_full_task_();

    // Stop workers: mark producers done and push poison pills
    {
        std::lock_guard lk(m_);
        producers_done_ = true;
        for (size_t i = 0; i < workers_.size(); ++i) {
            auto t = std::make_unique<BgzfTask>();
            t->poison = true;
            t->seq = UINT64_MAX;
            q_in_.push_back(std::move(t));
        }
    }
    cv_.notify_all();
    for (auto& th : workers_) th.join();
    workers_.clear();

    // Drain *all* remaining completed tasks strictly in-order
    drain_and_write_ready_();

    // Sanity: nothing left pending that could be written later
    {
        std::lock_guard lk(m_);
        if (!q_in_.empty() || !out_map_.empty()) {
            throw std::runtime_error("BGZFOutputStream::close: pending tasks remain (logic error)");
        }
    }

    // Write BGZF EOF sentinel and close
    fd_write_all_(std::span<const uint8_t>(BGZF_EOF, sizeof(BGZF_EOF)));
    ::fsync(fd_);
    ::close(fd_);
    fd_ = -1;
    closed_ = true;
}

size_t BGZFOutputStream::fd_write_all_(std::span<const uint8_t> data) {
    const uint8_t* p = data.data();
    size_t n = data.size();
    size_t total = 0;

    while (n > 0) {
        ssize_t w = ::write(fd_, p, n);
        if (w < 0) {
            if (errno == EINTR) continue;
            if (errno == EAGAIN || errno == EWOULDBLOCK) continue; // shouldn't happen on regular files
            throw std::system_error(errno, std::generic_category(), "write");
        }
        if (w == 0) { // should not happen on a regular file
            throw std::runtime_error("write returned 0");
        }
        p += static_cast<size_t>(w);
        n -= static_cast<size_t>(w);
        total += static_cast<size_t>(w);
    }
    return total;
}