#include "BGZFOutputStream.h"

#include <algorithm>
#include <cstring>
#include <fcntl.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/uio.h>
#include "../../timer/TraceTimeline.h"
#include "../utils/Constant.h"

// Canonical 28-byte BGZF EOF block
static const uint8_t BGZF_EOF[28] = {
    31,139,8,4,0,0,0,0,0,255,  6,0,  'B','C', 2,0,  27,0,  3,0,  0,0,0,0,  0,0,0,0
};

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
    if (threads == 0) {
        unsigned hw = std::thread::hardware_concurrency();
        threads = (hw > 16) ? (hw - 16) : 1;
    }
    const unsigned nworkers = std::clamp<unsigned>(threads, 1, 64);

    fd_ = open_write_trunc_(path_);

    // Start workers
    workers_.reserve(nworkers);
    for (unsigned i = 0; i < nworkers; ++i) {
        workers_.emplace_back([this, i] { this->worker_loop_(i); });
    }

    // Start writer
    writer_thread_ = std::jthread([this] { this->writer_loop_(); });

    // Initialize first task
    cur_ = std::make_unique<BgzfTask>();
    cur_->seq = next_seq_++;
    cur_->un.resize(UN_TARGET_);
    cur_->un_len = 0;
}

BGZFOutputStream::~BGZFOutputStream() {
    try { close(); } catch (...) {}
}

void BGZFOutputStream::worker_loop_(size_t id) {
    BQS_TRACE_SET_LANE(trace_tl::Role::BgzfWorker, id);
    try {
        std::unique_ptr<libdeflate_compressor, CompDeleter> comp(
            libdeflate_alloc_compressor(level_));
        if (!comp) throw std::runtime_error("libdeflate_alloc_compressor failed");

        while (true) {
            std::unique_ptr<BgzfTask> task;
            {
                std::unique_lock lk(m_);
                cv_workers_.wait(lk, [&] { return !q_in_.empty(); });
                task = std::move(q_in_.front());
                q_in_.pop_front();
            }

            if (task->poison) break;

            const size_t bound = libdeflate_deflate_compress_bound(comp.get(), task->un_len);
            task->cz.resize(bound);
            task->crc = libdeflate_crc32(0, task->un.data(), task->un_len);

            size_t got = libdeflate_deflate_compress(comp.get(),
                                                     task->un.data(), task->un_len,
                                                     task->cz.data(), task->cz.size());
            if (got == 0) throw std::runtime_error("libdeflate_deflate_compress failed");
            task->cz_len = got;

            {
                std::lock_guard lk(m_);
                out_map_.emplace(task->seq, std::move(task));
            }
            cv_writer_.notify_one();
        }
    } catch (...) {
        // Simple error handling for now
    }
}

void BGZFOutputStream::writer_loop_() {
    BQS_TRACE_SET_LANE(trace_tl::Role::BgzfWriter, 0);
    try {
        while (true) {
            std::unique_ptr<BgzfTask> next;
            {
                std::unique_lock lk(m_);
                cv_writer_.wait(lk, [&] {
                    return out_map_.contains(next_to_write_) || workers_done_;
                });

                if (workers_done_ && out_map_.empty()) break;

                auto it = out_map_.find(next_to_write_);
                if (it == out_map_.end()) continue;

                next = std::move(it->second);
                out_map_.erase(it);
            }

            build_and_write_block_(*next);
            next_to_write_++;
        }
    } catch (...) {}
}

void BGZFOutputStream::build_and_write_block_(const BgzfTask& t) {
    const size_t total = BGZF_HDR_SIZE + t.cz_len + BGZF_TRL_SIZE;
    uint8_t hdr[BGZF_HDR_SIZE];
    uint8_t trl[BGZF_TRL_SIZE];

    hdr[0] = 0x1f; hdr[1] = 0x8b; hdr[2] = 0x08; hdr[3] = 0x04;
    hdr[4] = 0; hdr[5] = 0; hdr[6] = 0; hdr[7] = 0; hdr[8] = 0; hdr[9] = 0xff;
    put_u16le(&hdr[10], 6);
    hdr[12] = 'B'; hdr[13] = 'C';
    put_u16le(&hdr[14], 2);
    put_u16le(&hdr[16], static_cast<uint16_t>(total - 1));

    put_u32le(&trl[0], t.crc);
    put_u32le(&trl[4], static_cast<uint32_t>(t.un_len));

    fd_write_all_(std::span<const uint8_t>(hdr, BGZF_HDR_SIZE));
    fd_write_all_(std::span<const uint8_t>(t.cz.data(), t.cz_len));
    fd_write_all_(std::span<const uint8_t>(trl, BGZF_TRL_SIZE));
}

void BGZFOutputStream::enqueue_full_task_() {
    {
        std::lock_guard lk(m_);
        q_in_.push_back(std::move(cur_));
    }
    cv_workers_.notify_one();

    cur_ = std::make_unique<BgzfTask>();
    cur_->seq = next_seq_++;
    cur_->un.resize(UN_TARGET_);
    cur_->un_len = 0;
}

size_t BGZFOutputStream::write(std::span<const uint8_t> data) {
    if (closed_) throw std::runtime_error("BGZFOutputStream::write on closed stream");

    const uint8_t* p = data.data();
    size_t n = data.size();
    size_t written = 0;

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
    if (cur_ && cur_->un_len > 0) enqueue_full_task_();

    std::unique_lock lk(m_);
    cv_writer_.wait(lk, [&] {
        return q_in_.empty() && out_map_.empty();
    });
}

void BGZFOutputStream::close() {
    if (closed_) return;

    if (cur_ && cur_->un_len > 0) enqueue_full_task_();
    cur_.reset();

    {
        std::lock_guard lk(m_);
        for (size_t i = 0; i < workers_.size(); ++i) {
            auto t = std::make_unique<BgzfTask>();
            t->poison = true;
            q_in_.push_back(std::move(t));
        }
    }
    cv_workers_.notify_all();

    for (auto& th : workers_) th.join();
    workers_.clear();

    {
        std::lock_guard lk(m_);
        workers_done_ = true;
    }
    cv_writer_.notify_one();

    if (writer_thread_.joinable()) writer_thread_.join();

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
            throw std::system_error(errno, std::generic_category(), "write");
        }
        p += w; n -= w; total += w;
    }
    return total;
}