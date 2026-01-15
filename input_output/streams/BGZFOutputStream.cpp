// BGZFOutputStream.cpp
#include "BGZFOutputStream.h"

#include <algorithm>
#include <cstring>
#include <fcntl.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/uio.h>
#include "../../timer/TraceTimeline.h"

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
    if (threads == 0) threads = 1;
    const unsigned nworkers = std::min<unsigned>(threads, 100);

    fd_ = open_write_trunc_(path_);

    out_ring_.resize(OUT_WINDOW_);
    for (size_t i = 0; i < OUT_WINDOW_; ++i) {
        out_ring_[i].seq = static_cast<uint64_t>(i);
        out_ring_[i].ready = false;
    }

    cur_ = std::make_unique<BgzfTask>();
    cur_->seq = next_seq_.fetch_add(1, std::memory_order_relaxed);
    cur_->un.resize(UN_TARGET_);
    cur_->un_len = 0;

    start_threads_(nworkers);
}

BGZFOutputStream::~BGZFOutputStream() {
    try { close(); } catch (...) {}
}

void BGZFOutputStream::set_exception_(std::exception_ptr e) {
    {
        std::lock_guard lk(ex_m_);
        if (!ex_) ex_ = e;
    }
    stop_.store(true, std::memory_order_release);
    in_cv_.notify_all();
    out_cv_.notify_all();
    done_cv_.notify_all();
}

void BGZFOutputStream::rethrow_if_exception_() {
    std::exception_ptr e;
    {
        std::lock_guard lk(ex_m_);
        e = ex_;
    }
    if (e) std::rethrow_exception(e);
}

void BGZFOutputStream::start_threads_(unsigned nworkers) {
    writer_ = std::jthread(&BGZFOutputStream::writer_loop_, this);
    workers_.reserve(nworkers);
    for (unsigned i = 0; i < nworkers; ++i) {
        workers_.emplace_back([this, i]{
            this->worker_loop_(i);
        });
    }
}

void BGZFOutputStream::enqueue_poison_pills_(size_t n) {
    for (size_t i = 0; i < n; ++i) {
        auto t = std::make_unique<BgzfTask>();
        t->poison = true;
        t->seq = UINT64_MAX;
        q_in_.push_back(std::move(t));
    }
}

void BGZFOutputStream::stop_threads_() {
    if (writer_.joinable()) writer_.request_stop();
    for (auto& th : workers_) th.request_stop();

    in_cv_.notify_all();
    out_cv_.notify_all();
    done_cv_.notify_all();

    for (auto& th : workers_) {
        if (th.joinable()) th.join();
    }
    workers_.clear();

    if (writer_.joinable()) writer_.join();
}

void BGZFOutputStream::worker_loop_(size_t worker_id) {
    BQS_TRACE_SET_LANE(trace_tl::Role::BgzfWorker, worker_id);
    try {
        std::unique_ptr<libdeflate_compressor, CompDeleter> comp(
            libdeflate_alloc_compressor(level_));
        if (!comp) throw std::runtime_error("libdeflate_alloc_compressor failed");

        std::unique_ptr<BgzfTask> task;

        for (;;) {
            {
                BQS_TRACE_SCOPE("bgzf_worker:wait_in_q");
                std::unique_lock lk(in_m_);
                in_cv_.wait(lk, [&]{
                    return stop_.load(std::memory_order_acquire) || !q_in_.empty();
                });
                if (stop_.load(std::memory_order_acquire)) return;

                task = std::move(q_in_.front());
                q_in_.pop_front();
            }

            if (task->poison) return;
            {
                BQS_TRACE_SCOPE("bgzf_worker:compress_total");

                const size_t bound = libdeflate_deflate_compress_bound(comp.get(), task->un_len);
                task->cz.resize(bound);

                task->crc = libdeflate_crc32(0, task->un.data(), task->un_len);

                const size_t got = libdeflate_deflate_compress(
                    comp.get(),
                    task->un.data(), task->un_len,
                    task->cz.data(), task->cz.size());
                if (got == 0) throw std::runtime_error("libdeflate_deflate_compress failed");
                task->cz_len = got;

            }

            {
                BQS_TRACE_SCOPE("bgzf_worker:publish_to_out_ring");

                std::unique_lock lk(out_m_);
                if (stop_.load(std::memory_order_acquire)) return;

                const uint64_t s = task->seq;
                const size_t idx = static_cast<size_t>(s % OUT_WINDOW_);
                auto& slot = out_ring_[idx];

                out_cv_.wait(lk, [&]{
                    return stop_.load(std::memory_order_acquire) || (slot.seq == s && !slot.ready);
                });

                if (stop_.load(std::memory_order_acquire)) return;

                slot.task = std::move(task);
                slot.ready = true;
            }

            out_cv_.notify_one();
        }
    } catch (...) {
        set_exception_(std::current_exception());
    }
}

void BGZFOutputStream::write_bgzf_task_(const BgzfTask& t) {
    BQS_TRACE_SCOPE("bgzf_writer:write_bgzf_task_total");

    const size_t total = BGZF_HDR_SIZE + t.cz_len + BGZF_TRL_SIZE;
    if (total > BGZF_MAX_TOTAL) {
        throw std::runtime_error("BGZF block > 64KiB (unexpected). Reduce UN_TARGET_.");
    }

    uint8_t hdr[BGZF_HDR_SIZE];
    uint8_t trl[BGZF_TRL_SIZE];

    hdr[0] = 0x1f; hdr[1] = 0x8b; hdr[2] = 0x08;
    hdr[3] = 0x04;
    hdr[4] = 0x00; hdr[5] = 0x00; hdr[6] = 0x00; hdr[7] = 0x00;
    hdr[8] = 0x00;
    hdr[9] = 0xff;

    put_u16le(&hdr[10], 6);
    hdr[12] = 'B'; hdr[13] = 'C';
    put_u16le(&hdr[14], 2);
    put_u16le(&hdr[16], static_cast<uint16_t>(total - 1));

    put_u32le(&trl[0], t.crc);
    put_u32le(&trl[4], static_cast<uint32_t>(t.un_len & 0xFFFFFFFFu));

    iovec iov[3];
    iov[0].iov_base = hdr;
    iov[0].iov_len  = BGZF_HDR_SIZE;
    iov[1].iov_base = const_cast<uint8_t*>(t.cz.data());
    iov[1].iov_len  = t.cz_len;
    iov[2].iov_base = trl;
    iov[2].iov_len  = BGZF_TRL_SIZE;

    fd_writev_all_(iov, 3);

}

void BGZFOutputStream::writer_loop_() {
    BQS_TRACE_SET_LANE(trace_tl::Role::BgzfWriter, 0);
    try {
        for (;;) {
            std::unique_ptr<BgzfTask> next;

            {
                BQS_TRACE_SCOPE("bgzf_writer:wait_next_block");

                std::unique_lock lk(out_m_);
                out_cv_.wait(lk, [&]{
                    if (stop_.load(std::memory_order_acquire)) return true;

                    const uint64_t want = next_to_write_.load(std::memory_order_relaxed);
                    const size_t idx = static_cast<size_t>(want % OUT_WINDOW_);
                    auto& slot = out_ring_[idx];
                    if (slot.seq == want && slot.ready) return true;

                    if (producers_done_.load(std::memory_order_acquire) &&
                        in_flight_.load(std::memory_order_acquire) == 0) return true;

                    return false;
                });

                if (stop_.load(std::memory_order_acquire)) return;
                if (producers_done_.load(std::memory_order_acquire) &&
                    in_flight_.load(std::memory_order_acquire) == 0) return;

                const uint64_t want = next_to_write_.load(std::memory_order_relaxed);
                const size_t idx = static_cast<size_t>(want % OUT_WINDOW_);
                auto& slot = out_ring_[idx];

                if (!(slot.seq == want && slot.ready)) continue;

                next = std::move(slot.task);
                slot.task.reset();
                slot.ready = false;
                slot.seq += OUT_WINDOW_;

                next_to_write_.store(want + 1, std::memory_order_relaxed);
            }

            write_bgzf_task_(*next);


            {
                BQS_TRACE_SCOPE("bgzf_writer:in_flight--");
                const size_t prev = in_flight_.fetch_sub(1, std::memory_order_acq_rel);
                if (prev == 1) {
                    std::lock_guard lk(done_m_);
                    done_cv_.notify_all();
                }
            }

            out_cv_.notify_all();
        }
    } catch (...) {
        set_exception_(std::current_exception());
    }
}

void BGZFOutputStream::enqueue_full_task_() {
    BQS_TRACE_SCOPE("bgzf:enqueue_full_task_total");

    if (stop_.load(std::memory_order_acquire)) return;

    in_flight_.fetch_add(1, std::memory_order_acq_rel);

    {
        BQS_TRACE_SCOPE("bgzf:push_q_in_lock");
        std::lock_guard lk(in_m_);
        q_in_.push_back(std::move(cur_));
    }

    in_cv_.notify_one();

    {
        BQS_TRACE_SCOPE("bgzf:alloc_new_task");
        cur_ = std::make_unique<BgzfTask>();
        cur_->seq = next_seq_.fetch_add(1, std::memory_order_relaxed);
        cur_->un.resize(UN_TARGET_);
        cur_->un_len = 0;
    }
}

size_t BGZFOutputStream::write(std::span<const uint8_t> data) {
    BQS_TRACE_SCOPE("bgzf:write_total");

    if (closed_) throw std::runtime_error("BGZFOutputStream::write on closed stream");
    rethrow_if_exception_();

    size_t written = 0;
    const uint8_t* p = data.data();
    size_t n = data.size();

    while (n > 0) {
        size_t space = UN_TARGET_ - cur_->un_len;
        if (space == 0) {
            enqueue_full_task_();
            space = UN_TARGET_;
        }
        const size_t take = std::min(space, n);
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

    if (cur_ && cur_->un_len > 0) {
        enqueue_full_task_();
    }

    std::unique_lock lk(done_m_);
    done_cv_.wait(lk, [&]{
        return stop_.load(std::memory_order_acquire) ||
               in_flight_.load(std::memory_order_acquire) == 0;
    });

    rethrow_if_exception_();
}

void BGZFOutputStream::close() {
    if (closed_) return;

    if (cur_ && cur_->un_len > 0) {
        enqueue_full_task_();
    }
    cur_.reset();

    producers_done_.store(true, std::memory_order_release);
    {
        std::lock_guard lk(in_m_);
        enqueue_poison_pills_(workers_.size());
    }
    in_cv_.notify_all();
    out_cv_.notify_all();

    {
        std::unique_lock lk(done_m_);
        done_cv_.wait(lk, [&]{
            return stop_.load(std::memory_order_acquire) ||
                   in_flight_.load(std::memory_order_acquire) == 0;
        });
    }

    stop_threads_();
    rethrow_if_exception_();

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
        const ssize_t w = ::write(fd_, p, n);
        if (w < 0) {
            if (errno == EINTR) continue;
            if (errno == EAGAIN || errno == EWOULDBLOCK) continue;
            throw std::system_error(errno, std::generic_category(), "write");
        }
        if (w == 0) throw std::runtime_error("write returned 0");
        p += static_cast<size_t>(w);
        n -= static_cast<size_t>(w);
        total += static_cast<size_t>(w);
    }
    return total;
}

void BGZFOutputStream::fd_writev_all_(struct iovec* iov, int iovcnt) {
    BQS_TRACE_SCOPE("bgzf_writer:fd_writev_all_total");

    int i = 0;
    while (i < iovcnt) {
        while (i < iovcnt && iov[i].iov_len == 0) ++i;
        if (i >= iovcnt) return;

        const ssize_t w = ::writev(fd_, &iov[i], iovcnt - i);
        if (w < 0) {
            if (errno == EINTR) continue;
            if (errno == EAGAIN || errno == EWOULDBLOCK) continue;
            throw std::system_error(errno, std::generic_category(), "writev");
        }
        if (w == 0) throw std::runtime_error("writev returned 0");

        size_t remaining = static_cast<size_t>(w);
        while (remaining > 0 && i < iovcnt) {
            const size_t len = iov[i].iov_len;
            if (remaining >= len) {
                remaining -= len;
                iov[i].iov_len = 0;
                ++i;
            } else {
                iov[i].iov_base = static_cast<char*>(iov[i].iov_base) + remaining;
                iov[i].iov_len -= remaining;
                remaining = 0;
            }
        }
    }
}