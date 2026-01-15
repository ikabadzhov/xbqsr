// BGZFOutputStream.h
#pragma once
#include <cstdint>
#include <cstddef>
#include <span>
#include <string>
#include <vector>
#include <deque>
#include <memory>
#include <mutex>
#include <condition_variable>
#include <thread>
#include <system_error>
#include <exception>
#include <atomic>

#include "OutputStream.h"
#include <libdeflate.h>

class BGZFOutputStream final : public OutputStream {
public:
    // compression_level: 1..12 (libdeflate range). threads: 1..16 (clamped)
    explicit BGZFOutputStream(const std::string& path,
                              int compression_level = 6,
                              unsigned threads = std::thread::hardware_concurrency());
    ~BGZFOutputStream() override;

    size_t write(std::span<const uint8_t> data) override;
    void   flush() override;
    void   close() override;

private:
    // ======== BGZF constants ========
    static constexpr size_t BGZF_MAX_TOTAL   = 65536;
    static constexpr size_t BGZF_HDR_SIZE    = 18;
    static constexpr size_t BGZF_TRL_SIZE    = 8;
    static constexpr size_t BGZF_PAYLOAD_MAX = BGZF_MAX_TOTAL - BGZF_HDR_SIZE - BGZF_TRL_SIZE;

    static constexpr size_t UN_TARGET_ = 65500; // a bit under 65510 for safety

    // Big ring so modulo collisions are practically nonexistent.
    static constexpr size_t OUT_WINDOW_ = 1u << 16; // 65536

    struct BgzfTask {
        uint64_t seq = 0;
        bool poison = false;

        std::vector<uint8_t> un;
        size_t un_len = 0;

        std::vector<uint8_t> cz;
        size_t cz_len = 0;

        uint32_t crc = 0;
    };

    struct OutSlot {
        uint64_t seq = 0;
        bool ready = false;
        std::unique_ptr<BgzfTask> task;
    };

    struct CompDeleter {
        void operator()(libdeflate_compressor* p) const noexcept {
            if (p) libdeflate_free_compressor(p);
        }
    };

    // ======== file ========
    std::string path_;
    int fd_{-1};
    bool closed_{false};

    // ======== compression ========
    int level_{6};
    std::vector<std::jthread> workers_;
    std::jthread writer_;

    // ======== input queue ========
    std::mutex in_m_;
    std::condition_variable in_cv_;
    std::deque<std::unique_ptr<BgzfTask>> q_in_;

    // ======== output ring ========
    std::mutex out_m_;
    std::condition_variable out_cv_;
    std::vector<OutSlot> out_ring_;

    // ======== global state ========
    std::atomic<bool> stop_{false};
    std::atomic<bool> producers_done_{false};
    std::atomic<uint64_t> next_seq_{0};
    std::atomic<uint64_t> next_to_write_{0};
    std::atomic<size_t> in_flight_{0};

    std::mutex ex_m_;
    std::exception_ptr ex_{nullptr};

    std::mutex done_m_;
    std::condition_variable done_cv_;

    // current assembling (producer thread only)
    std::unique_ptr<BgzfTask> cur_;

private:
    static inline void put_u16le(uint8_t* p, uint16_t x) noexcept {
        p[0] = uint8_t(x & 0xFF);
        p[1] = uint8_t(x >> 8);
    }
    static inline void put_u32le(uint8_t* p, uint32_t x) noexcept {
        p[0] = uint8_t(x & 0xFF);
        p[1] = uint8_t((x >> 8) & 0xFF);
        p[2] = uint8_t((x >> 16) & 0xFF);
        p[3] = uint8_t((x >> 24) & 0xFF);
    }

    static int open_write_trunc_(const std::string& path);

    void start_threads_(unsigned nworkers);
    void stop_threads_();

    void worker_loop_(size_t worker_id);
    void writer_loop_();

    void enqueue_full_task_();
    void enqueue_poison_pills_(size_t n);

    void write_bgzf_task_(const BgzfTask& t);

    size_t fd_write_all_(std::span<const uint8_t> data);
    void   fd_writev_all_(struct iovec* iov, int iovcnt);

    void set_exception_(std::exception_ptr e);
    void rethrow_if_exception_();
};