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
#include "../../pipeline/UnboundedChannel.h"
#include <libdeflate.h>
#include <unordered_map>

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
    static constexpr size_t UN_TARGET_       = 65500;

    struct BgzfTask {
        uint64_t seq = 0;
        bool poison = false;
        std::vector<uint8_t> un;
        size_t un_len = 0;
        std::vector<uint8_t> cz;
        size_t cz_len = 0;
        uint32_t crc = 0;
    };

    struct CompDeleter {
        void operator()(libdeflate_compressor* p) const noexcept {
            if (p) libdeflate_free_compressor(p);
        }
    };

    // ======== Status and Files ========
    std::string path_;
    int fd_{-1};
    bool closed_{false};
    int level_{6};

    // ======== Threading ========
    std::vector<std::jthread> workers_;
    std::jthread writer_thread_;

    // ======== Synchronization ========
    std::mutex m_;
    std::condition_variable cv_workers_;
    std::condition_variable cv_writer_;

    std::deque<std::unique_ptr<BgzfTask>> q_in_;
    std::unordered_map<uint64_t, std::unique_ptr<BgzfTask>> out_map_;

    bool workers_done_ = false;
    uint64_t next_seq_{0};
    uint64_t next_to_write_{0};

    // Current assembling task (producer side only)
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
    void worker_loop_(size_t id);
    void writer_loop_();

    void enqueue_full_task_();
    void build_and_write_block_(const BgzfTask& t);
    size_t fd_write_all_(std::span<const uint8_t> data);
};