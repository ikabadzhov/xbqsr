#pragma once
#include <cstdint>
#include <cstddef>
#include <span>
#include <string>
#include <vector>
#include <deque>
#include <unordered_map>
#include <memory>
#include <mutex>
#include <condition_variable>
#include <thread>
#include <system_error>

#include "OutputStream.h"          // your interface
#include <libdeflate.h>            // libdeflate_compressor, libdeflate_crc32

class BGZFOutputStream final : public OutputStream {
public:
    // compression_level: 1..12 (libdeflate range). threads: 1..16 (clamped)
    explicit BGZFOutputStream(const std::string& path,
                              int compression_level = 6,
                              unsigned threads = std::thread::hardware_concurrency());
    ~BGZFOutputStream() override;

    // OutputStream API
    size_t write(std::span<const uint8_t> data) override;
    void   flush() override;
    void   close() override;

private:
    // ======== BGZF constants ========
    // BGZF block total size must be <= 65536 (BSIZE field = total-1 fits uint16)
    // Header (gzip+FEXTRA+BGZF extra) = 18 bytes, Trailer (CRC32+ISIZE) = 8 bytes
    // So the max compressed payload we can store is 65536 - 18 - 8 = 65510
    static constexpr size_t BGZF_MAX_TOTAL   = 65536;
    static constexpr size_t BGZF_HDR_SIZE    = 18;
    static constexpr size_t BGZF_TRL_SIZE    = 8;
    static constexpr size_t BGZF_PAYLOAD_MAX = BGZF_MAX_TOTAL - BGZF_HDR_SIZE - BGZF_TRL_SIZE;

    // Target uncompressed size per block. Keep below PAYLOAD_MAX to avoid worst-case expansions.
    static constexpr size_t UN_TARGET_ = 65500; // a bit under 65510 for safety

    // ======== Task ========
    struct BgzfTask {
        uint64_t seq = 0;              // submission order
        bool poison = false;           // poison pill to stop worker

        std::vector<uint8_t> un;       // uncompressed bytes
        size_t un_len = 0;             // actual used bytes in 'un'

        std::vector<uint8_t> cz;       // compressed bytes
        size_t cz_len = 0;             // actual used bytes
        uint32_t crc = 0;              // CRC32 of uncompressed
    };

    struct CompDeleter {
        void operator()(libdeflate_compressor* p) const noexcept {
            if (p) libdeflate_free_compressor(p);
        }
    };

    // ======== Internals ========
    // file
    std::string path_;
    int fd_{-1};
    bool closed_{false};

    // compression level and threads
    int level_{6};
    std::vector<std::thread> workers_;

    // queues
    std::mutex m_;
    std::condition_variable cv_;
    std::deque<std::unique_ptr<BgzfTask>> q_in_;        // tasks to compress
    std::unordered_map<uint64_t, std::unique_ptr<BgzfTask>> out_map_; // completed, keyed by seq
    size_t in_flight_{0};                                // submitted not yet written
    bool producers_done_ = false;  // set true at close() after final enqueue

    // write order
    uint64_t next_seq_{0};       // sequence assigned to next task
    uint64_t next_to_write_{0};  // next sequence expected to be written

    // current assembling task
    std::unique_ptr<BgzfTask> cur_;

private:
    // helpers
    static inline void put_u16le(std::vector<uint8_t>& v, uint16_t x);
    static inline void put_u32le(std::vector<uint8_t>& v, uint32_t x);

    void start_workers_(unsigned n);
    void stop_workers_();
    void worker_loop_();

    void enqueue_full_task_();

    // build BGZF block (header + compressed payload + trailer) into 'out'
    void build_bgzf_block(const BgzfTask& t, std::vector<uint8_t>& out);

    // attempt to write any ready blocks in order (next_to_write_)
    void drain_and_write_ready_();

    // write all bytes to fd (handling partial writes/EINTR)
    size_t fd_write_all_(std::span<const uint8_t> data);

    // POSIX open helper
    static int open_write_trunc_(const std::string& path);
};