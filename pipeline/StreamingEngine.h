#pragma once
#include <thread>
#include <vector>
#include <functional>
#include <exception>
#include <stop_token>
#include <atomic>
#include <chrono>
#include <cstdio>
#include <memory>
#include <string>

#include "StreamingChannel.h"
#include "../sam_data_structures/SamAlignments.h"

class StreamingEngine {
public:
    using Batch      = SamAlignments;
    using StageFn    = std::function<void(Batch&)>;
    using FinalizeFn = std::function<void()>;
    using clock      = std::chrono::steady_clock;

    struct Config {
        size_t channel_capacity{4}; // bounded work queue
        //size_t pool_capacity{0};    // 0 => auto = channel_capacity + workers + producers
        size_t workers{1};          // consumers
    };

    struct Stage {
        std::string name;
        StageFn     fn;
        std::atomic<long long> ns_total{0};
        std::atomic<long long> calls{0};

        Stage(std::string n, StageFn f) : name(std::move(n)), fn(std::move(f)) {}

        Stage(Stage&& other) noexcept
            : name(std::move(other.name)), fn(std::move(other.fn)) {
            ns_total.store(other.ns_total.load(std::memory_order_relaxed), std::memory_order_relaxed);
            calls.store(other.calls.load(std::memory_order_relaxed), std::memory_order_relaxed);
        }

        Stage& operator=(Stage&& other) noexcept {
            if (this != &other) {
                name = std::move(other.name);
                fn   = std::move(other.fn);
                ns_total.store(other.ns_total.load(std::memory_order_relaxed), std::memory_order_relaxed);
                calls.store(other.calls.load(std::memory_order_relaxed), std::memory_order_relaxed);
            }
            return *this;
        }

        Stage(const Stage&)            = delete;
        Stage& operator=(const Stage&) = delete;
    };

    explicit StreamingEngine(Config cfg, std::vector<Stage>&& stages)
        : cfg_(cfg),
          stages_(std::move(stages)),
          chan_(cfg_.channel_capacity) {}

    StreamingEngine(Config cfg, std::vector<Stage>&& stages, std::vector<FinalizeFn>&& finals)
        : cfg_(cfg),
          stages_(std::move(stages)),
          finalizers_(std::move(finals)),
          chan_(cfg_.channel_capacity) {}

    void add_finalizer(FinalizeFn f) { finalizers_.push_back(std::move(f)); }

    // ReaderT is assumed pointer-like: readers[pid]->ParseAlignmentStreamed(...)
    template<class ReaderT>
    void run_streaming(std::vector<ReaderT>& readers, std::shared_ptr<Header> hdr) {
        const size_t P = readers.size();
        const size_t W = cfg_.workers;
        const size_t C = cfg_.channel_capacity;
        const size_t pool_sz = C + P + W;

        // Local fixed-size pool queue for this run.
        BoundedChannel<Batch> free_chan(pool_sz);

        // Prefill exactly pool_sz buffers (won't block).
        for (size_t i = 0; i < pool_sz; ++i) {
            const bool ok = free_chan.push(Batch(hdr));
            if (!ok) throw std::runtime_error("free pool closed during prefill");
        }

        // Reset errors for this run
        prod_ex_ = nullptr;
        cons_ex_ = nullptr;

        std::atomic<bool> shutting_down{false};
        auto shutdown_all = [&] {
            bool expected = false;
            if (shutting_down.compare_exchange_strong(expected, true, std::memory_order_acq_rel)) {
                chan_.close();
                free_chan.close();
            }
        };

        std::atomic<int> producers_left{ (int)P };

        // Stats (optional)
        std::atomic<long long> producer_ns_total_atomic{0};
        std::atomic<long long> producer_batches_atomic{0};

        // --- Consumers first ---
        workers_.clear();
        workers_.reserve(W);

        for (size_t wi = 0; wi < W; ++wi) {
            workers_.emplace_back([&, wi](std::stop_token st) {
                try {
                    for (;;) {
                        if (st.stop_requested()) { shutdown_all(); break; }

                        auto item = chan_.pop();
                        if (!item) break; // closed & drained

                        Batch batch = std::move(*item);

                        for (auto& s : stages_) {
                            const auto t0 = clock::now();
                            s.fn(batch);
                            const auto t1 = clock::now();
                            const auto dt = std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count();
                            s.ns_total.fetch_add(dt, std::memory_order_relaxed);
                            s.calls.fetch_add(1, std::memory_order_relaxed);
                        }

                        // CRITICAL: never block returning to pool (prevents circular deadlock).
                        (void)free_chan.try_push(std::move(batch)); // drop if full/closed (shouldn't happen)
                    }
                } catch (...) {
                    cons_ex_ = std::current_exception();
                    shutdown_all();
                }
            });
        }

        // --- Producers ---
        std::vector<std::jthread> producers;
        producers.reserve(P);

        for (size_t pid = 0; pid < P; ++pid) {
            producers.emplace_back([&, pid](std::stop_token st) {
                try {
                    auto t_prev = clock::now();

                    readers[pid]->ParseAlignmentStreamed(
                        hdr,
                        // Acquire a buffer from pool (blocking).
                        [&]() -> Batch {
                            if (st.stop_requested()) { shutdown_all(); throw std::runtime_error("producer stopped"); }
                            auto b = free_chan.pop();  // nullopt if closed & drained
                            if (!b) { shutdown_all(); throw std::runtime_error("free pool closed"); }
                            return std::move(*b);
                        },
                        // Emit to bounded work queue (blocking backpressure).
                        [&](Batch&& b) {
                            if (st.stop_requested()) { shutdown_all(); return; }

                            const auto t_now = clock::now();
                            const auto dt = std::chrono::duration_cast<std::chrono::nanoseconds>(t_now - t_prev).count();
                            producer_ns_total_atomic.fetch_add(dt, std::memory_order_relaxed);
                            producer_batches_atomic.fetch_add(1, std::memory_order_relaxed);
                            t_prev = t_now;

                            if (!chan_.push(std::move(b))) {
                                shutdown_all();
                                return;
                            }
                        }
                    );
                } catch (...) {
                    if (!prod_ex_) prod_ex_ = std::current_exception();
                    shutdown_all();
                }

                // Last producer closes work queue so consumers exit.
                if (producers_left.fetch_sub(1, std::memory_order_acq_rel) == 1) {
                    chan_.close();
                    free_chan.close();
                }
            });
        }

        // --- join ---
        for (auto& p : producers) p.join();
        for (auto& w : workers_) w.join();
        workers_.clear();

        // --- propagate errors ---
        if (prod_ex_) std::rethrow_exception(prod_ex_);
        if (cons_ex_) std::rethrow_exception(cons_ex_);

        // --- finalizers ---
        for (auto& f : finalizers_) f();

        // --- print stats ---
        const double total_s = (producer_ns_total_atomic.load(std::memory_order_relaxed) / 1e9) / producers.size();
        const auto   n_batches = producer_batches_atomic.load(std::memory_order_relaxed);
        const double avg_s = n_batches ? (total_s / n_batches) : 0.0;

        std::fprintf(stderr,
                     "[Streaming] %-20s: %10.3f sec total  (%8lld batches, avg %7.3f sec)\n",
                     "Read BAM (Avg Across Producers)", total_s, (long long)n_batches, avg_s);

        for (auto& s : stages_) {
            const auto ns  = s.ns_total.load(std::memory_order_relaxed);
            const auto n   = s.calls.load(std::memory_order_relaxed);
            const double sec = ns / 1e9;
            const double avg_sec = n ? (sec / n) : 0.0;
            std::fprintf(stderr,
                         "[Streaming] %-20s (Avg Across Workers): %10.3f sec total  (%8lld calls, avg %7.3f sec)\n",
                         s.name.c_str(), sec, (long long)n, avg_sec);
        }
    }

    template<class ReaderT>
    void run_streaming_serial(ReaderT& reader, std::shared_ptr<Header> hdr) {
        reader.ParseAlignmentStreamed(hdr, [&](Batch&& b){
            Batch batch = std::move(b);
            for (auto& s : stages_) s.fn(batch);
        });
        for (auto& f : finalizers_) f();
    }

private:
    Config cfg_;
    std::vector<Stage> stages_;
    std::vector<FinalizeFn> finalizers_;

    BoundedChannel<Batch> chan_;          // bounded work queue
    std::vector<std::jthread> workers_;  // consumers

    std::exception_ptr prod_ex_{}, cons_ex_{};
};