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

#include "IExecutor.h"
#include "StreamingChannel.h"
#include "../sam_data_structures/SamAlignments.h"


class StreamingEngine {
public:
    using Batch      = SamAlignments;
    using StageFn    = std::function<void(Batch&, size_t worker_id)>;
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
    void run_streaming(std::vector<ReaderT>& readers,
                       std::vector<std::unique_ptr<IExecutor>>& executors,
                       std::shared_ptr<Header> hdr)
    {
        const size_t P = readers.size();
        const size_t W = executors.size();
        if (W == 0) throw std::runtime_error("run_streaming: executors is empty");
        if (P == 0) throw std::runtime_error("run_streaming: readers is empty");

        const size_t C = cfg_.channel_capacity;
        const size_t pool_sz = C + P + W;

        // Local fixed-size pool queue for this run.
        BoundedChannel<Batch> free_chan(pool_sz);

        for (size_t i = 0; i < pool_sz; ++i) {
            if (!free_chan.push(Batch(hdr)))
                throw std::runtime_error("free pool closed during prefill");
        }

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

        // ---- wall-clock timestamps (ns since steady epoch) ----
        std::vector<std::atomic<long long>> prod_start_ns(P), prod_end_ns(P);
        std::vector<std::atomic<long long>> work_start_ns(W), work_end_ns(W);
        for (size_t i = 0; i < P; ++i) { prod_start_ns[i] = 0; prod_end_ns[i] = 0; }
        for (size_t i = 0; i < W; ++i) { work_start_ns[i] = 0; work_end_ns[i] = 0; }

        auto now_ns = []() -> long long {
            return std::chrono::duration_cast<std::chrono::nanoseconds>(
                clock::now().time_since_epoch()).count();
        };

        // ---- Consumers (one per executor) ----
        workers_.clear();
        workers_.reserve(W);

        for (size_t wi = 0; wi < W; ++wi) {
            workers_.emplace_back([&, wi](std::stop_token st) {
                try {
                    work_start_ns[wi].store(now_ns(), std::memory_order_relaxed);

                    for (;;) {
                        if (st.stop_requested()) { shutdown_all(); break; }

                        auto item = chan_.pop();
                        if (!item) break; // closed & drained

                        Batch batch = std::move(*item);

                        // Generic device-local processing
                        executors[wi]->Process(batch);

                        // Never block returning to pool (avoid deadlocks)
                        (void)free_chan.try_push(std::move(batch));
                    }

                    // Make sure executor flushes device work before we stop timing
                    executors[wi]->Finalize();

                    work_end_ns[wi].store(now_ns(), std::memory_order_relaxed);
                } catch (...) {
                    cons_ex_ = std::current_exception();
                    shutdown_all();
                }
            });
        }

        // ---- Producers ----
        std::vector<std::jthread> producers;
        producers.reserve(P);

        for (size_t pid = 0; pid < P; ++pid) {
            producers.emplace_back([&, pid](std::stop_token st) {
                try {
                    prod_start_ns[pid].store(now_ns(), std::memory_order_relaxed);

                    readers[pid]->ParseAlignmentStreamed(
                        hdr,
                        [&]() -> Batch {
                            if (st.stop_requested()) { shutdown_all(); throw std::runtime_error("producer stopped"); }
                            auto b = free_chan.pop();
                            if (!b) { shutdown_all(); throw std::runtime_error("free pool closed"); }
                            return std::move(*b);
                        },
                        [&](Batch&& b) {
                            if (st.stop_requested()) { shutdown_all(); return; }
                            if (!chan_.push(std::move(b))) { shutdown_all(); return; }
                        }
                    );

                    prod_end_ns[pid].store(now_ns(), std::memory_order_relaxed);
                } catch (...) {
                    if (!prod_ex_) prod_ex_ = std::current_exception();
                    shutdown_all();
                }

                if (producers_left.fetch_sub(1, std::memory_order_acq_rel) == 1) {
                    chan_.close();
                    free_chan.close();
                }
            });
        }

        for (auto& p : producers) p.join();
        for (auto& w : workers_) w.join();
        workers_.clear();

        if (prod_ex_) std::rethrow_exception(prod_ex_);
        if (cons_ex_) std::rethrow_exception(cons_ex_);

        for (auto& f : finalizers_) f();

        // ---- wall-clock totals ----
        auto min_nonzero = [](const auto& v) -> long long {
            long long m = 0;
            for (auto& a : v) {
                long long x = a.load(std::memory_order_relaxed);
                if (x == 0) continue;
                if (m == 0 || x < m) m = x;
            }
            return m;
        };
        auto max_val = [](const auto& v) -> long long {
            long long m = 0;
            for (auto& a : v) {
                long long x = a.load(std::memory_order_relaxed);
                if (x > m) m = x;
            }
            return m;
        };

        const long long read_t0  = min_nonzero(prod_start_ns);
        const long long read_t1  = max_val(prod_end_ns);
        const long long comp_t0  = min_nonzero(work_start_ns);
        const long long comp_t1  = max_val(work_end_ns);

        const double read_wall_s = (read_t1 > read_t0) ? double(read_t1 - read_t0) / 1e9 : 0.0;
        const double comp_wall_s = (comp_t1 > comp_t0) ? double(comp_t1 - comp_t0) / 1e9 : 0.0;

        std::fprintf(stderr, "[Streaming] Read wall-clock   : %10.3f sec (min start -> max end)\n", read_wall_s);
        std::fprintf(stderr, "[Streaming] Compute wall-clock : %10.3f sec (min start -> max end)\n", comp_wall_s);
    }
    template<class ReaderT>
    void run_streaming_serial(ReaderT& reader, std::shared_ptr<Header> hdr) {
        reader.ParseAlignmentStreamed(hdr, [&](Batch&& b){
            Batch batch = std::move(b);
            for (auto& s : stages_); //s.fn(batch);
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