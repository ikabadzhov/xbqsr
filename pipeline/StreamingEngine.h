#pragma once
#include <thread>
#include <vector>
#include <functional>
#include <exception>
#include <stop_token>
#include "StreamingChannel.h"
#include "../sam_data_structures/SamAlignments.h"

class StreamingEngine {
public:
    using Batch   = SamAlignments;
    using StageFn = std::function<void(Batch&)>;
    using FinalizeFn = std::function<void()>;
    using clock      = std::chrono::steady_clock;

    struct Config {
        size_t channel_capacity{4};
        size_t workers{1}; // keep 1 if stages share GPU/global state
    };

    // A named stage with its own counters
    struct Stage {
        std::string name;
        StageFn     fn;
        std::atomic<long long> ns_total{0};
        std::atomic<long long> calls{0};

        // construct from name and function
        Stage(std::string n, StageFn f)
            : name(std::move(n)), fn(std::move(f)) {}

        // move-only semantics
        Stage(Stage&& other) noexcept
            : name(std::move(other.name)),
              fn(std::move(other.fn)) {
            ns_total.store(other.ns_total.load(std::memory_order_relaxed),
                           std::memory_order_relaxed);
            calls.store(other.calls.load(std::memory_order_relaxed),
                        std::memory_order_relaxed);
        }

        Stage& operator=(Stage&& other) noexcept {
            if (this != &other) {
                name = std::move(other.name);
                fn   = std::move(other.fn);
                ns_total.store(other.ns_total.load(std::memory_order_relaxed),
                               std::memory_order_relaxed);
                calls.store(other.calls.load(std::memory_order_relaxed),
                            std::memory_order_relaxed);
            }
            return *this;
        }

        Stage(const Stage&)            = delete;
        Stage& operator=(const Stage&) = delete;
    };

    size_t ns_producer = 0;
    size_t producer_ns_total = 0;
    size_t producer_ns_last = 0;
    size_t producer_batches = 0;

    StreamingEngine(Config cfg, std::vector<Stage>&& stages) : cfg_(cfg), stages_(std::move(stages)), chan_(cfg_.channel_capacity) {}

    StreamingEngine(Config cfg, std::vector<Stage>&& stages, std::vector<FinalizeFn>&& finals) : cfg_(cfg), stages_(std::move(stages)), finalizers_(std::move(finals)), chan_(cfg_.channel_capacity) {}

    void add_finalizer(FinalizeFn f) { finalizers_.push_back(std::move(f)); }

    template<class ReaderT>
    void run_streaming(ReaderT& reader, std::shared_ptr<Header> hdr) {
        // Producer: push batches from reader
        std::jthread producer([&](std::stop_token st){
            try {
                // Reader must call the callback once per parsed batch:
                //   void ParseAlignmentStreamed(hdr, [&](Batch&& b){ chan_.push(std::move(b)); });
                auto t_total0 = clock::now();
                auto t_prev   = t_total0;

                reader.ParseAlignmentStreamed(hdr, [&](Batch&& b){
                    if (st.stop_requested()) return;

                    auto t_now = clock::now();
                    auto dt    = std::chrono::duration_cast<std::chrono::nanoseconds>(t_now - t_prev).count();

                    producer_ns_total += dt; //.fetch_add(dt, std::memory_order_relaxed);
                    producer_ns_last += dt; //.store(dt, std::memory_order_relaxed);
                    producer_batches += 1; //.fetch_add(1, std::memory_order_relaxed);

                    t_prev = t_now;               // baseline for the next iteration
                    chan_.push(std::move(b));
                });
                auto t_end = clock::now();
                const auto dt = std::chrono::duration_cast<std::chrono::nanoseconds>(t_end - t_prev).count();
                ns_producer = dt;
            } catch (...) {
                prod_ex_ = std::current_exception();
                for (auto& t : workers_) t.request_stop();
            }
            chan_.close();
        });

        // Consumers: apply all stages to each batch
        for (size_t i=0; i<cfg_.workers; ++i) {
            workers_.emplace_back([&] (std::stop_token st){
                try {
                    for (;;) {
                        if (st.stop_requested()) break;
                        auto item = chan_.pop();
                        if (!item) break;  // closed & drained
                        Batch& batch = *item;
                        auto batch_t0 = clock::now();
                        for (auto& s : stages_) {
                            const auto t0 = clock::now();
                            s.fn(batch);
                            const auto t1 = clock::now();
                            const auto dt = std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count();
                            s.ns_total.fetch_add(dt, std::memory_order_relaxed);
                            s.calls.fetch_add(1, std::memory_order_relaxed);
                        }
                    }
                } catch (...) {
                    cons_ex_ = std::current_exception();
                    for (auto& t : workers_) t.request_stop();
                    chan_.close();
                }
            });
        }

        producer.join();
        for (auto& t : workers_) t.join();
        workers_.clear(); // destroys jthreads (already joined) and frees memory

        // If any stage/producer failed, rethrow first (skip finalizers by default)
        if (prod_ex_) std::rethrow_exception(prod_ex_);
        if (cons_ex_) std::rethrow_exception(cons_ex_);

        // Run finalizers once, on the caller thread, in order
        for (auto& f : finalizers_) {
            f();
        }
        const double total_ms = producer_ns_total / 1e9;     // nanoseconds → milliseconds
        const double avg_ms   = total_ms / std::max<size_t>(1, producer_batches);
        const double n        = producer_batches;

        std::fprintf(stderr,
            "[Streaming] %-20s : %10.3f sec total  (%8lld batches, avg %7.3f sec)\n",
            "Read BAM", total_ms, (long long)n, avg_ms);

        for (auto& s : stages_) {
            const auto ns  = s.ns_total.load(std::memory_order_relaxed);
            const auto n   = s.calls.load(std::memory_order_relaxed);
            const double ms = ns / 1e9;
            const double avg_ms = n ? (ms / n) : 0.0;
            // Logger::getLogger().log(LOG_INFO, s.name, ": ", ms, " ms total (", n, " calls, avg ", avg_ms, " ms)");
            std::fprintf(stderr, "[Streaming] %-20s : %10.3f sec total  (%8lld calls, avg %7.3f sec)\n",
                         s.name.c_str(), ms, (long long)n, avg_ms);
        }
    }

    template<class ReaderT>
    void run_streaming_serial(ReaderT& reader, std::shared_ptr<Header> hdr) {
        // no channel, no threads: process batch -> stages inline
        try {
            reader.ParseAlignmentStreamed(hdr, [&](Batch&& b){
                Batch batch = std::move(b);     // materialize the batch locally
                for (auto& s : stages_) {
                    s.fn(batch);
                }
            });
        } catch (...) {
            // propagate immediately; finalizers won’t run if producer fails
            throw;
        }

        // run finalizers after all batches
        for (auto& f : finalizers_) f();
    }

private:
    Config cfg_;
    std::vector<Stage> stages_;
    std::vector<FinalizeFn> finalizers_;
    BoundedChannel<Batch> chan_;
    std::vector<std::jthread> workers_;
    std::exception_ptr prod_ex_{}, cons_ex_{};
};