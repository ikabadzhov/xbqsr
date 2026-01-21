// StreamingEngine.h
#pragma once

#include <thread>
#include <vector>
#include <functional>
#include <exception>
#include <stop_token>
#include <atomic>
#include <cstdio>
#include <memory>
#include <string>
#include <optional>
#include <stdexcept>
#include <utility>
#include <cstdint>
#include <mutex>

#include "BoundedChannel.h"
#include "UnboundedChannel.h"
#include "AssociativeChannel.h"
#include "BoundedPool.h"

#include "../sam_data_structures/SamAlignments.h"
#include "../input_output/bam_io/BamWriter.h"
#include "../input_output/formatter/BamBatchFormatter.h"
#include "../timer/TraceTimeline.h"


/*
 - Producers parse and push batches to in_chan_
 - Workers pop, run stages, return batch to free_in pool
 - OPTIONAL ordered path:
      workers format into OutBatch from an unbounded output pool,
      publish OutBatch by id into AssociativeChannel,
      writer takes next_id and writes in-order, returns OutBatch to pool.
 - Clean shutdown on error: abort result channel, close input channel/pools.
*/
class StreamingEngine {
public:
    using Batch    = SamAlignments;
    using OutBatch = OutputBatchBam;

    using StageFn     = std::function<void(Batch&, size_t worker_id)>;
    using TransformFn = std::function<void(Batch&, OutBatch&)>;
    using FinalizeFn  = std::function<void()>;

    struct Config {
        size_t channel_capacity{4}; // producers -> workers
        size_t workers{1};

        // Pools (0 => auto)
        size_t in_pool_size{0};

        // Output pool sizing (only if ordered_writer != nullptr)
        size_t out_pool_initial{0};          // 0 => auto
        size_t out_initial_capacity{DEFAULT_BUFFER_SIZE};
    };

    struct Stage {
        std::string name;
        StageFn fn;

        Stage(std::string n, StageFn f) : name(std::move(n)), fn(std::move(f)) {}

        Stage(Stage&&) noexcept = default;
        Stage& operator=(Stage&&) noexcept = default;

        Stage(const Stage&)            = delete;
        Stage& operator=(const Stage&) = delete;
    };

    StreamingEngine(Config cfg,
                    std::vector<Stage>&& stages,
                    std::vector<FinalizeFn>&& finals,
                    std::shared_ptr<AlignmentWriter> ordered_writer = nullptr)
        : cfg_(cfg),
          stages_(std::move(stages)),
          finalizers_(std::move(finals)),
          in_chan_(cfg_.channel_capacity),
          ordered_writer_(std::move(ordered_writer)) {}

    void add_finalizer(FinalizeFn f) { finalizers_.push_back(std::move(f)); }

    template<class ReaderT>
    void run_streaming(std::vector<ReaderT>& readers, std::shared_ptr<Header> hdr)
    {
        const size_t P = readers.size();
        const size_t W = cfg_.workers;

        if (P == 0) throw std::runtime_error("run_streaming: readers is empty");
        if (W == 0) throw std::runtime_error("run_streaming: workers is 0");

        const bool ordered = (ordered_writer_ != nullptr);
        const size_t C     = cfg_.channel_capacity;

        // ---------- Pool sizing ----------
        const size_t in_pool_sz =
            (cfg_.in_pool_size != 0) ? cfg_.in_pool_size : (2*P + C);

        // ---------- Batch pool ----------
        BoundedPool<Batch> free_in(in_pool_sz, [&](){return Batch(hdr);}, P+C);

        // ---------- Ordered writer plumbing ----------
        // Output pool: unbounded pool of reusable OutBatch buffers
        // Results: associative channel id->OutBatch
        const size_t out_init =
                (cfg_.out_pool_initial != 0) ? cfg_.out_pool_initial : (W);

        UnboundedChannel<OutBatch> out_pool(ordered ? out_init : 0, [&]() {
                    OutBatch ob(0);
                    ob.reset();
                    return ob;
                });

        AssociativeChannel<OutBatch, std::uint64_t> results;

        if (ordered) {
            transforms_ = [this](Batch& batch, OutBatch& out) {
                formatter_.FormatAlignmentBatchTo(batch, out);
            };
        }

        // ---------- Shutdown coordination ----------
        std::atomic<bool> shutting_down{false};
        auto shutdown_all = [&] {
            bool expected = false;
            if (shutting_down.compare_exchange_strong(expected, true, std::memory_order_acq_rel)) {
                in_chan_.close();   // unblock workers
                free_in.close();    // unblock producers waiting for buffers
                if (ordered) {
                    out_pool.close(); // prevent returning buffers during teardown mistakes
                    results.close();  // unblock writer if it is waiting on missing ids
                }
            }
        };

        std::atomic<int> producers_left{ (int)P };
        std::atomic<int> workers_left{ (int)W };

        std::exception_ptr prod_ex;
        std::exception_ptr work_ex;
        std::exception_ptr writer_ex;

        std::mutex ex_mutex;

        // ---------- Writer thread (ordered only) ----------
        std::jthread writer_thread;
        if (ordered) {
            writer_thread = std::jthread([&](std::stop_token st) {
                BQS_TRACE_SET_LANE(trace_tl::Role::Writer, 0);
                try {

                    std::uint64_t next_id = 0;

                    for (;;) {
                        if (st.stop_requested()) break;
                        if (shutting_down.load(std::memory_order_acquire)) break;

                        auto item = [&]{
                            BQS_TRACE_SCOPE("writer:results.take(next)");
                            return results.take(next_id);
                        }();

                        if (!item) {
                            // If we did not shut down due to an error, then this is a real gap.
                            break;
                        }

                        OutBatch b = std::move(*item);

                        {
                            BQS_TRACE_SCOPE("writer:write");
                            ordered_writer_->WriteFormattedBatchDirect(b);
                        }

                        {
                            BQS_TRACE_SCOPE("writer:return_out_pool");
                            b.reset();
                            (void)out_pool.push(std::move(b));
                        }

                        ++next_id;
                    }
                } catch (...) {
                    {
                        std::lock_guard<std::mutex> lk(ex_mutex);
                        if (!writer_ex) writer_ex = std::current_exception();
                    }
                    results.abort(writer_ex);
                    shutdown_all();
                }
            });
        }

        // ---------- Workers ----------
        workers_.clear();
        workers_.reserve(W);

        for (size_t wi = 0; wi < W; ++wi) {
            workers_.emplace_back([&, wi](std::stop_token st) {
                BQS_TRACE_SET_LANE(trace_tl::Role::Worker, (uint32_t)wi);
                try {
                    for (;;) {
                        if (st.stop_requested()) {
                            shutdown_all();
                            break;
                        }

                        auto item = [&]{
                            BQS_TRACE_SCOPE("worker:in_chan.pop");
                            return in_chan_.pop();
                        }();

                        if (!item) break; // closed & drained

                        Batch batch = std::move(*item);

                        // Run stages
                        for (auto& s : stages_) {
                            BQS_TRACE_SCOPE(s.name.c_str());
                            s.fn(batch, wi);
                        }

                        if (!ordered) {
                            // Return batch to pool
                            BQS_TRACE_SCOPE("worker:free_in.push");
                            (void)free_in.push(std::move(batch));
                            continue;
                        }

                        // Ordered path: format + publish
                        const std::uint64_t id = (std::uint64_t)batch.BatchId();

                        // Acquire out buffer from unbounded pool
                        OutBatch out = [&]{
                            BQS_TRACE_SCOPE("worker:out_pool.pop");
                            return out_pool.pop();
                        }();

                        {
                            BQS_TRACE_SCOPE("worker:transform/format");
                            out.reset();
                            out.set_batch_id(id);
                            transforms_(batch, out);
                        }

                        // Return input batch to pool as soon as possible
                        {
                            BQS_TRACE_SCOPE("worker:return_batch_free_in");
                            (void)free_in.push(std::move(batch));
                        }

                        // Publish result by id
                        {
                            BQS_TRACE_SCOPE("worker:results.put");
                            if (!results.put(id, std::move(out))) {
                              shutdown_all(); break;
                            }
                        }
                    }
                } catch (...) {
                    {
                        std::lock_guard<std::mutex> lk(ex_mutex);
                        if (!work_ex) work_ex = std::current_exception();
                    }
                    results.abort(work_ex);
                    shutdown_all();
                }

                // last worker closes results so writer can stop once it reaches the end
                if (workers_left.fetch_sub(1, std::memory_order_acq_rel) == 1) {
                    if (ordered) results.close();
                    // close pool to unblock any producer stuck in free_in.pop()
                    free_in.close();
                }
            });
        }

        // ---------- Producers ----------
        std::vector<std::jthread> producers;
        producers.reserve(P);

        for (size_t pid = 0; pid < P; ++pid) {
            producers.emplace_back([&, pid](std::stop_token st) {
                BQS_TRACE_SET_LANE(trace_tl::Role::Producer, (uint32_t)pid);
                try {
                    readers[pid]->ParseAlignmentStreamed(
                        hdr,
                        [&]() -> Batch {
                            if (st.stop_requested()) { shutdown_all(); throw std::runtime_error("producer stopped"); }
                            BQS_TRACE_SCOPE("producer:free_in.pop");
                            auto b = free_in.pop();
                            if (!b) { shutdown_all(); throw std::runtime_error("free_in closed"); }
                            return std::move(*b);
                        },
                        [&](Batch&& b) {
                            if (st.stop_requested()) { shutdown_all(); return; }
                            BQS_TRACE_SCOPE("producer:in_chan.push");
                            if (!in_chan_.push(std::move(b))) { shutdown_all(); return; }
                        }
                    );
                } catch (...) {
                    {
                        std::lock_guard<std::mutex> lk(ex_mutex);
                        if (!prod_ex) prod_ex = std::current_exception();
                    }
                    if (ordered) results.abort(prod_ex);
                    shutdown_all();
                }

                if (producers_left.fetch_sub(1, std::memory_order_acq_rel) == 1) {
                    // last producer closes input channel so workers can drain and exit
                    in_chan_.close();
                }
            });
        }

        // ---------- Join ----------
        for (auto& p : producers) p.join();
        // Ensure workers can exit even if producers ended without closing (defensive)
        in_chan_.close();

        for (auto& w : workers_) w.join();
        workers_.clear();

        if (ordered) {
            // Ensure writer unblocks even if workers already ended
            results.close();
            //writer_thread.request_stop();
            writer_thread.join();
        }

        // ---------- Tear down ----------
        free_in.close();

        // Propagate exceptions (writer first if ordered)
        if (writer_ex) std::rethrow_exception(writer_ex);
        if (prod_ex)   std::rethrow_exception(prod_ex);
        if (work_ex)   std::rethrow_exception(work_ex);

        for (auto& f : finalizers_) f();

        BQS_TRACE_DUMP();
    }

private:
    Config cfg_;

    std::vector<Stage> stages_;
    TransformFn transforms_;
    std::vector<FinalizeFn> finalizers_;

    BoundedChannel<Batch> in_chan_;
    std::vector<std::jthread> workers_;
    std::shared_ptr<AlignmentWriter> ordered_writer_;

    BamBatchFormatter formatter_;
};