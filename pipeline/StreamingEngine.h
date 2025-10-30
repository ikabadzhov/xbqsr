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

    struct Config {
        size_t channel_capacity{4};
        size_t workers{1}; // keep 1 if stages share GPU/global state
    };

    StreamingEngine(Config cfg, std::vector<StageFn> stages)
        : cfg_(cfg), stages_(std::move(stages)), chan_(cfg_.channel_capacity) {}

    // ctor with finalizers
    StreamingEngine(Config cfg, std::vector<StageFn> stages, std::vector<FinalizeFn> finalizers)
        : cfg_(cfg), stages_(std::move(stages)), finalizers_(std::move(finalizers)), chan_(cfg_.channel_capacity) {}

    void add_finalizer(FinalizeFn f) { finalizers_.push_back(std::move(f)); }

    template<class ReaderT>
    void run_streaming(ReaderT& reader, std::shared_ptr<Header> hdr) {
        // Producer: push batches from reader
        std::jthread producer([&](std::stop_token st){
            try {
                // Reader must call the callback once per parsed batch:
                //   void ParseAlignmentStreamed(hdr, [&](Batch&& b){ chan_.push(std::move(b)); });
                reader.ParseAlignmentStreamed(hdr, [&](Batch&& b){
                    if (st.stop_requested()) return;
                    chan_.push(std::move(b));
                });
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
                        for (auto& s : stages_) s(batch);
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
        for (auto& f : finalizers_) f();
    }

    template<class ReaderT>
    void run_streaming_serial(ReaderT& reader, std::shared_ptr<Header> hdr) {
        // no channel, no threads: process batch -> stages inline
        try {
            reader.ParseAlignmentStreamed(hdr, [&](Batch&& b){
                Batch batch = std::move(b);     // materialize the batch locally
                for (auto& s : stages_) s(batch);
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
    std::vector<StageFn> stages_;
    std::vector<FinalizeFn> finalizers_;
    BoundedChannel<Batch> chan_;
    std::vector<std::jthread> workers_;
    std::exception_ptr prod_ex_{}, cons_ex_{};
};