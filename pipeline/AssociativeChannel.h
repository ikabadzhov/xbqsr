// AssociativeChannel.h
#pragma once

#include <cstdint>
#include <condition_variable>
#include <mutex>
#include <optional>
#include <stdexcept>
#include <unordered_map>
#include <utility>

// A multi-producer / single-consumer (or even multi-consumer) channel where
// producers publish (id -> value) out-of-order, and the consumer can block
// waiting for a specific id.
//
// Semantics:
// - put(id, v): inserts exactly once. Duplicate id => throws (logic error).
// - take(id): blocks until either:
//     a) the value for id is available -> returns it
//     b) the channel is closed AND id will never arrive -> returns nullopt
// - close(): indicates no more items will be produced. Unblocks waiters.
// - abort(ex): like close(), but also stores an exception to be rethrown by take().
//
// Notes:
// - Memory grows with the number of out-of-order items retained.
//   If you want bounded memory, pair this with an explicit window/credit scheme
//   at the workers, or add an optional max_pending limit.
//
// Thread safety: all methods are thread-safe.
template <class T, class IdT = std::uint64_t>
class AssociativeChannel {
public:
    using id_type = IdT;
    using value_type = T;

    AssociativeChannel() = default;

    AssociativeChannel(const AssociativeChannel&) = delete;
    AssociativeChannel& operator=(const AssociativeChannel&) = delete;

    // Publish value for an id. Returns false if channel already closed/aborted.
    // Duplicate id => throws std::runtime_error.
    bool put(id_type id, T&& v) {
        std::unique_lock<std::mutex> lk(m_);

        if (closed_) return false;

        auto [it, inserted] = pending_.emplace(id, std::move(v));
        if (!inserted) {
            // Duplicate IDs are a hard bug under your assumptions.
            throw std::runtime_error("AssociativeChannel: duplicate id " + std::to_string((uint64_t)id));
        }

        lk.unlock();
        cv_.notify_all();
        return true;
    }

    bool put(id_type id, const T& v) {
        std::unique_lock<std::mutex> lk(m_);

        if (closed_) return false;

        auto [it, inserted] = pending_.emplace(id, v);
        if (!inserted) {
            throw std::runtime_error("AssociativeChannel: duplicate id " + std::to_string((uint64_t)id));
        }

        lk.unlock();
        cv_.notify_all();
        return true;
    }

    // Blocking take-by-id.
    // Returns std::nullopt if closed/aborted and the id is not present.
    std::optional<T> take(id_type id) {
        std::unique_lock<std::mutex> lk(m_);

        cv_.wait(lk, [&]{
            if (aborted_) return true;
            if (auto it = pending_.find(id); it != pending_.end()) return true;
            return closed_;
        });

        if (aborted_) {
            // If you want: return nullopt instead. But abort should fail fast.
            std::rethrow_exception(abort_ex_);
        }

        auto it = pending_.find(id);
        if (it == pending_.end()) {
            // closed_ and id not available => will never arrive
            return std::nullopt;
        }

        T out = std::move(it->second);
        pending_.erase(it);
        return out;
    }

    // Non-blocking check for a specific id.
    std::optional<T> try_take(id_type id) {
        std::lock_guard<std::mutex> lk(m_);
        if (aborted_) std::rethrow_exception(abort_ex_);

        auto it = pending_.find(id);
        if (it == pending_.end()) return std::nullopt;

        T out = std::move(it->second);
        pending_.erase(it);
        return out;
    }

    // Close: no more puts expected. Unblocks waiters.
    void close() {
        std::lock_guard<std::mutex> lk(m_);
        closed_ = true;
        cv_.notify_all();
    }

    // Abort: store an exception and unblock. take() will rethrow.
    void abort(std::exception_ptr ex) {
        if (!ex) ex = std::make_exception_ptr(std::runtime_error("AssociativeChannel: aborted"));
        std::lock_guard<std::mutex> lk(m_);
        abort_ex_ = ex;
        aborted_ = true;
        closed_ = true;
        cv_.notify_all();
    }

    // Introspection (optional)
    bool is_closed() const {
        std::lock_guard<std::mutex> lk(m_);
        return closed_;
    }

    std::size_t pending_size() const {
        std::lock_guard<std::mutex> lk(m_);
        return pending_.size();
    }

private:
    mutable std::mutex m_;
    std::condition_variable cv_;

    // Out-of-order storage
    std::unordered_map<id_type, T> pending_;

    bool closed_{false};

    bool aborted_{false};
    std::exception_ptr abort_ex_;
};