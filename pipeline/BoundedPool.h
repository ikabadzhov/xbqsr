// BoundedPool.h
#pragma once

#include <cstddef>
#include <condition_variable>
#include <functional>
#include <mutex>
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>

// BoundedPool<T> with BoundedChannel-like API.
//
// Think of it as a "bounded stack of reusable objects":
// - pop()/try_pop(): acquire an object.
//   * if one is free -> returns it
//   * else if created_ < capacity_ -> creates one via factory_ and returns it
//   * else -> blocks (pop) / returns nullopt (try_pop)
// - push()/try_push(): return an object to the pool.
//   * push blocks only if you choose to enforce "no more than created_ items returned"
//     (we don't need to block: pool can always accept returned objects).
//   * returns false if closed.
//
// Close semantics mirror BoundedChannel:
// - pop() returns nullopt when closed AND no free objects exist.
// - try_pop() returns nullopt when no free object and cannot grow (or closed).
template <class T>
class BoundedPool {
public:
    using Factory = std::function<T()>;

    BoundedPool() = default;

    explicit BoundedPool(std::size_t capacity, Factory factory, std::size_t prefill = 0)
        : capacity_(capacity), factory_(std::move(factory))
    {
        if (capacity_ == 0) throw std::runtime_error("BoundedPool: capacity must be > 0");
        if (!factory_)      throw std::runtime_error("BoundedPool: factory is null");
        if (prefill > capacity_) prefill = capacity_;

        free_.reserve(capacity_);
        for (std::size_t i = 0; i < prefill; ++i) {
            free_.push_back(factory_());
            ++created_;
        }
    }

    BoundedPool(const BoundedPool&)            = delete;
    BoundedPool& operator=(const BoundedPool&) = delete;

    // ------------------------------------------------------------------
    // Channel-like API
    // ------------------------------------------------------------------

    // Return an object to the pool. Never blocks. Returns false if closed.
    bool push(T&& v) {
        std::lock_guard<std::mutex> lk(m_);
        if (closed_) return false;
        free_.push_back(std::move(v));
        cv_.notify_one();
        return true;
    }

    // Same as push() (kept for symmetry with BoundedChannel).
    bool try_push(T&& v) { return push(std::move(v)); }

    // Acquire an object. Blocks only if:
    // - pool is open
    // - no free objects
    // - AND created_ == capacity_
    //
    // Returns nullopt when closed and no free objects exist.
    std::optional<T> pop() {
        // 1) Try fast path under lock: take free
        {
            std::unique_lock<std::mutex> lk(m_);
            if (!free_.empty()) {
                T v = std::move(free_.back());
                free_.pop_back();
                return v;
            }

            // 2) Grow if allowed (do factory outside lock)
            if (!closed_ && created_ < capacity_) {
                ++created_;
                lk.unlock();
                return factory_();
            }

            // 3) If closed and nothing free -> drained
            if (closed_) return std::nullopt;

            // 4) Otherwise wait for someone to push() a freed object, or close()
            cv_.wait(lk, [&] { return closed_ || !free_.empty(); });

            if (!free_.empty()) {
                T v = std::move(free_.back());
                free_.pop_back();
                return v;
            }
            // closed + drained
            return std::nullopt;
        }
    }

    // Non-blocking acquire:
    // - returns a free object if available
    // - else grows if created_ < capacity_
    // - else returns nullopt
    std::optional<T> try_pop() {
        std::unique_lock<std::mutex> lk(m_);

        if (!free_.empty()) {
            T v = std::move(free_.back());
            free_.pop_back();
            return v;
        }

        if (!closed_ && created_ < capacity_) {
            ++created_;
            lk.unlock();
            return factory_();
        }

        return std::nullopt;
    }

    void close() {
        std::lock_guard<std::mutex> lk(m_);
        closed_ = true;
        cv_.notify_all();
    }

    // Optional: drops all currently free objects (does not change created_).
    // Useful on abort paths if you want to release memory early.
    void clear() {
        std::lock_guard<std::mutex> lk(m_);
        free_.clear();
        // no notify needed; pop() waiters care about free_ becoming non-empty
    }

    // ------------------------------------------------------------------
    // Introspection (same style as BoundedChannel)
    // ------------------------------------------------------------------
    bool is_closed() const {
        std::lock_guard<std::mutex> lk(m_);
        return closed_;
    }

    std::size_t size() const {
        std::lock_guard<std::mutex> lk(m_);
        return free_.size(); // number of FREE objects currently available
    }

    std::size_t capacity() const noexcept { return capacity_; }

    // Total number of objects ever created (<= capacity_).
    std::size_t created() const {
        std::lock_guard<std::mutex> lk(m_);
        return created_;
    }

    // Convenience name matching your previous API
    std::size_t free_count() const { return size(); }

private:
    std::size_t capacity_{0};
    Factory factory_;

    mutable std::mutex m_;
    std::condition_variable cv_;

    std::vector<T> free_;     // LIFO stack of free objects
    std::size_t created_{0};  // ever-created count (in-use = created_ - free_.size())
    bool closed_{false};
};