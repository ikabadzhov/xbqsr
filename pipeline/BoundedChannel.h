// StreamingChannel.h
#pragma once
#include <mutex>
#include <condition_variable>
#include <deque>
#include <optional>
#include <cstddef>

template<class T>
class BoundedChannel {
public:
    explicit BoundedChannel(size_t capacity) : capacity_(capacity) {}

    // Blocking push. Returns false if channel is closed.
    bool push(T&& v) {
        std::unique_lock lk(m_);
        cv_not_full_.wait(lk, [&]{ return closed_ || q_.size() < capacity_; });
        if (closed_) return false;
        q_.emplace_back(std::move(v));
        lk.unlock();
        cv_not_empty_.notify_one();
        return true;
    }

    // Non-blocking push. Returns false if closed or full.
    bool try_push(T&& v) {
        std::lock_guard lk(m_);
        if (closed_) return false;
        if (q_.size() >= capacity_) return false;
        q_.emplace_back(std::move(v));
        cv_not_empty_.notify_one();
        return true;
    }

    // Blocking pop. Returns nullopt when channel is closed and drained.
    std::optional<T> pop() {
        std::unique_lock lk(m_);
        cv_not_empty_.wait(lk, [&]{ return closed_ || !q_.empty(); });
        if (q_.empty()) return std::nullopt; // closed & drained
        T v = std::move(q_.front());
        q_.pop_front();
        lk.unlock();
        cv_not_full_.notify_one();
        return v;
    }

    // Non-blocking pop. Returns nullopt if empty (even if open),
    // and also returns nullopt if closed+empty.
    std::optional<T> try_pop() {
        std::lock_guard lk(m_);
        if (q_.empty()) return std::nullopt;
        T v = std::move(q_.front());
        q_.pop_front();
        cv_not_full_.notify_one();
        return v;
    }

    void close() {
        std::lock_guard lk(m_);
        closed_ = true;
        cv_not_full_.notify_all();
        cv_not_empty_.notify_all();
    }

    // Introspection helpers (safe; take the mutex)
    bool is_closed() const {
        std::lock_guard lk(m_);
        return closed_;
    }

    bool empty() const {
        std::lock_guard lk(m_);
        return q_.empty();
    }

    size_t size() const {
        std::lock_guard lk(m_);
        return q_.size();
    }

    size_t capacity() const noexcept { return capacity_; }

    // Optional: clear buffered items (only use when you own the channel and
    // are aborting; this drops items on the floor).
    void clear() {
        std::lock_guard lk(m_);
        q_.clear();
        cv_not_full_.notify_all();
    }

private:
    size_t capacity_;
    mutable std::mutex m_;
    mutable std::condition_variable cv_not_full_, cv_not_empty_;
    std::deque<T> q_;
    bool closed_{false};
};