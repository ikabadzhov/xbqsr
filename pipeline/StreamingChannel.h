// StreamingChannel.h
#pragma once
#include <mutex>
#include <condition_variable>
#include <deque>
#include <optional>

template<class T>
class BoundedChannel {
public:
    explicit BoundedChannel(size_t capacity) : capacity_(capacity) {}

    bool push(T&& v) {
        std::unique_lock lk(m_);
        cv_not_full_.wait(lk, [&]{ return closed_ || q_.size() < capacity_; });
        if (closed_) return false;
        q_.emplace_back(std::move(v));
        cv_not_empty_.notify_one();
        return true;
    }

    bool try_push(T&& v) {
        std::lock_guard lk(m_);
        if (closed_) return false;
        if (q_.size() >= capacity_) return false;
        q_.emplace_back(std::move(v));
        cv_not_empty_.notify_one();
        return true;
    }

    std::optional<T> pop() {
        std::unique_lock lk(m_);
        cv_not_empty_.wait(lk, [&]{ return closed_ || !q_.empty(); });
        if (q_.empty()) return std::nullopt; // closed & drained
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

private:
    size_t capacity_;
    std::mutex m_;
    std::condition_variable cv_not_full_, cv_not_empty_;
    std::deque<T> q_;
    bool closed_{false};
};