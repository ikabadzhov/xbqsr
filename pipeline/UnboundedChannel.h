// UnboundedChannel.h
#pragma once

#include <deque>
#include <functional>
#include <mutex>
#include <stdexcept>
#include <utility>
#include <cstddef>

// Unbounded "pool" channel:
// - Pre-fills with N reusable elements built by factory()
// - pop() NEVER blocks:
//    * if queue has an element -> returns it
//    * if empty and open        -> returns factory() (grows the pool)
//    * if empty and closed      -> throws (logic error / shutdown path)
// - push() returns false if closed
template <class T>
class UnboundedChannel {
public:
    using Factory = std::function<T()>;

    UnboundedChannel() = default;

    UnboundedChannel(std::size_t initial_pool, Factory factory)
        : factory_(std::move(factory)) {
        if (!factory_) {
            throw std::runtime_error("UnboundedChannel: factory is null");
        }
        q_.clear();
        q_.shrink_to_fit();
        for (std::size_t i = 0; i < initial_pool; ++i) {
            q_.emplace_back(factory_());
        }
    }

    UnboundedChannel(const UnboundedChannel&)            = delete;
    UnboundedChannel& operator=(const UnboundedChannel&) = delete;

    bool push(T&& v) {
        std::lock_guard<std::mutex> lk(m_);
        if (closed_) return false;
        q_.emplace_back(std::move(v));
        return true;
    }

    // Never blocks. See class comment for behavior.
    T pop() {
        std::lock_guard<std::mutex> lk(m_);
        if (!q_.empty()) {
            T v = std::move(q_.front());
            q_.pop_front();
            return v;
        }
        if (closed_) {
            throw std::runtime_error("UnboundedChannel: pop() on closed channel");
        }
        // grow
        return factory_();
    }

    void close() {
        std::lock_guard<std::mutex> lk(m_);
        closed_ = true;
    }

    bool is_closed() const {
        std::lock_guard<std::mutex> lk(m_);
        return closed_;
    }

    std::size_t size() const {
        std::lock_guard<std::mutex> lk(m_);
        return q_.size();
    }

private:
    mutable std::mutex m_;
    std::deque<T> q_;
    Factory factory_;
    bool closed_{false};
};