// UnboundedChannel.h
#pragma once

#include <vector>
#include <functional>
#include <mutex>
#include <stdexcept>
#include <utility>
#include <cstddef>

// Unbounded "pool" channel:
// - Pre-fills with N reusable elements built by factory()
// - pop() NEVER blocks:
//    * if pool has an element -> returns it
//    * if empty and open      -> returns factory() (grows the pool)
//    * if empty and closed    -> throws (logic error / shutdown path)
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
        pool_.reserve(initial_pool);
        for (std::size_t i = 0; i < initial_pool; ++i) {
            pool_.emplace_back(factory_());
        }
    }

    UnboundedChannel(const UnboundedChannel&)            = delete;
    UnboundedChannel& operator=(const UnboundedChannel&) = delete;

    bool push(T&& v) {
        std::lock_guard<std::mutex> lk(m_);
        if (closed_) return false;
        pool_.emplace_back(std::move(v));
        return true;
    }

    // Never blocks. See class comment for behavior.
    T pop() {
        std::lock_guard<std::mutex> lk(m_);
        if (!pool_.empty()) {
            T v = std::move(pool_.back());
            pool_.pop_back();
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
        return pool_.size();
    }

private:
    mutable std::mutex m_;
    std::vector<T> pool_;
    Factory factory_;
    bool closed_{false};
};