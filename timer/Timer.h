#ifndef GENOMICPIPELINE_TIMER_H
#define GENOMICPIPELINE_TIMER_H


#include <chrono>
#include <iostream>

#include <chrono>

class Timer {
    using clock_t = std::chrono::high_resolution_clock;
    using second_t = std::chrono::duration<double, std::ratio<1>>;

    std::chrono::time_point<clock_t> m_beg;

public:
    Timer() : m_beg(clock_t::now()) {}

    void reset() {
        m_beg = clock_t::now();
    }

    double elapsed() const {
        return (std::chrono::duration_cast<second_t>(clock_t::now() - m_beg).count() < 0.0009) ? 0 : std::chrono::duration_cast<second_t>(clock_t::now() - m_beg).count();
    }

    long long elapsed_us() {
        return std::chrono::duration_cast<std::chrono::microseconds>(clock_t::now() - m_beg).count();
    }
};



#endif //GENOMICPIPELINE_TIMER_H
