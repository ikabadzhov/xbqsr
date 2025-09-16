#pragma once
#include <chrono>
#include <array>
#include <iostream>
#include <iomanip>
#include <string_view>

enum class StageId : size_t {
    AddReadGroup = 0,
    Transform,
    ComputeSkipVector,
    ComputeSNPSVector,
    ComputeContextVector,
    UpdateIndexes,
    Count // keep last → number of stages
};

static constexpr std::array<std::string_view, (size_t)StageId::Count> StageNames = {
    "Add Read Group",
    "Transform",
    "Compute Skip Vector",
    "Compute SNPS Vector",
    "Compute Context Vector",
    "Update Indexes"
};

class ProfilerBQSR {
    using Clock = std::chrono::high_resolution_clock;

    struct Stage {
        uint64_t acc_ns = 0;
        Clock::time_point start_time;
        bool running = false;
    };

    std::array<Stage, (size_t)StageId::Count> stages_{};

    // Private constructor → only accessible through instance()
    ProfilerBQSR() = default;

public:
    // Delete copy/move
    ProfilerBQSR(const ProfilerBQSR&) = delete;
    ProfilerBQSR& operator=(const ProfilerBQSR&) = delete;
    ProfilerBQSR(ProfilerBQSR&&) = delete;
    ProfilerBQSR& operator=(ProfilerBQSR&&) = delete;

    static ProfilerBQSR& instance() {
        static ProfilerBQSR instance_; // guaranteed thread-safe in C++11+
        return instance_;
    }

    void start(StageId id) {
        Stage& s = stages_[(size_t)id];
        if (!s.running) {
            s.start_time = Clock::now();
            s.running = true;
        }
    }

    void stop(StageId id) {
        Stage& s = stages_[(size_t)id];
        if (s.running) {
            auto end = Clock::now();
            s.acc_ns += std::chrono::duration_cast<std::chrono::nanoseconds>(end - s.start_time).count();
            s.running = false;
        }
    }

    void reset() {
        for (auto& s : stages_) {
            s.acc_ns = 0;
            s.running = false;
        }
    }

    void report(std::ostream& os = std::cout) const {
        os << "\n=== BQSR Profiler Report ===\n";
        os << std::left << std::setw(24) << "Stage" << std::right << std::setw(14) << "Time (sec)\n";
        os << "----------------------------------------------\n";
        for (size_t i = 0; i < (size_t)StageId::Count; i++) {
            os << std::left << std::setw(24) << StageNames[i]
               << std::right << std::setw(14) << std::fixed << std::setprecision(3)
               << (double) stages_[i].acc_ns / 1e9 << "\n";
        }
        os << "==============================================\n";
    }
};