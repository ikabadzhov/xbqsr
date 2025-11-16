#pragma once
#include <chrono>
#include <array>
#include <cstdint>
#include <iostream>
#include <iomanip>
#include <string_view>

// Turn profiling on/off globally
#ifndef ENABLE_IO_PROFILER
#define ENABLE_IO_PROFILER 1
#endif

enum class IOStageId : std::size_t {
    DataLoading = 0,
    ProcessBatch,
    EnsureCapacity,
    GetBoundaries,
    ParseSimpleTypes,
    ParseComplexTypes,
    Scan,
    Allocation,
    ParallelFor,
    Count
};

static constexpr std::array<std::string_view, (std::size_t)IOStageId::Count> IOStageNames = {
    "Data Loading",
    "Process Batch",
    "\tEnsure Capacity",
    "\tGet Boundaries",
    "\tParse Simple Types",
    "\tParse Complex Types",
    "\t\tScan",
    "\t\tAllocation",
    "\t\tParallelFor"
};

//=========================================================
// PROFILER CORE
//=========================================================
class ProfilerIO {
public:
    using Clock = std::chrono::steady_clock;

private:
    struct Stage {
        std::uint64_t acc_ns = 0;
        Clock::time_point start_time{};
        bool running = false;
    };

    std::array<Stage, (std::size_t)IOStageId::Count> stages_{};

    ProfilerIO() = default;

public:
    ProfilerIO(const ProfilerIO&) = delete;
    ProfilerIO& operator=(const ProfilerIO&) = delete;

    static ProfilerIO& instance() {
        static ProfilerIO inst;
        return inst;
    }

#if ENABLE_IO_PROFILER
    //-------------- EXPLICIT START/STOP -------------------
    inline void start(IOStageId id) noexcept {
        Stage& s = stages_[(std::size_t)id];
        if (!s.running) {
            s.start_time = Clock::now();
            s.running = true;
        }
    }

    inline void stop(IOStageId id) noexcept {
        Stage& s = stages_[(std::size_t)id];
        if (s.running) {
            using namespace std::chrono;
            auto ns = duration_cast<nanoseconds>(Clock::now() - s.start_time).count();
            s.acc_ns += static_cast<std::uint64_t>(ns);
            s.running = false;
        }
    }

    //-------------- FAST ADD-TIME (RAII USES THIS) --------
    inline void add_time(IOStageId id, std::uint64_t ns) noexcept {
        stages_[(std::size_t)id].acc_ns += ns;
    }

#else
    // Everything compiles-out to no-ops
    inline void start(IOStageId) noexcept {}
    inline void stop(IOStageId) noexcept {}
    inline void add_time(IOStageId, std::uint64_t) noexcept {}
#endif

    //-------------- RESET ----------------------------------
    inline void reset() noexcept {
        for (auto& s : stages_) {
            s.acc_ns = 0;
            s.running = false;
        }
    }

    //-------------- REPORT ---------------------------------
    void report(std::ostream& os = std::cout) const {
        os << "\n=== IO Profiler Report ===\n";
        os << std::left << std::setw(24) << "Stage"
           << std::right << std::setw(14) << "Time (sec)\n";
        os << "----------------------------------------------\n";

        for (std::size_t i = 0; i < (std::size_t)IOStageId::Count; ++i) {
            double sec = double(stages_[i].acc_ns) * 1e-9;
            os << std::left << std::setw(24) << IOStageNames[i]
               << std::right << std::setw(14)
               << std::fixed << std::setprecision(3)
               << sec << "\n";
        }

        os << "==============================================\n";
    }
};

//=========================================================
// START/STOP MACROS
//=========================================================

#if ENABLE_IO_PROFILER

#define IO_PROFILE_START(StageEnum) \
ProfilerIO::instance().start(IOStageId::StageEnum)

#define IO_PROFILE_STOP(StageEnum) \
ProfilerIO::instance().stop(IOStageId::StageEnum)

#else

#define IO_PROFILE_START(StageEnum)  ((void)0)
#define IO_PROFILE_STOP(StageEnum)   ((void)0)

#endif

//=========================================================
// RAII SCOPED TIMER
//=========================================================
#if ENABLE_IO_PROFILER
struct IOStageScope {
    IOStageId stage;
    ProfilerIO::Clock::time_point start;

    explicit IOStageScope(IOStageId id) noexcept
        : stage(id), start(ProfilerIO::Clock::now()) {}

    ~IOStageScope() noexcept {
        using namespace std::chrono;
        auto end = ProfilerIO::Clock::now();
        auto ns  = duration_cast<nanoseconds>(end - start).count();
        ProfilerIO::instance().add_time(stage,
                                        static_cast<std::uint64_t>(ns));
    }
};

#define IO_PROFILE_SCOPE(StageEnum) \
    IOStageScope _io_scope_##__LINE__{ IOStageId::StageEnum }

#else
// Profiling disabled → macro does nothing
#define IO_PROFILE_SCOPE(StageEnum) ((void)0)


#endif

