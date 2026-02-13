#pragma once
#include <cstdint>
#include <chrono>
#include <cstdlib>
#include <fstream>
#include <mutex>
#include <string>
#include <string_view>
#include <vector>
#include <algorithm>

#ifndef BQS_TRACE_COMPILED
#define BQS_TRACE_COMPILED 1
#endif

namespace trace_tl {

enum class Role : uint32_t {
    Producer=1, Worker=2, Writer=3,
    BgzfWorker=4, BgzfWriter=5,
    Reader=6, Loader=7, Decompress=8
};

#if BQS_TRACE_COMPILED

// ================= REAL IMPLEMENTATION =================


inline bool enabled() {
    static int v = -1;
    if (v == -1) {
        const char* e = std::getenv("BQS_TRACE");
        v = (e && *e && std::string_view(e) != "0") ? 1 : 0;
    }
    return v == 1;
}

inline std::string output_path() {
    const char* p = std::getenv("BQS_TRACE_PATH");
    return (p && *p) ? std::string(p) : std::string("trace.json");
}

inline uint64_t steady_us() {
    using namespace std::chrono;
    return (uint64_t)duration_cast<microseconds>(steady_clock::now().time_since_epoch()).count();
}
inline uint64_t now_us_rel() {
    static const uint64_t t0 = steady_us();
    return steady_us() - t0;
}

struct Event        { const char* name; uint32_t pid; uint32_t tid; uint64_t ts_us; uint64_t dur_us; };
struct InstantEvent { const char* name; uint32_t pid; uint32_t tid; uint64_t ts_us; };
struct CounterEvent { const char* name; uint32_t pid; uint32_t tid; uint64_t ts_us; int64_t value; };

class Timeline {
public:
    static Timeline& instance() { static Timeline t; return t; }

    void set_group(uint32_t pid_group, uint32_t tid_in_group) {
        if (!enabled()) return;
        tls_pid() = pid_group;
        tls_tid() = tid_in_group;
        ensure_registered();
    }

    void set_lane(Role role, uint32_t index) {
        if (!enabled()) return;
        uint32_t pid = 1;
        uint32_t tid = make_lane(role, index, pid);
        tls_pid() = pid;
        tls_tid() = tid;
        ensure_registered();
    }

    uint32_t pid()  const { return tls_pid(); }
    uint32_t lane() const { return tls_tid(); }

    void add(const char* name, uint64_t ts_us, uint64_t dur_us) {
        if (!enabled()) return;
        local_events().push_back({name, tls_pid(), tls_tid(), ts_us, dur_us});
    }
    void instant(const char* name, uint64_t ts_us) {
        if (!enabled()) return;
        local_inst().push_back({name, tls_pid(), tls_tid(), ts_us});
    }
    void counter(const char* name, uint64_t ts_us, int64_t v) {
        if (!enabled()) return;
        local_ctr().push_back({name, tls_pid(), tls_tid(), ts_us, v});
    }

    void dump() {
        if (!enabled()) return;

        std::lock_guard<std::mutex> lk(mu_);
        std::ofstream out(output_path(), std::ios::out | std::ios::trunc);

        out << "[\n";
        bool first = true;

        // collect pids
        std::vector<uint32_t> pids;
        auto add_pid = [&](uint32_t pid) {
            if (std::find(pids.begin(), pids.end(), pid) == pids.end())
                pids.push_back(pid);
        };
        for (auto& vec : all_)  for (auto& e : vec) add_pid(e.pid);
        for (auto& vec : inst_) for (auto& e : vec) add_pid(e.pid);
        for (auto& vec : ctr_)  for (auto& e : vec) add_pid(e.pid);

        // metadata first
        for (uint32_t pid : pids) {
            std::string proc_name;
            if (pid >= 100 && pid < 200) proc_name = "Reader " + std::to_string(pid - 100);
            else if (pid == 200)         proc_name = "Workers";
            else if (pid == 300)         proc_name = "Writer";
            else                         proc_name = "PID " + std::to_string(pid);

            int sort_index = 0;
            if (pid >= 100 && pid < 200) sort_index = (int)(pid - 100);
            else if (pid == 200)         sort_index = 1000;
            else if (pid == 300)         sort_index = 1100;
            else                         sort_index = 2000;

            if (!first) out << ",\n";
            first = false;
            out << "{\"ph\":\"M\",\"name\":\"process_name\",\"pid\":" << pid
                << ",\"tid\":0,\"args\":{\"name\":\"" << escape(proc_name.c_str()) << "\"}}";
            out << ",\n";
            out << "{\"ph\":\"M\",\"name\":\"process_sort_index\",\"pid\":" << pid
                << ",\"tid\":0,\"args\":{\"sort_index\":" << sort_index << "}}";
        }

        // events
        for (auto& vec : all_) for (auto& e : vec) {
            out << ",\n{\"name\":\"" << escape(e.name) << "\",\"ph\":\"X\",\"pid\":" << e.pid
                << ",\"tid\":" << e.tid << ",\"ts\":" << e.ts_us << ",\"dur\":" << e.dur_us << "}";
        }
        for (auto& vec : inst_) for (auto& e : vec) {
            out << ",\n{\"name\":\"" << escape(e.name) << "\",\"ph\":\"i\",\"s\":\"t\",\"pid\":" << e.pid
                << ",\"tid\":" << e.tid << ",\"ts\":" << e.ts_us << "}";
        }
        for (auto& vec : ctr_) for (auto& c : vec) {
            out << ",\n{\"name\":\"" << escape(c.name) << "\",\"ph\":\"C\",\"pid\":" << c.pid
                << ",\"tid\":" << c.tid << ",\"ts\":" << c.ts_us
                << ",\"args\":{\"v\":" << c.value << "}}";
        }

        out << "\n]\n";
    }

private:
    Timeline() = default;

    static uint32_t make_lane(Role r, uint32_t idx, uint32_t& pid_out) {
        switch (r) {
            case Role::Reader:     pid_out = 100u + idx; return 1u;
            case Role::Loader:     pid_out = 100u + idx; return 2u;
            case Role::Decompress: pid_out = 100u + idx; return 3u;

            case Role::Producer:   pid_out = 200u; return 1000u + idx;
            case Role::Worker:     pid_out = 200u; return 2000u + idx;

            case Role::Writer:     pid_out = 300u; return 1u;
            case Role::BgzfWriter: pid_out = 300u; return 10u;
            case Role::BgzfWorker: pid_out = 300u; return 100u + idx;
        }
        pid_out = 200u;
        return 900000u + idx;
    }

    static uint32_t& tls_pid() { thread_local uint32_t pid = 1; return pid; }
    static uint32_t& tls_tid() { thread_local uint32_t tid = 0; return tid; }
    static int& tls_index()    { thread_local int idx = -1; return idx; }

    void ensure_registered() {
        int& idx = tls_index();
        if (idx >= 0) return;
        std::lock_guard<std::mutex> lk(mu_);
        all_.push_back({});
        inst_.push_back({});
        ctr_.push_back({});
        idx = (int)all_.size() - 1;
    }

    std::vector<Event>&        local_events() { ensure_registered(); return all_[(size_t)tls_index()]; }
    std::vector<InstantEvent>& local_inst()   { ensure_registered(); return inst_[(size_t)tls_index()]; }
    std::vector<CounterEvent>& local_ctr()    { ensure_registered(); return ctr_[(size_t)tls_index()]; }

    static std::string escape(const char* s) {
        std::string r;
        for (; *s; ++s) {
            if (*s == '"' || *s == '\\') r.push_back('\\');
            r.push_back(*s);
        }
        return r;
    }

    std::mutex mu_;
    std::vector<std::vector<Event>>        all_;
    std::vector<std::vector<InstantEvent>> inst_;
    std::vector<std::vector<CounterEvent>> ctr_;
};

class Scope {
public:
    explicit Scope(const char* name) : name_(name), start_(now_us_rel()) {}
    ~Scope() { Timeline::instance().add(name_, start_, now_us_rel() - start_); }
private:
    const char* name_;
    uint64_t start_;
};

#else

// ================= NO-OP IMPLEMENTATION (compiled out) =================
inline constexpr bool enabled() { return false; }

class Timeline {
public:
    static Timeline& instance() { static Timeline t; return t; }
    void set_group(uint32_t, uint32_t) {}
    void set_lane(Role, uint32_t) {}
    uint32_t pid()  const { return 0; }
    uint32_t lane() const { return 0; }
    void add(const char*, uint64_t, uint64_t) {}
    void instant(const char*, uint64_t) {}
    void counter(const char*, uint64_t, int64_t) {}
    void dump() {}
};

class Scope {
public:
    explicit Scope(const char*) {}
};

inline uint64_t now_us_rel() { return 0; } // keeps macro compatibility

#endif

} // namespace trace_tl

// Macros: compile-time gated so nothing is emitted when disabled
#if BQS_TRACE_COMPILED
#define BQS_TRACE_SET_LANE(role, idx) \
    do { if (trace_tl::enabled()) trace_tl::Timeline::instance().set_lane((role), (uint32_t)(idx)); } while(0)
#define BQS_TRACE_SCOPE(name_literal_or_cstr) \
    trace_tl::Scope _bqs_trace_scope_##__LINE__(name_literal_or_cstr)
#define BQS_TRACE_INSTANT(name_literal_or_cstr) \
    do { if (trace_tl::enabled()) trace_tl::Timeline::instance().instant((name_literal_or_cstr), trace_tl::now_us_rel()); } while(0)
#define BQS_TRACE_COUNTER(name_literal_or_cstr, value) \
    do { if (trace_tl::enabled()) trace_tl::Timeline::instance().counter((name_literal_or_cstr), trace_tl::now_us_rel(), (int64_t)(value)); } while(0)
#define BQS_TRACE_DUMP() \
    do { trace_tl::Timeline::instance().dump(); } while(0)
#else
#define BQS_TRACE_SET_LANE(role, idx) do {} while(0)
#define BQS_TRACE_SCOPE(name_literal_or_cstr) do {} while(0)
#define BQS_TRACE_INSTANT(name_literal_or_cstr) do {} while(0)
#define BQS_TRACE_COUNTER(name_literal_or_cstr, value) do {} while(0)
#define BQS_TRACE_DUMP() do {} while(0)
#endif