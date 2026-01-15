#pragma once
#include <atomic>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <mutex>
#include <string>
#include <string_view>
#include <thread>
#include <vector>

namespace trace_tl {

enum class Role : uint32_t {
    Producer=1, Worker=2, Writer=3,
    BgzfWorker=4, BgzfWriter=5,
    Reader=6, Loader=7, Decompress=8
};

// --------- config (runtime) ----------
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

// --------- time ----------
inline uint64_t steady_us() {
    using namespace std::chrono;
    return (uint64_t)duration_cast<microseconds>(steady_clock::now().time_since_epoch()).count();
}
inline uint64_t now_us_rel() {
    static const uint64_t t0 = steady_us();
    return steady_us() - t0;
}

// --------- events ----------
struct Event        { const char* name; uint32_t pid; uint32_t tid; uint64_t ts_us; uint64_t dur_us; };          // ph:"X"
struct InstantEvent { const char* name; uint32_t pid; uint32_t tid; uint64_t ts_us; };                          // ph:"i"
struct CounterEvent { const char* name; uint32_t pid; uint32_t tid; uint64_t ts_us; int64_t value; };           // ph:"C"

// --------- per-thread buffer + global sink ----------
class Timeline {
public:
    static Timeline& instance() { static Timeline t; return t; }

    // Optional helper if you ever want explicit grouping:
    // - Use pid_group to group lanes under one "process" in Chrome.
    // - Use tid_in_group as a small lane id inside that group.
    void set_group(uint32_t pid_group, uint32_t tid_in_group) {
        if (!enabled()) return;
        tls_pid() = pid_group;
        tls_tid() = tid_in_group;
        ensure_registered();
    }

    // Main API: keep call sites unchanged (set_lane(role, idx))
    // - Readers (Reader/Loader/Decompress) are grouped FIRST in the UI by using pid in [100..]
    // - Others stay in pid=1.
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

        // ---- 1) Collect PIDs present in the trace ----
        std::vector<uint32_t> pids;
        auto add_pid = [&](uint32_t pid) {
            if (std::find(pids.begin(), pids.end(), pid) == pids.end())
                pids.push_back(pid);
        };

        for (auto& vec : all_)  for (auto& e : vec) add_pid(e.pid);
        for (auto& vec : inst_) for (auto& e : vec) add_pid(e.pid);
        for (auto& vec : ctr_)  for (auto& e : vec) add_pid(e.pid);

        // ---- 2) Emit process metadata FIRST (ph:"M") ----
        for (uint32_t pid : pids) {
            // Name
            std::string proc_name;
            if (pid >= 100 && pid < 200) proc_name = "Reader " + std::to_string(pid - 100);
            else if (pid == 200)         proc_name = "Workers";
            else if (pid == 300)         proc_name = "Writer";
            else                         proc_name = "PID " + std::to_string(pid);

            // Sort: smaller = shown earlier/top.
            // Put readers first, then pipeline.
            int sort_index = 0;
            if (pid >= 100 && pid < 200) sort_index = (int)(pid - 100); // Reader 0..N
            else if (pid == 200)         sort_index = 1000;             // after readers
            else if (pid == 300)         sort_index = 1100;             // after workers
            else                         sort_index = 2000;

            if (!first) out << ",\n";
            first = false;
            out << "{\"ph\":\"M\",\"name\":\"process_name\","
                << "\"pid\":" << pid << ",\"tid\":0,"
                << "\"args\":{\"name\":\"" << escape(proc_name.c_str()) << "\"}}";

            out << ",\n";
            out << "{\"ph\":\"M\",\"name\":\"process_sort_index\","
                << "\"pid\":" << pid << ",\"tid\":0,"
                << "\"args\":{\"sort_index\":" << sort_index << "}}";
        }

        // ---- 3) Your existing events (X / i / C) ----

        // ---- X events ----
        for (auto& vec : all_) {
            for (auto& e : vec) {
                out << ",\n";
                out << "{\"name\":\"" << escape(e.name) << "\","
                    << "\"ph\":\"X\","
                    << "\"pid\":" << e.pid << ","
                    << "\"tid\":" << e.tid << ","
                    << "\"ts\":"  << e.ts_us << ","
                    << "\"dur\":" << e.dur_us
                    << "}";
            }
        }

        // ---- Instant events ----
        for (auto& vec : inst_) {
            for (auto& e : vec) {
                out << ",\n";
                out << "{\"name\":\"" << escape(e.name) << "\","
                    << "\"ph\":\"i\",\"s\":\"t\","
                    << "\"pid\":" << e.pid << ","
                    << "\"tid\":" << e.tid << ","
                    << "\"ts\":"  << e.ts_us
                    << "}";
            }
        }

        // ---- Counter events ----
        for (auto& vec : ctr_) {
            for (auto& c : vec) {
                out << ",\n";
                out << "{\"name\":\"" << escape(c.name) << "\","
                    << "\"ph\":\"C\","
                    << "\"pid\":" << c.pid << ","
                    << "\"tid\":" << c.tid << ","
                    << "\"ts\":"  << c.ts_us << ","
                    << "\"args\":{\"v\":" << c.value << "}"
                    << "}";
            }
        }

        out << "\n]\n";
    }

private:
    Timeline() = default;

    // IMPORTANT:
    // - Reader-related roles are grouped under pid = 100 + index, so they show up BEFORE pid=1 in Chrome.
    // - Inside each reader pid, use small tids so lanes are compact and stable:
    //     tid=1 reader, tid=2 loader, tid=3 decompress.
    // - Everything else stays pid=1 and uses big tid ranges like before.
    static uint32_t make_lane(Role r, uint32_t idx, uint32_t& pid_out) {
        switch (r) {
            // ---- One process per reader ----
            case Role::Reader:     pid_out = 100u + idx; return 1u;
            case Role::Loader:     pid_out = 100u + idx; return 2u;
            case Role::Decompress: pid_out = 100u + idx; return 3u;

                // ---- Workers process (NO writer here) ----
            case Role::Producer:   pid_out = 200u; return 1000u + idx;
            case Role::Worker:     pid_out = 200u; return 2000u + idx;

                // ---- Writer process (writer + bgzf threads) ----
            case Role::Writer:     pid_out = 300u; return 1u;            // pipeline writer lane
            case Role::BgzfWriter: pid_out = 300u; return 10u;           // bgzf writer lane
            case Role::BgzfWorker: pid_out = 300u; return 100u + idx;    // bgzf worker lanes
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

// --------- scope helper ----------
class Scope {
public:
    explicit Scope(const char* name)
        : name_(name), start_(now_us_rel()) {}
    ~Scope() {
        Timeline::instance().add(name_, start_, now_us_rel() - start_);
    }
private:
    const char* name_;
    uint64_t start_;
};

} // namespace trace_tl

// Macros
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