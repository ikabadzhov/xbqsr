#pragma once
#include <cstddef>
#include <cstdint>

//TODO: This parameter must be determined at runtime
const int warp_size = 32;

template< bool B, class T = void >
using enable_if_t = typename std::enable_if<B,T>::type;

template<typename real_t>
constexpr real_t min(real_t a, real_t b) {
    return (b < a) ? b : a;
}

template<typename real_t>
constexpr real_t max(real_t a, real_t b) {
    return (a < b) ? b : a;
}

struct empty_t { };

template<typename... args_t>
void swallow(args_t...) { }

constexpr int clz(int x) {
    for(int i = 31; i >= 0; --i)
        if((1<< i) & x) return 31 - i;
    return 32;
}

constexpr int popc(unsigned x) {
    int c;
    for (c = 0; x; ++c)
        x &= x - 1;
    return c;
}

constexpr unsigned bfe(unsigned x, unsigned bit, unsigned num_bits) {
    unsigned result;
    result = ((1<< num_bits) - 1) & (x>> bit);
    return result;
}

constexpr bool is_pow2(int x) {
    return 0 == (x & (x - 1));
}

constexpr int find_log2(int x, bool round_up = false) {
    int a = 31 - clz(x);
    if(round_up) a += !is_pow2(x);
    return a;
}

SYCL_EXTERNAL constexpr int div_up(int x, int y) {
    return (x + y - 1) / y;
}

SYCL_EXTERNAL constexpr int64_t div_up(int64_t x, int64_t y) {
    return (x + y - 1) / y;
}
SYCL_EXTERNAL constexpr size_t div_up(size_t x, size_t y) {
    return (x + y - 1) / y;
}


SYCL_EXTERNAL constexpr int s_log2(int x, int p = 0) {
    int log = 0;
    while (x >>= 1) { // Equivalent to x = x / 2
        log++;
    }
    return log;
}

SYCL_EXTERNAL constexpr size_t s_log2(size_t x, size_t p = 0) {
    size_t log = 0;
    while (x >>= 1) { // Equivalent to x = x / 2
        log++;
    }
    return log;
}

////////////////////////////////////////////////////////////////////////////////
// Conditional typedefs.

// Typedef type_a if type_a is not empty_t.
// Otherwise typedef type_b.
template<typename type_a, typename type_b>
struct conditional_typedef_t {
    typedef typename std::conditional<
            !std::is_same<type_a, empty_t>::value,
            type_a,
            type_b
    >::type type_t;
};


////////////////////////////////////////////////////////////////////////////////
// Template unrolled looping construct.

template<int i, int count, bool valid = (i < count)>
struct iterate_t {
    template<typename func_t>
    static void eval(func_t f) {
        f(i);
        iterate_t<i + 1, count>::eval(f);
    }
};

template<int i, int count>
struct iterate_t<i, count, false> {
    template<typename func_t>
    static void eval(func_t f) { }
};

template<int begin, int end, typename func_t>
void iterate(func_t f) {
    iterate_t<begin, end>::eval(f);
}

template<int count, typename func_t>
void iterate(func_t f) {
    iterate<0, count>(f);
}

template<int count, typename type_t>
type_t reduce(const type_t(&x)[count]) {
    type_t y;
    iterate<count>([&](int i) { y = i ? x[i] + y : x[i]; });
    return y;
}

template<int count, typename type_t>
void fill(type_t(&x)[count], type_t val) {
    iterate<count>([&](int i) { x[i] = val; });
}


// Invoke unconditionally.
template<int nt, int vt, typename func_t>
void strided_iterate(func_t f, int tid) {
    iterate<vt>([=](int i) { f(i, nt * i + tid); });
}

// Check range.
template<int nt, int vt, int vt0 = vt, typename func_t>
void strided_iterate(func_t f, int tid, int count) {
    // Unroll the first vt0 elements of each thread.
    if(vt0 > 1 && count >= nt * vt0) {
        strided_iterate<nt, vt0>(f, tid);    // No checking
    } else {
        iterate<vt0>([=](int i) {
            int j = nt * i + tid;
            if(j < count) {
                f(i, j);
            }
        });
    }

    iterate<vt0, vt>([=](int i) {
        int j = nt * i + tid;
        if(j < count) {
            f(i, j);
        }
    });
}
template<int vt, typename func_t>
void thread_iterate(func_t f, int tid) {
    iterate<vt>([=](int i) { f(i, vt * tid + i); });
}

