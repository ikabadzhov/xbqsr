#ifndef TESTONEAPI_TYPES_H
#define TESTONEAPI_TYPES_H

#include "operator.h"


#define ALIGN(x) __attribute__((aligned(x)))

template<typename type_t, int size>
struct array_t {
    type_t data[size];

    type_t operator[](int i) const { return data[i]; }
    type_t& operator[](int i) { return data[i]; }

    array_t() = default;
    array_t(const array_t&) = default;
    array_t& operator=(const array_t&) = default;

    // Fill the array with x.
    array_t(type_t x) {
        iterate<size>([&](int i) { data[i] = x; });
    }
};

template<typename type_t>
struct array_t<type_t, 0> {
    type_t operator[](int i) const { return type_t(); }
    type_t& operator[](int i) { return *(type_t*)nullptr; }
};

// Reduce on components of array_t.
template<typename type_t, int size, typename op_t = plus_t<type_t> >
type_t reduce(array_t<type_t, size> x, op_t op = op_t()) {
    type_t a;
    iterate<size>([&](int i) {
        a = i ? op(a, x[i]) : x[i];
    });
    return a;
}

// Call the operator component-wise on all components.
template<typename type_t, int size, typename op_t>
array_t<type_t, size> combine(array_t<type_t, size> x,
                              array_t<type_t, size> y, op_t op) {

    array_t<type_t, size> z;
    iterate<size>([&](int i) { z[i] = op(x[i], y[i]); });
    return z;
}

template<typename type_t, int size>
array_t<type_t, size> operator+(
        array_t<type_t, size> a, array_t<type_t, size> b) {
    return combine(a, b, plus_t<type_t>());
}

template<typename type_t, int size>
array_t<type_t, size> operator-(
        array_t<type_t, size> a, array_t<type_t, size> b) {
    return combine(a, b, minus_t<type_t>());
}


template<typename key_t, typename val_t, int size>
struct kv_array_t {
    array_t<key_t, size> keys;
    array_t<val_t, size> vals;
};

enum bounds_t {
    bounds_lower,
    bounds_upper
};

struct ALIGN(8) range_t {
    int begin, end;
    int size() const { return end - begin; }
    int count() const { return size(); }
    bool valid() const { return end > begin; }
};

struct ALIGN(8) range_64_t {
    int64_t begin, end;
    int64_t size() const { return end - begin; }
    int64_t count() const { return size(); }
    bool valid() const { return end > begin; }
};


SYCL_EXTERNAL range_t get_tile(int cta, int nv, int count);

struct ALIGN(16) merge_range_t {
    int a_begin, a_end, b_begin, b_end;

    int a_count() const { return a_end - a_begin; }
    int b_count() const { return b_end - b_begin; }
    int total() const { return a_count() + b_count(); }

    range_t a_range() const {
        return range_t { a_begin, a_end };
    }

    range_t b_range() const {
        return range_t { b_begin, b_end };
    }

    merge_range_t to_local() const {
        return merge_range_t { 0, a_count(), a_count(), total() };
    }

    // Partition from mp to the end.
    merge_range_t partition(int mp0, int diag) const {
        return merge_range_t { a_begin + mp0, a_end, b_begin + diag - mp0, b_end };
    }

    // Partition from mp0 to mp1.
    merge_range_t partition(int mp0, int diag0,
                            int mp1, int diag1) const {
        return merge_range_t {
                a_begin + mp0,
                a_begin + mp1,
                b_begin + diag0 - mp0,
                b_begin + diag1 - mp1
        };
    }

    bool a_valid() const {
        return a_begin < a_end;
    }

    bool b_valid() const {
        return b_begin < b_end;
    }
};


struct ALIGN(16) merge_range_64_t {
    int64_t a_begin, a_end, b_begin, b_end;

    int64_t a_count() const { return a_end - a_begin; }
    int64_t b_count() const { return b_end - b_begin; }
    int64_t total() const { return a_count() + b_count(); }

    range_64_t a_range() const {
        return range_64_t { a_begin, a_end };
    }

    range_64_t b_range() const {
        return range_64_t { b_begin, b_end };
    }

    merge_range_64_t to_local() const {
        return merge_range_64_t { 0, a_count(), a_count(), total() };
    }

    // Partition from mp to the end.
    merge_range_64_t partition(int mp0, int diag) const {
        return merge_range_64_t { a_begin + mp0, a_end, b_begin + diag - mp0, b_end };
    }

    // Partition from mp0 to mp1.
    merge_range_64_t partition(int mp0, int diag0,
                               int mp1, int diag1) const {
        return merge_range_64_t {
                a_begin + mp0,
                a_begin + mp1,
                b_begin + diag0 - mp0,
                b_begin + diag1 - mp1
        };
    }

    bool a_valid() const {
        return a_begin < a_end;
    }

    bool b_valid() const {
        return b_begin < b_end;
    }
};


template<typename type_t, int size>
struct merge_pair_t {
    array_t<type_t, size> keys;
    array_t<int, size> indices;
};



#endif //TESTONEAPI_TYPES_H
