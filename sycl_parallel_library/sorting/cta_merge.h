#ifndef TESTONEAPI_CTA_MERGE_H
#define TESTONEAPI_CTA_MERGE_H

#include "../device_utils/types.h"
#include "../device_utils/meta.h"

template<bounds_t bounds = bounds_lower, typename a_keys_it, typename b_keys_it, typename int_t, typename comp_t>
int_t merge_path(a_keys_it a_keys, int_t a_count, b_keys_it b_keys, int_t b_count, int_t diag, comp_t comp) {

    typedef typename std::iterator_traits<a_keys_it>::value_type type_t;
    int_t begin = max(int_t(0), (int_t)(diag - b_count));
    int_t end = min(diag, a_count);

    while(begin < end) {
        int_t mid = (begin + end) / 2;
        type_t a_key = a_keys[mid];
        type_t b_key = b_keys[diag - 1 - mid];
        bool pred = (bounds_upper == bounds) ? comp(a_key, b_key) : !comp(b_key, a_key);

        if(pred) begin = mid + 1;
        else end = mid;
    }
    return begin;
}

template<bounds_t bounds, typename keys_it, typename comp_t>
int merge_path(keys_it keys, merge_range_t range, int diag, comp_t comp) {
    return merge_path<bounds>(
            keys + range.a_begin, range.a_count(),
            keys + range.b_begin, range.b_count(),
            diag, comp);
}

template<bounds_t bounds, bool range_check, typename type_t, typename comp_t>
bool merge_predicate(type_t a_key, type_t b_key, merge_range_t range, comp_t comp) {
    bool p;
    if(range_check && !range.a_valid()) p = false;
    else if(range_check && !range.b_valid()) p = true;
    else p = (bounds_upper == bounds) ? comp(a_key, b_key) : !comp(b_key, a_key);
    return p;
}

SYCL_EXTERNAL merge_range_t compute_merge_range(int a_count, int b_count, int partition, int spacing, int mp0, int mp1);
SYCL_EXTERNAL merge_range_64_t compute_merge_range_long(int64_t a_count, int64_t b_count, int64_t partition, int64_t spacing, int64_t mp0, int64_t mp1);


// Specialization that emits just one LD instruction. Can only reliably used
// with raw pointer types. Fixed not to use pointer arithmetic so that
// we don't get undefined behaviors with unaligned types.
template<int nt, int vt, typename type_t>
array_t<type_t, vt> load_two_streams_reg(const type_t* a, int a_count, const type_t* b, int b_count, int tid) {

    b -= a_count;
    array_t<type_t, vt> x;
    strided_iterate<nt, vt>([&](int i, int index) {
        const type_t* p = (index >= a_count) ? b : a;
        x[i] = p[index];
    }, tid, a_count + b_count);

    return x;
}

template<int nt, int vt, typename type_t, typename a_it, typename b_it>
enable_if_t<!(std::is_pointer<a_it>::value && std::is_pointer<b_it>::value), array_t<type_t, vt>> load_two_streams_reg(a_it a, int a_count, b_it b, int b_count, int tid) {
    b -= a_count;
    array_t<type_t, vt> x;
    strided_iterate<nt, vt>([&](int i, int index) {
        x[i] = (index < a_count) ? a[index] : b[index];
    }, tid, a_count + b_count);
    return x;
}

template<int nt, int vt, typename a_it, typename b_it, typename type_t, int shared_size>
void load_two_streams_shared(a_it a, int a_count, b_it b, int b_count, int tid, type_t (&shared)[shared_size], sycl::nd_item<1>& item, bool sync = true) {

    // Load into register then make an unconditional strided store into memory.
    array_t<type_t, vt> x = load_two_streams_reg<nt, vt, type_t>(a, a_count, b, b_count, tid);
    reg_to_shared_strided<nt>(x, tid, shared, item, sync);
}

template<int nt, int vt, typename type_t>
array_t<type_t, vt> gather_two_streams_strided(const type_t* a, int a_count, const type_t* b, int b_count, array_t<int, vt> indices, int tid) {

    ptrdiff_t b_offset = b - a - a_count;
    int count = a_count + b_count;

    array_t<type_t, vt> x;
    strided_iterate<nt, vt>([&](int i, int j) {
        ptrdiff_t gather = indices[i];
        if(gather >= a_count) gather += b_offset;
        x[i] = a[gather];
    }, tid, count);

    return x;
}

template<int nt, int vt, typename type_t, typename a_it, typename b_it>
enable_if_t<!(std::is_pointer<a_it>::value && std::is_pointer<b_it>::value), array_t<type_t, vt>> gather_two_streams_strided(a_it a, int a_count, b_it b, int b_count, array_t<int, vt> indices, int tid) {

    b -= a_count;
    array_t<type_t, vt> x;
    strided_iterate<nt, vt>([&](int i, int j) {
        x[i] = (indices[i] < a_count) ? a[indices[i]] : b[indices[i]];
        }, tid, a_count + b_count
    );
    return x;
}

template<int nt, int vt, typename a_it, typename b_it, typename c_it>
void transfer_two_streams_strided(a_it a, int a_count, b_it b,
                                              int b_count, array_t<int, vt> indices, int tid, c_it c) {
    typedef typename std::iterator_traits<a_it>::value_type type_t;
    array_t<type_t, vt> x = gather_two_streams_strided<nt, vt, type_t>(a, a_count, b, b_count, indices, tid);
    reg_to_mem_strided<nt>(x, tid, a_count + b_count, c);
}


// This function must be able to dereference keys[a_begin] and keys[b_begin],
// no matter the indices for each. The caller should allocate at least
// nt * vt + 1 elements for
template<bounds_t bounds, int vt, typename type_t, typename comp_t>
merge_pair_t<type_t, vt> serial_merge(const type_t* keys_shared, merge_range_t range, comp_t comp, sycl::nd_item<1>& item, bool sync = true) {

    type_t a_key = keys_shared[range.a_begin];
    type_t b_key = keys_shared[range.b_begin];

    merge_pair_t<type_t, vt> merge_pair;
    iterate<vt>([&](int i) {
        bool p = merge_predicate<bounds, true>(a_key, b_key, range, comp);
        int index = p ? range.a_begin : range.b_begin;

        merge_pair.keys[i] = p ? a_key : b_key;
        merge_pair.indices[i] = index;

        type_t c_key = keys_shared[++index];
        if(p) a_key = c_key, range.a_begin = index;
        else b_key = c_key, range.b_begin = index;
    });

    if(sync) item.barrier(sycl::access::fence_space::local_space);
    return merge_pair;
}

// Load arrays a and b from global memory and merge into register.
template<bounds_t bounds, int nt, int vt, typename a_it, typename b_it, typename type_t, typename comp_t, int shared_size>
merge_pair_t<type_t, vt> cta_merge_from_mem(a_it a, b_it b, merge_range_t range_mem, int tid, comp_t comp, type_t (&keys_shared)[shared_size], sycl::nd_item<1>& item) {

    static_assert(shared_size >= nt * vt + 1, "cta_merge_from_mem requires temporary storage of at least nt * vt + 1 items");

    // Load the data into shared memory.
    load_two_streams_shared<nt, vt>(a + range_mem.a_begin, range_mem.a_count(), b + range_mem.b_begin, range_mem.b_count(), tid, keys_shared, item, true);

    // Run a merge path to find the start of the serial merge for each thread.
    merge_range_t range_local = range_mem.to_local();
    int diag = vt * tid;
    int mp = merge_path<bounds>(keys_shared, range_local, diag, comp);

    // Compute the ranges of the sources in shared memory. The end iterators
    // of the range are inaccurate, but still facilitate exact merging, because
    // only vt elements will be merged.
    merge_pair_t<type_t, vt> merged = serial_merge<bounds, vt>(keys_shared, range_local.partition(mp, diag), comp, item);

    return merged;
};

#endif //TESTONEAPI_CTA_MERGE_H
