#ifndef TESTONEAPI_MERGESORT_H
#define TESTONEAPI_MERGESORT_H

#include <vector>
#include <sycl/queue.hpp>
#include "../memory/context.h"
#include "../device_utils/meta.h"
#include "cta_mergesort.h"

namespace syclgal_mergesort {
    const size_t MAX_NUM_ELEMENTS = 2000000000;
}
template<typename keys_it, typename comp_t>
mem_t<int> merge_sort_partitions(keys_it keys, int count, int coop,
                                 int spacing, comp_t comp, sycl::queue& q) {

    int num_partitions = div_up(count, spacing) + 1;

    mem_t<int> mem(num_partitions, q);

    int *p = mem.data();
    q.parallel_for(sycl::range<1>(num_partitions), [=](sycl::id<1> idx) {
        int index = idx[0];

        // Compute the merge range for this partition
        merge_range_t range = compute_mergesort_range(count, index, coop, spacing);

        // Calculate diag
        int diag = std::min(spacing * index, count) - range.a_begin;

        // Perform the merge path operation
        p[index] = merge_path<bounds_lower>(
                keys + range.a_begin, range.a_count(),
                keys + range.b_begin, range.b_count(),
                diag, comp);
    }).wait();
    return mem;
}

// Key-value mergesort.
template<typename launch_arg_t = empty_t, typename key_t, typename val_t,
        typename comp_t>
void mergesort(key_t* keys_input, val_t* vals_input, int count,
               comp_t comp, sycl::queue& q) {

    enum { has_values = !std::is_same<val_t, empty_t>::value };

    enum { nt = 128, vt = 11, nv = nt*vt };
    int num_ctas = div_up(count, nv);
    int num_passes = find_log2(num_ctas, true);

    mem_t<key_t> keys_temp(num_passes ? count : 0, q);
    key_t* keys_output = keys_temp.data();

    mem_t<val_t> vals_temp(has_values && num_passes ? count : 0, q);
    val_t* vals_output = vals_temp.data();

    key_t* keys_blocksort = (1 & num_passes) ? keys_output : keys_input;
    val_t* vals_blocksort = (1 & num_passes) ? vals_output : vals_input;

    typedef cta_sort_t<nt, vt, key_t, val_t> sort_t;

    union SharedUnion {
        typename sort_t::storage_t sort;
        key_t keys[nv];
        val_t vals[nv];
    };

    q.submit([&](sycl::handler &h) {
        // Define the local memory (shared memory)
        sycl::local_accessor<SharedUnion, 1> shared(sycl::range<1>(1), h);

        h.parallel_for(
                sycl::nd_range<1>(sycl::range<1>(num_ctas * nt),
                                  sycl::range<1>(nt)),
                [=](sycl::nd_item<1> item) {
                    size_t tid = item.get_local_id(); // Current work-item index
                    size_t cta = item.get_group(0); // Current work-group index

                    range_t tile = get_tile(cta, nv, count);

                    // Load the keys and values.
                    kv_array_t<key_t, val_t, vt> unsorted;
                    unsorted.keys = mem_to_reg_thread<nt, vt>(keys_input + tile.begin, tid, tile.count(), shared[0].keys, item);

                    if (has_values)
                        unsorted.vals = mem_to_reg_thread<nt, vt>(
                                vals_input + tile.begin, tid, tile.count(), shared[0].vals, item);

                    // Blocksort.
                    kv_array_t<key_t, val_t, vt> sorted = sort_t().block_sort(unsorted,
                                                                              tid, tile.count(), comp, shared[0].sort, item);

                    // Store the keys and values.
                    reg_to_mem_thread<nt, vt>(sorted.keys, tid, tile.count(),
                                              keys_blocksort + tile.begin, shared[0].keys,
                                              item);
                    if (has_values) {
                        reg_to_mem_thread<nt, vt>(sorted.vals, tid, tile.count(),
                                                  vals_blocksort + tile.begin, shared[0].vals,
                                                  item);
                    }
                });
    });

    q.wait();

    if(1 & num_passes) {
        std::swap(keys_input, keys_output);
        std::swap(vals_input, vals_output);
    }

    union SharedData {
        key_t keys[nv + 1];
        int indices[nv];
    };

    for(int pass = 0; pass < num_passes; ++pass) {
        int coop = 2 << pass;
        mem_t<int> partitions = merge_sort_partitions(keys_input, count, coop, nv, comp, q);
        q.wait();

        int* mp_data = partitions.data();

        q.submit([&](sycl::handler &h) {
            // Define the local memory (shared memory)
            sycl::local_accessor<SharedData, 1> shared(sycl::range<1>(1), h);

            h.parallel_for(
                sycl::nd_range<1>(sycl::range<1>(num_ctas * nt),
                                  sycl::range<1>(nt)),
                [=](sycl::nd_item<1> item) {
                    size_t tid = item.get_local_id(); // Current work-item index
                    size_t cta = item.get_group(0); // Current work-group index

                    range_t tile = get_tile(cta, nv, count);

                    // Load the range for this CTA and merge the values into register.
                    merge_range_t range = compute_mergesort_range(count, cta, coop, nv,
                                                                  mp_data[cta + 0], mp_data[cta + 1]);

                    merge_pair_t<key_t, vt> merge = cta_merge_from_mem<bounds_lower, nt, vt>(
                            keys_input, keys_input, range, tid, comp, shared[0].keys, item);

                    // Store merged values back out.
                    reg_to_mem_thread<nt>(merge.keys, tid, tile.count(),
                                          keys_output + tile.begin, shared[0].keys, item);

                    if(has_values) {
                        // Transpose the indices from thread order to strided order.
                        array_t<int, vt> indices = reg_thread_to_strided<nt>(
                                merge.indices, tid, shared[0].indices, item);

                        // Gather the input values and merge into the output values.
                        transfer_two_streams_strided<nt>(vals_input + range.a_begin,
                                                         range.a_count(), vals_input + range.b_begin, range.b_count(),
                                                         indices, tid, vals_output + tile.begin);
                    }
                });
        });
        q.wait();
        std::swap(keys_input, keys_output);
        std::swap(vals_input, vals_output);
    }
}

template<typename launch_arg_t = empty_t, typename key_t, typename val_t,
        typename comp_t>
void mergesort_developer_version(key_t* keys_input, val_t* vals_input, int count,
               comp_t comp, sycl::queue& q, key_t* keys_temp, val_t* vals_temp) {

    enum { has_values = !std::is_same<val_t, empty_t>::value };

    enum { nt = 128, vt = 11, nv = nt*vt };
    int num_ctas = div_up(count, nv);
    int num_passes = find_log2(num_ctas, true);

    //mem_t<key_t> keys_temp(num_passes ? count : 0, q);
    //keys_temp = (num_passes ? count : nullptr);
    key_t* keys_output = keys_temp;//.data();

    //mem_t<val_t> vals_temp(has_values && num_passes ? count : 0, q);
    val_t* vals_output = vals_temp;//.data();

    key_t* keys_blocksort = (1 & num_passes) ? keys_output : keys_input;
    val_t* vals_blocksort = (1 & num_passes) ? vals_output : vals_input;

    typedef cta_sort_t<nt, vt, key_t, val_t> sort_t;

    union SharedUnion {
        typename sort_t::storage_t sort;
        key_t keys[nv];
        val_t vals[nv];
    };

    q.submit([&](sycl::handler &h) {
        // Define the local memory (shared memory)
        sycl::local_accessor<SharedUnion, 1> shared(sycl::range<1>(1), h);

        h.parallel_for(
                sycl::nd_range<1>(sycl::range<1>(num_ctas * nt),
                                  sycl::range<1>(nt)),
                [=](sycl::nd_item<1> item) {
                    size_t tid = item.get_local_id(); // Current work-item index
                    size_t cta = item.get_group(0); // Current work-group index

                    range_t tile = get_tile(cta, nv, count);

                    // Load the keys and values.
                    kv_array_t<key_t, val_t, vt> unsorted;
                    unsorted.keys = mem_to_reg_thread<nt, vt>(keys_input + tile.begin, tid, tile.count(), shared[0].keys, item);

                    if (has_values)
                        unsorted.vals = mem_to_reg_thread<nt, vt>(
                                vals_input + tile.begin, tid, tile.count(), shared[0].vals, item);

                    // Blocksort.
                    kv_array_t<key_t, val_t, vt> sorted = sort_t().block_sort(unsorted,
                                                                              tid, tile.count(), comp, shared[0].sort, item);

                    // Store the keys and values.
                    reg_to_mem_thread<nt, vt>(sorted.keys, tid, tile.count(),
                                              keys_blocksort + tile.begin, shared[0].keys,
                                              item);
                    if (has_values) {
                        reg_to_mem_thread<nt, vt>(sorted.vals, tid, tile.count(),
                                                  vals_blocksort + tile.begin, shared[0].vals,
                                                  item);
                    }
                });
    });

    q.wait();

    if(1 & num_passes) {
        std::swap(keys_input, keys_output);
        std::swap(vals_input, vals_output);
    }

    union SharedData {
        key_t keys[nv + 1];
        int indices[nv];
    };

    for(int pass = 0; pass < num_passes; ++pass) {
        int coop = 2 << pass;
        mem_t<int> partitions = merge_sort_partitions(keys_input, count, coop, nv, comp, q);
        q.wait();

        int* mp_data = partitions.data();

        q.submit([&](sycl::handler &h) {
            // Define the local memory (shared memory)
            sycl::local_accessor<SharedData, 1> shared(sycl::range<1>(1), h);

            h.parallel_for(
                    sycl::nd_range<1>(sycl::range<1>(num_ctas * nt),
                                      sycl::range<1>(nt)),
                    [=](sycl::nd_item<1> item) {
                        size_t tid = item.get_local_id(); // Current work-item index
                        size_t cta = item.get_group(0); // Current work-group index

                        range_t tile = get_tile(cta, nv, count);

                        // Load the range for this CTA and merge the values into register.
                        merge_range_t range = compute_mergesort_range(count, cta, coop, nv,
                                                                      mp_data[cta + 0], mp_data[cta + 1]);

                        merge_pair_t<key_t, vt> merge = cta_merge_from_mem<bounds_lower, nt, vt>(
                                keys_input, keys_input, range, tid, comp, shared[0].keys, item);

                        // Store merged values back out.
                        reg_to_mem_thread<nt>(merge.keys, tid, tile.count(),
                                              keys_output + tile.begin, shared[0].keys, item);

                        if(has_values) {
                            // Transpose the indices from thread order to strided order.
                            array_t<int, vt> indices = reg_thread_to_strided<nt>(
                                    merge.indices, tid, shared[0].indices, item);

                            // Gather the input values and merge into the output values.
                            transfer_two_streams_strided<nt>(vals_input + range.a_begin,
                                                             range.a_count(), vals_input + range.b_begin, range.b_count(),
                                                             indices, tid, vals_output + tile.begin);
                        }
                    });
        });
        q.wait();
        std::swap(keys_input, keys_output);
        std::swap(vals_input, vals_output);
    }
}


// Key-only mergesort
template<typename launch_arg_t = empty_t, typename key_t, typename comp_t>
void mergesort(key_t* keys_input, int count, comp_t comp,
               sycl::queue& q) {
    mergesort<launch_arg_t>(keys_input, (empty_t*)nullptr, count, comp, q);
}

#endif //TESTONEAPI_MERGESORT_H
