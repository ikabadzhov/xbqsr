#ifndef TESTONEAPI_EXTERNAL_MERGE_H
#define TESTONEAPI_EXTERNAL_MERGE_H

#include <tbb/parallel_for.h>
#include "../device_utils/types.h"
#include "../device_utils/meta.h"
#include "../device_utils/operator.h"

#include "../memory/context.h"
#include "../sorting/cta_merge.h"
#include "merge.h"


template<bounds_t bounds = bounds_lower, typename a_keys_it, typename b_keys_it, typename int_t, typename comp_t>
int_t merge_path_cpu(a_keys_it a_keys, int_t a_count, b_keys_it b_keys, int_t b_count, int_t diag, comp_t comp) {

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



template<bounds_t bounds, typename a_keys_it, typename b_keys_it, typename comp_t>
std::vector<int64_t> merge_path_partitions_cpu(a_keys_it a, int64_t a_count, b_keys_it b, int64_t b_count, int64_t spacing, comp_t comp) {
    int64_t num_partitions = (int)div_up(a_count + b_count, spacing) + 1;
    std::vector<int64_t> p(num_partitions);

    tbb::parallel_for(0UL, static_cast<size_t>(num_partitions), [&](size_t index) {
        auto diag = static_cast<int64_t>(min(spacing * (int64_t)index, a_count + b_count));
        p[index] = merge_path_cpu<bounds>(a, (int64_t) a_count, b, (int64_t) b_count,
                                          (int64_t)diag, comp);
    });
    return p;
}

// Key-value merge.
template<typename launch_arg_t = empty_t,
        typename a_keys_it, typename a_vals_it,
        typename b_keys_it, typename b_vals_it,
        typename c_keys_it, typename c_vals_it,
        typename comp_t>
void external_merge(a_keys_it a_keys, a_vals_it a_vals, int64_t a_count,
                    b_keys_it b_keys, b_vals_it b_vals, int64_t b_count,
                    c_keys_it c_keys, c_vals_it c_vals, comp_t comp, sycl::queue& q, size_t max_sz) {

    typedef typename std::iterator_traits<a_keys_it>::value_type type_t;
    typedef typename std::iterator_traits<a_vals_it>::value_type val_t;

    enum { has_values = !std::is_same<val_t, empty_t>::value };

    size_t sz_val = has_values ? sizeof(val_t) : 0;
//    sycl::device device = q.get_device();
//    std::string device_name = device.get_info<sycl::info::device::name>();
//    size_t global_mem_sz = device.get_info<sycl::info::device::global_mem_size>();
    size_t key_val_sz = (sizeof(type_t) + sz_val);

    if ((a_count + b_count) * 2 * key_val_sz < max_sz) { // Everything fit in memory
        mem_t<type_t> input_key_1(a_count, q);
        mem_t<type_t> input_key_2(b_count, q);
        mem_t<type_t> out_key(a_count + b_count, q);

        if (has_values) {
            mem_t<val_t> input_val_1(a_count, q);
            mem_t<val_t> input_val_2(b_count, q);
            mem_t<val_t> out_val(a_count + b_count, q);
            q.memcpy(input_key_1.data(), a_keys, sizeof(type_t) * a_count);
            q.memcpy(input_val_1.data(), a_vals, sizeof(val_t) * a_count);
            q.memcpy(input_key_2.data(), b_keys, sizeof(type_t) * b_count);
            q.memcpy(input_val_2.data(), b_vals, sizeof(val_t) * b_count);
            q.wait();
            merge(input_key_1.data(), input_val_1.data(), a_count, input_key_2.data(), input_val_2.data(), b_count, out_key.data(), out_val.data(), comp, q);
            q.wait();
            q.memcpy(c_keys, out_key.data(), sizeof(type_t) * (a_count + b_count));
            q.memcpy(c_vals, out_val.data(), sizeof(val_t) * (a_count + b_count));
            q.wait();
        } else {
            q.memcpy(input_key_1.data(), a_keys, sizeof(type_t) * a_count);
            q.memcpy(input_key_2.data(), b_keys, sizeof(type_t) * b_count);
            q.wait();
            merge(input_key_1.data(), a_count, input_key_2.data(), b_count, out_key.data(), comp, q);
            q.wait();
            q.memcpy(c_keys, out_key.data(), sizeof(type_t) * (a_count + b_count));
            q.wait();
        }
    }
    else { // Not fit in memory
        size_t max_elements = max_sz / (sizeof(type_t) + sz_val);
        size_t mem_out = max_elements / 2;
        size_t mem_in = mem_out / 2;

        size_t spacing = mem_in;

        mem_t<type_t> input_key_1(mem_in, q);
        mem_t<type_t> input_key_2(mem_in, q);
        mem_t<type_t> out_key(mem_in * 2, q);

        std::vector<int64_t> partitions = merge_path_partitions_cpu<bounds_lower>(a_keys, a_count, b_keys, b_count,
                                                                                  spacing, comp);

        //std::cout << "Number of partitions: " << partitions.size() << std::endl;
        if (has_values) {
            mem_t<val_t> input_val_1(mem_in, q);
            mem_t<val_t> input_val_2(mem_in, q);
            mem_t<val_t> out_val(mem_in * 2, q);
            size_t offset = 0;
            for (size_t pIdx = 0; pIdx < partitions.size() - 1; pIdx++) {
                auto mp0 = partitions[pIdx + 0];
                auto mp1 = partitions[pIdx + 1];
                auto range = compute_merge_range_long(a_count, b_count, pIdx, spacing, mp0, mp1);
                q.memcpy(input_key_1.data(), a_keys + range.a_begin, sizeof(type_t) * range.a_count());
                q.memcpy(input_val_1.data(), a_vals + range.a_begin, sizeof(val_t) * range.a_count());
                q.memcpy(input_key_2.data(), b_keys + range.b_begin, sizeof(type_t) * range.b_count());
                q.memcpy(input_val_2.data(), b_vals + range.b_begin, sizeof(val_t) * range.b_count());
                q.wait();
                merge(input_key_1.data(), input_val_1.data(), range.a_count(), input_key_2.data(), input_val_2.data(),
                      range.b_count(), out_key.data(), out_val.data(), comp, q);
                q.wait();
                q.memcpy(c_keys + offset, out_key.data(), sizeof(type_t) * (range.a_count() + range.b_count()));
                q.memcpy(c_vals + offset, out_val.data(), sizeof(val_t) * (range.a_count() + range.b_count()));
                offset += (range.a_count() + range.b_count());
            }
            q.wait();
        } else {
            size_t offset = 0;
            for (size_t pIdx = 0; pIdx < partitions.size() - 1; pIdx++) {
                auto mp0 = partitions[pIdx + 0];
                auto mp1 = partitions[pIdx + 1];
                auto range = compute_merge_range_long(a_count, b_count, pIdx, spacing, mp0, mp1);
                q.memcpy(input_key_1.data(), a_keys + range.a_begin, sizeof(type_t) * range.a_count());
                q.memcpy(input_key_2.data(), b_keys + range.b_begin, sizeof(type_t) * range.b_count());
                q.wait();
                merge(input_key_1.data(), range.a_count(), input_key_2.data(), range.b_count(), out_key.data(), comp,
                      q);
                q.wait();
                q.memcpy(c_keys + offset, out_key.data(), sizeof(type_t) * (range.a_count() + range.b_count()));
                offset += (range.a_count() + range.b_count());
            }
            q.wait();
        }
    }
}


/**
 * TODO
 * Possible first solution: we have two large arrays on CPU
 * We split the arrays in partitions with the merge_path algorithm.
 * We copy only the part of the arrays to be merged.
 *
 * If GMPU MEM can accomodate M elements, the arrays size is N each (assuming same size for both), with N > M
 * the merge operator should work giving in output (M / 2) elements, half for the output and the other half split among
 * the two input arrays.
 *
 * In each iteration we copy the input, merge on GPU, writing out the output, and asynchronously we start copying the
 * next batch.
 *
 * */
template<typename launch_t = empty_t, typename a_keys_it, typename b_keys_it, typename c_keys_it, typename comp_t>
void external_merge(a_keys_it a_keys, int64_t a_count, b_keys_it b_keys, int64_t b_count, c_keys_it c_keys, comp_t comp, sycl::queue& q, size_t max_sz) {
    external_merge<launch_t>(a_keys, (const empty_t*)nullptr, a_count, b_keys, (const empty_t*)nullptr, b_count, c_keys, (empty_t*)nullptr, comp, q, max_sz);
}


#endif //TESTONEAPI_EXTERNAL_MERGE_H
