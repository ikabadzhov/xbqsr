#pragma once
#include <iterator>
#include "../device_utils/types.h"
#include "../memory/context.h"
#include "../device_utils/meta.h"
#include "../sorting/cta_merge.h"
#include "../device_utils/iterators.h"

template<bounds_t bounds, typename a_keys_it, typename b_keys_it, typename comp_t>
mem_t<int> merge_path_partitions(a_keys_it a, int64_t a_count, b_keys_it b, int64_t b_count, int64_t spacing, comp_t comp, sycl::queue& q) {

    //typedef typename std::iterator_traits<a_keys_it>::value_type int_t;
    int num_partitions = (int)div_up(a_count + b_count, spacing) + 1;
    mem_t<int> mem(num_partitions, q);
    int* p = mem.data();

    q.submit([&](sycl::handler &h) {
        h.parallel_for(
                sycl::range<1>(num_partitions), [=](sycl::id<1> idx) {
                    int index = idx[0];
                    int diag = (int) min(spacing * index, a_count + b_count);
                    p[index] = merge_path<bounds>(a, (int) a_count, b, (int) b_count,
                                                  diag, comp);
                });
    });
    q.wait();
    return mem;
}


template<bounds_t bounds, typename a_keys_it, typename b_keys_it, typename comp_t>
mem_t<size_t> merge_path_partitions_long(a_keys_it a, int64_t a_count, b_keys_it b, int64_t b_count, int64_t spacing, comp_t comp, sycl::queue& q) {

    typedef typename std::iterator_traits<a_keys_it>::value_type int_t;
    int num_partitions = (int)div_up(a_count + b_count, spacing) + 1;
    mem_t<size_t> mem(num_partitions, q);
    size_t* p = mem.data();

    q.submit([&](sycl::handler &h) {
        h.parallel_for(
                sycl::range<1>(num_partitions), [=](sycl::id<1> idx) {
                    int index = idx[0];
                    int diag = (int) min(spacing * index, a_count + b_count);
                    p[index] = merge_path<bounds>(a, (int) a_count, b, (int) b_count,
                                                  diag, comp);
                });
    });
    q.wait();
    return mem;
}

template<typename segments_it>
auto load_balance_partitions(int64_t dest_count, segments_it segments,
                             int num_segments, int spacing, sycl::queue& q) ->
mem_t<typename std::iterator_traits<segments_it>::value_type> {

    typedef typename std::iterator_traits<segments_it>::value_type int_t;
    return merge_path_partitions<bounds_upper>(counting_iterator_t<int_t>(0),
                                               dest_count, segments, num_segments, spacing, less_t<int_t>(), q);
}

template<typename segments_it>
auto load_balance_partitions_long(int64_t dest_count, segments_it segments,
                             int num_segments, int spacing, sycl::queue& q) ->
mem_t<typename std::iterator_traits<segments_it>::value_type> {

    typedef typename std::iterator_traits<segments_it>::value_type int_t;
    return merge_path_partitions_long<bounds_upper>(counting_iterator_t<int_t>(0),
                                               dest_count, segments, num_segments, spacing, less_t<int_t>(), q);
}