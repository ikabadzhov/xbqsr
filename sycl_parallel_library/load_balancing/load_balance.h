#pragma once

#include <sycl/sycl.hpp>
#include <iterator>
#include "cta_load_balance.h"
#include "../merge/search.h"


template<typename launch_arg_t = empty_t, typename func_t, typename segments_it, typename pointers_t, typename... args_t>
void transform_lbs(func_t f, int count, segments_it segments, int num_segments, pointers_t caching_iterators, sycl::queue& q, args_t... args) {
  enum {nt = 128, vt = 7, vt0 = vt, nv = nt * vt};

  typedef typename std::iterator_traits<segments_it>::value_type int_t;
  typedef syclgal::tuple_iterator_value_t<pointers_t> value_t;

  int num_ctas = div_up(count + num_segments, nv);
  mem_t<int_t> mp = load_balance_partitions(count, segments, num_segments, nv, q);

  const int_t* mp_data = mp.data();

    typedef cta_load_balance_t<nt, vt> load_balance_t;
    typedef detail::cached_segment_load_t<nt, pointers_t> cached_load_t;


  union SharedUnion {
    typename load_balance_t::storage_t lbs;
    typename cached_load_t::storage_t cached;
  };

  q.submit([&](sycl::handler &h) {
    // Define the local memory (shared memory)
    sycl::local_accessor<SharedUnion, 1> shared(sycl::range<1>(1), h);

    h.parallel_for(sycl::nd_range<1>(sycl::range<1>(num_ctas * nt), sycl::range<1>(nt)),
        [=](sycl::nd_item<1> item) {
          size_t tid = item.get_local_id(); // Current work-item index

          size_t cta = item.get_group(0); // Current work-group index


          // Compute the load-balancing search and materialize (index, seg, rank)
             // arrays.
          auto lbs = load_balance_t().load_balance(count, segments, num_segments,
                                                   tid, cta, mp_data, shared[0].lbs, item);

             // Load from the cached iterators. Use the placement range, not the
             // merge-path range for situating the segments.
             array_t<value_t, vt> cached_values = cached_load_t::template load<vt0>(
               tid, lbs.merge_range.a_count(), lbs.placement.range.b_range(),
               lbs.segments, shared[0].cached, caching_iterators, item);

             // Call the user-supplied functor f.
             strided_iterate<nt, vt, vt0>([=](int i, int j) {
               int index = lbs.merge_range.a_begin + j;
               int seg = lbs.segments[i];
               int rank = lbs.ranks[i];

               f(index, seg, rank, cached_values[i], args...);
             }, tid, lbs.merge_range.a_count());
    });
  });
  q.wait_and_throw();
}

/*// load-balancing search without caching.
template<typename launch_arg_t = empty_t, typename func_t,
  typename segments_it, typename... args_t>
void transform_lbs(func_t f, int count, segments_it segments,
  int num_segments, sycl::queue& q, args_t... args) {

  transform_lbs<launch_arg_t>(
    [=]MGPU_DEVICE(int index, int seg, int rank, tuple<>, args_t... args) {
      f(index, seg, rank, args...);    // drop the cached values.
    },
    count, segments, num_segments, tuple<>(), context, args...
  );
}

template<typename launch_arg_t = empty_t, typename segments_it,
  typename output_it>
void load_balance_search(int count, segments_it segments,
  int num_segments, output_it output, sycl::queue& q) {

  transform_lbs<launch_arg_t>([=]MGPU_DEVICE(int index, int seg, int rank) {
    output[index] = seg;
  }, count, segments, num_segments, context);
}*/
