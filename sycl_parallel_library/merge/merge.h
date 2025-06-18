#ifndef TESTONEAPI_MERGE_H
#define TESTONEAPI_MERGE_H

#include <sycl/sycl.hpp>
#include "search.h"
#include "../device_utils/meta.h"
#include "../device_utils/types.h"
#include "../memory/loadstore.h"
#include "../memory/memory.h"
#include "../memory/context.h"
#include "../sorting/cta_merge.h"


// Key-value merge.
template<typename launch_arg_t = empty_t,
        typename a_keys_it, typename a_vals_it,
        typename b_keys_it, typename b_vals_it,
        typename c_keys_it, typename c_vals_it,
        typename comp_t>
void merge(a_keys_it a_keys, a_vals_it a_vals, int a_count,
           b_keys_it b_keys, b_vals_it b_vals, int b_count,
           c_keys_it c_keys, c_vals_it c_vals, comp_t comp, sycl::queue& q) {

//    typedef typename conditional_typedef_t<launch_arg_t,
//            launch_box_t<
//            arch_20_cta<128, 15>,
//    arch_35_cta<128, 11>,
//            arch_52_cta<128, 15>
//                             >
//                             >::type_t launch_t;

    enum { nt = 128, vt = 11, nv = nt*vt };

    typedef typename std::iterator_traits<a_keys_it>::value_type type_t;
    typedef typename std::iterator_traits<a_vals_it>::value_type val_t;
    enum { has_values = !std::is_same<val_t, empty_t>::value };

    mem_t<int> partitions = merge_path_partitions<bounds_lower>(a_keys, a_count,
                                                                b_keys, b_count, nv, comp, q);
    int* mp_data = partitions.data();

    // if (a_count == 39719176){
    //     mem_t<int> h_partition(partitions.size(), q, memory_space_t::memory_space_shared);
    //     q.memcpy(h_partition.data(), partitions.data(), sizeof(int)*partitions.size()).wait();
    //     for (size_t idx = 0; idx < partitions.size(); idx++) {
    //         std::cout << h_partition.data()[idx] << std::endl;
    //     }
    // }

    int num_ctas = div_up(a_count + b_count, nv);

    union SharedUnion {
        type_t keys[nv + 1];
        int indices[nv];
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

                    int mp0 = mp_data[cta + 0];
                    int mp1 = mp_data[cta + 1];
                    merge_range_t range = compute_merge_range(a_count, b_count, cta, nv,
                                                              mp0, mp1);

                    merge_pair_t<type_t, vt> merge = cta_merge_from_mem<bounds_lower, nt, vt>(
                            a_keys, b_keys, range, tid, comp, shared[0].keys, item);

                    int dest_offset = nv * cta;
                    reg_to_mem_thread<nt>(merge.keys, tid, range.total(), c_keys + dest_offset,
                                          shared[0].keys, item);

                    if(has_values) {
                        // Transpose the indices from thread order to strided order.
                        array_t<int, vt> indices = reg_thread_to_strided<nt>(merge.indices, tid,
                                                                             shared[0].indices, item);

                        // Gather the input values and merge into the output values.
                        transfer_two_streams_strided<nt>(a_vals + range.a_begin, range.a_count(),
                                                         b_vals + range.b_begin, range.b_count(), indices, tid,
                                                         c_vals + dest_offset);
                    }
                });
    });
    q.wait();
}

template<typename launch_t = empty_t, typename a_keys_it, typename b_keys_it, typename c_keys_it, typename comp_t>
void merge(a_keys_it a_keys, int a_count, b_keys_it b_keys, int b_count, c_keys_it c_keys, comp_t comp, sycl::queue& q) {
    merge<launch_t>(a_keys, (const empty_t*)nullptr, a_count, b_keys, (const empty_t*)nullptr, b_count, c_keys, (empty_t*)nullptr, comp, q);
}


#endif //TESTONEAPI_MERGE_H
