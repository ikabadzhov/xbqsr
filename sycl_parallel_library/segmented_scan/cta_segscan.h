#pragma once

#include "../device_utils/meta.h"
#include "../device_utils/operator.h"
#include "../device_utils/intrinsics.h"

#include <sycl/sycl.hpp>

template<typename type_t>
struct segscan_result_t {
    type_t scan;
    type_t reduction;
    bool has_carry_in;
    int left_lane;
};

template<int nt, typename type_t>
struct cta_segscan_t {
    enum { num_warps = nt / warp_size };

    union storage_t {
        int delta[num_warps + nt] = 0;
        struct { type_t values[2 * nt]; int packed[nt]; };
    };

    inline int find_left_lane(int tid, bool has_head_flag,
                                   storage_t& storage, const sycl::nd_item<1> &item) const {
        int32_t warp = tid / warp_size;
        int32_t lane = (warp_size - 1) & tid;
        int warp_mask = 0xffffffff >> (31 - lane);   // inclusive search.
        int cta_mask = 0x7fffffff >> (31 - lane);    // exclusive search.

        // Build a head flag bitfield and store it into shared memory.
        // In CUDA would be:
        // int warp_bits = ballot(has_head_flag);
        auto sg = item.get_sub_group();
        int warp_bits = static_cast<int>(sycl::reduce_over_group(sg, has_head_flag ? (1UL << sg.get_local_linear_id()) : 0, sycl::plus<int>()));

        storage.delta[warp] = warp_bits;

        /*
        DPCT1065:73: Consider replacing sycl::nd_item::barrier() with
        sycl::nd_item::barrier(sycl::access::fence_space::local_space) for better
        performance if there is no access to global memory.
        */
        item.barrier(sycl::access::fence_space::local_space);
        /*
         * TODO: Check with people in SYCLOPS. Why if this condition is moved below it works while if kept here it hang forever?
         * */
        //if (tid < num_warps) {
//            /*
//            DPCT1086:24: __activemask() is migrated to 0xffffffff. You may need to
//            adjust the code.
//            */
//            //unsigned mask = 0xffffffff;
            // In CUDA would be:
            //  int cta_bits = ballot(0 != storage.delta[tid], item);
            //auto sg = item.get_sub_group();
            bool predicate = false;
            if (tid < num_warps) {
                 predicate = (0 != storage.delta[tid]);
            }

            size_t cta_bits = sycl::reduce_over_group(sg, predicate ? (1UL << sg.get_local_linear_id()) : 0, sycl::plus<size_t>());


            int warp_segment = 31 - clz(cta_mask & static_cast<int>(cta_bits));
            int start = (-1 != warp_segment) ?
                        (31 - clz(storage.delta[warp_segment]) + 32 * warp_segment) : 0;

        if (tid < num_warps) {
            storage.delta[num_warps + tid] = start;
        }
//        /*
//        DPCT1065:74: Consider replacing sycl::nd_item::barrier() with
//        sycl::nd_item::barrier(sycl::access::fence_space::local_space) for better
//        performance if there is no access to global memory.
//        */
        item.barrier(sycl::access::fence_space::local_space);
//
//        // Find the closest flag to the left of this thread within the warp.
//        // Include the flag for this thread.
        start = 31 - clz(warp_mask & warp_bits);
        if(-1 != start) start += ~31 & tid;
        else start = storage.delta[num_warps + warp];
        /*
//        DPCT1065:75: Consider replacing sycl::nd_item::barrier() with
//        sycl::nd_item::barrier(sycl::access::fence_space::local_space) for better
//        performance if there is no access to global memory.
//        */
        item.barrier(sycl::access::fence_space::local_space);

        return start;
    }

    template<typename op_t = plus_t<type_t> >
    segscan_result_t<type_t> segscan(int tid, bool has_head_flag,
                                                 bool has_carry_out, type_t x, storage_t& storage,
                                                 const sycl::nd_item<1> &item, type_t init = type_t(),
                                                 op_t op = op_t()) const {

        if(!has_carry_out) {
            x = init;
        }

        int left_lane = find_left_lane(tid, has_head_flag, storage, item);
        int tid_delta = tid - left_lane;

        // Store the has_carry_out flag.
        storage.packed[tid] = (int)has_carry_out | (left_lane<< 1);

        // Run an inclusive scan.
        int first = 0;
        storage.values[first + tid] = x;
        /*
        DPCT1065:76: Consider replacing sycl::nd_item::barrier() with
        sycl::nd_item::barrier(sycl::access::fence_space::local_space) for better
        performance if there is no access to global memory.
        */
        item.barrier(sycl::access::fence_space::local_space);

        int packed = storage.packed[left_lane];
        left_lane = packed>> 1;
        tid_delta = tid - left_lane;

        if(0 == (1 & packed)) {
            --tid_delta;
        }

        /*
        DPCT1109:25: Recursive functions cannot be called in SYCL device code. You
        need to adjust the code.
        */
        iterate<s_log2(nt)>([&](int pass) {
            int offset = 1<< pass;
            if(tid_delta >= offset)
                x = op(x, storage.values[first + tid - offset]);
            first = nt - first;
            storage.values[first + tid] = x;
            /*
            DPCT1065:78: Consider replacing sycl::nd_item::barrier() with
            sycl::nd_item::barrier(sycl::access::fence_space::local_space) for better
            performance if there is no access to global memory.
            */
            item.barrier(sycl::access::fence_space::local_space);
        });

        // Get the exclusive scan by fetching the preceding element. Also return
        // the carry-out value as the total.
        bool has_carry_in = tid ? (0 != (1 & storage.packed[tid - 1])) : false;

        segscan_result_t<type_t> result {
                (has_carry_in && tid) ? storage.values[first + tid - 1] : init,
                storage.values[first + nt - 1],
                has_carry_in,
                left_lane
        };
        /*
        DPCT1065:77: Consider replacing sycl::nd_item::barrier() with
        sycl::nd_item::barrier(sycl::access::fence_space::local_space) for better
        performance if there is no access to global memory.
        */
        item.barrier(sycl::access::fence_space::local_space);

        return result;
    }
};