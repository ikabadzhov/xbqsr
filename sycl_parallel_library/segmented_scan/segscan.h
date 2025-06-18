#pragma once

#include <sycl/sycl.hpp>
#include "../memory/memory.h"
#include "cta_segscan.h"
/*
namespace syclgal_segscan {
    ////////////////////////////////////////////////////////////////////////////////
    // Adds the carry-out for each segreduce CTA into the outputs.

    template<typename output_it, typename type_t, typename op_t>
    void segreduce_fixup(output_it output, const type_t* values,
                         const int* codes, int count, op_t op, type_t init,
                         sycl::queue& q) {

        enum { nt = 512 };
        int num_ctas = div_up(count, nt);

        mem_t<type_t> carry_out(num_ctas, q);
        mem_t<int> codes_out(num_ctas, q);
        type_t* carry_out_data = carry_out.data();
        int* codes_data = codes_out.data();

        typedef cta_segscan_t<nt, type_t> segscan_t;

        struct SharedMemory{
            bool head_flags[nt];
            typename segscan_t::storage_t segscan;
        } ;

        q.submit([&](sycl::handler &h) {
            // Define the local memory (shared memory)
            sycl::local_accessor<SharedMemory, 1> shared(sycl::range<1>(1), h);

            h.parallel_for(
                    sycl::nd_range<1>(sycl::range<1>(num_ctas * nt),
                                      sycl::range<1>(nt)),
                    [=](sycl::nd_item<1> item) {
                        size_t tid = item.get_local_id(); // Current work-item index
                        size_t cta = item.get_group(0); // Current work-group index
                        //typedef cta_segscan_t<nt, type_t> segscan_t;


                        range_t tile = get_tile(cta, nt, count);
                        int gid = tile.begin + tid;

                        ////////////////////////////////////////////////////////////////////////////
                        // As in the outer segmented reduce kernel, update the reductions for all
                        // segments that *start* in this CTA. That is, the first carry-out code
                        // for a segment must be mapped into this CTA to actually apply the
                        // accumulate. This CTA will return a partial reduction for the segment
                        // that overlaps this CTA but starts in a preceding CTA.

                        // We don't need to worry about storing new overwrite bits as this kernel
                        // will always add carry-in values to empty segments.

                        int code0 = (gid - 1 >= 0 && gid - 1 < count) ? codes[gid - 1] : -1;
                        int code1 = (gid < count) ? codes[gid] : -1;
                        int code2 = (gid + 1 < count) ? codes[gid + 1] : -1;
                        type_t value = (gid < count) ? values[gid] : init;

                        int seg0 = code0 >> 1;
                        int seg1 = code1 >> 1;
                        int seg2 = code2 >> 1;
                        bool has_head_flag = seg0 != seg1 || -1 == seg1;
                        bool has_carry_out = -1 != seg1 && seg1 == seg2;
                        bool has_end_flag = seg1 != seg2;

                        // Put the head flag in shared memory, because the last thread
                        // participating in a reduction in the CTA needs to check the head flag
                        // for the first thread in the reduction.
                        shared[0].head_flags[tid] = has_head_flag;

                        segscan_result_t<type_t> result = segscan_t().segscan(tid, has_head_flag,
                                                                              has_carry_out, value, shared[0].segscan, item,
                                                                              init, op);

                        bool carry_out_written = false;
                        if (-1 != seg1 && (has_end_flag || nt - 1 == tid)) {
                            // This is a valid reduction.
                            if (result.has_carry_in)
                                value = op(value, result.scan);

                            if (0 == result.left_lane && !shared[0].head_flags[result.left_lane]) {
                                carry_out_data[cta] = value;
                                codes_data[cta] = seg1 << 1;
                                carry_out_written = true;
                            } else {
                                int left_code = codes[tile.begin + result.left_lane - 1];
                                if (0 == (1 & left_code))     // Add in the value already stored.
                                    value = op(value, output[seg1]);
                                output[seg1] = value;
                            }
                        }
                        item.barrier(sycl::access::fence_space::local_space);
                        carry_out_written = sycl::any_of_group(item.get_group(), carry_out_written);
                        //__syncthreads_or(carry_out_written);
                        if (!carry_out_written && !tid)
                            codes_data[cta] = -1 << 1;
                    });
        });

        q.wait_and_throw();

        if(num_ctas > 1) {
            segreduce_fixup(output, carry_out_data, codes_data,
                            num_ctas, op, init, q);
        }
    }
} // namespace detail

////////////////////////////////////////////////////////////////////////////////
// Segmented reduction with loading from an input iterator. This does not
// require explicit materialization of the load-balancing search.

template<typename launch_arg_t = empty_t, typename input_it,
        typename segments_it, typename output_it, typename op_t, typename type_t>
void segscan(input_it input, int count, segments_it segments,
               int num_segments, output_it output, op_t op, type_t init,
               sycl::queue& q) {
    //enum { nt = 128, vt = 11, nv = nt*vt };
    enum { nt = 128, vt = 7, vt0 = 5, nv = nt*vt };

    int num_ctas = div_up(count + num_segments, nv);

    mem_t<type_t> carry_out(num_ctas, q);
    mem_t<int> codes(num_ctas, q);
    type_t* carry_out_data = carry_out.data();
    int* codes_data = codes.data();

    mem_t<int> mp = load_balance_partitions(count, segments, num_segments, nv, q);
    const int* mp_data = mp.data();

//    std::vector<int> tmp_mp(mp.size(), 0);
//    q.memcpy(tmp_mp.data(), mp.data(), sizeof(int)*mp.size()).wait();
//    for (size_t i = 0; i < mp.size(); i++) {
//        std::cout << tmp_mp[i] << std::endl;
//    }

    typedef detail::cta_segreduce_t<nt, vt, type_t> segreduce_t;
    typedef cta_segscan_t<nt, type_t> segscan_t;

    union SharedMemory {
        typename segreduce_t::storage_t segreduce;
        type_t values[nt * vt + 1];
        type_t indices[nt * vt + 2];
    };

    q.submit([&](sycl::handler &h) {
        // Define the local memory (shared memory)
        sycl::local_accessor<SharedMemory, 1> shared(sycl::range<1>(1), h);
        //sycl::stream out(1024*1024*1024, 64, h);

        h.parallel_for(
                sycl::nd_range<1>(sycl::range<1>(num_ctas * nt),
                                  sycl::range<1>(nt)),
                [=](sycl::nd_item<1> item) {
                    size_t tid = item.get_local_id(); // Current work-item index
                    size_t cta = item.get_group(0); // Current work-group index

                    merge_range_t merge_range = compute_merge_range(count, num_segments,
                                                                    cta, nt * vt, mp_data[cta], mp_data[cta + 1]);

                    // Cooperatively load values from input into shared.
                    mem_to_shared<nt, vt, vt0>(input + merge_range.a_begin, tid,
                                               merge_range.a_count(), shared[0].segreduce.values, item);

                    // Load segment data into the B region of shared. Search for the starting
                    // index of each thread for a merge.
                    int* b_shared = sizeof(type_t) > sizeof(int) ?
                                    (int*)(shared[0].segreduce.values + merge_range.a_count()) :
                                    ((int*)shared[0].segreduce.values + merge_range.a_count());
                    lbs_placement_t placement = cta_load_balance_place<nt, vt>(tid,
                                                                               merge_range, count, segments, num_segments, b_shared, item);

                    // Adjust the pointer so that dereferencing at the segment ID returns the
                    // offset of that segment.
                    b_shared -= placement.range.b_begin;
                    int cur_item = placement.a_index;
                    int cur_segment = placement.b_index;
                    array_t<bool, vt + 1> merge_bits;
                    iterate<vt + 1>([&](int i) {
                        bool p = cur_item < b_shared[cur_segment + 1];
                        if(p) ++cur_item;
                        else ++cur_segment;
                        merge_bits[i] = p;
                    });

                    // Compute the segmented reduction.
                    segscan_t().segscan(tid, has_head_flag,
                                                                  has_carry_out, x[vt - 1], storage.segscan, item, init, op);
                });
    }).wait_and_throw();

    if(num_ctas > 1) {
        syclgal_segscan::segreduce_fixup(output, carry_out_data, codes_data, num_ctas,
                                op, init, q);
    }
}*/