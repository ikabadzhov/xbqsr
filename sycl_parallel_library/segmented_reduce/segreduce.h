#pragma once
#include "search.h"
#include "../load_balancing/cta_load_balance.h"
#include "../segmented_scan/cta_segscan.h"
#include "../scan/cta_scan.h"
#include "../scan/scan.h"
#include "../memory/memory.h"
#include "../segmented_scan/cta_segscan.h"
#include "../device_utils/types.h"
#include "../memory/memory.h"
#include "../device_utils/meta.h"
#include "../merge/search.h"
#include <sycl/sycl.hpp>

namespace detail {

////////////////////////////////////////////////////////////////////////////////
// cta_segreduce_t is common intra-warp segmented reduction code for
// these kernels. Should clean up and move to cta_segreduce.hxx.

    template<int nt, int vt, typename type_t>
    struct cta_segreduce_t {
        typedef cta_segscan_t<nt, type_t> segscan_t;

        union storage_t {
            typename segscan_t::storage_t segscan;
            type_t values[nt * vt + 1];
        };

        // Values must be stored in storage.values on entry.
        template<typename op_t, typename output_it>
        void segreduce(merge_range_t merge_range,
                       lbs_placement_t placement, array_t<bool, vt + 1> p, int tid,
                       int cta, type_t init, op_t op, output_it output,
                       type_t* carry_out_values, int* carry_out_codes, storage_t& storage, const sycl::nd_item<1> &item) {

            int cur_item = placement.a_index;
            int begin_segment = placement.b_index;
            int cur_segment = begin_segment;
            bool carry_in = false;

            const type_t* a_shared = storage.values - merge_range.a_begin;
            type_t x[vt];
            int segments[vt + 1];
            iterate<vt>([&](int i) {
                if(p[i]) {
                    // This is a data node, so accumulate and advance the data ID.
                    x[i] = a_shared[cur_item++];
                    if(carry_in) x[i] = op(x[i - 1], x[i]);
                    carry_in = true;
                } else {
                    // This is a segment node, so advance the segment ID.
                    x[i] = init;
                    ++cur_segment;
                    carry_in = false;
                }
                segments[i] = cur_segment;
            });

            // Always flush at the end of the last thread.
            bool overwrite = (nt - 1 == tid) && (!p[vt - 1] && p[vt]);
            if(nt - 1 == tid) p[vt] = false;
            if(!p[vt]) ++cur_segment;
            segments[vt] = cur_segment;
            //overwrite = reduce_over_group(item.get_group(), overwrite, sycl::logical_or<>());//
            item.barrier(sycl::access::fence_space::local_space);
            overwrite = sycl::any_of_group(item.get_group(), overwrite);
            //overwrite = __syncthreads_or(overwrite);

            // Get the segment ID for the next item. This lets us find an end flag
            // for the last value in this thread.
            bool has_head_flag = begin_segment < segments[vt - 1];
            bool has_carry_out = p[vt - 1];

            // Compute the carry-in for each thread.
            segscan_result_t<type_t> result = segscan_t().segscan(tid, has_head_flag,
                                                                  has_carry_out, x[vt - 1], storage.segscan, item, init, op);

            // Add the carry-in back into each value and recompute the reductions.
            type_t* x_shared = storage.values - placement.range.b_begin;
            carry_in = result.has_carry_in && p[0];
            iterate<vt>([&](int i) {
                if(segments[i] < segments[i + 1]) {
                    // We've hit the end of this segment. Store the reduction to shared
                    // memory.
                    if(carry_in) x[i] = op(result.scan, x[i]);
                    x_shared[segments[i]] = x[i];
                    carry_in = false;
                }
            });
            item.barrier(sycl::access::fence_space::local_space);

            // Store the reductions for segments which begin in this tile.
            for(int i = merge_range.b_begin + tid; i < merge_range.b_end; i += nt)
                output[i] = x_shared[i];

            // Store the partial reduction for the segment which begins in the
            // preceding tile, if there is one.
            if(!tid) {
                if(segments[0] == merge_range.b_begin) {
                        segments[0] = -1;
                }

                int code = (segments[0]<< 1) | (int)overwrite;
                carry_out_values[cta] = (segments[0] != -1) ? x_shared[segments[0]] : init;
                carry_out_codes[cta] = code;
            }
        }
    };

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
void segreduce(input_it input, int count, segments_it segments,
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
    union SharedMemory {
        typename segreduce_t::storage_t segreduce;
        type_t values[nt * vt + 1];
        type_t indices[nt * vt + 2];
    } ;

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
                    //if (merge_range.b_begin==1 || merge_range.b_begin == 0)
                   // out << merge_range.a_begin << ", " << merge_range.a_begin + merge_range.a_count() << ", " << merge_range.b_begin << ", " << merge_range.b_begin + merge_range.b_count() << "\n";
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
                    segreduce_t().segreduce(merge_range, placement, merge_bits, tid, cta,
                                            init, op, output, carry_out_data, codes_data, shared[0].segreduce, item);
                });
    }).wait_and_throw();

    if(num_ctas > 1) {
        detail::segreduce_fixup(output, carry_out_data, codes_data, num_ctas,
                                op, init, q);
    }
}


////////////////////////////////////////////////////////////////////////////////
// lbs_segreduce

template<typename launch_arg_t = empty_t,
        typename func_t, typename segments_it, typename pointers_t,
        typename output_it, typename op_t, typename type_t, typename... args_t>
void lbs_segreduce(func_t f, int count, segments_it segments,
                   int num_segments, pointers_t caching_iterators, output_it output, op_t op,
                   type_t init, sycl::queue& q, args_t... args) {

    enum { nt = 128, vt = 11, vt0 = 8, nv = nt*vt };

    typedef syclgal::tuple_iterator_value_t<pointers_t> value_t;

    //cta_dim_t cta_dim = launch_t::cta_dim(context);
    int num_ctas = div_up(count + num_segments, nv);//cta_dim.num_ctas(count + num_segments);

    mem_t<type_t> carry_out(num_ctas, q);
    mem_t<int> codes(num_ctas, q);
    type_t* carry_out_data = carry_out.data();
    int* codes_data = codes.data();

    mem_t<int> mp = load_balance_partitions(count, segments, num_segments, nv, q);
    const int* mp_data = mp.data();

    typedef cta_load_balance_t<nt, vt> load_balance_t;
    typedef detail::cached_segment_load_t<nt, pointers_t> cached_load_t;
    typedef detail::cta_segreduce_t<nt, vt, type_t> segreduce_t;

    union SharedMemory{
        typename load_balance_t::storage_t lbs;
        typename cached_load_t::storage_t cached;
        typename segreduce_t::storage_t segreduce;
        type_t values[nt * vt + 1];
    } ;

    q.submit([&](sycl::handler &h) {
        // Define the local memory (shared memory)
        sycl::local_accessor<SharedMemory, 1> shared(sycl::range<1>(1), h);

        h.parallel_for<class BlockReduce>(
                sycl::nd_range<1>(sycl::range<1>(num_ctas * nt),
                                  sycl::range<1>(nt)),
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
                    array_t<type_t, vt> strided_values;
                    strided_iterate<nt, vt, vt0>([&](int i, int j) {
                        int index = lbs.merge_range.a_begin + j;
                        int seg = lbs.segments[i];
                        int rank = lbs.ranks[i];

                        strided_values[i] = f(index, seg, rank, cached_values[i], args...);
                        }, tid, lbs.merge_range.a_count()
                    );

                    // Store the values back to shared memory for segmented reduction.
                    reg_to_shared_strided<nt, vt>(strided_values, tid,
                                                  shared[0].segreduce.values, item);

                    // Split the flags.
                    array_t<bool, vt + 1> merge_bits;
                    iterate<vt + 1>([&](int i) {
                        merge_bits[i] = 0 != ((1<< i) & lbs.merge_flags);
                    });

                    // Compute the segmented reduction.
                    segreduce_t().segreduce(lbs.merge_range, lbs.placement, merge_bits,
                                            tid, cta, init, op, output, carry_out_data, codes_data,
                                            shared[0].segreduce, item);
                });
    });

    if(num_ctas > 1) {
        detail::segreduce_fixup(output, carry_out_data, codes_data, num_ctas,
                                op, init, q);
    }
}

// reduce_by_key front-end
// Not available in moderngpu new version. Refactored from version 1.1
template<typename launch_arg_t = empty_t,
        typename key_it, typename output_it, typename op_t, typename type_t, typename comp_t, typename... args_t>
void reduce_by_key_preprocessing(key_it d_keys, int count, output_it d_segments, key_it d_output_keys, int& num_segments, op_t op, comp_t comp, type_t init,
                   sycl::queue& q, args_t... args) {
    enum {
        nt = 128, vt = 7, vt0 = 5, nv = nt * vt
    };

    typedef cta_scan_t<nt, type_t> scan_t;
    typedef cta_segscan_t<nt, type_t> segscan_t;

    union SharedMemory {
        type_t keys[nt * (vt + 1)];
        typename segscan_t::storage_t indices;
        typename scan_t::storage_t scan;
    };

    int num_ctas = div_up(count, nv);

    mem_t<int> d_limits(num_ctas + 1, q);
    mem_t<int> d_limits_scan(num_ctas + 1, q);
    mem_t<int> d_thread_codes(num_ctas * nt, q);

    int* d_limits_data = d_limits.data();
    int* d_thread_codes_data = d_thread_codes.data();

    q.submit([&](sycl::handler &h) {
        // Define the local memory (shared memory)
        sycl::local_accessor<SharedMemory, 1> shared(sycl::range<1>(1), h);

        h.parallel_for(sycl::nd_range<1>(sycl::range<1>(num_ctas * nt),
                                  sycl::range<1>(nt)),
                [=](sycl::nd_item<1> item) {
                    size_t tid = item.get_local_id(); // Current work-item index
                    size_t cta = item.get_group(0); // Current work-group index
                    range_t tile = get_tile(cta, nv, count);
                    int count2 = tile.count();

                    mem_to_shared<nt, vt>(d_keys + tile.begin, tid, count2, shared[0].keys, item);

                    if(tid == 0 && tile.end < count)
                        shared[0].keys[nv] = d_keys[tile.end];

                    item.barrier(sycl::access::fence_space::local_space);

                    int endFlags = 0;
                    type_t key = shared[0].keys[vt * tid];
                    iterate<vt>([&](int i) {
                        int index = vt * tid + 1 + i;
                        // A position is valid if it lies within the tile,
                        // or if it's exactly at the tile's end (count2) and the tile is full (count2 == nv)
                        bool valid = (index < count2) || ((index == count2) && (count2 == nv));
                        // Use the extra key if we're at the tile boundary; otherwise load normally.
                        type_t next = (index < count2) ? shared[0].keys[index] : shared[0].keys[nv];
                        // If valid and the keys mismatch, set the corresponding bit.
                        endFlags |= (valid && (!comp(key, next))) ? (1 << i) : 0;
                        // Update key only if we have a valid next value.
                        key = valid ? next : key;
                    });

                    item.barrier(sycl::access::fence_space::local_space);

                    // Count the number of encountered end flags.

                    auto scan_result = scan_t().scan(tid, popc(endFlags), shared[0].scan, item);

                    int scan = scan_result.scan;
                    int total   = scan_result.reduction;

                    if(!tid) {
                        d_limits_data[cta] = total;
                    }

                    if(total) {
                        // Find the segmented scan start for this thread.
                        int left_lane = segscan_t().find_left_lane(tid, 0 != endFlags, shared[0].indices, item);
                        int tid_delta = tid - left_lane;
                        // threadCodes:
                        // 12:0 - end flags for up to 13 values per thread.
                        // 19:13 - tid delta for up to 128 threads.
                        // 30:20 - scan offset for streaming partials.
                        int threadCodes = endFlags | (tid_delta<< 13) | (scan<< 20);
                        d_thread_codes_data[nt * cta + tid] = threadCodes;
                    }
                });
    });
    q.wait_and_throw();

//    std::vector<int32_t> limit(d_limits.size());
//    q.memcpy(limit.data(), d_limits_data, (num_ctas+1) * sizeof(limit[0]));
//    q.wait_and_throw();

    // Scan the output counts.
    scan(d_limits_data, d_limits.size(), d_limits.data(), plus_t<int>(), discard_iterator_t<type_t>(), q);
    q.wait_and_throw();

    q.memcpy(&num_segments, d_limits.data() + num_ctas, 1 * sizeof(num_segments));
    q.wait_and_throw();

    num_segments++;

    union Shared {
        int indices[nv];
    };

    q.submit([&](sycl::handler &h) {
        // Define the local memory (shared memory)
        sycl::local_accessor<Shared, 1> shared(sycl::range<1>(1), h);

        h.parallel_for(sycl::nd_range<1>(sycl::range<1>(num_ctas * nt),
                                         sycl::range<1>(nt)),
                       [=](sycl::nd_item<1> item) {
                           int tid = item.get_local_id();
                           int block = item.get_group(0);
                           int gid = nv * block;

                           int limit0 = d_limits_data[block];
                           int limit1 = d_limits_data[block + 1];
                           int threadCodes = d_thread_codes_data[nt * block + tid];

                           int total = limit1 - limit0;
                           if (total) {
                               // Reconstruct row IDs from thread codes and the starting row offset.
                               int rows[vt + 1];

                               //DeviceExpandFlagsToRows<vt>(threadCodes >> 20, threadCodes, rows);
                               int first = threadCodes >> 20;
                               int endFlags = threadCodes;

                               rows[0] = first;

                               iterate<vt>([&](int i){
                                   if((1<< i) & endFlags) ++first;
                                   rows[i + 1] = first;
                               });

                               // Compact the location of the last element in each segment.
                               int index = gid + vt * tid;
                               iterate<vt>([&](int i){
                                   if (rows[i] != rows[i + 1]) {
                                       shared[0].indices[rows[i]] = index + i;
                                   }
                               });
                               item.barrier(sycl::access::fence_space::local_space);

                               // Copy a key from the start of each segment.
                               for (int i = tid; i < total; i += nt) {
                                   int pos = shared[0].indices[i] + 1;
                                   int seg = limit0 + i + 1;
                                   if (pos >= count) {
                                       pos = 0;
                                       seg = 0;
                                   }
                                   d_output_keys[seg] = d_keys[pos];
                                   d_segments[seg] = pos;//keys[pos];
                               }
                           }
                       });
    });
    q.memcpy(d_output_keys, d_keys, 1 * sizeof(d_keys[0]));
    q.wait_and_throw();
}

//template<typename launch_arg_t = empty_t,
//        typename key_it, typename val_it, typename output_vals_t, typename op_t, typename type_t, typename comp_t, typename... args_t>
//void reduce_by_key(key_it keys, val_it vals, int count, key_it output_keys, output_vals_t output_vals, int& output_sz, op_t op, comp_t comp, type_t init,
//                                 sycl::queue& q, args_t... args) {
//    mem_t<int> segments(count, q);
//    int num_segments = 0;
//
//    reduce_by_key_preprocessing(keys, count, segments.data(), output_keys, num_segments, op, comp, init, q);
//
//    // num_segments is to be computed
//    segreduce(vals, count, segments.data(), num_segments, output_vals, op, init, q);
//
//    output_sz = num_segments;
//}


template<typename launch_arg_t = empty_t,
        typename key_t, typename val_t, typename output_key_t, typename type_t, typename op_t, typename comp_t, typename... args_t>
void reduce_by_key(key_t* keys, val_t* vals, int count, std::vector<output_key_t>& output_keys, std::vector<type_t>& output_vals, int& output_sz, op_t op, comp_t comp, type_t init,
                   sycl::queue& q, memory_space_t memory_space_input_data = memory_space_t::memory_space_device , args_t... args) {

    if (memory_space_input_data == memory_space_device) {
        mem_t<key_t> d_keys(count, q);
        mem_t<val_t> d_vals(count, q);

        data_transfer(d_keys.data(), keys, count, q, false);
        data_transfer(d_vals.data(), vals, count, q);

        mem_t<int32_t> segments(count, q);
        int num_segments = 0;

        mem_t<output_key_t> d_output_keys(count, q);

        reduce_by_key_preprocessing(d_keys.data(), count, segments.data(), d_output_keys.data(), num_segments, op, comp, init, q);

        mem_t<type_t> d_output_vals(num_segments, q);

        segreduce(d_vals.data(), count, segments.data(), num_segments, d_output_vals.data(), op, init, q);

        output_keys.resize(num_segments);
        output_vals.resize(num_segments);

        data_transfer(output_keys.data(), d_output_keys.data(), num_segments, q, false);
        data_transfer(output_vals.data(), d_output_vals.data(), num_segments, q);

        output_sz = num_segments;
    }
}

template<typename launch_arg_t = empty_t,
        typename key_t, typename val_t, typename output_key_t, typename type_t, typename op_t, typename comp_t, typename... args_t>
void reduce_by_key(key_t* keys, val_t* vals, int count, std::vector<output_key_t>& output_keys, std::vector<type_t>& output_vals, int& output_sz, op_t op, comp_t comp, type_t init,
                   sycl::queue& q, std::vector<int32_t>& h_segments, memory_space_t memory_space_input_data = memory_space_t::memory_space_device , args_t... args) {

    if (memory_space_input_data == memory_space_device) {
        mem_t<key_t> d_keys(count, q);
        mem_t<val_t> d_vals(count, q);

        data_transfer(d_keys.data(), keys, count, q, false);
        data_transfer(d_vals.data(), vals, count, q);

        mem_t<int32_t> d_segments(count, q);
        int num_segments = 0;

        mem_t<output_key_t> d_output_keys(count, q);

        reduce_by_key_preprocessing(d_keys.data(), count, d_segments.data(), d_output_keys.data(), num_segments, op, comp, init, q);

        mem_t<type_t> d_output_vals(num_segments, q);

        segreduce(d_vals.data(), count, d_segments.data(), num_segments, d_output_vals.data(), op, init, q);

        output_keys.resize(num_segments);
        output_vals.resize(num_segments);

        data_transfer(output_keys.data(), d_output_keys.data(), num_segments, q, false);
        data_transfer(output_vals.data(), d_output_vals.data(), num_segments, q);

        output_sz = num_segments;
        h_segments.resize(num_segments);
        q.memcpy(h_segments.data(), d_segments.data(), sizeof(h_segments[0])*num_segments);
    }
    q.wait();
}
