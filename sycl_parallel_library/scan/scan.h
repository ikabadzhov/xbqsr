#pragma once

#include "cta_scan.h"
#include "../reduce/cta_reduce.h"
#include "../memory/memory.h"
#include "../device_utils/iterators.h"

template<scan_type_t scan_type = scan_type_exc,
        typename launch_arg_t = empty_t, typename input_it,
        typename output_it, typename op_t, typename reduction_it>
void scan_event(input_it input, int count, output_it output, op_t op,
                reduction_it reduction, sycl::queue& q) {

    enum {nt = 128, vt = 7, nv = nt*vt};

    typedef typename std::iterator_traits<input_it>::value_type type_t;
    typedef cta_reduce_t<nt, type_t> reduce_t;
    typedef cta_scan_t<nt, type_t> scan_t;


    int num_ctas = div_up(count, nv);

    if(num_ctas > 8) {
        mem_t<type_t> partials(num_ctas, q);
        type_t* partials_data = partials.data();

        union SharedMemory{
                typename reduce_t::storage_t reduce;
        };

        ////////////////////////////////////////////////////////////////////////////
        // Upsweep phase. Reduce each tile to a scalar and store to partials.

        q.submit([&](sycl::handler &h) {
            // Define the local memory (shared memory)
            sycl::local_accessor<SharedMemory, 1> shared(sycl::range<1>(1), h);

            h.parallel_for(
                    sycl::nd_range<1>(sycl::range<1>(num_ctas * nt),
                                      sycl::range<1>(nt)),
                    [=](sycl::nd_item<1> item) {
                        size_t tid = item.get_local_id(); // Current work-item index
                        size_t cta = item.get_group(0); // Current work-group index

                        // Load the tile's data into register.
                        range_t tile = get_tile(cta, nv, count);
                        array_t<type_t, vt> x = mem_to_reg_strided<nt, vt>(input + tile.begin,
                                                                           tid, tile.count());

                        // Reduce the thread's values into a scalar.
                        type_t scalar = type_t();
                        strided_iterate<nt, vt>([&](int i, int j) {
                            scalar = i ? op(scalar, x[i]) : x[0];
                        }, tid, tile.count());

                        // Reduce across all threads.
                        type_t all_reduce = reduce_t().reduce(tid, scalar, shared[0].reduce, item,
                                                              tile.count(), op);

                        // Store the final reduction to the partials.
                        if(!tid)
                            partials_data[cta] = all_reduce;
                    });
        });
        q.wait_and_throw();
//        std::vector<type_t> h_partials(partials.size());
//        q.memcpy(h_partials.data(), partials.data(), sizeof (type_t) * partials.size()).wait();
//        for (size_t i = 0; i < h_partials.size(); i++) {
//            std::cout << "cta: "<< i << ", partial: " <<  h_partials[i] << std::endl;
//        }
        ////////////////////////////////////////////////////////////////////////////
        // Spine phase. Recursively call scan on the CTA partials.

        scan_event<scan_type_exc>(partials_data, num_ctas, partials_data,
                                  op, reduction, q);

        q.wait_and_throw();

        // Record the event. This lets the caller wait on just the reduction
        // part of the operation. It's useful when writing the reduction to
        // host-side paged-locked memory; the caller can read out the value more
        // quickly to allocate memory and launch the next kernel.
//        if(event)
//            cudaEventRecord(event, context.stream());

        union SharedMemory2 {
            typename scan_t::storage_t scan;
            type_t values[nv];
        };
        ////////////////////////////////////////////////////////////////////////////
        // Downsweep phase. Perform an intra-tile scan and add the scan of the
        // partials as carry-in.
        q.submit([&](sycl::handler &h) {
            // Define the local memory (shared memory)
            sycl::local_accessor<SharedMemory2, 1> shared(sycl::range<1>(1), h);

            h.parallel_for(
                    sycl::nd_range<1>(sycl::range<1>(num_ctas * nt),
                                      sycl::range<1>(nt)),
                    [=](sycl::nd_item<1> item) {
                        size_t tid = item.get_local_id(); // Current work-item index
                        size_t cta = item.get_group(0); // Current work-group index

                        // Load a tile to register in thread order.
                        range_t tile = get_tile(cta, nv, count);
                        array_t<type_t, vt> x = mem_to_reg_thread<nt, vt>(input + tile.begin,
                                                                          tid, tile.count(), shared[0].values, item);

                        // Scan the array with carry-in from the partials.
                        array_t<type_t, vt> y = scan_t().scan(tid, x, shared[0].scan, item,
                                                              partials_data[cta], cta > 0, tile.count(), op, type_t(),
                                                              scan_type).scan;

                        // Store the scanned values to the output.
                        reg_to_mem_thread<nt, vt>(y, tid, tile.count(), output + tile.begin,
                                                  shared[0].values, item);
                    });
        });
        q.wait_and_throw();

    } else {
//        typedef launch_params_t<512, 3> spine_params_t;

        enum { nt = 512, vt = 3, nv = nt * vt };
        typedef cta_scan_t<nt, type_t> scan_t;

        int num_ctas = 1;

        union SharedMemory{
            typename scan_t::storage_t scan;
            type_t values[nv];
        } ;

        ////////////////////////////////////////////////////////////////////////////
        // Small input specialization. This is the non-recursive branch.
        q.submit([&](sycl::handler &h) {
            // Define the local memory (shared memory)
            sycl::local_accessor<SharedMemory, 1> shared(sycl::range<1>(1), h);

            h.parallel_for(
                    sycl::nd_range<1>(sycl::range<1>(num_ctas * nt),
                                      sycl::range<1>(nt)),
                    [=](sycl::nd_item<1> item) {
                        size_t tid = item.get_local_id();
                        type_t carry_in = type_t();
                        for(int cur = 0; cur < count; cur += nv) {
                            // Cooperatively load values into register.
                            int count2 = min<int>(count - cur, nv);
                            array_t<type_t, vt> x = mem_to_reg_thread<nt, vt>(input + cur,
                                                                              tid, count2, shared[0].values, item);

                            scan_result_t<type_t, vt> result = scan_t().scan(tid, x, shared[0].scan, item,
                                                                             carry_in, cur > 0, count2, op, type_t(), scan_type);

                            // Store the scanned values back to global memory.
                            reg_to_mem_thread<nt, vt>(result.scan, tid, count2,
                                                      output + cur, shared[0].values, item);

                            // Roll the reduction into carry_in.
                            carry_in = result.reduction;
                        }

                        // Store the carry-out to the reduction pointer. This may be a
                        // discard_iterator_t if no reduction is wanted.
                        if(!tid)
                            *reduction = carry_in;
                    });
        });
       q.wait_and_throw();

        // Record the event. This lets the caller wait on just the reduction
        // part of the operation. It's useful when writing the reduction to
        // host-side paged-locked memory; the caller can read out the value more
        // quickly to allocate memory and launch the next kernel.
//        if(event)
//            cudaEventRecord(event, context.stream());
    }
}

template<scan_type_t scan_type = scan_type_exc,
        typename launch_arg_t = empty_t, typename input_it,
        typename output_it, typename op_t, typename reduction_it>
void scan(input_it input, int count, output_it output, op_t op,
          reduction_it reduction, sycl::queue& q) {
    return scan_event<scan_type, launch_arg_t>(input, count, output, op,
                                               reduction, q);
}

template<scan_type_t scan_type = scan_type_exc,
        typename launch_arg_t = empty_t,
        typename input_it, typename output_it>
void scan(input_it input, int count, output_it output, sycl::queue& q) {

    typedef typename std::iterator_traits<input_it>::value_type type_t;
    scan<scan_type, launch_arg_t>(input, count, output, plus_t<type_t>(),
                                  discard_iterator_t<type_t>(), q);
}

template<typename type_t, scan_type_t scan_type = scan_type_exc,
        typename launch_arg_t = empty_t, typename func_t, typename output_it,
        typename op_t, typename reduction_it>
void transform_scan_event(func_t f, int count, output_it output, op_t op,
                          reduction_it reduction, sycl::queue& q) {

    scan_event<scan_type, launch_arg_t>(make_load_iterator<type_t>(f),
                                        count, output, op, reduction, q);
}

template<typename type_t, scan_type_t scan_type = scan_type_exc,
        typename launch_arg_t = empty_t, typename func_t, typename output_it,
        typename op_t, typename reduction_it>
void transform_scan(func_t f, int count, output_it output, op_t op,
                    reduction_it reduction, sycl::queue& q) {

    transform_scan_event<type_t, scan_type, launch_arg_t>(f, count, output, op,
                                                          reduction, q);
}