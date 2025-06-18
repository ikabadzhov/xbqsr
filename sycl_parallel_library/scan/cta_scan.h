#pragma once
#include "../memory/loadstore.h"
#include "../device_utils/intrinsics.h"
#include "../device_utils/types.h"
#include "../device_utils/meta.h"

enum scan_type_t {
    scan_type_exc,
    scan_type_inc
};

template<typename type_t, int vt = 0, bool is_array = (vt > 0)>
struct scan_result_t {
    type_t scan;
    type_t reduction;
};

template<typename type_t, int vt>
struct scan_result_t<type_t, vt, true> {
    array_t<type_t, vt> scan;
    type_t reduction;
};

////////////////////////////////////////////////////////////////////////////////

template<int nt, typename type_t>
struct cta_scan_t {
    enum { num_warps = nt / warp_size, capacity = nt + num_warps };
    union storage_t {
        type_t data[2 * nt];
        struct { type_t threads[nt], warps[num_warps]; };
    };

    //////////////////////////////////////////////////////////////////////////////
    // Standard CTA scan code that does not use shfl intrinsics.

    template<typename op_t = plus_t<type_t> >
    SYCL_EXTERNAL scan_result_t<type_t> scan(int tid, type_t x, storage_t& storage, sycl::nd_item<1>& item, int count = nt, op_t op = op_t(),
         type_t init = type_t(), scan_type_t type = scan_type_exc) const {

        int first = 0;
        storage.data[first + tid] = x;

        item.barrier(sycl::access::fence_space::local_space);

        iterate<s_log2(nt)>([&](int pass) {
            int offset = 1<< pass;
            if(tid >= offset)
                x = op(storage.data[first + tid - offset], x);
            first = nt - first;
            storage.data[first + tid] = x;
            item.barrier(sycl::access::fence_space::local_space);
        });

        scan_result_t<type_t> result;
        result.reduction = storage.data[first + count - 1];
        result.scan = (tid < count) ?
                      (scan_type_inc == type ? x :
                       (tid ? storage.data[first + tid - 1] : init)) :
                      result.reduction;

        item.barrier(sycl::access::fence_space::local_space);

        return result;
    }

    //////////////////////////////////////////////////////////////////////////////
    // CTA vectorized scan. Accepts multiple values per thread and adds in
    // optional global carry-in.

    template<int vt, typename op_t = plus_t<type_t> >
    SYCL_EXTERNAL scan_result_t<type_t, vt> scan(int tid, array_t<type_t, vt> x, storage_t& storage, sycl::nd_item<1>& item,
         type_t carry_in = type_t(), bool use_carry_in = false,
         int count = nt, op_t op = op_t(), type_t init = type_t(),
         scan_type_t type = scan_type_exc) const {

        // Start with an inclusive scan of the in-range elements.
        if(count >= nt * vt) {
            iterate<vt>([&](int i) {
                x[i] = i ? op(x[i], x[i - 1]) : x[i];
            });
        } else {
            iterate<vt>([&](int i) {
                int index = vt * tid + i;
                x[i] = i ?
                       ((index < count) ? op(x[i], x[i - 1]) : x[i - 1]) :
                       (x[i] = (index < count) ? x[i] : init);
            });
        }

        // Scan the thread-local reductions for a carry-in for each thread.
        scan_result_t<type_t> result = scan(tid, x[vt - 1], storage, item,
                                            div_up(count, vt), op, init, scan_type_exc);

        // Perform the scan downsweep and add both the global carry-in and the
        // thread carry-in to the values.
        if(use_carry_in) {
            result.reduction = op(carry_in, result.reduction);
            result.scan = tid ? op(carry_in, result.scan) : carry_in;
        } else
            use_carry_in = tid > 0;

        array_t<type_t, vt> y;
        iterate<vt>([&](int i) {
            if(scan_type_exc == type) {
                y[i] = i ? x[i - 1] : result.scan;
                if(use_carry_in && i > 0) y[i] = op(result.scan, y[i]);
            } else {
                y[i] = use_carry_in ? op(x[i], result.scan) : x[i];
            }
        });

        return scan_result_t<type_t, vt> { y, result.reduction };
    }
};

////////////////////////////////////////////////////////////////////////////////
// Overload for scan of bools.

template<int nt>
struct cta_scan_t<nt, bool> {
    enum { num_warps = nt / warp_size };
    struct storage_t {
        int warps[num_warps];
    };

    scan_result_t<int> scan(int tid, bool x, storage_t& storage, sycl::nd_item<1>& item) const {

        // Store the bit totals for each warp.
        int lane = (warp_size - 1) & tid;
        int warp = tid / warp_size;

        //int bits = ballot(x);
        auto sg = item.get_sub_group();
        int bits = static_cast<int>(sycl::reduce_over_group(sg, x, sycl::ext::oneapi::plus<>()));

        storage.warps[warp] = popc(bits);
        item.barrier(sycl::access::fence_space::local_space);

#if defined(__CUDA_ARCH__) && __CUDA_ARCH__ >= 300
        if(tid < num_warps) {
      // Cooperative warp scan of partial reductions.
      int scan = storage.warps[tid];
      iterate<s_log2(num_warps)>([&](int i) {
        scan = shfl_up_op(scan, 1<< i, plus_t<int>(), num_warps);
      });
      storage.warps[tid] = scan;
    }
    __syncthreads();
#else

        if(0 == tid) {
            // Inclusive scan of partial reductions..
            int scan = 0;
            iterate<num_warps>([&](int i) {
                storage.warps[i] = scan += storage.warps[i];
            });
        }
        item.barrier(sycl::access::fence_space::local_space);

#endif

        int scan = ((warp > 0) ? storage.warps[warp - 1] : 0) +
                   popc(bfe(bits, 0, lane));
        int reduction = storage.warps[num_warps - 1];
        item.barrier(sycl::access::fence_space::local_space);


        return scan_result_t<int> { scan, reduction };
    }
};

