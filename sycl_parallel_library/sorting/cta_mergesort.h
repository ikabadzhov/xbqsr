#ifndef TESTONEAPI_CTA_MERGESORT_H
#define TESTONEAPI_CTA_MERGESORT_H

#include <sycl/sycl.hpp>
#include "../device_utils/meta.h"
#include "../device_utils/types.h"
#include "sort_network.h"
#include "../memory/loadstore.h"
#include "cta_merge.h"

SYCL_EXTERNAL int out_of_range_flags(int first, int vt, int count);

SYCL_EXTERNAL merge_range_t compute_mergesort_frame(int partition, int coop, int spacing);

SYCL_EXTERNAL merge_range_t compute_mergesort_range(int count, int partition, int coop, int spacing);

SYCL_EXTERNAL merge_range_t compute_mergesort_range(int count, int partition, int coop, int spacing, int mp0, int mp1);

template<int nt, int vt, typename key_t, typename val_t>
struct cta_sort_t {
    enum {
        has_values = !std::is_same<val_t, empty_t>::value,
        num_passes = s_log2(nt)
    };

    union storage_t {
        key_t keys[nt * vt + 1];
        val_t vals[nt * vt];
    };

    static_assert(is_pow2(nt), "cta_sort_t requires pow2 number of threads");

    template<typename comp_t>
    kv_array_t<key_t, val_t, vt> merge_pass(kv_array_t<key_t, val_t, vt> x, int tid, int count, int pass, comp_t comp, storage_t& storage, sycl::nd_item<1>& item) const {

        // Divide the CTA's keys into lists.
        int coop = 2<< pass;
        merge_range_t range = compute_mergesort_range(count, tid, coop, vt);
        int diag = vt * tid - range.a_begin;

        // Store the keys into shared memory for searching.
        reg_to_shared_thread<nt, vt>(x.keys, tid, storage.keys, item);

        // Search for the merge path for this thread within its list.
        int mp = merge_path<bounds_lower>(storage.keys, range, diag, comp);

        // Run a serial merge and return.
        merge_pair_t<key_t, vt> merge = serial_merge<bounds_lower, vt>(
                storage.keys, range.partition(mp, diag), comp, item);
        x.keys = merge.keys;

        if(has_values) {
            // Reorder values through shared memory.
            reg_to_shared_thread<nt, vt>(x.vals, tid, storage.vals, item);
            x.vals = shared_gather<nt, vt>(storage.vals, merge.indices, item);
        }

        return x;
    }

    template<typename comp_t>
    kv_array_t<key_t, val_t, vt>
    block_sort(kv_array_t<key_t, val_t, vt> x, int tid, int count, comp_t comp, storage_t& storage, sycl::nd_item<1>& item) const {

        // Sort the inputs within each thread. If any threads have fewer than
        // vt items, use the segmented sort network to prevent out-of-range
        // elements from contaminating the sort.
        if(count < nt * vt) {
            int head_flags = out_of_range_flags(vt * tid, vt, count);
            x = odd_even_sort(x, comp, head_flags);
        } else
            x = odd_even_sort(x, comp);

        // Merge threads starting with a pair until all values are merged.
        for(int pass = 0; pass < num_passes; ++pass)
            x = merge_pass(x, tid, count, pass, comp, storage, item);

        return x;
    }
};

#endif //TESTONEAPI_CTA_MERGESORT_H
