#pragma once
#include <tbb/parallel_for.h>
#include <sycl/sycl.hpp>
#include "../device_utils/types.h"
#include "../device_utils/meta.h"
#include "../device_utils/operator.h"
#include "../memory/memory.h"
#include "../merge/external_merge.h"
#include "mergesort.h"

// Key-value mergesort.
template<typename launch_arg_t = empty_t, typename key_t, typename val_t,
        typename comp_t>
void external_mergesort(key_t* keys_input, val_t* vals_input, int64_t count,
                        comp_t comp, sycl::queue& q) {

    enum { has_values = !std::is_same<val_t, empty_t>::value };
    size_t sz_val = has_values ? sizeof(val_t) : 0;

    size_t global_mem_size = q.get_device().get_info<sycl::info::device::global_mem_size>();
    size_t mem_size_available = global_mem_size - pow(2,30);

    size_t sorting_mem_sz = mem_size_available / 2;
    size_t max_num_elements = sorting_mem_sz / (sizeof(key_t) + sz_val);

    size_t n_subarrays = (count + max_num_elements - 1) / max_num_elements;

    //auto start = std::chrono::high_resolution_clock::now();

    // Sorting phase
    {
        if (has_values) {
            mem_t<key_t> input_key_d(max_num_elements, q);
            mem_t<val_t> input_val_d(max_num_elements, q);
            size_t offset = 0;
            for (size_t i = 0; i < n_subarrays; i++) {
                size_t n_elements = min(max_num_elements, count - offset);
                q.memcpy(input_key_d.data(), keys_input + offset, n_elements * sizeof(key_t));
                q.memcpy(input_val_d.data(), vals_input + offset, n_elements * sizeof(val_t));
                q.wait();
                mergesort(input_key_d.data(), input_val_d.data(), n_elements, comp, q);
                q.memcpy(keys_input + offset, input_key_d.data(), n_elements * sizeof(key_t));
                q.memcpy(vals_input + offset, input_val_d.data(), n_elements * sizeof(val_t));
                q.wait();
                offset += n_elements;
            }
        }else {
            mem_t<key_t> input_d(max_num_elements, q);
            size_t offset = 0;
            for (size_t i = 0; i < n_subarrays; i++) {
                size_t n_elements = min(max_num_elements, count - offset);
                q.memcpy(input_d.data(), keys_input + offset, n_elements * sizeof(key_t));
                q.wait();
                mergesort(input_d.data(), n_elements, comp, q);
                q.memcpy(keys_input + offset, input_d.data(), n_elements * sizeof(key_t));
                q.wait();
                offset += n_elements;
            }
        }
    }

    //auto end = std::chrono::high_resolution_clock::now();
    //auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    //std::cout << "\tTime to sort chunks: " << (double)duration / 1000.0 << std::endl;

    //start = std::chrono::high_resolution_clock::now();

    // Merging phase
    if (n_subarrays > 1) {
        if (has_values) {
            std::vector<key_t> tmp_output_key(keys_input, keys_input + count);
            std::vector<val_t> tmp_output_val(vals_input, vals_input + count);

            int num_passes = find_log2(n_subarrays, true);

            key_t *key_output = tmp_output_key.data();
            val_t *val_output = tmp_output_val.data();

            if (1 & num_passes) {
                std::swap(keys_input, key_output);
                std::swap(vals_input, val_output);
            }

            size_t offset_in1 = 0;
            size_t offset_in2 = 0;
            size_t offset_out = 0;
            size_t chunk_sz1 = 0;
            size_t chunk_sz2 = 0;
            size_t tmp_n_subarrays = n_subarrays;

            size_t n_elements = max_num_elements;

            for (size_t pass = 0; pass < num_passes; pass++) {
                offset_out = 0;

                size_t merged_subarrays = 0; // Count how many subarrays are processed this pass
                for (size_t n = 0; n + 1 < tmp_n_subarrays; n += 2) {
                    offset_in1 = n * n_elements;
                    // Ensure we don't exceed the total count
                    chunk_sz1 = std::min(n_elements, count - offset_in1);

                    offset_in2 = (n + 1) * n_elements;
                    if (offset_in2 >= count) {
                        // If the second chunk doesn't exist, break and handle leftover after the loop
                        break;
                    }

                    chunk_sz2 = std::min(n_elements, count - offset_in2);

                    external_merge(keys_input + offset_in1, vals_input + offset_in1, chunk_sz1,
                                   keys_input + offset_in2, vals_input + offset_in2, chunk_sz2,
                                   key_output + offset_out, val_output + offset_out, comp, q, mem_size_available);

                    offset_out += (chunk_sz1 + chunk_sz2);
                    merged_subarrays += 2;
                }

                // Check if there's a leftover subarray (in case tmp_n_subarrays is odd)
                if (merged_subarrays < tmp_n_subarrays) {
                    size_t leftover_index = merged_subarrays;
                    size_t leftover_offset = leftover_index * n_elements;
                    if (leftover_offset < count) {
                        // Copy the leftover subarray as-is
                        size_t leftover_count = std::min(n_elements, count - leftover_offset);
                        std::memcpy(key_output + offset_out, keys_input + leftover_offset,
                                    leftover_count * sizeof(key_t));
                        std::memcpy(val_output + offset_out, vals_input + leftover_offset,
                                    leftover_count * sizeof(val_t));
                        offset_out += leftover_count;
                        merged_subarrays += 1;
                    }
                }

                // Each pass doubles the size of the subarrays we merge
                n_elements *= 2;
                // Half the number of subarrays (rounding up if odd)
                tmp_n_subarrays = (merged_subarrays + 1) / 2;

                std::swap(keys_input, key_output);
                std::swap(vals_input, val_output);
            }
        } else {
            std::vector<key_t> tmp_output(keys_input, keys_input + count);
            int num_passes = find_log2(n_subarrays, true);

            key_t *key_output = tmp_output.data();
            if (1 & num_passes) {
                std::swap(keys_input, key_output);
            }

            size_t offset_in1 = 0;
            size_t offset_in2 = 0;
            size_t offset_out = 0;
            size_t chunk_sz1 = 0;
            size_t chunk_sz2 = 0;
            size_t tmp_n_subarrays = n_subarrays;

            size_t n_elements = max_num_elements;

            for (size_t pass = 0; pass < num_passes; pass++) {
                offset_out = 0;

                size_t merged_subarrays = 0; // Count how many subarrays are processed this pass
                for (size_t n = 0; n + 1 < tmp_n_subarrays; n += 2) {
                    offset_in1 = n * n_elements;
                    // Ensure we don't exceed the total count
                    chunk_sz1 = std::min(n_elements, count - offset_in1);

                    offset_in2 = (n + 1) * n_elements;
                    if (offset_in2 >= count) {
                        // If the second chunk doesn't exist, break and handle leftover after the loop
                        break;
                    }

                    chunk_sz2 = std::min(n_elements, count - offset_in2);

                    external_merge(keys_input + offset_in1, chunk_sz1,
                                   keys_input + offset_in2, chunk_sz2,
                                   key_output + offset_out, comp, q, mem_size_available);

                    offset_out += (chunk_sz1 + chunk_sz2);
                    merged_subarrays += 2;
                }

                // Check if there's a leftover subarray (in case tmp_n_subarrays is odd)
                if (merged_subarrays < tmp_n_subarrays) {
                    size_t leftover_index = merged_subarrays;
                    size_t leftover_offset = leftover_index * n_elements;
                    if (leftover_offset < count) {
                        // Copy the leftover subarray as-is
                        size_t leftover_count = std::min(n_elements, count - leftover_offset);
                        std::memcpy(key_output + offset_out, keys_input + leftover_offset,
                                    leftover_count * sizeof(key_t));
                        offset_out += leftover_count;
                        merged_subarrays += 1;
                    }
                }

                // Each pass doubles the size of the subarrays we merge
                n_elements *= 2;
                // Half the number of subarrays (rounding up if odd)
                tmp_n_subarrays = (merged_subarrays + 1) / 2;

                std::swap(keys_input, key_output);
            }
        }
    }

    //end = std::chrono::high_resolution_clock::now();
    //duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    //std::cout << "\tTime to merge chunks: " << (double)duration / 1000.0 << std::endl;
}


// Key-only mergesort
template<typename launch_arg_t = empty_t, typename key_t, typename comp_t>
void external_mergesort(key_t* keys_input, int64_t count, comp_t comp, sycl::queue& q) {
    external_mergesort<launch_arg_t>(keys_input, (empty_t*)nullptr, count, comp, q);
}
