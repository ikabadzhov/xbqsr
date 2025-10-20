#include "ApplyBqsrKernelData.h"

#include "../BaseRecalibratorTablesAccessor.h"
#include "../../timer/Timer.h"
#include "../../logger/Logger.h"


static size_t global_chunk_offset = 0;


size_t compute_chunk_sz(const ApplyBqsrSequencesInput& in, const AppyBqsrKernelData& kernel_data, size_t start_idx, size_t n_elements);

ApplyBqsrSequencesInput::ApplyBqsrSequencesInput(bqsrEntry* quality_table, bqsrEntry* cycle_table, bqsrEntry* context_table, size_t quality_table_sz, size_t cycle_table_sz, size_t context_table_sz, BaseRecalibratorTablesAccessor& table_acc,
                            char* sequences, SpanField* sequence_ptrs, uint8_t* qualities, uint16_t* flags, int32_t* read_group, combinedBqsrEntry* combined_bqsr_table, size_t combined_bqsr_table_size,
                            uint8_t* quantized_quality_scores, size_t quantized_quality_scores_size, uint8_t* static_quantized_quality_scores, size_t static_quantized_quality_scores_size, int* recalibrated_score_cache, size_t n_elements) : InputData(n_elements) {

    d_recalibrated_score_cache = recalibrated_score_cache;
    h_quality_table = quality_table;
    h_cycle_table = cycle_table;
    h_context_table = context_table;

    sz_qual_score_cov = quality_table_sz;


    sz_cycle_cov = cycle_table_sz;


    sz_context_cov = context_table_sz;

    h_table_acc = table_acc;

    h_sequences = sequences;
    h_sequences_ptr = sequence_ptrs;
    h_qualities = qualities;
    h_flags = flags;
    h_read_group = read_group;
    h_combined_bqsr_table = combined_bqsr_table;
    h_quantized_quality_scores = quantized_quality_scores;
    h_static_quantized_quality_scores = static_quantized_quality_scores;
    this->static_quantized_quality_scores_size = static_quantized_quality_scores_size;
    this->quantized_quality_scores_size = quantized_quality_scores_size;
    this->combined_bqsr_table_size = combined_bqsr_table_size;

}

void AppyBqsrKernelData::swap(AppyBqsrKernelData& rhs){
    std::swap(_q, rhs._q);
    std::swap(memory_buffer, rhs.memory_buffer);
    std::swap(max_memory, rhs.max_memory);
    std::swap(offset_of_constant_allocation, rhs.offset_of_constant_allocation);
    //std::swap(chunk_offset, rhs.chunk_offset);

    std::swap(d_sequences, rhs.d_sequences);
    std::swap(d_sequences_ptr, rhs.d_sequences_ptr);

    std::swap(d_qualities, rhs.d_qualities);

    std::swap(d_flag, rhs.d_flag);
    std::swap(d_reads_group_covariate, rhs.d_reads_group_covariate);

    std::swap(d_stranded_seq, rhs.d_stranded_seq);
    std::swap(d_stranded_seq_ptr, rhs.d_stranded_seq_ptr);
    std::swap(d_base_context, rhs.d_base_context);
    std::swap(d_base_context_ptr, rhs.d_base_context_ptr);

    std::swap(d_table_acc, rhs.d_table_acc);

    std::swap(chunk_sz, rhs.chunk_sz);
    std::swap(prev_batch_max_index, rhs.prev_batch_max_index);

    std::swap(d_quality_table, rhs.d_quality_table);
    std::swap(d_cycle_table, rhs.d_cycle_table);
    std::swap(d_context_table, rhs.d_context_table);

    std::swap(d_combined_bqsr_table, rhs.d_combined_bqsr_table);
    std::swap(combined_bqsr_table_size, rhs.combined_bqsr_table_size);

    std::swap(d_quantized_quality_scores, rhs.d_quantized_quality_scores);
    std::swap(quantized_quality_scores_size, rhs.quantized_quality_scores_size);

    std::swap(d_static_quantized_scores, rhs.d_static_quantized_scores);
    std::swap(static_quantized_scores_size, rhs.static_quantized_scores_size);

    std::swap(d_temporary_segments_sequences_offsets, rhs.d_temporary_segments_sequences_offsets);
    std::swap(d_recalibrated_score_cache, rhs.d_recalibrated_score_cache);

}


AppyBqsrKernelData::~AppyBqsrKernelData() {
    sycl::free(memory_buffer, _q);
    memory_buffer = nullptr;
    d_sequences = nullptr;
    d_sequences_ptr = nullptr;
    d_qualities = nullptr;
    d_flag = nullptr;
    d_reads_group_covariate = nullptr;

    d_stranded_seq = nullptr;
    d_stranded_seq_ptr = nullptr;
    d_base_context = nullptr;
    d_base_context_ptr = nullptr;

    d_quality_table = nullptr;
    d_cycle_table = nullptr;
    d_context_table = nullptr;
    d_combined_bqsr_table = nullptr;

    d_quantized_quality_scores = nullptr;
    d_static_quantized_scores = nullptr;
    d_temporary_segments_sequences_offsets = nullptr;

    max_memory = 0;
    chunk_sz = 0;
    prev_batch_max_index = 0;
    static_quantized_scores_size = 0;
}

void AppyBqsrKernelData::memcpy_host_to_dev(sycl::queue& transfer_queue, AppyBqsrKernelData& dest, size_t dest_offset, const ApplyBqsrSequencesInput& src,
                                            size_t src_offset, size_t n_elements) {
    Timer t;
    t.reset();
    try {
        char *base = static_cast<char *>(dest.memory_buffer);
        size_t offset = 0;

        if (dest.is_first_copy) {
            dest.d_recalibrated_score_cache = src.d_recalibrated_score_cache;

            dest.d_quality_table = reinterpret_cast<bqsrEntry*>(base + offset);
            size_t sz = src.sz_qual_score_cov * sizeof(src.h_quality_table[0]);
            transfer_queue.memcpy(dest.d_quality_table, src.h_quality_table, sz);
            offset += sz;

            size_t alignment = alignof(decltype(dest.d_cycle_table[0]));
            size_t padding = (alignment - (offset % alignment)) % alignment;
            offset += padding;
            dest.d_cycle_table = reinterpret_cast<bqsrEntry*>(base + offset);
            sz = src.sz_cycle_cov * sizeof(src.h_cycle_table[0]);
            transfer_queue.memcpy(dest.d_cycle_table, src.h_cycle_table, sz);
            offset += sz;

            alignment = alignof(decltype(dest.d_context_table[0]));
            padding = (alignment - (offset % alignment)) % alignment;
            offset += padding;
            dest.d_context_table = reinterpret_cast<bqsrEntry*>(base + offset);
            sz = src.sz_context_cov * sizeof(src.h_context_table[0]);
            transfer_queue.memcpy(dest.d_context_table, src.h_context_table, src.sz_context_cov * sizeof(src.h_context_table[0]));
            offset += sz;

            alignment = alignof(decltype(dest.d_combined_bqsr_table[0]));
            padding = (alignment - (offset % alignment)) % alignment;
            offset += padding;
            dest.d_combined_bqsr_table = reinterpret_cast<combinedBqsrEntry*>(base + offset);
            sz = src.combined_bqsr_table_size * sizeof(src.h_combined_bqsr_table[0]);
            transfer_queue.memcpy(dest.d_combined_bqsr_table, src.h_combined_bqsr_table, sz);
            offset += sz;

            alignment = alignof(decltype(dest.d_static_quantized_scores[0]));
            padding = (alignment - (offset % alignment)) % alignment;
            offset += padding;
            dest.d_static_quantized_scores = reinterpret_cast<uint8_t*>(base + offset);
            sz = src.static_quantized_quality_scores_size * sizeof(src.h_static_quantized_quality_scores[0]);
            transfer_queue.memcpy(dest.d_static_quantized_scores, src.h_static_quantized_quality_scores, sz);
            offset += sz;

            alignment = alignof(decltype(dest.d_quantized_quality_scores[0]));
            padding = (alignment - (offset % alignment)) % alignment;
            offset += padding;
            dest.d_quantized_quality_scores = reinterpret_cast<uint8_t*>(base + offset);
            sz = src.quantized_quality_scores_size * sizeof(src.h_quantized_quality_scores[0]);
            transfer_queue.memcpy(dest.d_quantized_quality_scores, src.h_quantized_quality_scores, sz);
            offset += sz;

            transfer_queue.wait();

            dest.d_table_acc = src.h_table_acc;
            dest.combined_bqsr_table_size = src.combined_bqsr_table_size;
            dest.quantized_quality_scores_size = src.quantized_quality_scores_size;
            dest.static_quantized_scores_size = src.static_quantized_quality_scores_size;

            dest.offset_of_constant_allocation = offset;
            dest.is_first_copy = false;
        }

        offset = dest.offset_of_constant_allocation;
        if (n_elements > 0) {
            size_t chunk_sz = 0;

            // 1. d_sequences
            {
                dest.d_sequences = reinterpret_cast<char *>(base + offset);
                chunk_sz = (src.h_sequences_ptr[src_offset + n_elements - 1].start +
                                   src.h_sequences_ptr[src_offset + n_elements - 1].len) -
                                  (src.h_sequences_ptr[src_offset].start);
                size_t chunk_offset = src.h_sequences_ptr[src_offset].start;
                size_t sequences_data_size = sizeof(src.h_sequences[0]) * chunk_sz;
                transfer_queue.memcpy(dest.d_sequences, src.h_sequences + chunk_offset, sequences_data_size).wait();
                offset += sequences_data_size;
                dest.sequence_span_offset = chunk_offset;
            }

            dest.chunk_sz = chunk_sz;
            dest.chunk_offset = global_chunk_offset;
            global_chunk_offset += chunk_sz;
            // 2. Allocate d_sequences_ptr (an array of SpanField)
            //    Align the pointer properly.
            {
                size_t alignment = alignof(decltype(dest.d_sequences_ptr[0]));
                size_t padding = (alignment - (offset % alignment)) % alignment;
                offset += padding;
                dest.d_sequences_ptr = reinterpret_cast<SpanField *>(base + offset);
                size_t sequences_ptr_size = sizeof(dest.d_sequences_ptr[0]) * n_elements;
                transfer_queue.memcpy(dest.d_sequences_ptr, src.h_sequences_ptr + src_offset,
                                      sequences_ptr_size);
                offset += sequences_ptr_size;
            }

            // 5. Allocate d_qualities (an array of char)
            {
                size_t alignment = alignof(decltype(dest.d_qualities[0]));
                size_t padding = (alignment - (offset % alignment)) % alignment;
                offset += padding;
                dest.d_qualities = reinterpret_cast<char *>(base + offset);
                chunk_sz = (src.h_sequences_ptr[src_offset + n_elements - 1].start +
                            src.h_sequences_ptr[src_offset + n_elements - 1].len) -
                           (src.h_sequences_ptr[src_offset].start);
                size_t chunk_offset = src.h_sequences_ptr[src_offset].start;
                size_t qualities_data_size = sizeof(src.h_qualities[0]) * chunk_sz;
                transfer_queue.memcpy(dest.d_qualities, src.h_qualities + chunk_offset, qualities_data_size);
                offset += qualities_data_size;
            }

            // 10. Allocate d_flag (an array of int8_t)
            {
                size_t alignment = alignof(decltype(dest.d_flag[0]));
                size_t padding = (alignment - (offset % alignment)) % alignment;
                offset += padding;
                dest.d_flag = reinterpret_cast<uint16_t *>(base + offset);
                size_t flag_data_size = sizeof(dest.d_flag[0]) * n_elements;
                transfer_queue.memcpy(dest.d_flag, src.h_flags + src_offset, flag_data_size);
                offset += flag_data_size;
            }
            // 11. Allocate d_reads_group_covariate (an array of int8_t)
            {
                size_t alignment = alignof(decltype(dest.d_reads_group_covariate[0]));
                size_t padding = (alignment - (offset % alignment)) % alignment;
                offset += padding;
                dest.d_reads_group_covariate = reinterpret_cast<int32_t *>(base + offset);
                size_t reads_group_covariate_data_size = sizeof(dest.d_reads_group_covariate[0]) * n_elements;
                transfer_queue.memcpy(dest.d_reads_group_covariate, src.h_read_group + src_offset,
                                      reads_group_covariate_data_size);
                offset += reads_group_covariate_data_size;
            }

            // 16. Allocate d_stranded_seq (an array of char)
            {
                size_t alignment = alignof(decltype(dest.d_stranded_seq[0]));
                size_t padding = (alignment - (offset % alignment)) % alignment;
                offset += padding;
                dest.d_stranded_seq = reinterpret_cast<char *>(base + offset);
                chunk_sz = (src.h_sequences_ptr[src_offset + n_elements - 1].start +
                            src.h_sequences_ptr[src_offset + n_elements - 1].len) -
                           src.h_sequences_ptr[src_offset].start;
                size_t stranded_seq_size = sizeof(dest.d_stranded_seq[0]) * chunk_sz;
                size_t chunk_offset = src.h_sequences_ptr[src_offset].start;
                transfer_queue.memcpy(dest.d_stranded_seq, src.h_sequences + chunk_offset, stranded_seq_size);
                offset += stranded_seq_size;
            }

            // 17. Allocate d_stranded_seq_ptr (an array of SpanField)
            {
                size_t alignment = alignof(decltype(dest.d_stranded_seq_ptr[0]));
                size_t padding = (alignment - (offset % alignment)) % alignment;
                offset += padding;
                dest.d_stranded_seq_ptr = reinterpret_cast<SpanField*>(base + offset);
                size_t stranded_seq_ptr_size = sizeof(dest.d_stranded_seq_ptr[0]) * n_elements;
                transfer_queue.memcpy(dest.d_stranded_seq_ptr, src.h_sequences_ptr + src_offset,
                                      stranded_seq_ptr_size);
                offset += stranded_seq_ptr_size;
            }

            // 18. Allocate d_base_context (an array of int32_t)
            {
                size_t alignment = alignof(decltype(dest.d_base_context[0]));
                size_t padding = (alignment - (offset % alignment)) % alignment;
                offset += padding;
                dest.d_base_context = reinterpret_cast<int32_t *>(base + offset);
                chunk_sz = (src.h_sequences_ptr[src_offset + n_elements - 1].start +
                            src.h_sequences_ptr[src_offset + n_elements - 1].len) -
                           src.h_sequences_ptr[src_offset].start;
                size_t base_context_size = sizeof(dest.d_base_context[0]) * chunk_sz;
                offset += base_context_size;
            }

            // 19. Allocate d_base_context_ptr (an array of SpanField)
            {
                size_t alignment = alignof(decltype(dest.d_base_context_ptr[0]));
                size_t padding = (alignment - (offset % alignment)) % alignment;
                offset += padding;
                dest.d_base_context_ptr = reinterpret_cast<SpanField *>(base + offset);
                size_t base_context_ptr_size = sizeof(dest.d_base_context_ptr[0]) * n_elements;
                transfer_queue.memcpy(dest.d_base_context_ptr, src.h_sequences_ptr + src_offset, base_context_ptr_size);
                offset += base_context_ptr_size;
            }
            // 28. Allocate d_temporary_segments_offsets, capped at 2 BILLIONS
            {
                size_t alignment = alignof(decltype(dest.d_temporary_segments_sequences_offsets[0]));
                size_t padding = (alignment - (offset % alignment)) % alignment;
                offset += padding;
                dest.d_temporary_segments_sequences_offsets = reinterpret_cast<int32_t*>(base + offset);
                chunk_sz = n_elements;
                size_t segments_size = sizeof(dest.d_temporary_segments_sequences_offsets[0]) * chunk_sz;
                offset += segments_size;
            }
        }
        if (offset >= dest.max_memory) {
            exit(-1);
        }
        dest.max_num_elements = apply_bqsr::update_batch_num_elements(src, dest, src_offset + n_elements);
        transfer_queue.wait_and_throw();
    } catch (sycl::exception& e) {
        std::cout << "Error during memory copy from host to device: " << e.what() << std::endl;
    }
    Logger::getLogger().log(LOG_DEBUG, "    [Load batch] Time: "+
                                      std::to_string(t.elapsed()) + " sec");
}
//[[clang::optnone]]
void AppyBqsrKernelData::memcpy_dev_to_host(sycl::queue &transfer_queue, ApplyBqsrOutput &dest, size_t dest_offset,
                                                    AppyBqsrKernelData &src, size_t src_offset,
                                                    size_t n_elements) {
    Timer t;
    t.reset();

    transfer_queue.memcpy(dest.h_qual + src.chunk_offset, src.d_qualities, sizeof(dest.h_qual[0]) * src.chunk_sz);
    transfer_queue.wait();
    Logger::getLogger().log(LOG_DEBUG, "    [Moving out flags] Time: " + std::to_string(t.elapsed()) + " sec");
}



size_t apply_bqsr::update_batch_num_elements(const ApplyBqsrSequencesInput &in, const AppyBqsrKernelData &kernel_data,
                                 size_t start_idx) {
    size_t n_elements = 1;
    size_t idx = start_idx + n_elements;
    //size_t base = 0;

    if (idx > in._n) {
        return 0;
    }
    size_t offset = compute_chunk_sz(in, kernel_data, start_idx, n_elements);

    while (offset < kernel_data.max_memory && idx < in._n) {
        n_elements = std::min(n_elements * 2, in._n-start_idx);
        offset = compute_chunk_sz(in, kernel_data, start_idx, n_elements);
        idx = start_idx + n_elements;
    }

    // At this point, n_elements is an upper bound that may exceed memory.
    // Set up binary search between [low, high], where:
    //   low = n_elements/2 (last count that fit) and high = n_elements.
    size_t low = n_elements / 2;
    size_t high = n_elements;
    size_t best = low;  // best so far that fits.
    size_t usage = 0;

    while (low <= high) {
        size_t mid = low + (high - low) / 2;
        if (start_idx + mid > in._n) {  // avoid going out of bounds
            if (mid==0) continue;
            high = mid - 1;
            continue;
        }
        usage = compute_chunk_sz(in, kernel_data, start_idx, mid);
        if (usage < kernel_data.max_memory) {
            best = mid;    // mid fits within the memory
            low = mid + 1; // try for more elements
        } else {
            if (mid==0) continue;
            high = mid - 1;
        }
    }
    usage = compute_chunk_sz(in, kernel_data, start_idx, best);
    n_elements = best;
    return n_elements;
}

//[[clang::optnone]]
size_t compute_chunk_sz(const ApplyBqsrSequencesInput& in, const AppyBqsrKernelData& kernel_data, size_t start_idx, size_t n_elements) {

    size_t offset = 0;

    /*
     * Constant size
     */
    size_t table_sz = in.sz_qual_score_cov * sizeof(in.h_quality_table[0]);
    offset += table_sz;

    size_t alignment = alignof(decltype(in.h_cycle_table[0]));
    size_t padding = (alignment - (offset % alignment)) % alignment;
    offset += padding;
    table_sz = in.sz_cycle_cov * sizeof(in.h_cycle_table[0]);
    offset += table_sz;

    alignment = alignof(decltype(in.h_context_table[0]));
    padding = (alignment - (offset % alignment)) % alignment;
    offset += padding;
    table_sz = in.sz_context_cov * sizeof(in.h_context_table);
    offset += table_sz;

    alignment = alignof(decltype(in.h_combined_bqsr_table[0]));
    padding = (alignment - (offset % alignment)) % alignment;
    offset += padding;
    table_sz = in.combined_bqsr_table_size * sizeof(in.h_combined_bqsr_table[0]);
    offset += table_sz;

    alignment = alignof(decltype(in.h_static_quantized_quality_scores[0]));
    padding = (alignment - (offset % alignment)) % alignment;
    offset += padding;
    table_sz = in.static_quantized_quality_scores_size * sizeof(in.h_static_quantized_quality_scores[0]);
    offset += table_sz;

    alignment = alignof(decltype(in.h_quantized_quality_scores[0]));
    padding = (alignment - (offset % alignment)) % alignment;
    offset += padding;
    table_sz = in.quantized_quality_scores_size * sizeof(in.h_quantized_quality_scores[0]);
    offset += table_sz;


    /*
     * Variable sizes
     *
     */
    size_t chunk_sz = (in.h_sequences_ptr[start_idx + n_elements - 1].start +
                       in.h_sequences_ptr[start_idx + n_elements - 1].len) -
                      (in.h_sequences_ptr[start_idx].start);
    {
        alignment = alignof(decltype(kernel_data.d_sequences[0]));
        padding = (alignment - (offset % alignment)) % alignment;
        offset += padding;
        size_t sequences_data_size = sizeof(kernel_data.d_sequences[0]) * chunk_sz;
        offset += sequences_data_size;
    }

    {
        alignment = alignof(decltype(kernel_data.d_sequences_ptr[0]));
        padding = (alignment - (offset % alignment)) % alignment;
        offset += padding;
        size_t sequences_ptr_size = sizeof(kernel_data.d_sequences_ptr[0]) * n_elements;
        offset += sequences_ptr_size;
    }
    {
        alignment = alignof(decltype(kernel_data.d_qualities[0]));
        padding = (alignment - (offset % alignment)) % alignment;
        offset += padding;
        size_t qualities_data_size = sizeof(kernel_data.d_qualities[0]) * chunk_sz;
        offset += qualities_data_size;
    }
    {
        alignment = alignof(decltype(kernel_data.d_flag[0]));
        padding = (alignment - (offset % alignment)) % alignment;
        offset += padding;
        size_t flag_data_size = sizeof(kernel_data.d_flag[0]) * n_elements;
        offset += flag_data_size;
    }
    {
        alignment = alignof(decltype(kernel_data.d_reads_group_covariate[0]));
        padding = (alignment - (offset % alignment)) % alignment;
        offset += padding;
        size_t reads_group_covariate_data_size = sizeof(kernel_data.d_reads_group_covariate[0]) * n_elements;
        offset += reads_group_covariate_data_size;
    }
    {
        alignment = alignof(decltype(kernel_data.d_reads_group_covariate[0]));
        padding = (alignment - (offset % alignment)) % alignment;
        offset += padding;
        size_t reads_group_covariate_data_size = sizeof(kernel_data.d_reads_group_covariate[0]) * n_elements;
        offset += reads_group_covariate_data_size;
    }
    {
        alignment = alignof(decltype(kernel_data.d_stranded_seq[0]));
        padding = (alignment - (offset % alignment)) % alignment;
        offset += padding;
        size_t stranded_seq_size = sizeof(kernel_data.d_stranded_seq[0]) * chunk_sz;
        offset += stranded_seq_size;
    }
    {
        alignment = alignof(decltype(kernel_data.d_stranded_seq_ptr[0]));
        padding = (alignment - (offset % alignment)) % alignment;
        offset += padding;
        size_t stranded_seq_ptr_size = sizeof(kernel_data.d_stranded_seq_ptr[0]) * n_elements;
        offset += stranded_seq_ptr_size;
    }
    {
        alignment = alignof(decltype(kernel_data.d_base_context[0]));
        padding = (alignment - (offset % alignment)) % alignment;
        offset += padding;
        size_t base_context_size = sizeof(kernel_data.d_base_context[0]) * chunk_sz;
        offset += base_context_size;
    }
    {
        alignment = alignof(decltype(kernel_data.d_base_context_ptr[0]));
        padding = (alignment - (offset % alignment)) % alignment;
        offset += padding;
        size_t base_context_ptr_size = sizeof(kernel_data.d_base_context_ptr[0]) * n_elements;
        offset += base_context_ptr_size;
    }
    // 28. Allocate d_temporary_segments_offsets, capped at 2 BILLIONS
    {
        alignment = alignof(decltype(kernel_data.d_temporary_segments_sequences_offsets[0]));
        padding = (alignment - (offset % alignment)) % alignment;
        offset += padding;
        chunk_sz = n_elements;
        size_t segments_size = sizeof(kernel_data.d_temporary_segments_sequences_offsets[0]) * chunk_sz;
        offset += segments_size;
    }
    return offset;
}
