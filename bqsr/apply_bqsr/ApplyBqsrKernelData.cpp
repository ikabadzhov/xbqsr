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
    char*  base   = static_cast<char*>(dest.memory_buffer);
    size_t offset = 0;

    try {
        // ------------------------------------------------------------------
        // Constant region (copied once per run)
        // ------------------------------------------------------------------
        if (dest.is_first_copy) {
            dest.d_recalibrated_score_cache = src.d_recalibrated_score_cache;

            // Quality table
            dest.d_quality_table =
                alloc_and_copy<bqsrEntry>(transfer_queue, base, offset,
                                          src.h_quality_table, src.sz_qual_score_cov);

            // Cycle table
            dest.d_cycle_table =
                alloc_and_copy<bqsrEntry>(transfer_queue, base, offset,
                                          src.h_cycle_table, src.sz_cycle_cov);

            // Context table
            dest.d_context_table =
                alloc_and_copy<bqsrEntry>(transfer_queue, base, offset,
                                          src.h_context_table, src.sz_context_cov);

            // Combined BQSR table
            dest.d_combined_bqsr_table =
                alloc_and_copy<combinedBqsrEntry>(transfer_queue, base, offset,
                                                  src.h_combined_bqsr_table, src.combined_bqsr_table_size);

            // Static quantized scores (optional)
            if (src.static_quantized_quality_scores_size > 0) {
                dest.d_static_quantized_scores =
                    alloc_and_copy<uint8_t>(transfer_queue, base, offset,
                                            src.h_static_quantized_quality_scores,
                                            src.static_quantized_quality_scores_size);
            } else {
                // still advance offset by 0 so layout is consistent
                dest.d_static_quantized_scores =
                    alloc_simple<uint8_t>(base, offset, /*count=*/0);
            }

            // Quantized quality scores (optional)
            if (src.quantized_quality_scores_size > 0) {
                dest.d_quantized_quality_scores =
                    alloc_and_copy<uint8_t>(transfer_queue, base, offset,
                                            src.h_quantized_quality_scores,
                                            src.quantized_quality_scores_size);
            } else {
                dest.d_quantized_quality_scores =
                    alloc_simple<uint8_t>(base, offset, /*count=*/0);
            }

            transfer_queue.wait_and_throw();

            dest.d_table_acc                   = src.h_table_acc;
            dest.combined_bqsr_table_size      = src.combined_bqsr_table_size;
            dest.quantized_quality_scores_size = src.quantized_quality_scores_size;
            dest.static_quantized_scores_size  = src.static_quantized_quality_scores_size;

            dest.offset_of_constant_allocation = offset;
            dest.is_first_copy                 = false;
        }

        // ------------------------------------------------------------------
        // Per-batch region (changes every call)
        // ------------------------------------------------------------------
        offset = dest.offset_of_constant_allocation;

        if (n_elements > 0) {
            // Figure the contiguous span in the flat buffers
            const size_t first_start = src.h_sequences_ptr[src_offset].start;
            const size_t last_end    =
                src.h_sequences_ptr[src_offset + n_elements - 1].start +
                src.h_sequences_ptr[src_offset + n_elements - 1].len;
            const size_t chunk_sz    = last_end - first_start;

            // 1) sequences (bytes)
            dest.d_sequences = alloc_simple<char>(base, offset, chunk_sz);
            transfer_queue.memcpy(dest.d_sequences,
                                  src.h_sequences + first_start,
                                  sizeof(char) * chunk_sz);

            dest.sequence_span_offset = first_start;
            dest.chunk_sz             = chunk_sz;
            dest.chunk_offset         = global_chunk_offset;
            global_chunk_offset      += chunk_sz;

            // 2) sequences_ptr (SpanField[n_elements])
            dest.d_sequences_ptr =
                alloc_and_copy<SpanField>(transfer_queue, base, offset,
                                          src.h_sequences_ptr + src_offset, n_elements);

            // 3) qualities (bytes)
            dest.d_qualities = alloc_simple<char>(base, offset, chunk_sz);
            transfer_queue.memcpy(dest.d_qualities,
                                  src.h_qualities + first_start,
                                  sizeof(char) * chunk_sz);

            // 4) flags (uint16_t[n_elements])  NOTE: type must match kernel and host
            dest.d_flag =
                alloc_and_copy<uint16_t>(transfer_queue, base, offset,
                                         src.h_flags + src_offset, n_elements);

            // 5) read-group covariate (int32_t[n_elements])
            dest.d_reads_group_covariate =
                alloc_and_copy<int32_t>(transfer_queue, base, offset,
                                        src.h_read_group + src_offset, n_elements);

            // 6) stranded_seq (bytes) — currently copied from h_sequences
            dest.d_stranded_seq = alloc_simple<char>(base, offset, chunk_sz);
            transfer_queue.memcpy(dest.d_stranded_seq,
                                  src.h_sequences + first_start,
                                  sizeof(char) * chunk_sz);

            // 7) stranded_seq_ptr (SpanField[n_elements])
            dest.d_stranded_seq_ptr =
                alloc_and_copy<SpanField>(transfer_queue, base, offset,
                                          src.h_sequences_ptr + src_offset, n_elements);

            // 8) base_context (int32_t[chunk_sz]) — reserved; kernel fills later
            dest.d_base_context =
                alloc_simple<int32_t>(base, offset, chunk_sz);

            // 9) base_context_ptr (SpanField[n_elements])
            dest.d_base_context_ptr =
                alloc_and_copy<SpanField>(transfer_queue, base, offset,
                                          src.h_sequences_ptr + src_offset, n_elements);

            // 10) temporary segments offsets (int32_t[n_elements])
            dest.d_temporary_segments_sequences_offsets =
                alloc_simple<int32_t>(base, offset, n_elements);
        }

        if (offset > dest.max_memory) {
            // hard guard
            std::abort();
        }

        transfer_queue.wait_and_throw();
    }
    catch (sycl::exception& e) {
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

    transfer_queue.memcpy(dest.h_qual + dest_offset, src.d_qualities, sizeof(dest.h_qual[0]) * n_elements);
    transfer_queue.wait();
    Logger::getLogger().log(LOG_DEBUG, "    [Moving out flags] Time: " + std::to_string(t.elapsed()) + " sec");
}



size_t apply_bqsr::update_batch_num_elements(const ApplyBqsrSequencesInput &in, const AppyBqsrKernelData &kd,
                                 size_t start_idx) {
    if (start_idx >= in._n) return 0;

    const size_t Nmax = in._n - start_idx;

    // 1) Probe: can 1 element fit at all?
    size_t usage = compute_chunk_sz(in, kd, start_idx, 1);
    if (usage > kd.max_memory) {
        return 0;  // caller MUST handle n==0 (no indexing with n-1!)
    }

    // 2) Exponential search for an upper bound that DOESN'T fit
    size_t best = 1;   // last known good (fits)
    size_t high = 1;
    while (true) {
        size_t next = std::min(high * 2, Nmax);
        usage = compute_chunk_sz(in, kd, start_idx, next);
        if (usage <= kd.max_memory) {
            best = next;
            if (next == Nmax) { high = next; break; } // exhausted
            high = next;
        } else {
            // first non-fitting count is our upper bound
            high = next;
            break;
        }
    }

    // 3) Binary search in [best, high] for the max that fits.
    // Invariant: best fits; high may or may not fit.
    size_t low = best;
    while (low < high) {
        size_t mid = low + (high - low + 1) / 2; // upper mid to avoid infinite loop
        usage = compute_chunk_sz(in, kd, start_idx, mid);
        if (usage <= kd.max_memory) {
            low = mid;   // mid fits → move up
        } else {
            high = mid - 1; // mid doesn't fit → move down
        }
    }
    return low;  // >= 1 here
}

//[[clang::optnone]]
size_t compute_chunk_sz(const ApplyBqsrSequencesInput& in,
                        const AppyBqsrKernelData&     kd,
                        size_t start_idx,
                        size_t n_elements)
{
    size_t offset = 0;

    // -----------------------------
    // Constant region (copied once)
    // -----------------------------
    simulate_alloc<bqsrEntry>(offset, in.sz_qual_score_cov);                 // quality table
    simulate_alloc<bqsrEntry>(offset, in.sz_cycle_cov);                      // cycle table
    simulate_alloc<bqsrEntry>(offset, in.sz_context_cov);                    // context table
    simulate_alloc<combinedBqsrEntry>(offset, in.combined_bqsr_table_size);  // combined table
    simulate_alloc<uint8_t>(offset, in.static_quantized_quality_scores_size);// static quant map (may be 0)
    simulate_alloc<uint8_t>(offset, in.quantized_quality_scores_size);       // dynamic quant map (may be 0)

    // -----------------------------
    // Variable region (per-batch)
    // -----------------------------
    // Flat span size across the selected reads
    const size_t first_start = in.h_sequences_ptr[start_idx].start;
    const size_t last_end    =
        in.h_sequences_ptr[start_idx + n_elements - 1].start +
        in.h_sequences_ptr[start_idx + n_elements - 1].len;
    const size_t chunk_sz = last_end - first_start;

    // Sequences bytes
    simulate_alloc<std::remove_pointer_t<decltype(kd.d_sequences)>>(offset, chunk_sz); // char

    // Sequence pointers (SpanField per element)
    simulate_alloc<std::remove_pointer_t<decltype(kd.d_sequences_ptr)>>(offset, n_elements); // SpanField

    // Qualities bytes
    simulate_alloc<std::remove_pointer_t<decltype(kd.d_qualities)>>(offset, chunk_sz); // char

    // Flags
    simulate_alloc<std::remove_pointer_t<decltype(kd.d_flag)>>(offset, n_elements); // uint16_t

    // Read-group covariate  (appeared twice in your original — keep only once)
    simulate_alloc<std::remove_pointer_t<decltype(kd.d_reads_group_covariate)>>(offset, n_elements); // int32_t

    // Stranded seq bytes
    simulate_alloc<std::remove_pointer_t<decltype(kd.d_stranded_seq)>>(offset, chunk_sz); // char

    // Stranded seq ptrs
    simulate_alloc<std::remove_pointer_t<decltype(kd.d_stranded_seq_ptr)>>(offset, n_elements); // SpanField

    // Base context (kernel writes later)
    simulate_alloc<std::remove_pointer_t<decltype(kd.d_base_context)>>(offset, chunk_sz); // int32_t

    // Base context ptrs
    simulate_alloc<std::remove_pointer_t<decltype(kd.d_base_context_ptr)>>(offset, n_elements); // SpanField

    // Temporary per-element offsets
    simulate_alloc<std::remove_pointer_t<decltype(kd.d_temporary_segments_sequences_offsets)>>(offset, n_elements); // int32_t

    return offset;
}