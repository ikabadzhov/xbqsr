#ifndef APPLYBQSRKERNELDATA_H
#define APPLYBQSRKERNELDATA_H
#include "../BaseRecalibratorTablesAccessor.h"
#include "../TablesBQSR.h"
#include "../../device/batch.h"
#include "../../data_types/VectorCollection.h"

struct combinedBqsrEntry;
struct bqsrEntry;


class ApplyBqsrSequencesInput : public InputData {
public:
    int32_t* h_read_group;
    char* h_sequences;
    SpanField* h_sequences_ptr;
    uint8_t* h_qualities;
    uint16_t* h_flags;

    int32_t sz_qual_score_cov;
    int32_t sz_cycle_cov;
    int32_t sz_context_cov;

    BaseRecalibratorTablesAccessor h_table_acc;

    bqsrEntry* h_quality_table;
    bqsrEntry* h_cycle_table;
    bqsrEntry* h_context_table;

    combinedBqsrEntry* h_combined_bqsr_table;
    size_t combined_bqsr_table_size;

    uint8_t* h_static_quantized_quality_scores;
    size_t static_quantized_quality_scores_size;

    uint8_t* h_quantized_quality_scores;
    size_t quantized_quality_scores_size;

    int* d_recalibrated_score_cache;

    ApplyBqsrSequencesInput(bqsrEntry* quality_table, bqsrEntry* cycle_table, bqsrEntry* context_table, size_t quality_table_sz, size_t cycle_table_sz, size_t context_table_sz, BaseRecalibratorTablesAccessor& table_acc,
                            char* sequences, SpanField* sequence_ptrs, uint8_t* qualities, uint16_t* flags, int32_t* read_group, combinedBqsrEntry* combined_bqsr_table, size_t combined_bqsr_table_size,
                            uint8_t* quantized_quality_scores, size_t quantized_quality_scores_size, uint8_t* static_quantized_quality_scores, size_t static_quantized_quality_scores_size, int* recalibrated_score_cache, size_t n_elements);
};

class ApplyBqsrOutput : public OutputData {
public:
   uint8_t* h_qual;
    ApplyBqsrOutput(uint8_t* qual) : h_qual(qual) {};
};


class AppyBqsrKernelConfiguration {
public:
    size_t max_memory = 0;
    size_t max_n_elements_first_batch = 0;
};

class AppyBqsrKernelData {
public:
    sycl::queue _q;
    memory_space_t _space = memory_space_device;
    size_t max_num_elements {};
    bool is_first_copy = true;
    size_t chunk_offset = 0;
    size_t offset_of_constant_allocation = 0;

    size_t sequence_span_offset = 0;

    char* memory_buffer = nullptr;
    size_t max_memory = 0;

    char* d_sequences = nullptr;
    SpanField* d_sequences_ptr = nullptr;

    char* d_qualities = nullptr;

    uint16_t* d_flag = nullptr;
    int32_t* d_reads_group_covariate = nullptr;

    char* d_stranded_seq = nullptr;
    SpanField* d_stranded_seq_ptr = nullptr;
    int32_t* d_base_context = nullptr;
    SpanField* d_base_context_ptr = nullptr;

    BaseRecalibratorTablesAccessor d_table_acc;

    size_t chunk_sz = 0;
    int32_t prev_batch_max_index = 0;

    bqsrEntry* d_quality_table = nullptr;
    bqsrEntry* d_cycle_table = nullptr;
    bqsrEntry* d_context_table = nullptr;

    combinedBqsrEntry* d_combined_bqsr_table = nullptr;
    size_t combined_bqsr_table_size;

    uint8_t* d_quantized_quality_scores = nullptr;
    size_t quantized_quality_scores_size;

    uint8_t* d_static_quantized_scores = nullptr;
    int32_t static_quantized_scores_size = 0;

    int32_t* d_temporary_segments_sequences_offsets = nullptr;

    int32_t* d_recalibrated_score_cache = nullptr;

public:
    AppyBqsrKernelData() = default;

    AppyBqsrKernelData(sycl::queue &q, AppyBqsrKernelConfiguration &config,
                                                           memory_space_t space = memory_space_device) : _q(q), _space(space) {
        //assert(config.max_memory > (config.interval_memory + config.reference_memory));

        size_t max_memory = config.max_memory;// - config.interval_memory - config.reference_memory;

        if (space != memory_space_t::none) {

            memory_buffer = sycl::malloc_device<char>(max_memory, q);

            if (!memory_buffer) {
                std::cerr << "Error: Memory allocation failed for one or more pointers." << std::endl;
                exit(-1);
            }
        }
        this->max_memory = max_memory;
        this->max_num_elements = config.max_n_elements_first_batch;

    }

    void swap(AppyBqsrKernelData& rhs);

    // Destructor to manage dynamic memory if ownership is taken
    ~ AppyBqsrKernelData();

    // Disable copy semantics to prevent accidental copying
    AppyBqsrKernelData(const AppyBqsrKernelData&) = delete;
    AppyBqsrKernelData& operator=(const AppyBqsrKernelData&) = delete;

    // Enable move semantics
    AppyBqsrKernelData(AppyBqsrKernelData&& rhs) noexcept;

    AppyBqsrKernelData& operator=(AppyBqsrKernelData&& rhs) noexcept;

    /**
     * @brief Generic memcpy function to copy specified fields from src to dest.
     *
     * @tparam Members Variadic template parameters representing member pointers.
     * @param dest Destination DeviceData object.
     * @param src Source DeviceData object.
     * @param src_from Starting index in the source.
     * @param dest_to Starting index in the destination.
     * @param n Number of elements to copy.
     * @param members Member pointers specifying which fields to copy.
     */
    static void memcpy_host_to_dev(sycl::queue& transfer_queue, AppyBqsrKernelData& dest, size_t dest_offset, const ApplyBqsrSequencesInput& src,
                                   size_t src_offset, size_t n_elements);

    static void memcpy_dev_to_host(sycl::queue& transfer_queue, ApplyBqsrOutput& dest, size_t dest_offset, AppyBqsrKernelData& src,
                                   size_t src_offset, size_t n_elements);

    template <typename InputT>
    static size_t compute_max_elements(const InputT& in,
                                      const AppyBqsrKernelData& kd,
                                      size_t start_idx);
};

namespace apply_bqsr {
    size_t update_batch_num_elements(const ApplyBqsrSequencesInput &in, const AppyBqsrKernelData &kernel_data,
                                     size_t start_idx);
}
template<typename InputT>
size_t AppyBqsrKernelData::compute_max_elements(const InputT &in, const AppyBqsrKernelData &kd, size_t start_idx) {
    return apply_bqsr::update_batch_num_elements(in, kd, start_idx);
}



#endif //APPLYBQSRKERNELDATA_H
