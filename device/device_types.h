#ifndef TESTONEAPI_DEVICE_TYPES_H
#define TESTONEAPI_DEVICE_TYPES_H

#include <cstddef>
#include <cstdint>
#include <vector>
#include <sycl/queue.hpp>
#include <sycl/sycl.hpp>
#include <sycl/stream.hpp>
#include <math.h>
#include "../sycl_parallel_library/memory/context.h"

const size_t MAX_SUPERBATCH_SZ = 2*static_cast<size_t>(std::pow(2,30));

class InputData {
public:
    int32_t _n_loaded;
    ssize_t _n_to_process;
    ssize_t _n;
    explicit InputData(ssize_t n) : _n_loaded(0), _n_to_process(n), _n(n) {}
    ~InputData() = default;
};

class OutputData {
public:
    ssize_t _n_written = 0;
    OutputData() = default;
    ~OutputData() = default;
};

template <typename kernel_data_t>
concept HasQueueMember = requires(kernel_data_t t) {
    { t._q } -> std::same_as<sycl::queue&>;
};

template <typename kernel_data_t>
concept HasMaxNumElelementsMember = requires(kernel_data_t t) {
    { t.max_num_elements } -> std::same_as<size_t &>;
};


// Concept to check for a default constructor
template <typename kernel_data_t>
concept HasDefaultConstructor = std::default_initializable<kernel_data_t>;

// Concept to check for a specific constructor with parameters
template <typename kernel_data_t, typename kernel_data_config_t>
concept HasSpecificConstructor = requires(sycl::queue& q, kernel_data_config_t& config, memory_space_t mem_space) {
    kernel_data_t(q, config, mem_space);
};

// Concept to check for a swap method
template <typename kernel_data_t>
concept HasSwapMethod = requires(kernel_data_t t, kernel_data_t& other) {
    { t.swap(other) } -> std::same_as<void>;
};

// Concept to check for a non-default destructor
template <typename kernel_data_t>
concept HasNonDefaultDestructor = !std::is_trivially_destructible_v<kernel_data_t>;

// Concept to check that copy constructor and copy assignment are deleted
template <typename kernel_data_t>
concept DeletedCopySemantics = !std::copy_constructible<kernel_data_t> && !std::assignable_from<kernel_data_t, kernel_data_t>;

// Concept to check for move semantics
template <typename kernel_data_t>
concept HasMoveSemantics = std::movable<kernel_data_t> && std::assignable_from<kernel_data_t&, kernel_data_t>;

// Concept to check for the static method memcpy_host_to_dev
template <typename kernel_data_t, typename input_data_t>
concept HasMemcpyHostToDev = requires(sycl::queue& transfer_queue, kernel_data_t& dest, size_t dest_offset,
const input_data_t& src, size_t src_offset, size_t n_elements) {
{ kernel_data_t::memcpy_host_to_dev(transfer_queue, dest, dest_offset, src, src_offset, n_elements) } -> std::same_as<void>;
};

// Concept to check for the static method memcpy_dev_to_host
template <typename kernel_data_t, typename output_data_t>
concept HasMemcpyDevToHost = requires(sycl::queue& transfer_queue, output_data_t& dest, size_t dest_offset,
        kernel_data_t& src, size_t src_offset, size_t n_elements) {
{ kernel_data_t::memcpy_dev_to_host(transfer_queue, dest, dest_offset, src, src_offset, n_elements) } -> std::same_as<void>;
};

template <typename kernel_data_t, typename batch_config_t, typename input_data_t, typename output_data_t>
concept DeviceDataInterface =
    HasQueueMember<kernel_data_t> &&
    HasDefaultConstructor<kernel_data_t> &&
    HasSpecificConstructor<kernel_data_t, batch_config_t> &&
    HasSwapMethod<kernel_data_t> &&
    HasNonDefaultDestructor<kernel_data_t> &&
    DeletedCopySemantics<kernel_data_t> &&
    HasMoveSemantics<kernel_data_t> &&
    HasMemcpyHostToDev<kernel_data_t, input_data_t> &&
    HasMemcpyDevToHost<kernel_data_t, output_data_t>;

template<typename input_data_t>
concept DerivedFromInputDataT = std::is_base_of_v<InputData, input_data_t>;

template<typename output_data_t>
concept DerivedFromOutputDataT = std::is_base_of_v<OutputData, output_data_t>;

template <typename kernel_data_t, typename batch_config_t>
class Batch {
public:
    sycl::queue _q_transfer;
    sycl::queue _q_process;
    kernel_data_t _d_data;

    //ssize_t batch_offset = 0;
    size_t _max_num_elements{};
    size_t _n;
    size_t _n_loaded;
    ssize_t _n_processed{};


    explicit Batch(sycl::queue &q_transfer, sycl::queue &q_process, batch_config_t& config) : _q_transfer(q_transfer),
                                                                                 _q_process(q_process), _d_data(q_transfer, config, memory_space_t::memory_space_device), _n(0), _n_loaded(0) {
        _max_num_elements = _d_data.max_num_elements;
    }

    ~Batch() = default;

    template<typename InputT>
    requires DerivedFromInputDataT<InputT>
    void memcpy_host_to_dev(InputT& input) {
        size_t n_reads_loading = input._n - input._n_loaded;
        if (n_reads_loading <= 0){
            this->_n = 0;
            return;
        }

        if (n_reads_loading > _max_num_elements) {
            n_reads_loading = _max_num_elements;
        }
        this->_n = n_reads_loading;
        size_t first_read_id = input._n_loaded;
        // TODO: Maybe we should change the method so that the destination offset is 0, as we always fill from 0.
        kernel_data_t::memcpy_host_to_dev(_q_transfer, this->_d_data, 0, input, first_read_id, n_reads_loading);
        input._n_to_process -= n_reads_loading;
        input._n_loaded += n_reads_loading;
        // TODO: Check that this doesn't break the compatibility in the other algorithms, e.g. mark duplicates
        _max_num_elements = _d_data.max_num_elements;
    }

    template<typename OutputT>
    void memcpy_dev_to_host(OutputT& output) {
        size_t to_store = this->_n;
        if (to_store <= 0){
            return;
        }

        kernel_data_t::memcpy_dev_to_host(_q_transfer, output, output._n_written, this->_d_data, 0, to_store);
        output._n_written += this->_n;
    }

    void transfer_wait() {
        _q_transfer.wait();
    }

    void process_wait() {
        _q_process.wait();
    }
};

#endif //TESTONEAPI_DEVICE_TYPES_H
