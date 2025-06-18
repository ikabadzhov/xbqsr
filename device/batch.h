#ifndef TESTONEAPI_SUPERBATCH_H
#define TESTONEAPI_SUPERBATCH_H

#include <future>
#include "device_types.h"
#include "../timer/Timer.h"
#include "../logger/Logger.h"



template <typename kernel_data_t, typename batch_config_t, typename input_data_t>
void load_batch(input_data_t& input, std::shared_ptr<Batch<kernel_data_t, batch_config_t>> batch_transfer){
    Timer t;
    t.reset();
    batch_transfer->memcpy_host_to_dev(input);
    batch_transfer->transfer_wait();
   // std::cout << "Time to load batch : " << t.elapsed_us() / 1000000.0 << "sec" << std::endl;
}


template <typename kernel_data_t, typename batch_config_t, typename output_data_t>
void store_batch(output_data_t& output, std::shared_ptr<Batch<kernel_data_t, batch_config_t>> batch_transfer){
    Timer t;
    t.reset();
    batch_transfer->memcpy_dev_to_host(output);
    batch_transfer->transfer_wait();
    //std::cout << "Time to write out tables: " << t.elapsed_us() /1000000.0 << " sec" << std::endl;
}


template <typename kernel_data_t, typename batch_config_t, typename input_data_t, typename output_data_t, typename proc_func_t>
requires DerivedFromInputDataT<input_data_t> && DerivedFromOutputDataT<output_data_t> && DeviceDataInterface<kernel_data_t, batch_config_t, input_data_t, output_data_t>
void  batch_processing(sycl::device& device, input_data_t& input, output_data_t& output, proc_func_t processing_f, batch_config_t& batch_config) {


    sycl::queue q_processing(device, sycl::async_handler{},
                             sycl::property_list{sycl::property::queue::enable_profiling()});

    sycl::queue q_transfer(q_processing.get_context(), device, sycl::async_handler{},
                           sycl::property_list{sycl::property::queue::enable_profiling()});

    std::shared_ptr<Batch<kernel_data_t, batch_config_t>> batch_process = std::make_shared<Batch<kernel_data_t, batch_config_t>>(q_transfer, q_processing, batch_config);
    std::shared_ptr<Batch<kernel_data_t, batch_config_t>> batch_transfer = std::make_shared<Batch<kernel_data_t, batch_config_t>>(q_transfer, q_processing, batch_config);

    bool first_batch = true;
    int32_t ith_batch = 0;
    // process until there is no more reads to process
    // except first batch where there is not yet data
    while (first_batch || batch_process->_n != 0)
    {
        std::future<void> transferAsync;
        transferAsync = std::async(std::launch::async, store_batch<kernel_data_t, batch_config_t, output_data_t>, std::ref(output),
                                   batch_transfer);

        //store_batch<kernel_data_t, batch_config_t, output_data_t>(output, batch_transfer);

        auto processAsync = std::async(std::launch::async, processing_f, batch_process); // does nothing if process_data->n_seqs is 0
        //processing_f(batch_process);

        transferAsync.wait();

        auto inputAsync = std::async(std::launch::async, load_batch<kernel_data_t, batch_config_t, input_data_t>, std::ref(input),
                                     batch_transfer);
        //load_batch<kernel_data_t, batch_config_t, input_data_t>(input, batch_transfer);

        inputAsync.wait();
        processAsync.wait();
        Logger::getLogger().log(LOG_DEBUG, "    Processed batch ", ith_batch);
        ith_batch++;
        batch_process->_n_processed = batch_process->_n;

        std::swap(batch_process, batch_transfer);

        first_batch = false;
    }
    //store_batch<kernel_data_t, batch_config_t, output_data_t>(output, batch_transfer);

    auto transferAsync = std::async(std::launch::async, store_batch<kernel_data_t, batch_config_t, output_data_t>, std::ref(output),
                                    batch_transfer);
    transferAsync.wait();
}


template <typename kernel_data_t, typename batch_config_t, typename input_data_t, typename output_data_t, typename proc_func_t>
requires DerivedFromInputDataT<input_data_t> && DerivedFromOutputDataT<output_data_t> && DeviceDataInterface<kernel_data_t, batch_config_t, input_data_t, output_data_t>
void batch_processing_variable_length(sycl::device& device, input_data_t& input, output_data_t& output, proc_func_t processing_f, batch_config_t& batch_config) {


    sycl::queue q_processing(device, sycl::async_handler{},
                             sycl::property_list{sycl::property::queue::enable_profiling()});

    sycl::queue q_transfer(q_processing.get_context(), device, sycl::async_handler{},
                           sycl::property_list{sycl::property::queue::enable_profiling()});

    std::shared_ptr<Batch<kernel_data_t, batch_config_t>> batch_process = std::make_shared<Batch<kernel_data_t, batch_config_t>>(q_transfer, q_processing, batch_config);
    std::shared_ptr<Batch<kernel_data_t, batch_config_t>> batch_transfer = std::make_shared<Batch<kernel_data_t, batch_config_t>>(q_transfer, q_processing, batch_config);

    bool first_batch = true;
    int32_t ith_batch = 0;

    // process until there is no more reads to process
    // except first batch where there is not yet data
    while (first_batch || batch_process->_n != 0)
    {
        std::future<void> transferAsync;
        transferAsync = std::async(std::launch::async, store_batch<kernel_data_t, batch_config_t, output_data_t>, std::ref(output),
                                   batch_transfer);
//      store_batch<kernel_data_t, batch_config_t, output_data_t>(output, batch_transfer);

        auto processAsync = std::async(std::launch::async, processing_f, batch_process); // does nothing if process_data->n_seqs is 0
//      processing_f(batch_process);

        transferAsync.wait();

        auto inputAsync = std::async(std::launch::async, load_batch<kernel_data_t, batch_config_t, input_data_t>, std::ref(input),
                                     batch_transfer);
//      load_batch<kernel_data_t, batch_config_t, input_data_t>(input, batch_transfer);

        inputAsync.wait();
        processAsync.wait();

        Logger::getLogger().log(LOG_DEBUG, "    Processed batch ", ith_batch);
        ith_batch++;

        batch_process->_n_processed = batch_process->_n;

        auto tmp_max_sz = batch_transfer->_max_num_elements;
        std::swap(batch_process, batch_transfer);
        batch_transfer->_max_num_elements = tmp_max_sz;

        first_batch = false;
    }

    auto transferAsync = std::async(std::launch::async, store_batch<kernel_data_t, batch_config_t, output_data_t>, std::ref(output),
                                    batch_transfer);
    transferAsync.wait();
    //store_batch<kernel_data_t, batch_config_t, output_data_t>(output, batch_transfer);

}
#endif //TESTONEAPI_SUPERBATCH_H
