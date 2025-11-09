#ifndef TESTONEAPI_SUPERBATCH_H
#define TESTONEAPI_SUPERBATCH_H

#include <future>
#include <mutex>
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
void batch_processing_variable_length(sycl::queue &queue, input_data_t& input, output_data_t& output, proc_func_t processing_f, batch_config_t& batch_config) {


    sycl::queue q_processing(queue.get_context(), queue.get_device(), sycl::async_handler{},
                             sycl::property_list{sycl::property::queue::enable_profiling()});

    sycl::queue q_transfer(queue.get_context(), queue.get_device(), sycl::async_handler{},
                           sycl::property_list{sycl::property::queue::enable_profiling()});

    std::shared_ptr<Batch<kernel_data_t, batch_config_t>> batch_process = std::make_shared<Batch<kernel_data_t, batch_config_t>>(q_transfer, q_processing, batch_config);
    std::shared_ptr<Batch<kernel_data_t, batch_config_t>> batch_transfer = std::make_shared<Batch<kernel_data_t, batch_config_t>>(q_transfer, q_processing, batch_config);


    // --- loader sync ---
    std::mutex m_load;
    std::condition_variable cv_load;
    bool load_request = false;   // main -> loader
    bool load_done    = false;   // loader -> main
    std::exception_ptr loader_ex;

    // preload first batch synchronously so the loop can start processing immediately
    load_batch(input, batch_process);

    // --- loader thread ---
    // waits on 'load_request', fills batch_transfer, posts 'load_done'
    std::jthread loader([&](std::stop_token st) {
        try {
            for (;;) {
                std::unique_lock lk(m_load);
                cv_load.wait(lk, [&]{ return load_request || st.stop_requested(); });
                if (st.stop_requested()) break;
                load_request = false;  // consume request
                lk.unlock();

                load_batch(input, batch_transfer);

                lk.lock();
                load_done = true;
                lk.unlock();
                cv_load.notify_one();
            }
        } catch (...) {
            loader_ex = std::current_exception();
            std::lock_guard lk(m_load);
            load_done = true;
            cv_load.notify_one();
        }
    });

    // --------------- main pipeline ---------------
    bool   first_batch = true;
    //int32_t ith_batch  = 0;

    while (first_batch || batch_process->_n != 0) {
        // kick the loader to fetch the NEXT batch into batch_transfer
        {
            std::lock_guard lk(m_load);
            load_request = true;
            load_done    = false;
        }
        cv_load.notify_one();

        // process the CURRENT batch while the loader is working
        processing_f(batch_process);

        // (Optional) store the just-processed results here, synchronously.
        // This uses batch_process only; no races.
        store_batch(output, batch_process);

        // wait for loader to finish; propagate error if any
        {
            std::unique_lock lk(m_load);
            cv_load.wait(lk, [&]{ return load_done || loader_ex != nullptr; });
        }
        if (loader_ex) std::rethrow_exception(loader_ex);

        // swap roles: loaded buffer becomes the next to process
        batch_process->_n_processed = batch_process->_n;

        //auto tmp_max = batch_transfer->_max_num_elements; // if you need to carry this across
        std::swap(batch_process, batch_transfer);
        //batch_transfer->_max_num_elements = tmp_max;

        first_batch = false;
        //++ith_batch;
    }

    // graceful shutdown of loader
    loader.request_stop();   // signal
    {
        std::lock_guard lk(m_load);
        load_request = true; // in case it's waiting; triggers loop to see stop_requested()
    }
    cv_load.notify_one();
    loader.join();

}
#endif //TESTONEAPI_SUPERBATCH_H
