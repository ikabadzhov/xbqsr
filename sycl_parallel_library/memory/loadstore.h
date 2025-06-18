#ifndef TESTONEAPI_LOADSTORE_H
#define TESTONEAPI_LOADSTORE_H

#include "../device_utils/types.h"
#include <sycl/sycl.hpp>

////////////////////////////////////////////////////////////////////////////////
// reg<->shared

template<int nt, int vt, typename type_t, int shared_size>
void reg_to_shared_thread(array_t<type_t, vt> x, int tid,
                                      type_t (&shared)[shared_size], const sycl::nd_item<1> &item, bool sync = true) {

    static_assert(shared_size >= nt * vt,
                  "reg_to_shared_thread must have at least nt * vt storage");

    thread_iterate<vt>([&](int i, int j) {
        shared[j] = x[i];
    }, tid);
    /*
    DPCT1118:11: SYCL group functions and algorithms must be encountered in
    converged control flow. You may need to adjust the code.
    */
    /*
    DPCT1065:48: Consider replacing sycl::nd_item::barrier() with
    sycl::nd_item::barrier(sycl::access::fence_space::local_space) for better
    performance if there is no access to global memory.
    */
    if (sync) item.barrier(sycl::access::fence_space::local_space);
}

template<int nt, int vt, typename type_t, int shared_size>
array_t<type_t, vt> shared_to_reg_thread(
        const type_t (&shared)[shared_size], int tid,
        const sycl::nd_item<1> &item, bool sync = true) {

    static_assert(shared_size >= nt * vt,
                  "reg_to_shared_thread must have at least nt * vt storage");

    array_t<type_t, vt> x;
    thread_iterate<vt>([&](int i, int j) {
        x[i] = shared[j];
    }, tid);
    /*
    DPCT1118:12: SYCL group functions and algorithms must be encountered in
    converged control flow. You may need to adjust the code.
    */
    /*
    DPCT1065:49: Consider replacing sycl::nd_item::barrier() with
    sycl::nd_item::barrier(sycl::access::fence_space::local_space) for better
    performance if there is no access to global memory.
    */
    if (sync) item.barrier(sycl::access::fence_space::local_space);
    return x;
}

template<int nt, int vt, typename type_t, int shared_size>
void reg_to_shared_strided(array_t<type_t, vt> x, int tid,
                                       type_t (&shared)[shared_size], const sycl::nd_item<1> &item, bool sync = true) {

    static_assert(shared_size >= nt * vt,
                  "reg_to_shared_strided must have at least nt * vt storage");

    strided_iterate<nt, vt>([&](int i, int j) { shared[j] = x[i]; }, tid);
    /*
    DPCT1118:13: SYCL group functions and algorithms must be encountered in
    converged control flow. You may need to adjust the code.
    */
    /*
    DPCT1065:50: Consider replacing sycl::nd_item::barrier() with
    sycl::nd_item::barrier(sycl::access::fence_space::local_space) for better
    performance if there is no access to global memory.
    */
    if (sync) item.barrier(sycl::access::fence_space::local_space);
}

template<int nt, int vt, typename type_t, int shared_size>
array_t<type_t, vt> shared_to_reg_strided(
        const type_t (&shared)[shared_size], int tid,
        const sycl::nd_item<1> &item, bool sync = true) {

    static_assert(shared_size >= nt * vt,
                  "shared_to_reg_strided must have at least nt * vt storage");

    array_t<type_t, vt> x;
    strided_iterate<nt, vt>([&](int i, int j) { x[i] = shared[j]; }, tid);
    /*
    DPCT1118:14: SYCL group functions and algorithms must be encountered in
    converged control flow. You may need to adjust the code.
    */
    /*
    DPCT1065:51: Consider replacing sycl::nd_item::barrier() with
    sycl::nd_item::barrier(sycl::access::fence_space::local_space) for better
    performance if there is no access to global memory.
    */
    if (sync) item.barrier(sycl::access::fence_space::local_space);
    return x;
}

template<int nt, int vt, typename type_t, int shared_size>
array_t<type_t, vt> shared_gather(const type_t(&data)[shared_size],
                                              array_t<int, vt> indices, const sycl::nd_item<1> &item, bool sync = true) {

    static_assert(shared_size >= nt * vt,
                  "shared_gather must have at least nt * vt storage");

    array_t<type_t, vt> x;
    iterate<vt>([&](int i) { x[i] = data[indices[i]]; });
    /*
    DPCT1118:15: SYCL group functions and algorithms must be encountered in
    converged control flow. You may need to adjust the code.
    */
    /*
    DPCT1065:52: Consider replacing sycl::nd_item::barrier() with
    sycl::nd_item::barrier(sycl::access::fence_space::local_space) for better
    performance if there is no access to global memory.
    */
    if (sync) item.barrier();
    return x;
}

template<int nt, int vt, typename type_t, int shared_size>
array_t<type_t, vt> thread_to_strided(array_t<type_t, vt> x, int tid, type_t (&shared)[shared_size], const sycl::nd_item<3> &item_ct1) {

    reg_to_shared_thread<nt, vt>(x, tid, shared, item_ct1);
    return shared_to_reg_strided<nt, vt>(shared, tid, item_ct1);
}



////////////////////////////////////////////////////////////////////////////////
// reg<->memory

template<int nt, int vt, int vt0 = vt, typename type_t, typename it_t>
void reg_to_mem_strided(array_t<type_t, vt> x, int tid,
                                    int count, it_t mem) {

    strided_iterate<nt, vt, vt0>([=](int i, int j) {
        mem[j] = x[i];
    }, tid, count);
}

template<int nt, int vt, int vt0 = vt, typename it_t>
array_t<typename std::iterator_traits<it_t>::value_type, vt> mem_to_reg_strided(it_t mem, int tid, int count) {
    typedef typename std::iterator_traits<it_t>::value_type type_t;
    array_t<type_t, vt> x;
    strided_iterate<nt, vt, vt0>([&](int i, int j) {
        x[i] = mem[j];
        }, tid, count);
    return x;
}

template<int nt, int vt, int vt0 = vt, typename type_t, typename it_t, int shared_size>
void reg_to_mem_thread(array_t<type_t, vt> x, int tid, int count, it_t mem, type_t (&shared)[shared_size], const sycl::nd_item<1> &item) {
    reg_to_shared_thread<nt>(x, tid, shared, item);
    array_t<type_t, vt> y = shared_to_reg_strided<nt, vt>(shared, tid, item);
    reg_to_mem_strided<nt, vt, vt0>(y, tid, count, mem);
}

template<int nt, int vt, int vt0 = vt, typename type_t, typename it_t, int shared_size>
array_t<type_t, vt> mem_to_reg_thread(it_t mem, int tid, int count, type_t (&shared)[shared_size], sycl::nd_item<1>& item) {

    array_t<type_t, vt> x = mem_to_reg_strided<nt, vt, vt0>(mem, tid, count);
    reg_to_shared_strided<nt, vt>(x, tid, shared, item);
    array_t<type_t, vt> y = shared_to_reg_thread<nt, vt>(shared, tid, item);
    return y;
}

template<int nt, int vt, int vt0 = vt, typename input_it, typename output_it>
void mem_to_mem(input_it input, int tid, int count,
                            output_it output) {
    typedef typename std::iterator_traits<input_it>::value_type type_t;
    type_t x[vt];

    strided_iterate<nt, vt, vt0>([&](int i, int j) {
        x[i] = input[j];
    }, tid, count);
    strided_iterate<nt, vt, vt0>([&](int i, int j) {
        output[j] = x[i];
    }, tid, count);
}

////////////////////////////////////////////////////////////////////////////////
// memory<->memory

template<int nt, int vt, int vt0 = vt, typename type_t, typename it_t>
void mem_to_shared(it_t mem, int tid, int count, type_t* shared,
                               const sycl::nd_item<1> &item_ct1,
                               bool sync = true) {

    array_t<type_t, vt> x = mem_to_reg_strided<nt, vt, vt0>(mem, tid, count);
    strided_iterate<nt, vt, vt0>([&](int i, int j) {
        shared[j] = x[i];
    }, tid, count);
    /*
    DPCT1118:16: SYCL group functions and algorithms must be encountered in
    converged control flow. You may need to adjust the code.
    */
    /*
    DPCT1065:53: Consider replacing sycl::nd_item::barrier() with
    sycl::nd_item::barrier(sycl::access::fence_space::local_space) for better
    performance if there is no access to global memory.
    */
    if (sync) item_ct1.barrier();
}

template<int nt, int vt, typename type_t, typename it_t>
void shared_to_mem(const type_t* shared, int tid, int count,
                               it_t mem, const sycl::nd_item<3> &item_ct1, bool sync = true) {

    strided_iterate<nt, vt>([&](int i, int j) {
        mem[j] = shared[j];
    }, tid, count);
    /*
    DPCT1118:17: SYCL group functions and algorithms must be encountered in
    converged control flow. You may need to adjust the code.
    */
    /*
    DPCT1065:54: Consider replacing sycl::nd_item::barrier() with
    sycl::nd_item::barrier(sycl::access::fence_space::local_space) for better
    performance if there is no access to global memory.
    */
    if (sync) item_ct1.barrier();
}

////////////////////////////////////////////////////////////////////////////////
// reg<->reg

template<int nt, int vt, typename type_t, int shared_size>
array_t<type_t, vt> reg_thread_to_strided(array_t<type_t, vt> x, int tid, type_t (&shared)[shared_size], const sycl::nd_item<1> &item) {
    reg_to_shared_thread<nt>(x, tid, shared, item);
    return shared_to_reg_strided<nt, vt>(shared, tid, item);
}

template<int nt, int vt, typename type_t, int shared_size>
array_t<type_t, vt> reg_strided_to_thread(array_t<type_t, vt> x, int tid, type_t (&shared)[shared_size], const sycl::nd_item<1> &item) {
    reg_to_shared_strided<nt>(x, tid, shared, item);
    return shared_to_reg_thread<nt, vt>(shared, tid, item);
}



#endif //TESTONEAPI_LOADSTORE_H
