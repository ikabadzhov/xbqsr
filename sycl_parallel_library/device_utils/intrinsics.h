#pragma once

#include <sycl/sycl.hpp>
#include <assert.h>

#ifndef MEMBERMASK
    #define MEMBERMASK 0xffffffff
#endif


//unsigned int ballot(int predicate,
//                const sycl::nd_item<1> &item) {
//    //unsigned y = 0;
//   // assert(item.get_sub_group().get_group_range() == 32);
//
////    y = sycl::reduce_over_group(item.get_sub_group(),
////            (mask & (0x1 << item.get_sub_group().get_local_linear_id())) && predicate ? (0x1 << item.get_sub_group().get_local_linear_id()) : 0,
////            sycl::ext::oneapi::plus<>());
//
//    return 0;
//}