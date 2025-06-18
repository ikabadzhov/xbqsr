#ifndef TESTONEAPI_SORT_NETWORK_H
#define TESTONEAPI_SORT_NETWORK_H

////////////////////////////////////////////////////////////////////////////////
// Odd-even transposition sorting network. Sorts keys and values in-place in
// register.
// http://en.wikipedia.org/wiki/Odd%E2%80%93even_sort

#include "../device_utils/types.h"

template<typename type_t, int vt, typename comp_t>
array_t<type_t, vt> odd_even_sort(array_t<type_t, vt> x, comp_t comp, int flags = 0) {
    iterate<vt>([&](int I) {
        #pragma unroll
        for(int i = 1 & I; i < vt - 1; i += 2) {
            if((0 == ((2<< i) & flags)) && comp(x[i + 1], x[i])) {
                swap_t(x[i], x[i + 1]);
            }
        }
    });
    return x;
}

template<typename key_t, typename val_t, int vt, typename comp_t>
kv_array_t<key_t, val_t, vt> odd_even_sort(kv_array_t<key_t, val_t, vt> x, comp_t comp, int flags = 0) {
    iterate<vt>([&](int I) {
        #pragma unroll
        for(int i = 1 & I; i < vt - 1; i += 2) {
            if((0 == ((2<< i) & flags)) && comp(x.keys[i + 1], x.keys[i])) {
                swap_t(x.keys[i], x.keys[i + 1]);
                swap_t(x.vals[i], x.vals[i + 1]);
            }
        }
    });
    return x;
}

#endif //TESTONEAPI_SORT_NETWORK_H
