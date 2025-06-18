#ifndef TESTONEAPI_MEMORY_H
#define TESTONEAPI_MEMORY_H

#include "context.h"


template <typename type_t>
void data_transfer(type_t *dest, const type_t *source, size_t count, sycl::queue& q, bool sync = true) {
    if(count) {
        q.memcpy(dest, source, sizeof(type_t) * count);
        if (sync) {
            q.wait_and_throw();
        }
    }

}


template<typename type_t, typename func_t>
mem_t<type_t> fill_function(func_t f, size_t count, sycl::queue& context) {
    mem_t<type_t> mem(count, context);
    type_t* p = mem.data();
    transform([=](int index) {
        p[index] = f(index);
    }, count, context);
    return mem;
}

#endif //TESTONEAPI_MEMORY_H
