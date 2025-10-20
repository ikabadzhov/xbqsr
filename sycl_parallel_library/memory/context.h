#ifndef TESTONEAPI_CONTEXT_H
#define TESTONEAPI_CONTEXT_H

#include <cstddef>
#include <oneapi/dpl/algorithm>

enum memory_space_t {
    memory_space_device = 0,
    memory_space_host = 1,
    memory_space_shared = 2,
    none = 3
};


template<typename type_t>
class mem_t {
    sycl::queue _q;
    type_t* _pointer;
    size_t _size;
    memory_space_t _space;

public:
    void swap(mem_t& rhs) {
        std::swap(_pointer, rhs._pointer);
        std::swap(_size, rhs._size);
        std::swap(_space, rhs._space);
        std::swap(_q, rhs._q);
    }

    mem_t() :  _pointer(nullptr), _size(0),
              _space(memory_space_device) { }

    mem_t& operator=(const mem_t& rhs) = delete;
    mem_t(const mem_t& rhs) = delete;

    mem_t(size_t size, sycl::queue& q,
          memory_space_t space = memory_space_device) :
            _q(q), _pointer(nullptr), _size(size), _space(space) {
        if (space==memory_space_device) {
            _pointer = sycl::malloc_device<type_t>(size, q);
        } else if (space==memory_space_shared) {
            _pointer = sycl::malloc_shared<type_t>(size, q);
        } else {
            _pointer = sycl::malloc_host<type_t>(size, q);
        }
    }

    mem_t(mem_t&& rhs)  noexcept : mem_t() {
        swap(rhs);
    }
    mem_t& operator=(mem_t&& rhs)  noexcept {
        swap(rhs);
        return *this;
    }

    ~mem_t() {
        if (_pointer) {
            sycl::free(_pointer, _q);
        }
        _pointer = nullptr;
        _size = 0;
        _q.wait_and_throw();
    }

    sycl::queue& context() { return _q; }
    [[nodiscard]] size_t size() const { return _size; }
    type_t* data() const { return _pointer; }
    [[nodiscard]] memory_space_t space() const { return _space; }
};

#endif //TESTONEAPI_CONTEXT_H
