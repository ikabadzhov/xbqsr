#ifndef TESTONEAPI_OPERATOR_H
#define TESTONEAPI_OPERATOR_H

#include <sycl/sycl.hpp>



////////////////////////////////////////////////////////////////////////////////
// Device-side comparison operators.

template<typename type_t>
struct less_t /*: public std::binary_function<type_t, type_t, bool>*/ {
    SYCL_EXTERNAL bool operator()(type_t a, type_t b) const {
        return a < b;
    }
};
template<typename type_t>
struct less_equal_t /*: public std::binary_function<type_t, type_t, bool>*/ {
    SYCL_EXTERNAL bool operator()(type_t a, type_t b) const {
        return a <= b;
    }
};
template<typename type_t>
struct greater_t /*: public std::binary_function<type_t, type_t, bool>*/ {
    SYCL_EXTERNAL bool operator()(type_t a, type_t b) const {
        return a > b;
    }
};
template<typename type_t>
struct greater_equal_t /*: public std::binary_function<type_t, type_t, bool>*/ {
    SYCL_EXTERNAL bool operator()(type_t a, type_t b) const {
        return a >= b;
    }
};
template<typename type_t>
struct equal_to_t /*: public std::binary_function<type_t, type_t, bool>*/ {
    SYCL_EXTERNAL bool operator()(type_t a, type_t b) const {
        return a == b;
    }
};
template<typename type_t>
struct not_equal_to_t /*: public std::binary_function<type_t, type_t, bool>*/ {
    SYCL_EXTERNAL bool operator()(type_t a, type_t b) const {
        return a != b;
    }
};

////////////////////////////////////////////////////////////////////////////////
// Device-side arithmetic operators.

template<typename type_t>
struct plus_t /*: public std::binary_function<type_t, type_t, type_t>*/ {
    SYCL_EXTERNAL type_t operator()(type_t a, type_t b) const {
        return a + b;
    }
};

template<typename type_t>
struct minus_t /*: public std::binary_function<type_t, type_t, type_t>*/ {
    SYCL_EXTERNAL type_t operator()(type_t a, type_t b) const {
        return a - b;
    }
};

template<typename type_t>
struct multiplies_t /*: public std::binary_function<type_t, type_t, type_t>*/ {
    SYCL_EXTERNAL type_t operator()(type_t a, type_t b) const {
        return a * b;
    }
};

template<typename type_t>
struct maximum_t  /*: public std::binary_function<type_t, type_t, type_t>*/ {
    SYCL_EXTERNAL type_t operator()(type_t a, type_t b) const {
        return max(a, b);
    }
};

template<typename type_t>
struct minimum_t  /*: public std::binary_function<type_t, type_t, type_t>*/ {
    SYCL_EXTERNAL type_t operator()(type_t a, type_t b) const {
        return min(a, b);
    }
};


template<typename type_t>
SYCL_EXTERNAL void swap_t(type_t &a, type_t &b) {
    type_t c = a;
    a = b;
    b = c;
}


#endif //TESTONEAPI_OPERATOR_H
