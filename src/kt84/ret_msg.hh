#pragma once
#include <iostream>

namespace kt84 {
    template <typename T>
    inline T ret_msg(const T& ret_val, const char* msg) {
        std::cout << msg << std::endl;
        return ret_val;
    }
    inline void ret_msg(const char* msg) {
        std::cout << msg << std::endl;
    }
}
