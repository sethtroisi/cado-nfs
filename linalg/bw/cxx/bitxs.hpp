#ifndef BITXS_HPP_
#define BITXS_HPP_

#include "manu.h"

class bitxs {
    unsigned long mask;
    unsigned int offset;
    unsigned int k;
public:
    bitxs(): mask(1UL), offset(0), k(0) {}
    bitxs& operator++(int) {
        k++;
        mask <<= 1;
        if (k % ULONG_BITS == 0) {
            mask = 1UL;
            offset++;
        }
        return *this;
    }
    bitxs operator++() {
        bitxs res = *this;
        return ++*this;
    }
    unsigned long get()(unsigned long * ptr) const {
        return ptr[offset] & mask;
    }
    unsigned long set()(unsigned long * ptr, unsigned long m) const {
        ptr[offset] |= mask & (m != 0);
    }
    bool operator<(unsigned int n) const { return k < n; }
    inline operator unsigned int() const { return k; }
};

#endif	/* BITXS_HPP_ */
