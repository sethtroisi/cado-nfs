#ifndef MATRIX_HPP_
#define MATRIX_HPP_

#include <istream>
#include <vector>
#include <boost/cstdint.hpp>

// #include "matrix_repr_prime.hpp"

#if 0
inline void fill_matrix_data(std::istream& mtx,
                std::streampos pos, uint32_t nc,
                uint i0, uint i1,
                uint32_t * idx, int32_t * val)
__attribute__((deprecated));

/* Legacy wrappers */

void fill_matrix_data(std::istream& mtx,
                std::streampos pos, uint32_t nc,
                uint i0, uint i1,
                uint32_t * idx, int32_t * val)
{
    typedef matrix_repr_prime T;

    T::ptr p;
    p.idx = idx;
    p.val = val;
    T::fill(mtx, pos, nc, i0, i1);
}
#endif

extern void count_matrix_coeffs(std::istream&,
		unsigned int,
		std::vector<std::streampos>&,
		std::vector<boost::uint32_t>&);

#endif	/* MATRIX_HPP_ */
