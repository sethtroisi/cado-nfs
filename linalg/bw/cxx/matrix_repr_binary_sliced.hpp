#ifndef MATRIX_REPR_BINARY_SLICED_HPP_
#define MATRIX_REPR_BINARY_SLICED_HPP_

#include <iterator>
#include <istream>
#include <iostream>

#include "matrix.hpp"
#include "matrix_line.hpp"
#include "matrix_repr.hpp"
#include <vector>

/* This should be more cache-friendly than matrix_repr_binary.
 * Activate with USE_SLICED_MATRIX in cxx/binary_ulong_traits.hpp.
 * Results indicate no speedup yet, though.
 */

class matrix_repr_binary_sliced {
public:
    struct matrix_rowset {
        typedef std::vector<uint32_t> data_t;
        typedef data_t::iterator ptr;
        typedef data_t::const_iterator const_ptr;
        data_t data;
        unsigned int nrows_slice;
        public:

	void fill(std::istream& mtx, matrix_slice const & slice)
	{
	    mtx.seekg(slice.pos);
	    std::istream_iterator<matrix_line> mit(mtx);
            unsigned int npack = 2048;
            data.push_back(0);
            for(unsigned int i = slice.i0 ; i < slice.i1 ; i += npack) {
                unsigned int next = i + npack;
                if (next > slice.i1) {
                    npack = slice.i1 - i;
                    next = slice.i1;
                }
                std::cout << "Packing " << npack << " rows from " << i
                    << " to " << next << std::endl;
                typedef std::vector<std::pair<uint32_t, uint32_t> > L_t;
                typedef L_t::const_iterator Lci_t;
                L_t L;
                for(unsigned int di = 0 ; di < npack ; di++) {
                    matrix_line l = *mit++;
                    typedef matrix_line::const_iterator lit_t;
                    for(lit_t lit = l.begin() ; lit != l.end() ; lit++) {
                        L.push_back(std::make_pair(lit->first,di));
                    }
                }
                std::sort(L.begin(), L.end());
                /* Now L is sorted by column. Group together the blocks
                 * corresponding to one given column. */

                typedef std::vector<uint32_t> touchedrows_t;
                typedef std::vector<std::pair<uint32_t, touchedrows_t > > C_t;
                C_t C;
                Lci_t k0,k1;
                for(k0 = L.begin() ; k0 != L.end() ; k0 = k1) {
                    k1 = k0;
                    uint32_t j = k0->first;
                    touchedrows_t touched;
                    for( ; k1 != L.end() && k1->first == j ; ++k1) {
                        touched.push_back(k1->second);
                    }
                    C.push_back(make_pair(j, touched));
                }
                L.clear();

                data.push_back(npack);
                data.push_back(C.size());
                uint32_t j = 0;
                for(C_t::const_iterator cp = C.begin() ; cp != C.end() ; ++cp) {
                    data.push_back(cp->first - j); j = cp->first;
                    data.push_back(cp->second.size());
                    data.insert(data.end(),
                            cp->second.begin(),
                            cp->second.end());
                }
                data[0]++;
            }
	}

	template<typename traits>
	void mul(
		typename traits::wide_scalar_t * dst,
		const typename traits::scalar_t * src) const
	{
            asm("# multiplication code\n");
	    const_ptr q = data.begin();
            uint32_t nhstrips = *q++;
            uint32_t i = 0;
            for(uint32_t s = 0 ; s < nhstrips ; s++) {
                uint32_t j = 0;
                uint32_t nrows_packed = *q++;
                for(uint32_t di = 0 ; di < nrows_packed ; di++) {
                    traits::zero(dst[i+di]);
                }
                asm("# critical loop\n");
                uint32_t ncols_used = *q++;
                for(uint32_t c = 0 ; c < ncols_used ; c++) {
                    j += *q++;
                    typename traits::scalar_t x = src[j];
                    uint32_t nrows_touched_bycol = *q++;
                    for(uint32_t r = 0 ; r < nrows_touched_bycol ; r++) {
                        uint32_t di = *q++;
                        traits::addmul(dst[i + di], x);
                    }
                }
                i += nrows_packed;
                asm("# end of critical loop\n");
            }
            asm("# end of multiplication code\n");
	}
    };	/* matrix_rowset */
};

/* vim: set sw=4: */

#endif	/* MATRIX_REPR_BINARY_SLICED_HPP_ */
