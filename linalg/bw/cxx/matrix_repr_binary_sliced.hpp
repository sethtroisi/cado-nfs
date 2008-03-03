#ifndef MATRIX_REPR_BINARY_SLICED_HPP_
#define MATRIX_REPR_BINARY_SLICED_HPP_

#include <iterator>
#include <istream>
#include <iostream>

#include "matrix.hpp"
#include "matrix_line.hpp"
#include "matrix_repr.hpp"
#include "limits.h"
#include <vector>

/* This should be more cache-friendly than matrix_repr_binary.
 * Activate with USE_SLICED_MATRIX in cxx/binary_ulong_traits.hpp.
 * Results indicate no speedup yet, though.
 */

class matrix_repr_binary_sliced {
public:
    struct matrix_rowset {
        typedef std::vector<uint16_t> data_t;
        typedef data_t::iterator ptr;
        typedef data_t::const_iterator const_ptr;
        data_t data;
        unsigned int nrows_slice;
        inline void push(uint32_t x) {
            BUG_ON(x > UINT16_MAX);
            data.push_back(x);
        }
        public:

	void fill(std::istream& mtx, matrix_slice const & slice)
	{
            uint32_t datamax = 0;
            using namespace std;
	    mtx.seekg(slice.pos);
	    std::istream_iterator<matrix_line> mit(mtx);
            unsigned int npack = 3072;
            data.push_back(0);
            for(unsigned int i = slice.i0 ; i < slice.i1 ; i += npack) {
                unsigned int next = i + npack;
                if (next > slice.i1) {
                    npack = slice.i1 - i;
                    next = slice.i1;
                }
                /*
                std::cout << "Packing " << npack << " rows from " << i
                    << " to " << next << std::endl;
                    */
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
                /* L is the list of (column index, row index) of all
                 * coefficients in the current horizontal slice */
                std::sort(L.begin(), L.end());

                /* Tetracapillectomy is the fact of having a separate
                 * array for each column ; in most cases it's a bad idea,
                 * since the corresponding array will be rather small,
                 * and its size will be quite erratic. Hence we will have
                 * many mispredictions.
                 *
                 * So it's not expected to be a smart move.
                 */
#ifdef TETRACAPILLECTOMY
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

                push(npack);
                push(C.size() & ((1u << 16) - 1));
                push(C.size() >> 16);
                uint32_t j = 0;
                for(C_t::const_iterator cp = C.begin() ; cp != C.end() ; ++cp) {
                    push(cp->first - j);
                    j = cp->first;
                    push(cp->second.size());
                    /* Take it easy for this one -- we know that values
                     * are bounded by npack anyway */
                    data.insert(data.end(),
                            cp->second.begin(),
                            cp->second.end());
                }
#else
                push(npack);
                push(L.size() & ((1u << 16) - 1));
                push(L.size() >> 16);
                uint32_t j = 0;
                for(Lci_t lp = L.begin() ; lp != L.end() ; ++lp) {
                    push(lp->first - j); j = lp->first;
                    push(lp->second);
                }
                data[0]++;
#endif


            }
            cout << "Max coeff in data array is " << datamax << endl;
	}

	template<typename traits>
	void mul(
		typename traits::wide_scalar_t * dst,
		const typename traits::scalar_t * src) const
	{
            asm("# multiplication code\n");
	    const_ptr q = data.begin();
            uint16_t nhstrips = *q++;
            uint32_t i = 0;
            for(uint16_t s = 0 ; s < nhstrips ; s++) {
                uint32_t j = 0;
                uint16_t nrows_packed = *q++;
                for(uint32_t di = 0 ; di < nrows_packed ; di++) {
                    traits::zero(dst[i+di]);
                }
                asm("# critical loop\n");
#ifdef  TETRACAPILLECTOMY
                uint32_t ncols_used;
                ncols_used = *q++;
                ncols_used |= ((uint32_t) *q++) << 16;
                for(uint32_t c = 0 ; c < ncols_used ; c++) {
                    j += *q++;
                    typename traits::scalar_t x = src[j];
                    uint16_t nrows_touched_bycol = *q++;
                    for(uint16_t r = 0 ; r < nrows_touched_bycol ; r++) {
                        uint32_t di = *q++;
                        traits::addmul(dst[i + di], x);
                    }
                }
#else
                uint32_t ncoeffs_slice;
                ncoeffs_slice = *q++;
                ncoeffs_slice |= ((uint32_t) *q++) << 16;
                for(uint32_t c = 0 ; c < ncoeffs_slice ; c++) {
                    j += *q++;
                    typename traits::scalar_t x = src[j];
                    uint32_t di = *q++;
                    traits::addmul(dst[i + di], x);
                }
#endif
                i += nrows_packed;
                asm("# end of critical loop\n");
            }
            asm("# end of multiplication code\n");
	}
    };	/* matrix_rowset */
};

/* vim: set sw=4: */

#endif	/* MATRIX_REPR_BINARY_SLICED_HPP_ */
