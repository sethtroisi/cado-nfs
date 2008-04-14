#ifndef MATRIX_REPR_BINARY_SLICED_HPP_
#define MATRIX_REPR_BINARY_SLICED_HPP_

#include <iterator>
#include <istream>
#include <iostream>
#include <fstream>

#include "matrix.hpp"
#include "matrix_line.hpp"
#include "matrix_repr.hpp"
#include "limits.h"
#include <vector>
#include "ticks.hpp"

/* This should be more cache-friendly than matrix_repr_binary.
 * Activate with USE_SLICED_MATRIX in cxx/binary_ulong_traits.hpp.
 * Results indicate no speedup yet, though.
 */

/* Activate this flag to have a dj.stats file, containing some
 * detailed ns/coeff info.
 */
#define xxxSLICE_STATS

class matrix_repr_binary_sliced {
public:
    struct matrix_rowset {
        typedef std::vector<uint16_t> data_t;
        typedef data_t::iterator ptr;
        typedef data_t::const_iterator const_ptr;
        data_t data;
        inline void push(uint32_t x) {
            BUG_ON(x > UINT16_MAX);
            data.push_back(x);
        }

#ifdef  SLICE_STATS
        struct slice_info {
            unsigned int nrows;
            unsigned int ncoeffs;
            double spent_time;
            unsigned int npasses;
            double dj_avg;
        };

        mutable std::vector<slice_info> slices_info;
#endif
        unsigned int nslices0;
        unsigned int nslices1;
        unsigned int nslices;
        unsigned int packbase;

        public:

	void fill(std::istream& mtx, matrix_slice const & slice)
	{
            using namespace std;
	    mtx.seekg(slice.pos);
	    std::istream_iterator<matrix_line> mit(mtx);
            unsigned int npack = 3072;
            nslices = iceildiv(slice.i1-slice.i0, npack);
            data.push_back(0);
            packbase = (slice.i1-slice.i0) / nslices;
            /* How many slices of size packbase+1 */
            nslices1 = (slice.i1-slice.i0) % nslices;
            /* How many slices of size packbase */
            nslices0 = nslices - nslices1;

            unsigned int i;
            unsigned int next = slice.i0;

#ifdef  SLICE_STATS
            slices_info.assign(nslices, slice_info());
#endif

            for(unsigned int s = 0 ; s < nslices ; s++) {
                i = next;
                npack = packbase + (s < nslices1);
                next = i + npack;

#ifdef  SLICE_STATS
                memset(&(slices_info[s]),0,sizeof(slices_info));
                slices_info[s].nrows = npack;
#endif

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
#ifdef  SLICE_STATS
                double sumdj=0;
#endif
                for(Lci_t lp = L.begin() ; lp != L.end() ; ++lp) {
#ifdef  SLICE_STATS
                    slices_info[s].ncoeffs++;
                    sumdj+=lp->first - j;
#endif
                    push(lp->first - j); j = lp->first;
                    push(lp->second);
                }
                data[0]++;
#endif
#ifdef  SLICE_STATS
                slices_info[s].dj_avg = sumdj / slices_info[s].ncoeffs;
#endif
            }
        }

        void info() const {
            std::cout
                << "// " << nslices1
                << " sub-slices of " << (packbase+1) << " rows\n";
            if (nslices0) {
                std::cout
                    << "// " << nslices0
                    << " sub-slices of " << packbase << " rows\n";
            }
#ifdef  SLICE_STATS
            std::cout << "dj per sub-slice:";
            for(unsigned int s = 0 ; s < nslices ; s++) {
                std::cout << " " << slices_info[s].dj_avg;
            }
            std::cout << "\n";
#endif
            std::cout << std::flush;
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
#ifdef  SLICE_STATS
            std::vector<slice_info>::iterator sit = slices_info.begin();
            double tick = oncpu_ticks();
#endif
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
#ifdef  SLICE_STATS
                double ntick = oncpu_ticks();
                sit->spent_time += ntick - tick;
                tick = ntick;
                sit->npasses++;
                sit++;
#endif
            }
            asm("# end of multiplication code\n");
	}
        void report() const {
#ifdef  SLICE_STATS
            std::ofstream o("dj.stats");
            o << "// Report of timing per slice\n";
            o << "// <snum> <nrows> <ncoeffs> <dj> <npasses> <spent> <M0>\n";
            for(unsigned int s = 0 ; s < slices_info.size() ; s++) {
                const slice_info& t(slices_info[s]);
                o << s
                    << " " << t.nrows
                    << " " << t.ncoeffs
                    << " " << t.dj_avg
                    << " " << t.npasses
                    << " " << t.spent_time
                    << " " << (t.spent_time/t.npasses/t.ncoeffs * 1.0e9)
                    << "\n";
            }
            o << std::flush;
#endif
        }
    };	/* matrix_rowset */
};

/* vim: set sw=4: */

#endif	/* MATRIX_REPR_BINARY_SLICED_HPP_ */
