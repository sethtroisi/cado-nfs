#ifndef MATRIX_REPR_BINARY_SLICED_HPP_
#define MATRIX_REPR_BINARY_SLICED_HPP_

#include <iterator>
#include <istream>
#include <iostream>
#include <fstream>
#include <algorithm>

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

/* Use this whenever you encounter the dreaded "x > UINT16_MAX" bug. It
 * has a non-trivial negative impact on the program's speed and memory
 * requirements.
 */
#define xxxUSE_32BIT_MATRIX_INDICES

class matrix_repr_binary_sliced {
public:
    struct matrix_rowset {
#ifdef  USE_32BIT_MATRIX_INDICES
        typedef std::vector<uint32_t> data_t;
#else
        typedef std::vector<uint16_t> data_t;
#endif
        typedef data_t::iterator ptr;
        typedef data_t::const_iterator const_ptr;
        data_t data;
        inline void push(uint32_t x) {
#ifndef  USE_32BIT_MATRIX_INDICES
            BUG_ON(x > UINT16_MAX);
#endif
            data.push_back(x);
        }

#ifdef  SLICE_STATS
        struct slice_info {
            unsigned int nrows;
            unsigned int ncoeffs;
            double spent_time;
            unsigned int npasses;
            double dj_avg;
            double dj_sdev;
            double di_avg;
            double di_sdev;
            int32_t di_max;
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

            /* There we're handling the horizontal strips */
            for(unsigned int s = 0 ; s < nslices ; s++) {
                i = next;
                npack = packbase + (s < nslices1);
                next = i + npack;

#ifdef  SLICE_STATS
                slice_info si;
                memset(&si,0,sizeof(si));
                si.nrows = npack;
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

                push(npack);
#ifdef  USE_32BIT_MATRIX_INDICES
                push(L.size());
#else
                push(L.size() & ((1u << 16) - 1));
                push(L.size() >> 16);
#endif
                uint32_t j = 0;
                uint32_t i = 0;
#ifdef  SLICE_STATS
                int32_t di_max = 0;
                double sumdj=0;
                double sumdj2=0;
                double sumdi=0;
                double sumdi2=0;
#endif
                for(Lci_t lp = L.begin() ; lp != L.end() ; ++lp) {
                    uint32_t dj = lp->first - j; j += dj;
                    int32_t di = lp->second - i; i += di;
#ifdef  SLICE_STATS
                    si.ncoeffs++;

                    if (di<0) di = -di;
                    if (di > di_max) di_max = di;

                    double ddj = dj;
                    double ddi = di;
                    sumdj+=ddj;
                    sumdj2+=ddj*ddj;
                    sumdi+=ddi;
                    sumdi2+=ddi*ddi;
#endif
                    push(dj); push(i);
                }
                data[0]++;
#ifdef  SLICE_STATS
                double dj2_avg = sumdj2 / slices_info[s].ncoeffs;
                double dj_avg = sumdj / slices_info[s].ncoeffs;
                double dj_sdev = std::sqrt(dj2_avg-dj_avg*dj_avg);
                double di2_avg = sumdi2 / slices_info[s].ncoeffs;
                double di_avg = sumdi / slices_info[s].ncoeffs;
                double di_sdev = std::sqrt(di2_avg-di_avg*di_avg);
                si.dj_avg = dj_avg;
                si.dj_sdev = dj_sdev;
                si.di_avg = di_avg;
                si.di_sdev = di_sdev;

                slices_info.push_back(si);
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
            std::cout << "info per sub-slice:";
            typedef std::vector<slice_info>::const_iterator si_t;
            for(si_t s = slices_info.begin() ; s != slices_info.end() ; s++) {
                std::cout
                    << (s-slices_info.begin())
                    << " " << s->dj_avg << " " << s->dj_sdev
                    << " " << s->di_avg << " " << s->di_sdev
                    << " " << s->di_max
                    << "\n";
            }
            // std::cout << "\n";
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
                uint32_t ncoeffs_slice;
#ifdef  USE_32BIT_MATRIX_INDICES
                ncoeffs_slice = *q++;
#else
                ncoeffs_slice = *q++;
                ncoeffs_slice |= ((uint32_t) *q++) << 16;
#endif
                for(uint32_t c = 0 ; c < ncoeffs_slice ; c++) {
                    j += *q++;
                    typename traits::scalar_t x = src[j];
                    uint32_t di = *q++;
                    traits::addmul(dst[i + di], x);
                }
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
                    << " " << t.dj_sdev
                    << " " << t.di_avg
                    << " " << t.di_sdev
                    << " " << t.di_max
                    << "\n";
            }
            o << std::flush;
#endif
        }
    };	/* matrix_rowset */
};

/* vim: set sw=4: */

#endif	/* MATRIX_REPR_BINARY_SLICED_HPP_ */
