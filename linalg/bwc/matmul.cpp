#define __STDC_LIMIT_MACROS
/* Manage the in-memory data for the matrix */
/* It's in C++ because the STL is handy, but that's really all there is
 * to it... */

#include <cstdio>
#include <cstdlib>
#include <cerrno>

#include <climits>
#include <cmath>

#include <stdint.h>

// C++ headers.
#include <string>
#include <vector>
#include <algorithm>    // sort
#include <iostream>     // cout

using namespace std;

#include "matmul.h"
#include "readmat.h"
#include "abase.h"
#include "manu.h"

/* Here is how the matrix is stored in memory, in the "matmul_data_s" type.
 * The matrix is cut in slices, where each slice is a set of contiguous
 * rows. The size of a slice is hardcoded (?) as about 3072, with some
 * adjustment to handle non-exact divisibility: some slices will have
 * packbase rows and some others packbase+1 rows.
 * Within a slice, entries are stored as a list of pairs
 *   (column index, row index),
 * sorted according to column index. Then, only the difference between
 * two consecutive column indices is actually stored, so that it will fit
 * in 16 bits. Also the row index fits in 16 bits, so that a slice is
 * actually a list of 16-bit unsigned integers. 
 * All the slices are stored in a big array of uint16_t, in the data
 * field. Within this data field, the information is organized like so:
 *   data[0] : total number of slices
 *   data[1] : number of rows in slice number 1
 *   data[2..3]: number of entries in slice number 1, say k_1
 *   data[4..(4+2*k_1)]: list of entries of slice number 1
 *   data[(4+2k_1+1)]: number of rows in slice number 2
 *   etc...
 * Additionally, for each slice, a slice_info structure is stored, giving
 * statistics about the slice, and the offset in the data array where it
 * is stored.
 */


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
    size_t data_offset;
};

typedef vector<uint16_t> data_t;

struct matmul_data_s {
    abobj_t xab;
    data_t data;
    vector<slice_info> dslices_info;
    void push(uint32_t x)
    {
        BUG_ON(x > UINT16_MAX);
        data.push_back(x);
    }
};

#define MM      ((struct matmul_data_s *) mm)

void matmul_clear(matmul_ptr mm)
{
    MM->dslices_info.clear();
    MM->data.clear();
    free(MM);
}

matmul_ptr matmul_build(abobj_ptr xx, const char * filename)
{
    struct matmul_data_s * mm;
    mm = (struct matmul_data_s *) malloc(sizeof(struct matmul_data_s));

    MM->xab = xx;

    sparse_mat_t smat;
    FILE * f;

    sparse_mat_init(smat);

    f = fopen(filename, "r");
    if (f == NULL) {
        fprintf(stderr, "fopen(%s): %s\n", filename, strerror(errno));
        exit(1);
    }
    read_matrix_header(f, smat);

    uint32_t i0 = 0;
    uint32_t i1 = smat->nrows;

    unsigned int npack = 3072;
    unsigned int nslices = iceildiv(i1-i0, npack);
    MM->push(0);
    unsigned int packbase = (i1-i0) / nslices;
    /* How many slices of size packbase+1 */
    unsigned int nslices1 = (i1-i0) % nslices;
    /* How many slices of size packbase */
    unsigned int nslices0 = nslices - nslices1;

    unsigned int current;
    unsigned int next = i0;

    /* There we're handling the horizontal strips */
    unsigned int s;

    for(s = 0 ; s < nslices ; s++) {
        current = next;
        npack = packbase + (s < nslices1);
        next = current + npack;

        slice_info si;
        memset(&si,0,sizeof(si));
        si.nrows = npack;
        si.data_offset = MM->data.size();

        /*
           std::cout << "Packing " << npack << " rows from " << current
           << " to " << next << std::endl;
           */
        typedef std::vector<std::pair<uint32_t, uint32_t> > L_t;
        typedef L_t::const_iterator Lci_t;
        L_t L;
        for(unsigned int di = 0 ; di < npack ; di++) {
            /* Our dst pointer never changes, since we write matrix rows
             * over and over again in memory. */
            read_matrix_row(f,smat,smat->data,1);
            for(unsigned int j = 0 ; j < smat->data[0] ; j++) {
                L.push_back(std::make_pair(smat->data[1+j],di));
            }
        }
        /* L is the list of (column index, row index) of all
         * coefficients in the current horizontal slice */
        std::sort(L.begin(), L.end());

        MM->push(npack);
        MM->push(L.size() & ((1u << 16) - 1));
        MM->push(L.size() >> 16);

        uint32_t j = 0;
        uint32_t i = 0;
        int32_t di_max = 0;
        double sumdj=0;
        double sumdj2=0;
        double sumdi=0;
        double sumdi2=0;
        uint32_t weight=0;
        for(Lci_t lp = L.begin() ; lp != L.end() ; ++lp) {
            uint32_t dj = lp->first - j; j += dj;
            int32_t di = lp->second - i; i += di;
            weight++;

            if (di<0) di = -di;
            if (di > di_max) di_max = di;

            double ddj = dj;
            double ddi = di;
            sumdj+=ddj;
            sumdj2+=ddj*ddj;
            sumdi+=ddi;
            sumdi2+=ddi*ddi;

            /* If di exceeds +/- 128, then we can emulate the
             * wide displacement by something like:
             * (di_max,0)(0,0)(di_max,0)(0,0)...(di_remainder,dj)
             *
             * It's also really likely that one would gain by
             * moving _both_ up and down in the strip, so that
             * the negative offsets remain significant
             * (otherwise, we'd rather do a cmov based on the
             * sign, or something).
             */

            MM->push(dj); MM->push(i);
        }
        MM->data[0]++;
        double dj2_avg = sumdj2 / weight;
        double dj_avg = sumdj / weight;
        double dj_sdev = sqrt(dj2_avg-dj_avg*dj_avg);
        double di2_avg = sumdi2 / weight;
        double di_avg = sumdi / weight;
        double di_sdev = sqrt(di2_avg-di_avg*di_avg);
        si.dj_avg = dj_avg;
        si.dj_sdev = dj_sdev;
        si.di_avg = di_avg;
        si.di_sdev = di_sdev;
        si.ncoeffs = weight;

#ifdef  SPARSE_STRIPS
        if (si.dj_avg > 0.5) {
            std::cerr << "hstrip #" << s
                << " is too sparse ; switching to vstrips\n";
            // data.erase(data.begin()+si.data_offset, data.end());
            // break;
        }
#endif
        MM->dslices_info.push_back(si);
    }
    /* There's an other option for the dense strips. For a given
     * set of source values (fitting in L2), store delta_j's for
     * a given destination value. Could actually be eight bits.
     * When the direction changes, increase the destination
     * pointer
     */
#ifdef  SPARSE_STRIPS
    /* We haven't finished the set of horizontal strips. Treat
     * them as vertical sparse strips then. */
#endif
    sparse_mat_clear(smat);
    fclose(f);

    /* TODO gather info directly from the dslices_info vector
     * now that it has so much data */
    if (nslices1) {
        std::cout
            << "// " << nslices1
            << " sub-slices of " << (packbase+1) << " rows\n";
    }
    std::cout
        << "// " << nslices0
        << " sub-slices of " << packbase << " rows\n";
#ifdef  SLICE_STATS
    std::cout << "info per sub-slice:";
    typedef std::vector<slice_info>::const_iterator si_t;
    for(si_t s = dslices_info.begin() ; s != dslices_info.end() ; s++) {
        std::cout
            << (s-dslices_info.begin())
            << " " << s->dj_avg << " " << s->dj_sdev
            << " " << s->di_avg << " " << s->di_sdev
            << " " << s->di_max
            << "\n";
    }
    // std::cout << "\n";
#endif
    std::cout << std::flush;

    return (matmul_ptr)mm;
}

matmul_ptr matmul_reload_cache(abobj_ptr xx, const char * filename)
{
    string base(filename);
    FILE * f;

    base += ".bin";

    f = fopen(base.c_str(), "r");
    if (f == NULL) {
        // It is not an error.
        // fprintf(stderr, "fopen(%s): %s\n", base.c_str(), strerror(errno));
        return NULL;
    }

    size_t n;
    unsigned int rc;
    struct matmul_data_s * mm;
    mm = (struct matmul_data_s *) malloc(sizeof(struct matmul_data_s));

    MM->xab = xx;


    rc = fread(&n, sizeof(size_t), 1, f);
    FATAL_ERROR_CHECK(rc < 1, "No valid data in cached matrix file");

    MM->data.resize(n);
    rc = fread(&(MM->data.front()), sizeof(uint16_t), n, f);
    FATAL_ERROR_CHECK(rc < n, "Short read from cached matrix file");

    fclose(f);

    return (matmul_ptr) mm;
}

void matmul_save_cache(matmul_ptr mm, const char * filename)
{
    string base(filename);
    FILE * f;

    base += ".bin";

    f = fopen(base.c_str(), "w");
    if (f == NULL) {
        fprintf(stderr, "fopen(%s): %s\n", base.c_str(), strerror(errno));
        exit(1);
    }

    size_t n = MM->data.size();
    unsigned int rc;

    rc = fwrite(&n, sizeof(size_t), 1, f);
    FATAL_ERROR_CHECK(rc < 1, "Cannot write to cached matrix file");

    rc = fwrite(&(MM->data.front()), sizeof(uint16_t), n, f);
    FATAL_ERROR_CHECK(rc < n, "Short write to cached matrix file");

    fclose(f);
}

void matmul(matmul_ptr mm, abt * dst, abt const * src)
{
    asm("# multiplication code\n");
    data_t::const_iterator q = MM->data.begin();
    uint16_t nhstrips = *q++;
    uint32_t i = 0;
#ifdef  SLICE_STATS
    std::vector<slice_info>::iterator sit = dslices_info.begin();
    double tick = oncpu_ticks();
#endif
    abobj_ptr x = MM->xab;
    for(uint16_t s = 0 ; s < nhstrips ; s++) {
        uint32_t j = 0;
        uint16_t nrows_packed = *q++;
        abzero(x, dst + aboffset(x, i), nrows_packed);
        asm("# critical loop\n");
        uint32_t ncoeffs_slice;
        ncoeffs_slice = *q++;
        ncoeffs_slice |= ((uint32_t) *q++) << 16;
        for(uint32_t c = 0 ; c < ncoeffs_slice ; c++) {
            j += *q++;
            uint32_t di = *q++;
            abadd(x, dst + aboffset(x, i + di), src + aboffset(x, j));
        }
        i += nrows_packed;
        asm("# end of critical loop\n");
#ifdef  SLICE_STATS
        if (!MM->dslices_info.empty()) {
            double ntick = oncpu_ticks();
            sit->spent_time += ntick - tick;
            tick = ntick;
            sit->npasses++;
            sit++;
        }
#endif
    }
    asm("# end of multiplication code\n");
}

void matmul_report(matmul_ptr mm MAYBE_UNUSED) {
#ifdef  SLICE_STATS
    if (MM->dslices_info.empty())
        return;
    std::ofstream o("dj.stats");
    o << "// Report of timing per slice\n";
    o << "// <snum> <nrows> <ncoeffs> <dj> <npasses> <spent> <M0>\n";
    for(unsigned int s = 0 ; s < dslices_info.size() ; s++) {
        const slice_info& t(MM->dslices_info[s]);
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

/* vim: set sw=4: */
