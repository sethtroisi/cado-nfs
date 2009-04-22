#ifndef __cplusplus
#define _GNU_SOURCE         /* asprintf */
#endif
#define _DARWIN_C_SOURCE    /* for asprintf. _ANSI_SOURCE must be undefined */

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
#include <string.h>

// C++ headers.
// #include <string>
#include <vector>
#include <algorithm>    // sort
#include <iostream>     // cout

using namespace std;

#include "matmul-sliced.h"
#ifndef DISABLE_ASM     /* Better define this if we don't have GNU as */
#ifdef  __x86_64
#include "matmul-sliced-asm.h"
#define ENABLE_ASM
#endif
#endif

#ifdef  __GNUC__
#define ASM_COMMENT(x)  __asm__("#\t" x "\n")
#else
#define ASM_COMMENT(x)  /**/
#endif

#include "readmat.h"
#include "abase.h"

// #define L1_CACHE_SIZE   32768
// take only 3/4 of the L1 cache.
#define L1_CACHE_SIZE   24576
// for 8-bytes abt values, this gives 3072 items.

/* Here is how the matrix is stored in memory, in the
 * "matmul_sliced_data_s" type.  The matrix is cut in slices, where each
 * slice is a set of contiguous rows. The size of a slice is tune so as
 * to fill the L1 cache as per the macro above, with some adjustment to
 * handle non-exact divisibility: some slices will have packbase rows and
 * some others packbase+1 rows.  Within a slice, entries are stored as a
 * list of pairs
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

/* This extension is used to distinguish between several possible
 * implementations of the product. The upper word correspond to the
 * implementation itself, the lower one to the n-th binary incompatible
 * change (make sure to bump it) */
#define MM_EXTENSION   "-sliced"

#define MM_MAGIC_FAMILY        0xa000UL
#define MM_MAGIC_VERSION       0x1007UL
#define MM_MAGIC (MM_MAGIC_FAMILY << 16 | MM_MAGIC_VERSION)

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

struct matmul_sliced_data_s {
    /* repeat the fields from the public_ interface */
    struct matmul_public_s public_[1];
    /* now our private fields */
    abobj_t xab;
    data_t data;
    vector<slice_info> dslices_info;
    int iteration[2];
    void push(uint32_t x)
    {
        ASSERT_ALWAYS(x >> 16 == 0);
        data.push_back(x);
    }
    void push32(uint64_t x)
    {
        ASSERT_ALWAYS(x >> 32 == 0);
        data.push_back(x & ((1u << 16) - 1));
        x >>= 16;
        data.push_back(x & ((1u << 16) - 1));
    }
    static inline uint32_t read32(data_t::const_iterator & q) {
        uint32_t res;
        res = *q++;
        res |= ((uint32_t) *q++) << 16;
        return res;
    }
    static inline uint32_t read32(const uint16_t * & q) {
        uint32_t res;
        res = *q++;
        res |= ((uint32_t) *q++) << 16;
        return res;
    }
};

void matmul_sliced_clear(struct matmul_sliced_data_s * mm)
{
    mm->dslices_info.clear();
    mm->data.clear();
    delete mm;
}

static struct matmul_sliced_data_s * matmul_sliced_init()
{
    return new matmul_sliced_data_s;
}


struct matmul_sliced_data_s * matmul_sliced_build(abobj_ptr xx MAYBE_UNUSED, const char * filename, param_list pl)
{
    struct matmul_sliced_data_s * mm = matmul_sliced_init();

    abobj_init_set(mm->xab, xx);

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

    mm->public_->dim[0] = smat->nrows;
    mm->public_->dim[1] = smat->ncols;
    mm->public_->ncoeffs = 0;
    mm->iteration[0]=0;
    mm->iteration[1]=0;

    unsigned int npack = L1_CACHE_SIZE;
    if (pl) param_list_parse_uint(pl, "bwc_mm_l1_cache_size", &npack);
    npack /= abbytes(mm->xab,1);
    unsigned int nslices = iceildiv(i1-i0, npack);

    unsigned int nslices_index = mm->data.size();
    mm->push(0);
    mm->push(0);        // placeholder for alignment

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
        si.data_offset = mm->data.size();

        /*
           std::cout << "Packing " << npack << " rows from " << current
           << " to " << next << std::endl;
           */
        typedef std::vector<std::pair<uint32_t, uint32_t> > L_t;
        typedef L_t::const_iterator Lci_t;
        L_t L;
        for(unsigned int di = 0 ; di < npack ; di++) {
            /* Our dst pointer for read_matrix_row never changes, since
             * we write matrix rows over and over again in memory. */
            read_matrix_row(f,smat,smat->data,1);
            for(unsigned int j = 0 ; j < smat->data[0] ; j++) {
                L.push_back(std::make_pair(smat->data[1+j],di));
            }
            mm->public_->ncoeffs += smat->data[0];
        }
        /* L is the list of (column index, row index) of all
         * coefficients in the current horizontal slice */
        std::sort(L.begin(), L.end());

        mm->push32(npack);
        mm->push32(L.size());

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

            mm->push(dj); mm->push(i);
        }
        mm->data[nslices_index]++;
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
        mm->dslices_info.push_back(si);
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

    return mm;
}

struct matmul_sliced_data_s * matmul_sliced_reload_cache(abobj_ptr xx MAYBE_UNUSED, const char * filename, param_list pl MAYBE_UNUSED)
{
    char * base;
    FILE * f;

    {
        int rc = asprintf(&base, "%s" MM_EXTENSION ".bin", filename);
        FATAL_ERROR_CHECK(rc < 0, "out of memory");
    }

    f = fopen(base, "r");

    if (f == NULL) {
        // fprintf(stderr, "fopen(%s): %s\n", base, strerror(errno));
        fprintf(stderr, "no cache file %s\n", base);
        free(base);
        return NULL;
    }
    free(base);

    size_t n;
    size_t rc;
    struct matmul_sliced_data_s * mm = matmul_sliced_init();

    abobj_init_set(mm->xab, xx);

    mm->iteration[0]=0;
    mm->iteration[1]=0;

    uint32_t magic_check;
    rc = fread(&magic_check, sizeof(uint32_t), 1, f);
    FATAL_ERROR_CHECK(rc < 1, "No valid data in cached matrix file");
    
    if (magic_check != MM_MAGIC) {
        fprintf(stderr, "Wrong magic in cached matrix file\n");
        return NULL;
    }

    uint32_t nbytes_check;
    rc = fread(&nbytes_check, sizeof(uint32_t), 1, f);
    FATAL_ERROR_CHECK(rc < 1, "No valid data in cached matrix file");
    
    /* It's not fatal. It only deserves a warning */
    if (nbytes_check != abbytes(mm->xab, 1)) {
        fprintf(stderr, "Warning: cached matrix file fits data with different striding\n");
    }

    {
        uint32_t v;
        rc = fread(&v, sizeof(uint32_t), 1, f);
        FATAL_ERROR_CHECK(rc < 1, "No valid data in cached matrix file");
        mm->public_->dim[0] = v;
    }
    {
        uint32_t v;
        rc = fread(&v, sizeof(uint32_t), 1, f);
        FATAL_ERROR_CHECK(rc < 1, "No valid data in cached matrix file");
        mm->public_->dim[1] = v;
    }
    {
        uint64_t v;
        rc = fread(&v, sizeof(uint64_t), 1, f);
        FATAL_ERROR_CHECK(rc < 1, "No valid data in cached matrix file");
        mm->public_->ncoeffs = v;
    }
    {
        uint32_t v;
        rc = fread(&v, sizeof(uint32_t), 1, f);
        FATAL_ERROR_CHECK(rc < 1, "No valid data in cached matrix file");
        n = v;
    }

    mm->data.resize(n);
    rc = fread(&(mm->data.front()), sizeof(uint16_t), n, f);
    FATAL_ERROR_CHECK(rc < n, "Short read from cached matrix file");

    fclose(f);

    data_t::const_iterator q = mm->data.begin();

    return mm;
}

void matmul_sliced_save_cache(struct matmul_sliced_data_s * mm, const char * filename)
{
    char * base;
    FILE * f;

    {
        int rc = asprintf(&base, "%s" MM_EXTENSION ".bin", filename);
        FATAL_ERROR_CHECK(rc < 0, "out of memory");
    }

    f = fopen(base, "w");
    if (f == NULL) {
        fprintf(stderr, "fopen(%s): %s\n", base, strerror(errno));
        free(base);
        exit(1);
    }
    free(base);

    size_t n = mm->data.size();
    uint32_t magic = MM_MAGIC;
    size_t rc;

    rc = fwrite(&magic, sizeof(uint32_t), 1, f);
    FATAL_ERROR_CHECK(rc < 1, "Cannot write to cached matrix file");

    uint32_t nbytes = abbytes(mm->xab, 1);
    rc = fwrite(&nbytes, sizeof(uint32_t), 1, f);
    FATAL_ERROR_CHECK(rc < 1, "Cannot write to cached matrix file");

    {
        uint32_t v = mm->public_->dim[0];
        rc = fwrite(&v, sizeof(uint32_t), 1, f);
        FATAL_ERROR_CHECK(rc < 1, "Cannot write to cached matrix file");
    }
    {
        uint32_t v = mm->public_->dim[1];
        rc = fwrite(&v, sizeof(uint32_t), 1, f);
        FATAL_ERROR_CHECK(rc < 1, "Cannot write to cached matrix file");
    }
    {
        uint64_t v = mm->public_->ncoeffs;
        rc = fwrite(&v, sizeof(uint64_t), 1, f);
        FATAL_ERROR_CHECK(rc < 1, "Cannot write to cached matrix file");
    }
    {
        uint32_t v = n;
        rc = fwrite(&v, sizeof(uint32_t), 1, f);
        FATAL_ERROR_CHECK(rc < 1, "Cannot write to cached matrix file");
    }

    rc = fwrite(&(mm->data.front()), sizeof(uint16_t), n, f);
    FATAL_ERROR_CHECK(rc < n, "Short write to cached matrix file");

    fclose(f);
}

void matmul_sliced_mul(struct matmul_sliced_data_s * mm, abt * dst, abt const * src, int d)
{
    ASM_COMMENT("multiplication code");
    const uint16_t * q = &(mm->data.front());

    uint16_t nhstrips = *q++;
    q++;        // alignment.
    uint32_t i = 0;
    abobj_ptr x = mm->xab;

    if (d == 1) {
#ifdef  SLICE_STATS
        std::vector<slice_info>::iterator sit = dslices_info.begin();
        double tick = oncpu_ticks();
#endif
        for(uint16_t s = 0 ; s < nhstrips ; s++) {
            uint32_t nrows_packed = matmul_sliced_data_s::read32(q);
            abt * where = dst + aboffset(x, i);
            abzero(x, where, nrows_packed);
            ASM_COMMENT("critical loop");
            abt const * from = src;
            /* Make sure that the assembly function is only called if it
             * matches correctly the abase header !! */
#if defined(ENABLE_ASM) && defined(ABASE_U64_H_)
            q = matmul_sliced_asm(x, where, from, (uint16_t const *) q);
#else
            /* The external function must have the same semantics as this
             * code block */
            uint32_t ncoeffs_slice = matmul_sliced_data_s::read32(q);
            uint32_t j = 0;
            for(uint32_t c = 0 ; c < ncoeffs_slice ; c++) {
                j += *q++;
                uint32_t di = *q++;
                abadd(x, where + aboffset(x, di), from + aboffset(x, j));
            }
#endif
            ASM_COMMENT("end of critical loop");
            i += nrows_packed;
#ifdef  SLICE_STATS
            if (!mm->dslices_info.empty()) {
                double ntick = oncpu_ticks();
                sit->spent_time += ntick - tick;
                tick = ntick;
                sit->npasses++;
                sit++;
            }
#endif
        }
    } else {
        /* d == 0 */
        /* BEWARE, it's wildly sub-optimal ! */
        if (mm->iteration[d] == 10) {
            fprintf(stderr, "Warning: Doing many iterations with bad code\n");
        }

        abzero(x, dst, mm->public_->dim[1]);
        for(uint16_t s = 0 ; s < nhstrips ; s++) {
            uint32_t j = 0;
            uint32_t nrows_packed = matmul_sliced_data_s::read32(q);
            ASM_COMMENT("critical loop");
            uint32_t ncoeffs_slice = matmul_sliced_data_s::read32(q);
            for(uint32_t c = 0 ; c < ncoeffs_slice ; c++) {
                j += *q++;
                uint32_t di = *q++;
                abadd(x, dst + aboffset(x, j), src + aboffset(x, i + di));
            }
            i += nrows_packed;
            ASM_COMMENT("end of critical loop");
        }
    }
    ASM_COMMENT("end of multiplication code");
    mm->iteration[d]++;
}

void matmul_sliced_report(struct matmul_sliced_data_s * mm MAYBE_UNUSED) {
#ifdef  SLICE_STATS
    if (mm->dslices_info.empty())
        return;
    std::ofstream o("dj.stats");
    o << "// Report of timing per slice\n";
    o << "// <snum> <nrows> <ncoeffs> <dj> <npasses> <spent> <M0>\n";
    for(unsigned int s = 0 ; s < dslices_info.size() ; s++) {
        const slice_info& t(mm->dslices_info[s]);
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

void matmul_sliced_auxv(struct matmul_sliced_data_s * mm MAYBE_UNUSED, int op MAYBE_UNUSED, va_list ap MAYBE_UNUSED)
{
}

void matmul_sliced_aux(struct matmul_sliced_data_s * mm, int op, ...)
{
    va_list ap;
    va_start(ap, op);
    matmul_sliced_auxv (mm, op, ap);
    va_end(ap);
}
/* vim: set sw=4: */
