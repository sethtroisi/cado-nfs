#ifndef __cplusplus
#define _GNU_SOURCE         /* asprintf */
#endif
#define _DARWIN_C_SOURCE    /* for asprintf. _ANSI_SOURCE must be undefined */

#define __STDC_LIMIT_MACROS
#define __STDC_FORMAT_MACROS
/* Manage the in-memory data for the matrix */

/* It's in C++ because the STL is handy, but that's really all there is
 * to it ; a conversion to C would not be extremely difficult */

#include <cstdio>
#include <cstdlib>
#include <cerrno>

#include <climits>
#include <cmath>

#include <stdint.h>
#include <inttypes.h>
#include <string.h>

// C++ headers.
// #include <string>
#include <vector>
#include <list>
#include <algorithm>    // sort
#include <iostream>     // cout
#include "bwc_config.h"
using namespace std;

#include "matmul-bucket.h"
#include "abase.h"

/* Make sure that the assembly function is only called if it matches
 * correctly the abase header !! */
#if defined(HAVE_GCC_STYLE_AMD64_ASM) && defined(ABASE_U64_H_) && !defined(DISABLE_ASM)
#include "matmul-sub-small1.h"
#include "matmul-sub-small2.h"
#include "matmul-sub-large-fbd.h"
#include "matmul-sub-large-fbi.h"
#define ENABLE_ASM
#endif

#include "readmat-easy.h"
#include "matmul-common.h"

// take only a portion of the L1 cache.
#define L1_CACHE_SIZE   28000
// for 8-bytes abt values, this gives 3500 items.

// make sure that it's something each core can use (i.e. divide by two
// the L2 cache size for a dual-core w/ shared L2).
#define L2_CACHE_SIZE   1600000

/* To what should we compare the average dj value for determining when we
 * should switch from dense (small) to sparse (large) blocks */
#define DJ_CUTOFF1   2.0
/* Same, but for switching from large to huge blocks */
// In accordance with input data considered, 10.0 does not seem too
// small. And perhaps it's even a bit large.
#define DJ_CUTOFF2   4.0

/* How many buckets in large slices, and in large slices which are
 * contained in small slices. */
#define LSL_NBUCKETS_MAX      256

/* How many large slices MINIMUM must go in a huge slice ? Below this
 * value, we keep going producing large slices. */
#define HUGE_MPLEX_MIN        2

/* Here, I've seen 97 slices be handled best with 2 mplexed rounds. So
 * let's fix 64. */
#define HUGE_MPLEX_MAX        64

#define xxxDEBUG_BUCKETS

/* This implementation builds upon the ``sliced'' variant. 
 * The format used with the -sliced implementation is still used, but
 * there is also another format used for very dense strips. We refer to
 * matmul-sliced.cpp for the first format, which is called ``small1''
 * here (and there also).
 *
 * The ``small2'' format is specialised for denser strips.  Within an
 * horizontal strip of less than 4096 rows, when two non-zero columns are
 * never more than 16 positions apart, then the information (dj,i) is
 * stored in 16 bits instead of 32. This reduces the memory footprint
 * quite considerably, since many coefficients are in dense strips.
 * Speed-wise, small2 strips provide a nice speed-up sometimes (not
 * always).
 *
 * When we reach sparser areas, we use another method. Dense (small2 and
 * small1) slices having been discarded, the rest is split in ``large''
 * slices of LSL_NBUCKETS_MAX*256 rows at most (in fact we arrange for the
 * slices to have equal size, so it's a bit less).  For each such large
 * slice, we use scrap space equal to the L2 cache size for performing
 * matrix multiplication. While reading column vector data, we copy
 * coefficient to (up to) LSL_NBUCKETS_MAX ``buckets'', each representing
 * 256 possible row indices. The dispatching information, telling in
 * which bucket a given coefficient has to be copied, is an 8-bit integer
 * stored in a table called ``main''. Each entry in table ``main'' also
 * has a second 8-bit integer indicating the offset to the next useful
 * column entry (because we're considering very sparse blocks, it is not
 * at all uncommon to see average values ~ 10 here, even when
 * LSL_NBUCKETS_MAX==256, so that 65536 rows are considered
 * simultaneously).
 *
 * The number of columns considered while filling the scrap space is
 * limited by the scrap space size. Once it's filled, we ``apply the
 * buckets'', so as to be able to reuse our scrap space buffer for
 * another batch of columns. Thus a ``large slice'' is typically split
 * into several ``vertical blocks''. Extremely sparse large slices will
 * have only one vertical block.
 *
 * ``Applying the buckets'' means taking a list of column coefficients
 * which have been directed to the bucket because it is known that they
 * affect at least one of the 256 corresponding row indices on output. So
 * in addition the the N stored coefficients in the bucket, we have a
 * list of N 8-bit integers indicating which row index is to be modified.
 *
 *
 * That's it for the basic idea. Now on top of that, we use an extra
 * indirection level for very sparse blocks. Several ``large blocks'' are
 * considered together as part of a ``huge block''. Call this number
 * ``nlarge''. Source coefficients are first dispatched in nlarge ``big
 * buckets'' of appropriate size. The number of columns considered is
 * guided by the fact that we want none of these big buckets to exceed
 * the L2 cache size.  Afterwards, we do a second bispatching from each
 * of these big buckets (in turn) to (up to) LSL_NBUCKETS_MAX small
 * buckets in the same manner as above (except that we no longer need the
 * information on the link to the next useful column.
 */

/* This extension is used to distinguish between several possible
 * implementations of the product. The upper word correspond to the
 * implementation itself, the lower one to the n-th binary incompatible
 * change (make sure to bump it) */
#define MM_EXTENSION   "-bucket"

#define MM_MAGIC_FAMILY        0xa003UL

#define MM_MAGIC_VERSION       0x100fUL
#define MM_MAGIC (MM_MAGIC_FAMILY << 16 | MM_MAGIC_VERSION)

/* see matmul-basic.c */
#define MM_DIR0_PREFERS_TRANSP_MULT   1

template<typename T> inline T * ptrbegin(vector<T>& v) { return &v.front(); }
template<typename T> inline T const * ptrbegin(vector<T> const & v) { return &v.front(); }
template<typename T> inline T * ptrend(vector<T>& v) { return v.size() + &v.front(); }
template<typename T> inline T const * ptrend(vector<T> const & v) { return v.size() + &v.front(); }

struct matmul_bucket_data_s {
    /* repeat the fields from the public_ interface */
    struct matmul_public_s public_[1];
    /* now our private fields */
    abobj_t xab;
    size_t npack;
    size_t scrapsize;
    size_t bigscrapsize;
    abt * scrap;
    abt * bigscrap;
    vector<uint16_t> t16;       /* For small (dense) slices */
    vector<uint8_t> t8;         /* For large (sparse) slices */
    vector<unsigned int> aux;   /* Various descriptors -- fairly small */
    unsigned long mms1_ncoeffs;
    unsigned long mms2_ncoeffs;
    unsigned long fbi_ncoeffs;
    unsigned long fbd_ncoeffs;
    unsigned long asb_ncoeffs;
    clock_t mms1_time;
    clock_t mms2_time;
    clock_t fbi_time;
    clock_t fbd_time;
    clock_t asb_time;
};

void matmul_bucket_clear(struct matmul_bucket_data_s * mm)
{
    if (mm->scrap) abclear(mm->xab, mm->scrap, mm->scrapsize * HUGE_MPLEX_MAX);
    if (mm->bigscrap) abclear(mm->xab, mm->bigscrap, mm->bigscrapsize * HUGE_MPLEX_MAX);
    matmul_common_clear(mm->public_);
    // delete properly calls the destructor for members as well.
    delete mm;
}

static void mm_finish_init(struct matmul_bucket_data_s * mm);

struct matmul_bucket_data_s * matmul_bucket_init(abobj_ptr xx MAYBE_UNUSED, param_list pl, int optimized_direction)
{
    struct matmul_bucket_data_s * mm;
    mm = new matmul_bucket_data_s;
    memset(mm, 0, sizeof(struct matmul_bucket_data_s));
    abobj_init_set(mm->xab, xx);

    unsigned int npack = L1_CACHE_SIZE;
    if (pl) param_list_parse_uint(pl, "l1_cache_size", &npack);
    npack /= abbytes(mm->xab,1);
    mm->npack = npack;

    unsigned int scrapsize = L2_CACHE_SIZE/2;
    if (pl) param_list_parse_uint(pl, "l2_cache_size", &scrapsize);
    scrapsize /= abbytes(mm->xab,1);
    mm->scrapsize = scrapsize;

    unsigned int bigscrapsize = scrapsize; // (1 << 24)/abbytes(mm->xab,1));
    mm->bigscrapsize = bigscrapsize;

    int suggest = optimized_direction ^ MM_DIR0_PREFERS_TRANSP_MULT;
    mm->public_->store_transposed = suggest;
    if (pl) {
        param_list_parse_uint(pl, "mm_store_transposed",
                &mm->public_->store_transposed);
        if (mm->public_->store_transposed != (unsigned int) suggest) {
            fprintf(stderr, "Warning, mm_store_transposed"
                    " overrides suggested matrix storage ordering\n");
        }   
    }

    return mm;
}

/* This moves an element at the tail of a list with no copy, transferring
 * the ownership to the container argument */
template<typename T>
void transfer(list<T> * ctr, T * elem)
{
    ctr->push_back(T());
    swap(ctr->back(), *elem);
}

struct builder {
    uint32_t * data[2];
    uint32_t nrows_t;
    uint32_t ncols_t;
    vector<uint32_t *> colheads;
    vector<uint32_t *> colptrs;
    const char * rowname;
    const char * colname;
    /* Statistics on the prime parameter affecting performance of
     * subroutines. All are computed as (total sum, #samples) */
    uint32_t fbi_avg[2];
    uint32_t fbd_avg[2];
    uint32_t asb_avg[2]; /* For asb: reuse ratio for output cells */
};

void builder_init(builder * mb, const char * filename, struct matmul_bucket_data_s * mm)
{
    read_easy(filename,
            &mb->data[ mm->public_->store_transposed],
            &mb->data[!mm->public_->store_transposed],
            &mm->public_->dim[0],
            &mm->public_->dim[1]);

    mb->nrows_t = mm->public_->dim[ mm->public_->store_transposed];
    mb->ncols_t = mm->public_->dim[!mm->public_->store_transposed];

    mb->rowname = rowcol[ mm->public_->store_transposed];
    mb->colname = rowcol[!mm->public_->store_transposed];
    memset(mb->fbi_avg, 0, sizeof(mb->fbi_avg));
    memset(mb->fbd_avg, 0, sizeof(mb->fbd_avg));
    memset(mb->asb_avg, 0, sizeof(mb->asb_avg));
}

void builder_switch_to_vertical(builder * mb, unsigned int i)
{
    /*  column pointers */
    mb->colheads.resize(mb->ncols_t, NULL);
    mb->colptrs.resize(mb->ncols_t, NULL);
    {
        uint32_t * ptr = mb->data[1];
        for(uint32_t j = 0 ; j < mb->ncols_t ; j++) {
            uint32_t len = *ptr;
            uint32_t * v0 = ptr + 1;
            mb->colheads[j] = ptr;
            ptr += 1 + len;
            uint32_t * c = v0;
            for( ; c - v0 < (ptrdiff_t) len && *c < (uint32_t) i ; c++) ;
            mb->colptrs[j] = c;
        }
    }
}

/************************************/
/* Elementary (!) buliding routines */

/* there's room for a slice type where data is stored row-major. It would
 * make it possible to avoid most stores, for a nice performance gain. We
 * might even consider accumulating four rows at a time to four different
 * registers, so that the loads are compensated. The in-memory format
 * must account for that.
 */

#define SMALL_SLICES_I_BITS     12
#define SMALL_SLICES_DJ_BITS    (16 - SMALL_SLICES_I_BITS)

/* {{{ small slices */

struct small_slice_t {
    uint32_t i0;
    uint32_t i1;
    unsigned int ncoeffs;
    double dj_avg;
    uint32_t dj_max;
    int is_small2;

    typedef std::vector<std::pair<uint32_t, uint32_t> > L_t;
    typedef L_t::const_iterator Lci_t;
    typedef L_t::iterator Li_t;

    L_t L;
};

#define UPDATE_DIFFERENCE(i,di,mem) do { \
    di = mem - i; i = mem; mem = di;     \
} while (0)


int builder_do_small_slice(builder * mb, struct small_slice_t * S, uint32_t ** ptr, uint32_t i0, uint32_t i1)
{
    S->i0 = i0;
    S->i1 = i1;

    ASSERT_ALWAYS(i1-i0 <= (1 << SMALL_SLICES_I_BITS) );

    uint32_t * ptr0 = *ptr;
    /* We're doing a new slice */
    for(uint32_t i = i0 ; i < i1 ; i++) {
        for(unsigned int j = 0 ; j < *(*ptr) ; j++) {
            S->L.push_back(std::make_pair((*ptr)[1+j],i-i0));
        }
        (*ptr) += 1 + *(*ptr);
    }
    /* L is the list of (column index, row index) of all
     * coefficients in the current horizontal slice */
    std::sort(S->L.begin(), S->L.end());

    /* Convert all j indices to differences */
    S->dj_max = 0;
    S->ncoeffs = S->L.size();
    uint32_t j = 0;

    for(small_slice_t::Li_t lp = S->L.begin() ; lp != S->L.end() ; ++lp) {
        uint32_t dj = lp->first - j;
        j = lp->first;
        // UPDATE_DIFFERENCE(j,dj,lp->first);
        if (dj > S->dj_max) { S->dj_max = dj; }
    }

    S->dj_avg = mb->ncols_t / (double) S->ncoeffs;

    /* There are two possible reasons for this slice to be discarded.
     * Either the max dj is too large -- larger than 16 bits -- or the
     * average dj is too large -- larger than some guessed value.
     */
    int keep1 = S->dj_max < (1UL << 16) && S->dj_avg < DJ_CUTOFF1;
    int keep2 = S->dj_max < (1 << SMALL_SLICES_DJ_BITS) && S->dj_avg < DJ_CUTOFF1;
    S->is_small2 = keep2;

    if (!keep1 && !keep2) {
        *ptr = ptr0;
    }
    return keep1 || keep2;
}

/* }}} */

/* {{{ large slices */

struct large_slice_vblock_t {
    uint32_t j0;
    uint32_t j1;
    unsigned int n;             // number of coeffs.
    unsigned int pad;           // (half the) number of padding coeffs.
    vector<uint8_t> t8c;        // read: contribution to t8.
    vector<unsigned int> auxc;  // read: contribution to aux.
};

struct large_slice_t {
    uint32_t i0;
    uint32_t i1;
    unsigned int ncoeffs;
    double dj_avg;
    uint32_t dj_max;
    list<large_slice_vblock_t> vbl;
};

struct large_slice_raw_t {
    vector<uint8_t> ind[LSL_NBUCKETS_MAX];
    vector<uint8_t> main;
    vector<size_t> col_sizes;
    vector<size_t> pad_sizes;
};

void split_large_slice_in_vblocks(builder * mb, large_slice_t * L, large_slice_raw_t * R, unsigned int scrapsize)
{
    /* Now split into vslices */
    uint8_t * mp = ptrbegin(R->main);
    uint8_t * ip[LSL_NBUCKETS_MAX];
    for(int k = 0 ; k < LSL_NBUCKETS_MAX ; k++) {
        ip[k] = ptrbegin(R->ind[k]);
    }
    unsigned int vblocknum = 0;
    double vbl_ncols_variance = 0;
    for(uint32_t j = 0 ; j < mb->ncols_t ; vblocknum++) {
        uint32_t j0 = j;
        size_t n = 0;
        size_t np = 0;
        for( ; j < mb->ncols_t ; j++) {
            size_t dn = R->col_sizes[j];
            size_t dnp = R->pad_sizes[j];
            if (n + dn + 2 * (np + dnp) > scrapsize) {
                /* If dn doesn't fit, then we're in really, really big
                 * trouble ! */
                ASSERT_ALWAYS(n != 0);
                break;
            }
            n += dn;
            np += dnp;
        }
        uint32_t j1 = j;

        large_slice_vblock_t V;
        V.j0 = j0;
        V.j1 = j1;
        /* First the main block, then all the bucket blocks. These come
         * in order. */
        /* We have nf coeffs, counting padding. This makes 2nf bytes in
         * main, and nf bytes in the bucket blocks */
        V.n = n;
        V.pad = np;
        size_t nf = n + 2 * np;
        V.auxc.reserve(258);
        V.auxc.push_back(V.j1 - V.j0);
        V.t8c.resize(3 * nf);
        uint8_t * q = ptrbegin(V.t8c);
        /* XXX Notice that the reader process will consider as a j0 index
         * the farthest possible index from the previous vblock, which
         * means that we _must_ start with a zero in the main block here
         * no matter what -- even if with regard to the way data is set
         * up at the upper level (huge unique vblock), we would have had
         * a positive integer. This digression does not apply to the
         * first vblock. */
        memcpy(q, mp, 2 * nf * sizeof(uint8_t));
        if (vblocknum) {
            q[0] = 0;
        }
        unsigned int ind_sizes[LSL_NBUCKETS_MAX] = {0,};
        q += 2*nf;
        for(size_t k = 0 ; k < nf ; k++) {
            ind_sizes[mp[2 * k + 1]]++;
        }
        mp += 2 * nf;
        V.auxc.push_back(nf);
        for(int k = 0 ; k < LSL_NBUCKETS_MAX ; k++) {
            ASSERT(ip[k]+ind_sizes[k] <= ptrend(R->ind[k]));
            memcpy(q, ip[k], ind_sizes[k] * sizeof(uint8_t));
            q += ind_sizes[k];
            ip[k] += ind_sizes[k];
            V.auxc.push_back(ind_sizes[k]);

            mb->asb_avg[0] += ind_sizes[k];
            mb->asb_avg[1] += MIN(L->i1 - L->i0 - k * 256, 256);
        }
        ASSERT(q == ptrend(V.t8c));

        vbl_ncols_variance += ((double)(j1-j0)) * ((double)(j1-j0));
#if 0
        printf(" vbl%u", vblocknum);
        printf(": %ss %u..%u ", mb->colname, j0, j1);
        printf("; w %u", V.n);
        printf("; avg dj %.1f", (j1 - j0) / (double) V.n);
        if (V.pad) printf("; pad 6*%u", V.pad);
        printf("\n");
#endif
        transfer(&(L->vbl), &V);
    }
    double vbl_ncols_mean = mb->ncols_t;
    vbl_ncols_mean /= vblocknum;
    vbl_ncols_variance /= vblocknum;
    double vbl_ncols_sdev = sqrt(vbl_ncols_variance - vbl_ncols_mean * vbl_ncols_mean);
    printf(" %u vblocks, sdev/avg = %.2f\n", vblocknum, vbl_ncols_sdev / vbl_ncols_mean);
}

int builder_do_large_slice(builder * mb, struct large_slice_t * L, uint32_t i0, uint32_t i1, unsigned int scrapsize)
{
    L->i0 = i0;
    L->i1 = i1;
    L->ncoeffs = 0;
    L->dj_max = 0;

    large_slice_raw_t R[1];

    R->col_sizes.assign(mb->ncols_t, 0);
    R->pad_sizes.assign(mb->ncols_t, 0);

    uint32_t last_j = 0;

    /* First we create a huge unique vblock, and later on decide on
     * how to split it. */
    for(uint32_t j = 0 ; j < mb->ncols_t ; j++) {
        uint32_t len  = * mb->colheads[j];
        uint32_t * v0 = mb->colheads[j] + 1;
        uint32_t * c  = mb->colptrs[j];

        if (!(c - v0 < (ptrdiff_t) len && *c < i1))
            continue;

        uint32_t diff = j - last_j;
        if (diff > L->dj_max) {
            L->dj_max = diff;
        }
        for( ; diff > 255 ; ) {
            R->main.push_back(255);
            R->main.push_back(0);
            R->ind[0].push_back(0);
            R->main.push_back(0);
            R->main.push_back(0);
            R->ind[0].push_back(0);
            diff -= 255;
            R->pad_sizes[j] += 1;
        }

        /* Only the first entry in the column has a non-zero dj. So place
         * 0 instead inside the loop, and fix it in the end. */
        R->main.push_back(diff);
        last_j = j;
        for( ; c - v0 < (ptrdiff_t) len && *c < i1 ; c++) {
            uint32_t i = *c - i0;
            uint8_t w = i / 256;
            R->main.push_back(w);
            R->main.push_back(0);
            R->ind[w].push_back((uint8_t) i);
            L->ncoeffs ++;
        }
        R->main.pop_back();
        R->col_sizes[j] = c - mb->colptrs[j];
        mb->colptrs[j] = c;
    }

    /* The average dj is quite easy */
    L->dj_avg = mb->ncols_t / (double) L->ncoeffs;

    printf(" w=%u, avg dj=%.1f, max dj=%u, bucket hit=1/%.1f",
            L->ncoeffs, L->dj_avg, L->dj_max, LSL_NBUCKETS_MAX * L->dj_avg);

    if (L->dj_avg > DJ_CUTOFF2) {
        printf("-> too sparse");
        if (mb->nrows_t - i0 < HUGE_MPLEX_MIN * LSL_NBUCKETS_MAX * 256) {
            printf("; kept because of short tail\n");
        } else {
            /* We won't keep this slice */
            printf("\n");
            for(uint32_t j = 0 ; j < mb->ncols_t ; j++) {
                mb->colptrs[j] -= R->col_sizes[j];
            }
            return 0;
        }
    }

    printf("\n");

    split_large_slice_in_vblocks(mb, L, R, scrapsize);

    return 1;
}

/* }}} */

/* {{{ huge slices */

struct huge_slice_t {
    uint32_t i0;
    uint32_t i1;
    unsigned int ncoeffs;
    double dj_avg;
    uint32_t dj_max;
    unsigned int nlarge;
    list<large_slice_vblock_t> vbl;
};

struct huge_subslice_raw_t {
    vector<uint8_t> ind[LSL_NBUCKETS_MAX];
    vector<uint8_t> main;
};

struct huge_subslice_ptrblock_t {
    uint8_t * ip[LSL_NBUCKETS_MAX];
    uint8_t * mp;
};

struct huge_slice_raw_t {
    vector<uint8_t> super;
    vector<huge_subslice_raw_t> subs;
    vector<size_t> col_sizes;
    vector<size_t> pad_sizes;
};

void split_huge_slice_in_vblocks(builder * mb, huge_slice_t * H, huge_slice_raw_t * R, unsigned int bigscrapsize)/*{{{*/
{
    /* Now split into vslices */
    unsigned int lsize = iceildiv(H->i1 - H->i0, H->nlarge);
    uint8_t * sp = ptrbegin(R->super);
    vector<huge_subslice_ptrblock_t> ptrs(R->subs.size());
    for(unsigned int i = 0 ; i < R->subs.size() ; i++) {
        ptrs[i].mp = ptrbegin(R->subs[i].main);
        for(int k = 0 ; k < LSL_NBUCKETS_MAX ; k++) {
            ptrs[i].ip[k] = ptrbegin(R->subs[i].ind[k]);
        }
    }
    double vbl_ncols_variance = 0;
    unsigned int vblocknum = 0;
    for(uint32_t j = 0 ; j < mb->ncols_t ; vblocknum++) {
        uint32_t j0 = j;
#if 0
        /* First way to advance to the next j -- compare the total number
         * of coeffs to the bigscrapsize argument */
        size_t n = 0;
        size_t np = 0;
        for( ; j < mb->ncols_t ; j++) {
            size_t dn = R->col_sizes[j];
            size_t dnp = R->pad_sizes[j];
            if (n + dn + 2 * (np + dnp) > bigscrapsize) {
                /* If dn doesn't fit, then we're in really, really big
                 * trouble ! */
                ASSERT_ALWAYS(n != 0);
                break;
            }
            n += dn;
            np += dnp;
        }
        uint32_t j1 = j;
#else
        /* However, there's another way. We might as well arrange so that
         * no mega bucket exceeds the L2 cache. */
        uint32_t j1 = j;
        uint8_t * sp0 = sp;
        uint8_t * spc = sp0;
        /* XXX Notice that the reader process will consider as a j0 index
         * the farthest possible index from the previous vblock, which
         * means that we _must_ start with a zero in the main block here
         * no matter what -- even if with regard to the way data is set
         * up at the upper level (huge unique vblock), we would have had
         * a positive integer. This digression does not apply to the
         * first vblock. */
        if (vblocknum) {
            sp0[0]=0;
        }
        {
            vector<unsigned int> Lsizes(H->nlarge, 0);
            for( ; j < mb->ncols_t && sp != ptrend(R->super) ; ) {
                unsigned int dj = *sp;
                j += dj;
                if (dj) {
                    /* beginning of a column. */
                    spc = sp;
                    j1 = j;
                }
                sp++;
                unsigned int k = *sp++;
                Lsizes[k]++;
                // use the bigscrapsize argument in both cases, even
                // though it doesn't mean the same.
                if (Lsizes[k] > bigscrapsize) {
                    Lsizes[k]--;
                    sp = spc;
                    break;
                }
            }
            if (sp == ptrend(R->super)) {
                spc = sp;
                j1 = mb->ncols_t;
            }
        }
        /* TODO: abandon the n/np distinction. It's useful for debugging
         * at best, and it would not hurt to have it displayed only at the
         * upper level.
         */
        size_t n = 0;
        size_t np = 0;
        for(j = j0 ; j < j1 ; j++) {
            n += R->col_sizes[j];
            np += R->pad_sizes[j];
        }
        ASSERT_ALWAYS((n+2*np)*2 == (size_t) (spc - sp0));
        sp = sp0;
#endif

        /* That's the same type, although it contains more stuff. */
        large_slice_vblock_t V;
        V.j0 = j0;
        V.j1 = j1;
        V.auxc.reserve(1 + H->nlarge * (1 + LSL_NBUCKETS_MAX));
        V.auxc.push_back(V.j1 - V.j0);
        /* First the super block, and then for each contained large
         * slice, the main block, then all the bucket blocks. These come
         * in order. */
        /* We have nf coeffs, counting padding. This makes 2nf bytes in
         * super, nf in mains, and nf in the bucket blocks */
        V.n = n;
        V.pad = np;
        size_t nf = n + 2 * np;
        V.t8c.resize(4 * nf);
        uint8_t * q = ptrbegin(V.t8c);
        memcpy(q, sp, 2 * nf * sizeof(uint8_t));
        /* TODO: If the solution above is to be kept, then of course we'd
         * rather reuse this thing above, since it has already been
         * computed.
         */
        vector<unsigned int> L_sizes(H->nlarge, 0);
        q += 2*nf;
        for(size_t k = 0 ; k < nf ; k++) {
            ASSERT(sp[2 * k + 1] < H->nlarge);
            L_sizes[sp[2 * k + 1]]++;
        }
        sp += 2 * nf;
        V.auxc.push_back(nf);
        for(unsigned int l = 0 ; l < H->nlarge ; l++) {
            V.auxc.push_back(L_sizes[l]);
        }
        for(unsigned int l = 0 ; l < H->nlarge ; l++) {
            unsigned int ind_sizes[LSL_NBUCKETS_MAX] = {0,};
            memcpy(q, ptrs[l].mp, L_sizes[l] * sizeof(uint8_t));
            q += L_sizes[l];
            for(unsigned int k = 0 ; k < L_sizes[l] ; k++) {
                ind_sizes[ptrs[l].mp[k]]++;
            }
            ptrs[l].mp += L_sizes[l];

            unsigned int i_size = MIN(H->i1 - H->i0 - l * lsize, lsize);

            for(int k = 0 ; k < LSL_NBUCKETS_MAX ; k++) {
                ASSERT(ptrs[l].ip[k]-(ptrbegin(R->subs[l].ind[k]))+ind_sizes[k] <= (ptrdiff_t) R->subs[l].ind[k].size());
                memcpy(q, ptrs[l].ip[k], ind_sizes[k] * sizeof(uint8_t));
                q += ind_sizes[k];
                ptrs[l].ip[k] += ind_sizes[k];
                V.auxc.push_back(ind_sizes[k]);
                mb->asb_avg[0] += ind_sizes[k];
                mb->asb_avg[1] += MIN(i_size - k * 256, 256);
            }
        }
        ASSERT(q == ptrend(V.t8c));

        vbl_ncols_variance += ((double)(j1-j0)) * ((double)(j1-j0));

        /* TODO: have some sort of verbosity level */
#if 0
        printf(" vbl%u", vblocknum);
        printf(": %ss %u..%u ", mb->colname, j0, j1);
        printf("; w %u", V.n);
        printf("; avg dj %.1f", (j1 - j0) / (double) V.n);
        if (V.pad) printf("; pad 6*%u", V.pad);
        printf("\n");
#endif

        transfer(&(H->vbl), &V);
    }
    double vbl_ncols_mean = mb->ncols_t;
    vbl_ncols_mean /= vblocknum;
    vbl_ncols_variance /= vblocknum;
    double vbl_ncols_sdev = sqrt(vbl_ncols_variance - vbl_ncols_mean * vbl_ncols_mean);
    printf(" %u vblocks, sdev/avg = %.2f\n", vblocknum, vbl_ncols_sdev / vbl_ncols_mean);
}/*}}}*/

int builder_do_huge_slice(builder * mb, struct huge_slice_t * H, uint32_t i0, uint32_t i1, unsigned int bigscrapsize)
{
    H->i0 = i0;
    H->i1 = i1;
    H->ncoeffs = 0;
    H->dj_max = 0;

    huge_slice_raw_t R[1];

    R->col_sizes.assign(mb->ncols_t, 0);
    R->pad_sizes.assign(mb->ncols_t, 0);

    /* How many large slices in this huge slice ? */

    H->nlarge = iceildiv(i1 - i0, LSL_NBUCKETS_MAX * 256);
    ASSERT(H->nlarge >= HUGE_MPLEX_MIN);
    ASSERT(H->nlarge  < HUGE_MPLEX_MAX);
    unsigned int lsize = iceildiv(i1 - i0, H->nlarge);
    ASSERT(lsize <= LSL_NBUCKETS_MAX * 256);
    ASSERT(lsize * H->nlarge >= i1 - i0);
    R->subs.assign(H->nlarge, huge_subslice_raw_t());

    printf(" (%u*%u) ", H->nlarge, lsize);
    uint32_t last_j = 0;

    /* First we create a huge unique vblock, and later on decide on
     * how to split it. */
    uint32_t next_dot = mb->ncols_t / H->nlarge;
    for(uint32_t j = 0 ; j < mb->ncols_t ; j++) {
        uint32_t len  = * mb->colheads[j];
        uint32_t * v0 = mb->colheads[j] + 1;
        uint32_t * c  = mb->colptrs[j];

        if (j > next_dot) {
            printf(".");
            fflush(stdout);
            next_dot += mb->ncols_t / H->nlarge;
        }

        if (!(c - v0 < (ptrdiff_t) len && *c < i1))
            continue;

        uint32_t diff = j - last_j;
        if (diff > H->dj_max) {
            H->dj_max = diff;
        }
        for( ; diff > 255 ; ) {
            R->super.push_back(255);    // dj
            R->super.push_back(0);      // k
            R->subs[0].main.push_back(0);
            R->subs[0].ind[0].push_back(0);
            R->super.push_back(0);      // dj
            R->super.push_back(0);      // k
            R->subs[0].main.push_back(0);
            R->subs[0].ind[0].push_back(0);
            diff -= 255;
            R->pad_sizes[j] += 1;
        }

        /* Only the first entry in the column has a non-zero dj. So place
         * 0 instead inside the loop, and fix it in the end. */
        R->super.push_back(diff);
        last_j = j;
        for( ; c - v0 < (ptrdiff_t) len && *c < i1 ; c++) {
            uint32_t i = *c - i0;
            uint8_t w = i / lsize;
            uint8_t wq = (i % lsize) / 256;
            uint8_t wr = i % lsize;
            R->super.push_back(w);
            R->super.push_back(0);
            R->subs[w].main.push_back(wq);
            R->subs[w].ind[wq].push_back(wr);
            H->ncoeffs ++;
        }
        R->super.pop_back();
        R->col_sizes[j] = c - mb->colptrs[j];
        mb->colptrs[j] = c;
    }
    printf("\n");

    /* The average dj is quite easy */
    H->dj_avg = mb->ncols_t / (double) H->ncoeffs;

    printf(" w=%u, avg dj=%.1f, max dj=%u, bucket block hit=1/%.1f\n",
            H->ncoeffs, H->dj_avg, H->dj_max,
            H->nlarge * H->dj_avg);

    split_huge_slice_in_vblocks(mb, H, R, bigscrapsize);

    return 1;
}

/* }}} */

/******************************************/
/* Iteratively call the building routines */

/* {{{ small slices */
void builder_do_all_small_slices(builder * mb, uint32_t * p_i0, list<small_slice_t> * Sq, unsigned int npack)
{
    /* npack is a guess for the expected size of small slices ; they are
     * arranged later to all have approximately equal size.
     */
    unsigned int nslices = iceildiv(mb->nrows_t, npack);
    unsigned int s;
    uint32_t * ptr = mb->data[0];
    for(s = 0 ; s < nslices ; s++) {
        uint32_t i0 =  s    * mb->nrows_t / nslices;
        uint32_t i1 = (s+1) * mb->nrows_t / nslices;

        small_slice_t S;

        printf("Ssl%u: %ss %u+%u...", s, mb->rowname, i0, i1-i0);
        fflush(stdout);

        int keep = builder_do_small_slice(mb, &S, &ptr, i0, i1);

        printf(" w %" PRIu32 " ; avg dj %.1f ; max dj %u%s\n",
                S.ncoeffs, S.dj_avg, S.dj_max,
                S.is_small2 ? " [packed]" : "");
        fflush(stdout);

        if (!keep) {
            printf("Switching to large slices. Ssl%u to be redone\n", s);
            break;
        }
        transfer(Sq, &S);
    }
    *p_i0 = s * mb->nrows_t / nslices;
}
/* }}} */

/* {{{ large slices */
void builder_do_all_large_slices(builder * mb, uint32_t * p_i0, list<large_slice_t> * Lq, unsigned int scrapsize)
{
    unsigned int rem_nrows = mb->nrows_t - *p_i0;
    unsigned int nlarge_slices = iceildiv(rem_nrows,LSL_NBUCKETS_MAX * 256);
    uint32_t done = 0;
    for(unsigned int s = 0 ; s < nlarge_slices ; s++) {
        large_slice_t L[1];
        uint32_t i0 = * p_i0 +  s      * rem_nrows / nlarge_slices;
        uint32_t i1 = * p_i0 + (s + 1) * rem_nrows / nlarge_slices;

        printf("Lsl%u %ss %u+%u", s, mb->rowname, i0, i1-i0);
        fflush(stdout);

        int keep = builder_do_large_slice(mb, L, i0, i1, scrapsize);

        if (!keep) {
            printf("Switching to huge slices. Lsl%u to be redone\n", s);
            break;
        }

        transfer(Lq, L);
        done += i1 - i0;
    }
    *p_i0 += done;
}
/* }}} */

/* {{{ huge slices */
void builder_do_all_huge_slices(builder * mb, uint32_t * p_i0, list<huge_slice_t> * Hq, unsigned int bigscrapsize)
{
    unsigned int rem_nrows = mb->nrows_t - *p_i0;
    unsigned int nhuge_slices = iceildiv(rem_nrows, HUGE_MPLEX_MAX * LSL_NBUCKETS_MAX * 256);
    uint32_t done = 0;
    for(unsigned int s = 0 ; s < nhuge_slices ; s++) {
        huge_slice_t H[1];
        uint32_t i0 = * p_i0 +  s      * rem_nrows / nhuge_slices;
        uint32_t i1 = * p_i0 + (s + 1) * rem_nrows / nhuge_slices;

        printf("Hsl%u %ss %u+%u", s, mb->rowname, i0, i1-i0);
        fflush(stdout);
        builder_do_huge_slice(mb, H, i0, i1, bigscrapsize);
        transfer(Hq, H);
        done += i1 - i0;
    }
    *p_i0 += done;
}
/* }}} */

/*************************************************************************/
/* Pushing slices to mm ; all these routines clear the given slice stack */

/* {{{ small slices */
void builder_push_small_slices(struct matmul_bucket_data_s * mm, list<small_slice_t> * Sq)
{
    mm->t16.push_back(Sq->size());
    mm->t16.push_back(0);       // for alignment.

    printf("Flushing %zu small slices\n", Sq->size());

    for(unsigned int s = 0 ; !Sq->empty() ; s++, Sq->pop_front()) {
        small_slice_t * S = &Sq->front();

        mm->t16.push_back(S->i1-S->i0);
        /* Older versions had a padding zero here */
        mm->t16.push_back(S->is_small2);

        if (S->is_small2) {
            mm->t16.push_back(S->L.size() & 0xffff);
            mm->t16.push_back(S->L.size() >> 16);
            /* small2 slices are smaller, and faster -- but more
             * restrictive.  */

            /* Because a small2 slice has data in 16bits chunks, we must
             * ensure proper alignment in the end. It's easiest right here. */
            for(unsigned int pad = (2 - 1) & - S->L.size() ; pad-- ; ) {
                mm->t16.push_back(0);
            }

            unsigned int j = 0;
            small_slice_t::Lci_t lp;
            for(lp = S->L.begin() ; lp != S->L.end() ; ++lp) {
                unsigned int dj = lp->first - j;
                j = lp->first;
                ASSERT(dj < (1 << SMALL_SLICES_DJ_BITS) );
                ASSERT(lp->second < (1 << SMALL_SLICES_I_BITS) );
                mm->t16.push_back((dj << SMALL_SLICES_I_BITS) | lp->second);
            }
        } else {
            /* How many vertical parts -- this merely has to do with index
             * wraparound, since we're processing data column-major anyway.
             */
            unsigned int nvstrips = 0;
            unsigned int j0 = 0;
            small_slice_t::Lci_t lp;
            for(lp = S->L.begin() ; lp != S->L.end() ; nvstrips++) {
                unsigned long stripsize_pos = mm->t16.size();
                mm->t16.push_back(0);
                mm->t16.push_back(0);
                unsigned int j = j0;
                unsigned int stripsize = 0;
                for( ; lp != S->L.end() ; ++lp, stripsize++) {
                    if (lp->first - j0 >= (1UL << 16)) 
                        break;
                    // mm->t16.push_back(lp->first - j);
                    mm->t16.push_back(lp->first);
                    j = lp->first;
                    mm->t16.push_back(lp->second);
                }
                j0 += (1UL << 16);
                mm->t16[stripsize_pos] = stripsize & 0xffff;
                mm->t16[stripsize_pos + 1] = stripsize >> 16;
            }
            mm->aux.push_back(nvstrips);
        }
        mm->public_->ncoeffs += S->ncoeffs;
    }
}
/* }}} */

/* {{{ large slices */
void builder_push_large_slices(struct matmul_bucket_data_s * mm, list<large_slice_t> * Lq)
{
    printf("Flushing %zu large slices\n", Lq->size());

    mm->aux.push_back(Lq->size());

    for(unsigned int l = 0 ; !Lq->empty() ; l++, Lq->pop_front()) {
        large_slice_t * L = &Lq->front();
        mm->aux.push_back(L->i1 - L->i0);
        for( ; ! L->vbl.empty() ; L->vbl.pop_front()) {
            large_slice_vblock_t & V(L->vbl.front());
            mm->aux.insert(mm->aux.end(), V.auxc.begin(), V.auxc.end());
            unsigned t8_size =  V.t8c.size();
            mm->t8.insert(mm->t8.end(), V.t8c.begin(), V.t8c.end());
            // large slices, but not huge slices, may put an odd number
            // of coefficients in t8
            if (t8_size & 1) { mm->t8.push_back(0); }
        }
        mm->public_->ncoeffs += L->ncoeffs;
    }
}
/* }}} */

/* {{{ huge slices */
void builder_push_huge_slices(struct matmul_bucket_data_s * mm, list<huge_slice_t> * Hq)
{
    printf("Flushing %zu huge slices\n", Hq->size());

    mm->aux.push_back(Hq->size());
    for(unsigned int l = 0 ; !Hq->empty() ; l++, Hq->pop_front()) {
        huge_slice_t * H = &Hq->front();
        mm->aux.push_back(H->i1 - H->i0);
        mm->aux.push_back(H->nlarge);
        for( ; ! H->vbl.empty() ; H->vbl.pop_front()) {
            large_slice_vblock_t & V(H->vbl.front());
            mm->aux.insert(mm->aux.end(), V.auxc.begin(), V.auxc.end());
            unsigned t8_size =  V.t8c.size();
            ASSERT_ALWAYS((t8_size & 1) == 0);
            mm->t8.insert(mm->t8.end(), V.t8c.begin(), V.t8c.end());
        }
        mm->public_->ncoeffs += H->ncoeffs;
    }
}
/* }}} */


void matmul_bucket_build_cache(struct matmul_bucket_data_s * mm)
{
    builder mb[1];
    builder_init(mb, mm->public_->filename, mm);


    uint32_t main_i0 = 0;

    list<small_slice_t> Sq[1];
    builder_do_all_small_slices(mb, &main_i0, Sq, mm->npack);
    builder_push_small_slices(mm, Sq);

    free(mb->data[0]);
    mb->data[0] = NULL;

    builder_switch_to_vertical(mb, main_i0);

    list<large_slice_t> Lq[1];
    builder_do_all_large_slices(mb, &main_i0, Lq, mm->scrapsize);
    builder_push_large_slices(mm, Lq);

    /* Don't deallocate right now, since data is still to be used for
     * huge slices.  */

    list<huge_slice_t> Hq[1];
    builder_do_all_huge_slices(mb, &main_i0, Hq, mm->bigscrapsize);
    builder_push_huge_slices(mm, Hq);

    printf("avg for asb: %.1f\n", mb->asb_avg[0] / (double) mb->asb_avg[1]);

    /* done, at last ! */
    free(mb->data[1]);
    mb->data[1] = NULL;

    mm_finish_init(mm);
}

int matmul_bucket_reload_cache(struct matmul_bucket_data_s * mm)/* {{{ */
{
    FILE * f;

    f = matmul_common_reload_cache_fopen(abbytes(mm->xab,1), mm->public_, MM_MAGIC);
    if (f == NULL) { return 0; }

    MATMUL_COMMON_READ_ONE32(mm->scrapsize, f);
    size_t n16, n8, naux;
    MATMUL_COMMON_READ_ONE32(n16, f);
    MATMUL_COMMON_READ_ONE32(n8, f);
    MATMUL_COMMON_READ_ONE32(naux, f);
    mm->t16.resize(n16);
    mm->t8.resize(n8);
    mm->aux.resize(naux);
    MATMUL_COMMON_READ_MANY16(ptrbegin(mm->t16), n16, f);
    MATMUL_COMMON_READ_MANY8(ptrbegin(mm->t8), n8, f);
    MATMUL_COMMON_READ_MANY32(ptrbegin(mm->aux), naux, f);

    fclose(f);

    mm_finish_init(mm);

    return 1;
}/*}}}*/

void matmul_bucket_save_cache(struct matmul_bucket_data_s * mm)/*{{{*/
{
    FILE * f;

    f = matmul_common_save_cache_fopen(abbytes(mm->xab,1), mm->public_, MM_MAGIC);
    MATMUL_COMMON_WRITE_ONE32(mm->scrapsize, f);
    size_t n16 = mm->t16.size();
    size_t n8 = mm->t8.size();
    size_t naux = mm->aux.size();
    MATMUL_COMMON_WRITE_ONE32(n16, f);
    MATMUL_COMMON_WRITE_ONE32(n8, f);
    MATMUL_COMMON_WRITE_ONE32(naux, f);
    MATMUL_COMMON_WRITE_MANY16(ptrbegin(mm->t16), n16, f);
    MATMUL_COMMON_WRITE_MANY8(ptrbegin(mm->t8), n8, f);
    MATMUL_COMMON_WRITE_MANY32(ptrbegin(mm->aux), naux, f);

    fclose(f);
}/*}}}*/

static inline uint32_t read32(uint16_t const * & q) /* {{{ */
{
    uint32_t res;
    res = *q++;
    res |= ((uint32_t) *q++) << 16;
    return res;
} /* }}} */

struct pos_desc {
    const uint16_t * q16;
    const uint8_t * q8;
    const unsigned int * ql;
    uint32_t i;
    uint32_t nrows_t;
    uint32_t ncols_t;
};

static inline void matmul_bucket_mul_small(struct matmul_bucket_data_s * mm, abt * dst, abt const * src, int d, struct pos_desc * pos)
{
    abobj_ptr x = mm->xab;

    ASM_COMMENT("multiplication code -- small (dense) slices"); /* {{{ */
    uint16_t nhstrips = *pos->q16++;
    pos->q16++;        // alignment.
    if (d == !mm->public_->store_transposed) {
        for(uint16_t s = 0 ; s < nhstrips ; s++) {
            uint32_t nrows_packed = *pos->q16++;
            int is_small2 = *pos->q16++;
            ASSERT(pos->i + nrows_packed <= pos->nrows_t);
            abt * where = dst + aboffset(x, pos->i);
            abzero(x, where, nrows_packed);
            ASM_COMMENT("critical loop");
            abt const * from = src;
            if (is_small2) {
                mm->mms2_time -= clock();
#if defined(ENABLE_ASM)
                pos->q16 = matmul_sub_small2(x, where, from, (uint16_t const *) pos->q16);
#else
                uint32_t ncoeffs_slice = read32(pos->q16);
                unsigned int pad = (2 - 1) & - ncoeffs_slice;
                pos->q16 += pad;
                uint32_t j = 0;
                for(uint32_t c = 0 ; c < ncoeffs_slice ; c++) {
                    uint16_t h = *pos->q16++;
                    j += h >> SMALL_SLICES_I_BITS;
                    ASSERT(j < pos->ncols_t);
                    uint32_t di = h & ((1UL << SMALL_SLICES_I_BITS) - 1);
                    abadd(x, where + aboffset(x, di), from + aboffset(x, j));
                }
#endif
                mm->mms2_time += clock();
            } else {
                unsigned int nv = *pos->ql++;
                for(unsigned int v = 0 ; v < nv ; v++) {
                    mm->mms1_time -= clock();
#if defined(ENABLE_ASM)
                    pos->q16 = matmul_sub_small1(x, where, from, (uint16_t const *) pos->q16);
#else
                    uint32_t ncoeffs_slice = read32(pos->q16);
                    uint32_t j = 0;
                    for(uint32_t c = 0 ; c < ncoeffs_slice ; c++) {
                        j = *pos->q16++;
                        ASSERT(j < pos->ncols_t);
                        uint32_t di = *pos->q16++;
                        abadd(x, where + aboffset(x, di), from + aboffset(x, j));
                    }
#endif
                    mm->mms1_time += clock();
                    from += aboffset(x, 1UL << 16);
                }
            }
            pos->i += nrows_packed;
            ASM_COMMENT("end of critical loop");
        }
    } else {
        /* d == mm->public_->store_transposed */
        abt * dst0 = dst;
        for(uint16_t s = 0 ; s < nhstrips ; s++) {
            dst = dst0;
            uint32_t nrows_packed = *pos->q16++;
            int is_small2 = *pos->q16++;
            if (is_small2) {
                mm->mms2_time -= clock();
                uint32_t ncoeffs_slice = read32(pos->q16);
                unsigned int pad = (2 - 1) & - ncoeffs_slice;
                pos->q16 += pad;
                uint32_t j = 0;
                for(uint32_t c = 0 ; c < ncoeffs_slice ; c++) {
                    uint16_t h = *pos->q16++;
                    j += h >> SMALL_SLICES_I_BITS;
                    ASSERT(j < pos->ncols_t);
                    uint32_t di = h & ((1UL << SMALL_SLICES_I_BITS) - 1);
                    abadd(x, dst + aboffset(x, j), src + aboffset(x, pos->i + di));
                }
                mm->mms2_time += clock();
            } else {
                unsigned int nv = *pos->ql++;
                for(unsigned int v = 0 ; v < nv ; v++) {
                    mm->mms1_time -= clock();
                    uint32_t ncoeffs_slice = read32(pos->q16);
                    uint32_t j = 0;
                    for(uint32_t c = 0 ; c < ncoeffs_slice ; c++) {
                        j = *pos->q16++;
                        ASSERT(j < (1UL << 16));
                        uint32_t di = *pos->q16++;
                        abadd(x, dst + aboffset(x, j), src + aboffset(x, pos->i + di));
                    }
                    mm->mms1_time += clock();
                    dst += aboffset(x, 1UL << 16);
                }
            }
            pos->i += nrows_packed;
        }
    }
    ASM_COMMENT("end of small(dense) slices"); /* }}} */
}

static inline void prepare_buckets_and_fences(abobj_ptr x, abt ** b, abt ** f, abt * z, const unsigned int * ql, unsigned int n)
{
    for(unsigned int k = 0 ; k < n ; k++) {
        b[k] = z;
        f[k] = b[k] + ql[k];
        z += aboffset(x, ql[k]);
    }
}

static inline void prepare_buckets(abobj_ptr x, abt ** b, abt * z, const unsigned int * ql, unsigned int n)
{
    for(unsigned int k = 0 ; k < n ; k++) {
        b[k] = z;
        z += aboffset(x, ql[k]);
    }
}

#ifndef ENABLE_ASM
static inline void fill_buckets_indirect(abobj_ptr x, abt ** sb, const abt * z, const uint8_t * q, unsigned int n)
{
    /* Dispatch data found in z[0]...z[f(n-1)] such that z[f(i)] is in
     * array pointed to by sb[q[2*i+1]]. The function f(i) is given by
     * the sum q[0]+q[2]+...+q[2*(i-1)]. Exactly 2n coefficients are
     * expected in q[] All the sb[] pointers are increased */
    for(unsigned int c = 0 ; c < n ; c++) {
        z += aboffset(x, *q);
        q++;
        abcopy(x, sb[*q], z, 1);
        sb[*q]+= aboffset(x, 1);
        q++;
    }
}
#endif

static inline void
unfill_buckets_indirect(abobj_ptr x, abt ** sb, abt * z, const uint8_t * q, unsigned int n)
{
    /* Does the converse of the above */
    for(unsigned int c = 0 ; c < n ; c++) {
        z += aboffset(x, *q);
        q++;
        abadd(x, z, sb[*q]);
        sb[*q]+= aboffset(x, 1);
        q++;
    }
}

static void apply_small_buckets(abobj_ptr x, abt * dst, const abt * z, const uint8_t * q, const unsigned int * ql) __attribute__((__noinline__));
static void apply_small_buckets(abobj_ptr x, abt * dst, const abt * z, const uint8_t * q, const unsigned int * ql)
{
    /* This ``applies'' the LSL_NBUCKETS_MAX small buckets whose
     * respective lengths are given by ql[0] to ql[LSL_NBUCKETS_MAX-1].
     *
     * For ql[k] <= i < ql[k+1], z[i] is added to dst[k*256+qk[i]], with
     * qk = q + ql[0] + ... + ql[k-1].
     */
    for(int k = 0 ; k < LSL_NBUCKETS_MAX ; k++) {
        unsigned int l = ql[k];
        for( ; l-- ; ) {
            /* For padding coeffs, the assertion can fail if
             * we choose a row not equal to (0,0) -- first in
             * the first bucket.
             */
            abadd(x, dst + aboffset(x, *q), z);
            z+= aboffset(x, 1);
            q++;
        }
        dst += aboffset(x, 256);
    }
}

static inline void unapply_small_buckets(abobj_ptr x, const abt * src, abt * z, const uint8_t * q, const unsigned int * ql)
{
    /* converse of the above */
    for(int k = 0 ; k < LSL_NBUCKETS_MAX ; k++) {
        unsigned int l = ql[k];
        for( ; l-- ; ) {
            /* For padding coeffs, the assertion can fail if
             * we choose a row not equal to (0,0) -- first in
             * the first bucket.
             */
            abcopy(x, z, src + aboffset(x, *q), 1);
            z += aboffset(x, 1);
            q++;
        }
        src += aboffset(x, 256);
    }
}

#ifndef ENABLE_ASM
static void
fill_buckets_direct(abobj_ptr x, abt ** sb, const abt * z, const uint8_t * q, unsigned int n)
{
    /* Dispatch data found in z[0]...z[n] such that z[i] is in array
     * pointed to by sb[q[i]]. Exactly n coefficients are expected. All
     * the sb[] pointers are increased */
    for(unsigned int c = 0 ; c < n ; c++) {
        abcopy(x, sb[q[c]], z, 1);
        sb[q[c]]+= aboffset(x, 1);
        z += aboffset(x, 1);
    }
}
#endif

static inline void
unfill_buckets_direct(abobj_ptr x, abt ** sb, abt * z, const uint8_t * q, unsigned int n)
{
    /* Does the converse of the above */
    for(unsigned int c = 0 ; c < n ; c++) {
        abcopy(x, z, sb[q[c]], 1);
        sb[q[c]]+= aboffset(x, 1);
        z += aboffset(x, 1);
    }
}

static inline void matmul_bucket_mul_large(struct matmul_bucket_data_s * mm, abt * dst, abt const * src, int d, struct pos_desc * pos)
{
    abobj_ptr x = mm->xab;

    abt * scrap = mm->scrap;
    ASM_COMMENT("multiplication code -- large (sparse) slices"); /* {{{ */
    unsigned int nlarge = *pos->ql++;
    if (d == !mm->public_->store_transposed) {
        for(unsigned int s = 0 ; s < nlarge ; s++) {
            uint32_t di = *pos->ql++;
            ASSERT(pos->i + di <= pos->nrows_t);
            ASSERT(di <= LSL_NBUCKETS_MAX * 256);
            uint32_t j = 0;
            abzero(x, dst + aboffset(x, pos->i), di);
            for( ; j < pos->ncols_t ; ) {
                uint32_t j1 = j + *pos->ql++;
                uint32_t n = *pos->ql++;
                abt * bucket[LSL_NBUCKETS_MAX];
                abt const * inp = src + aboffset(x, j);
                abt * outp = dst + aboffset(x, pos->i);
                prepare_buckets(x,bucket,scrap,pos->ql,LSL_NBUCKETS_MAX);
                mm->fbi_time -= clock();
                ASSERT_ALWAYS((((unsigned long)pos->q8)&1)==0);
                fill_buckets_indirect(x, bucket, inp, pos->q8, n);
                mm->fbi_time += clock();
                mm->asb_time -= clock();
                apply_small_buckets(x, outp, scrap, pos->q8+2*n, pos->ql);
                mm->asb_time += clock();
                pos->q8 += 3*n;
                // fix alignment !
                pos->q8 += n & 1;
                pos->ql += LSL_NBUCKETS_MAX;
                j = j1;
            }
            pos->i += di;
        }
    } else {
        for(unsigned int s = 0 ; s < nlarge ; s++) {
            uint32_t di = *pos->ql++;
            ASSERT(pos->i + di <= pos->nrows_t);
            ASSERT(di <= LSL_NBUCKETS_MAX * 256);
            uint32_t j = 0;
            for( ; j < pos->ncols_t ; ) {
                uint32_t j1 = j + *pos->ql++;
                uint32_t n = *pos->ql++;
                abt * bucket[LSL_NBUCKETS_MAX];
                abt * outp = dst + aboffset(x, j);
                abt const *  inp = src + aboffset(x, pos->i);
                prepare_buckets(x,bucket,scrap,pos->ql,LSL_NBUCKETS_MAX);
                mm->asb_time -= clock();
                unapply_small_buckets(x, inp, scrap, pos->q8+2*n, pos->ql);
                mm->asb_time += clock();
                mm->fbi_time -= clock();
                ASSERT_ALWAYS((((unsigned long)pos->q8)&1)==0);
                unfill_buckets_indirect(x, bucket, outp, pos->q8, n);
                mm->fbi_time += clock();
                pos->q8 += 3 * n;
                // fix alignment !
                pos->q8 += n & 1;
                pos->ql += LSL_NBUCKETS_MAX;
                j = j1;
            }
            pos->i += di;
        }
    }
    ASM_COMMENT("end of large (sparse) slices"); /* }}} */
}

static inline void matmul_bucket_mul_huge(struct matmul_bucket_data_s * mm, abt * dst, abt const * src, int d, struct pos_desc * pos)
{
    abobj_ptr x = mm->xab;

    abt * scrap = mm->scrap;
    ASM_COMMENT("multiplication code -- huge (very sparse) slices"); /* {{{ */
    unsigned int nhuge = *pos->ql++;
    if (d == !mm->public_->store_transposed) {
        /* ok -- hold your breath a second. */
        for(unsigned int s = 0 ; s < nhuge ; s++) {
            uint32_t di = *pos->ql++;
            ASSERT(pos->i + di <= pos->nrows_t);
            uint32_t j = 0;
            abzero(x, dst + aboffset(x, pos->i), di);
            unsigned int nlarge = *pos->ql++;
            uint32_t di_sub = iceildiv(di, nlarge);
            for( ; j < pos->ncols_t ; ) {
                uint32_t j1 = j + *pos->ql++;
                unsigned int n = *pos->ql++;
                ASSERT_ALWAYS(n <= HUGE_MPLEX_MAX * mm->bigscrapsize);
                abt * bigscrap = mm->bigscrap;
                abt const * inp = src + aboffset(x, j);
                abt * bucket[HUGE_MPLEX_MAX];
                const unsigned int * Lsizes = pos->ql;
                prepare_buckets(x,bucket,bigscrap,pos->ql,nlarge);
                pos->ql += nlarge;
                mm->fbi_time -= clock();
                ASSERT_ALWAYS((((unsigned long)pos->q8)&1)==0);
                fill_buckets_indirect(x, bucket, inp, pos->q8, n);
                mm->fbi_time += clock();
                pos->q8 += 2 * n;
                for(unsigned int k = 0 ; k < nlarge ; k++) {
                    abt * sbucket[LSL_NBUCKETS_MAX];
                    prepare_buckets(x,sbucket,scrap,pos->ql,LSL_NBUCKETS_MAX);
                    bucket[k] -= aboffset(x, Lsizes[k]);
                    mm->fbd_time -= clock();
                    fill_buckets_direct(x, sbucket, bucket[k], pos->q8, Lsizes[k]);
                    mm->fbd_time += clock();
                    pos->q8 += Lsizes[k];
                    abt * outp = dst + aboffset(x, pos->i + k * di_sub);
                    mm->asb_time -= clock();
                    apply_small_buckets(x, outp, scrap, pos->q8, pos->ql);
                    mm->asb_time += clock();
                    pos->q8 += Lsizes[k];
                    pos->ql += LSL_NBUCKETS_MAX;
                }
                j = j1;
            }
            pos->i += di;
        }
    } else {
        for(unsigned int s = 0 ; s < nhuge ; s++) {
            uint32_t di = *pos->ql++;
            ASSERT(pos->i + di <= pos->nrows_t);
            uint32_t j = 0;
            unsigned int nlarge = *pos->ql++;
            uint32_t di_sub = iceildiv(di, nlarge);
            for( ; j < pos->ncols_t ; ) {
                uint32_t j1 = j + *pos->ql++;
                unsigned int n = *pos->ql++;
                abt * bigscrap = mm->bigscrap;
                abt * outp = dst + aboffset(x, j);
                abt * bucket[HUGE_MPLEX_MAX];
                const unsigned int * Lsizes = pos->ql;
                prepare_buckets(x,bucket,bigscrap,pos->ql,nlarge);
                pos->ql += nlarge;
                const uint8_t * q8_saved = pos->q8;
                pos->q8 += 2 * n;
                for(unsigned int k = 0 ; k < nlarge ; k++) {
                    abt * sbucket[LSL_NBUCKETS_MAX];
                    prepare_buckets(x,sbucket,scrap,pos->ql,LSL_NBUCKETS_MAX);
                    const uint8_t * fill = pos->q8;
                    const uint8_t * apply = pos->q8 + Lsizes[k];
                    const abt * inp = src + aboffset(x, pos->i + k * di_sub);
                    mm->asb_time -= clock();
                    unapply_small_buckets(x, inp, scrap, apply, pos->ql);
                    mm->asb_time += clock();

                    mm->fbd_time -= clock();
                    unfill_buckets_direct(x, sbucket, bucket[k], fill, Lsizes[k]);
                    mm->fbd_time += clock();

                    pos->q8 += 2 * Lsizes[k];
                    pos->ql += LSL_NBUCKETS_MAX;
                }
                swap(pos->q8, q8_saved);
                mm->fbi_time -= clock();
                ASSERT_ALWAYS((((unsigned long)pos->q8)&1)==0);
                unfill_buckets_indirect(x, bucket, outp, pos->q8, n);
                mm->fbi_time += clock();
                pos->q8 += 2 * n;
                swap(pos->q8, q8_saved);
                j = j1;
            }
            pos->i += di;
        }
    }
    ASM_COMMENT("end of huge (very sparse) slices"); /* }}} */
}


///////////////////////////////////////////////////////////////////////
// just count how many times an iterations schedules a coefficient in
// the fbi/fbd/asb routines.

static inline void mm_finish_init(struct matmul_bucket_data_s * mm)
{
    struct pos_desc pos[1];

    pos->q16 = ptrbegin(mm->t16);
    pos->q8 = ptrbegin(mm->t8);
    pos->ql = ptrbegin(mm->aux);
    pos->i = 0;
    pos->nrows_t = mm->public_->dim[ mm->public_->store_transposed];
    pos->ncols_t = mm->public_->dim[!mm->public_->store_transposed];


    uint16_t nhstrips = *pos->q16++;
    pos->q16++;        // alignment.
    for(uint16_t s = 0 ; s < nhstrips ; s++) {
        uint32_t nrows_packed = *pos->q16++;
        int is_small2 = *pos->q16++;
        ASM_COMMENT("critical loop");
        if (is_small2) {
            uint32_t ncoeffs_slice = read32(pos->q16);
            unsigned int pad = (2 - 1) & - ncoeffs_slice;
            pos->q16 += pad;
            mm->mms2_ncoeffs += ncoeffs_slice;
            pos->q16 += ncoeffs_slice;
        } else {
            unsigned int nv = *pos->ql++;
            for(unsigned int v = 0 ; v < nv ; v++) {
                uint32_t ncoeffs_slice = read32(pos->q16);
                mm->mms1_ncoeffs += ncoeffs_slice;
                pos->q16 += 2 * ncoeffs_slice;
            }
        }
        pos->i += nrows_packed;
    }
    // we can now forget about q8
    unsigned int nlarge = *pos->ql++;
    for(unsigned int s = 0 ; s < nlarge ; s++) {
        uint32_t di = *pos->ql++;
        uint32_t j = 0;
        for( ; j < pos->ncols_t ; ) {
            uint32_t j1 = j + *pos->ql++;
            uint32_t n = *pos->ql++;
            mm->fbi_ncoeffs += n;
            mm->asb_ncoeffs += n;
            // pos->q8 += 3*n;
            pos->ql += LSL_NBUCKETS_MAX;
            j = j1;
        }
        pos->i += di;
    }
    unsigned int nhuge = *pos->ql++;
    for(unsigned int s = 0 ; s < nhuge ; s++) {
        uint32_t di = *pos->ql++;
        uint32_t j = 0;
        unsigned int nlarge = *pos->ql++;
        for( ; j < pos->ncols_t ; ) {
            uint32_t j1 = j + *pos->ql++;
            unsigned int n = *pos->ql++;
            mm->fbi_ncoeffs += n;
            mm->fbd_ncoeffs += n;
            mm->asb_ncoeffs += n;
            pos->ql += nlarge;
            pos->ql += nlarge * LSL_NBUCKETS_MAX;
            // pos->q8 += 4*n;
            j = j1;
        }
        pos->i += di;
    }

    mm->scrap = (abt *) abinit(mm->xab, mm->scrapsize * HUGE_MPLEX_MAX);
    mm->bigscrap = (abt *) abinit(mm->xab, mm->bigscrapsize * HUGE_MPLEX_MAX);
}

void matmul_bucket_mul(struct matmul_bucket_data_s * mm, abt * dst, abt const * src, int d)
{
    struct pos_desc pos[1];

    pos->q16 = ptrbegin(mm->t16);
    pos->q8 = ptrbegin(mm->t8);
    pos->ql = ptrbegin(mm->aux);
    pos->i = 0;
    pos->nrows_t = mm->public_->dim[ mm->public_->store_transposed];
    pos->ncols_t = mm->public_->dim[!mm->public_->store_transposed];

    abobj_ptr x = mm->xab;

    if (d == !mm->public_->store_transposed) {
        /* That the ``normal'' case (matrix times vector). */
    } else {
        /* d == mm->public_->store_transposed */
        /* BEWARE, it's a priori sub-optimal ! In practice, the
         * difference isn't so striking though. */
        if (mm->public_->iteration[d] == 10) {
            fprintf(stderr, "Warning: Doing many iterations with bad code\n");
        }
        /* We zero out the dst area beforehand */
        abzero(x, dst, pos->ncols_t);
    }
    matmul_bucket_mul_small(mm, dst, src, d, pos);
    matmul_bucket_mul_large(mm, dst, src, d, pos);
    matmul_bucket_mul_huge(mm, dst, src, d, pos);

    mm->public_->iteration[d]++;
}

void matmul_bucket_report(struct matmul_bucket_data_s * mm MAYBE_UNUSED, double scale)
{
    double scale0;
    scale0 = (mm->public_->iteration[0] + mm->public_->iteration[1]);

    printf("n %lu\n", mm->public_->ncoeffs);
#define one_timing_slot(x) do {						\
    double t_ ## x = (double) mm-> x ## _time / CLOCKS_PER_SEC;         \
    double a_ ## x = 1.0e9 * t_ ## x / mm-> x ## _ncoeffs / scale0;     \
    printf(#x " %.2fs ; n=%lu c ; %.2f ns/c ; scaled*%.2f : %.2f/c\n",	\
            t_ ## x,                                                    \
            mm-> x ## _ncoeffs,				                \
            a_ ## x, scale, a_ ## x * scale);                           \
} while (0)

    one_timing_slot(mms1);
    one_timing_slot(mms2);
    one_timing_slot(fbi);
    one_timing_slot(fbd);
    one_timing_slot(asb);
}

void matmul_bucket_auxv(struct matmul_bucket_data_s * mm MAYBE_UNUSED, int op MAYBE_UNUSED, va_list ap MAYBE_UNUSED) { }

void matmul_bucket_aux(struct matmul_bucket_data_s * mm, int op, ...)
{
    va_list ap;
    va_start(ap, op);
    matmul_bucket_auxv (mm, op, ap);
    va_end(ap);
}
/* vim: set sw=4: */
