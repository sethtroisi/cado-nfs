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
/* TODO re-enable some sort of assembly.
*/

#include "abase.h"
#include "readmat-easy.h"
#include "matmul-common.h"

// take only 3/4 of the L1 cache.
#define L1_CACHE_SIZE   24576
// for 8-bytes abt values, this gives 3072 items.

// make sure that it's something each core can use (i.e. divide by two
// the L2 cache size for a dual-core w/ shared L2).
#define L2_CACHE_SIZE   (1UL << 21)

/* To what should we compare the average dj value for determining when we
 * should switch from dense (small) to sparse (large) blocks */
#define DJ_CUTOFF1   1.0
/* Same, but for switching from large to huge blocks */
#define DJ_CUTOFF2   10.0

/* How many large slices in a huge slice ? */
#define HUGE_MPLEX_MIN        8
#define HUGE_MPLEX_MAX        256

#define xxxDEBUG_BUCKETS

/* This implementation builds upon the ``sliced'' variant. Several strips
 * are stored in just the same way, and when we reach sparser areas, we
 * use another method. Dense slices having been discarded, the rest is
 * split in ``large'' slices of 65536 rows at most (in fact we arrange
 * for the slices to have equal size, so it's a bit less).  For each such
 * large slice, we use scrap space equal to the L2 cache size for
 * performing matrix multiplication. While reading column vector data, we
 * copy coefficient to (up to) 256 ``buckets'', each representing 256
 * possible row indices. The dispatching information, telling in which
 * bucket a given coefficient has to be copied, is an 8-bit integer
 * stored in a table called ``main''. Each entry in table ``main'' also
 * has a second 8-bit integer indicating the offset to the next useful
 * column entry (because we're considering very sparse blocks, it is not
 * at all uncommon to see average values ~ 10 here, even when 65536 rows
 * are considered simultaneously).
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
 * indirection level for very sparse blocks. To be documented once it's
 * finished.
 */

/* This extension is used to distinguish between several possible
 * implementations of the product. The upper word correspond to the
 * implementation itself, the lower one to the n-th binary incompatible
 * change (make sure to bump it) */
#define MM_EXTENSION   "-bucket"

#define MM_MAGIC_FAMILY        0xa003UL

#define MM_MAGIC_VERSION       0x1009UL
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
    size_t scrapsize;
    size_t bigscrapsize;
    vector<uint16_t> t16;       /* For small (dense) slices */
    vector<uint8_t> t8;         /* For large (sparse) slices */
    vector<unsigned int> aux;   /* Various descriptors -- fairly small */
};

void matmul_bucket_clear(struct matmul_bucket_data_s * mm)
{
    free(mm);
}

static struct matmul_bucket_data_s * matmul_bucket_init(abobj_ptr xx MAYBE_UNUSED, param_list pl, int optimized_direction)
{
    struct matmul_bucket_data_s * mm;
    mm = (struct matmul_bucket_data_s *) malloc(sizeof(struct matmul_bucket_data_s));
    memset(mm, 0, sizeof(struct matmul_bucket_data_s));
    abobj_init_set(mm->xab, xx);

    int suggest = optimized_direction ^ MM_DIR0_PREFERS_TRANSP_MULT;
    matmul_common_init_post(mm->public_, pl, suggest);

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

/* {{{ small slices */

struct small_slice_t {
    uint32_t i0;
    uint32_t i1;
    unsigned int ncoeffs;
    double dj_avg;
    uint32_t dj_max;

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
        uint32_t dj;
        UPDATE_DIFFERENCE(j,dj,lp->first);
        if (dj > S->dj_max) { S->dj_max = dj; }
    }

    S->dj_avg = mb->ncols_t / (double) S->ncoeffs;

    /* There are two possible reasons for this slice to be discarded.
     * Either the max dj is too large -- larger than 16 bits -- or the
     * average dj is too large -- larger than some guessed value.
     */
    int keep = S->dj_max < (1UL << 16) && S->dj_avg < DJ_CUTOFF1;
    if (!keep) {
        *ptr = ptr0;
    }
    return keep;
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
    vector<uint8_t> ind[256];
    vector<uint8_t> main;
    vector<size_t> col_sizes;
    vector<size_t> pad_sizes;
};

void split_large_slice_in_vblocks(builder * mb, large_slice_t * L, large_slice_raw_t * R, unsigned int scrapsize)
{
    /* Now split into vslices */
    uint8_t * mp = ptrbegin(R->main);
    uint8_t * ip[256];
    for(int k = 0 ; k < 256 ; k++) {
        ip[k] = ptrbegin(R->ind[k]);
    }
    unsigned int vblocknum = 0;
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
        unsigned int ind_sizes[256] = {0,};
        q += 2*nf;
        for(size_t k = 0 ; k < nf ; k++) {
            ind_sizes[mp[2 * k + 1]]++;
        }
        mp += 2 * nf;
        V.auxc.push_back(nf);
        for(int k = 0 ; k < 256 ; k++) {
            ASSERT(ip[k]+ind_sizes[k] <= ptrend(R->ind[k]));
            memcpy(q, ip[k], ind_sizes[k] * sizeof(uint8_t));
            q += ind_sizes[k];
            ip[k] += ind_sizes[k];
            V.auxc.push_back(ind_sizes[k]);
        }
        ASSERT(q == ptrend(V.t8c));

        printf(" vbl%u", vblocknum);
        printf(": %ss %u..%u ", mb->colname, j0, j1);
        printf("; w %u", V.n);
        printf("; avg dj %.1f", (j1 - j0) / (double) V.n);
        if (V.pad) printf("; pad 6*%u", V.pad);
        printf("\n");

        transfer(&(L->vbl), &V);
    }
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
            L->ncoeffs, L->dj_avg, L->dj_max, 256.0 * L->dj_avg);

    if (L->dj_avg > DJ_CUTOFF2) {
        printf("-> too sparse");
        if (mb->nrows_t - i0 < HUGE_MPLEX_MIN * 65536) {
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
    vector<uint8_t> ind[256];
    vector<uint8_t> main;
};

struct huge_subslice_ptrblock_t {
    uint8_t * ip[256];
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
    uint8_t * sp = ptrbegin(R->super);
    vector<huge_subslice_ptrblock_t> ptrs(R->subs.size());
    for(unsigned int i = 0 ; i < R->subs.size() ; i++) {
        ptrs[i].mp = ptrbegin(R->subs[i].main);
        for(int k = 0 ; k < 256 ; k++) {
            ptrs[i].ip[k] = ptrbegin(R->subs[i].ind[k]);
        }
    }
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
         * at best, at it would not hurt to have it displayed only at the
         * upper level.
         */
        size_t n = 0;
        size_t np = 0;
        for(j = j0 ; j < j1 ; j++) {
            n += R->col_sizes[j];
            np += R->pad_sizes[j];
        }
        ASSERT((n+np)*2 == (size_t) (spc - sp0));
        sp = sp0;
#endif

        /* That's the same type, although it contains more stuff. */
        large_slice_vblock_t V;
        V.j0 = j0;
        V.j1 = j1;
        V.auxc.reserve(1 + H->nlarge * (1 + 256));
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
        /* XXX Notice that the reader process will consider as a j0 index
         * the farthest possible index from the previous vblock, which
         * means that we _must_ start with a zero in the main block here
         * no matter what -- even if with regard to the way data is set
         * up at the upper level (huge unique vblock), we would have had
         * a positive integer. This digression does not apply to the
         * first vblock. */
        memcpy(q, sp, 2 * nf * sizeof(uint8_t));
        if (vblocknum) {
            q[0] = 0;
        }
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
            unsigned int ind_sizes[256] = {0,};
            memcpy(q, ptrs[l].mp, L_sizes[l] * sizeof(uint8_t));
            q += L_sizes[l];
            for(unsigned int k = 0 ; k < L_sizes[l] ; k++) {
                ind_sizes[ptrs[l].mp[k]]++;
            }
            ptrs[l].mp += L_sizes[l];

            for(int k = 0 ; k < 256 ; k++) {
                ASSERT(ptrs[l].ip[k]-(ptrbegin(R->subs[l].ind[k]))+ind_sizes[k] <= (ptrdiff_t) R->subs[l].ind[k].size());
                memcpy(q, ptrs[l].ip[k], ind_sizes[k] * sizeof(uint8_t));
                q += ind_sizes[k];
                ptrs[l].ip[k] += ind_sizes[k];
                V.auxc.push_back(ind_sizes[k]);
            }
        }
        ASSERT(q == ptrend(V.t8c));

        printf(" vbl%u", vblocknum);
        printf(": %ss %u..%u ", mb->colname, j0, j1);
        printf("; w %u", V.n);
        printf("; avg dj %.1f", (j1 - j0) / (double) V.n);
        if (V.pad) printf("; pad 6*%u", V.pad);
        printf("\n");

        transfer(&(H->vbl), &V);
    }
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

    H->nlarge = iceildiv(i1 - i0, 65536);
    ASSERT(H->nlarge >= HUGE_MPLEX_MIN);
    ASSERT(H->nlarge  < HUGE_MPLEX_MAX);
    unsigned int lsize = iceildiv(i1 - i0, H->nlarge);
    ASSERT(lsize <= 65536);
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

        printf("Ssl %u: %ss %u..%u...", s, mb->rowname, i0, i1);
        fflush(stdout);

        int keep = builder_do_small_slice(mb, &S, &ptr, i0, i1);

        printf(" w %" PRIu32 " ; avg dj %.1f ; max dj %u\n",
                S.ncoeffs, S.dj_avg, S.dj_max);
        fflush(stdout);

        if (!keep) {
            printf("Switching to large slices. Ssl %u to be redone\n", s);
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
    unsigned int nlarge_slices = iceildiv(rem_nrows,65536);
    uint32_t done = 0;
    for(unsigned int s = 0 ; s < nlarge_slices ; s++) {
        large_slice_t L[1];
        uint32_t i0 = * p_i0 +  s      * rem_nrows / nlarge_slices;
        uint32_t i1 = * p_i0 + (s + 1) * rem_nrows / nlarge_slices;

        printf("Lsl %u %ss %u..%u", s, mb->rowname, i0, i1);
        fflush(stdout);

        int keep = builder_do_large_slice(mb, L, i0, i1, scrapsize);

        if (!keep) {
            printf("Switching to huge slices. Lsl %u to be redone\n", s);
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
    unsigned int nhuge_slices = iceildiv(rem_nrows, HUGE_MPLEX_MAX * 65536);
    uint32_t done = 0;
    for(unsigned int s = 0 ; s < nhuge_slices ; s++) {
        huge_slice_t H[1];
        uint32_t i0 = * p_i0 +  s      * rem_nrows / nhuge_slices;
        uint32_t i1 = * p_i0 + (s + 1) * rem_nrows / nhuge_slices;

        printf("Hsl %u %ss %u..%u", s, mb->rowname, i0, i1);
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
        mm->t16.push_back(0);
        mm->t16.push_back(S->L.size() & 0xffff);
        mm->t16.push_back(S->L.size() >> 16);
        for(small_slice_t::Lci_t lp = S->L.begin() ; lp != S->L.end() ; ++lp) {
            mm->t16.push_back(lp->first);
            mm->t16.push_back(lp->second);
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
            mm->t8.insert(mm->t8.end(), V.t8c.begin(), V.t8c.end());
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
            mm->t8.insert(mm->t8.end(), V.t8c.begin(), V.t8c.end());
        }
        mm->public_->ncoeffs += H->ncoeffs;
    }
}
/* }}} */

struct matmul_bucket_data_s * matmul_bucket_build(abobj_ptr xx, const char * filename, param_list pl, int optimized_direction)
{
    struct matmul_bucket_data_s * mm;
    builder mb[1];
    mm = matmul_bucket_init(xx, pl, optimized_direction);
    builder_init(mb, filename, mm);


    uint32_t main_i0 = 0;

    unsigned int npack = L1_CACHE_SIZE;
    if (pl) param_list_parse_uint(pl, "l1_cache_size", &npack);
    npack /= abbytes(mm->xab,1);

    list<small_slice_t> Sq[1];
    builder_do_all_small_slices(mb, &main_i0, Sq, npack);
    builder_push_small_slices(mm, Sq);

    free(mb->data[0]);
    mb->data[0] = NULL;

    builder_switch_to_vertical(mb, main_i0);

    unsigned int scrapsize = L2_CACHE_SIZE;
    if (pl) param_list_parse_uint(pl, "l2_cache_size", &scrapsize);
    scrapsize /= abbytes(mm->xab,1);
    mm->scrapsize = scrapsize;

    list<large_slice_t> Lq[1];
    builder_do_all_large_slices(mb, &main_i0, Lq, scrapsize);
    builder_push_large_slices(mm, Lq);

    /* Don't deallocate right now, since data is still to be used for
     * huge slices.  */

    unsigned int bigscrapsize = scrapsize; // (1 << 24)/abbytes(mm->xab,1));
    list<huge_slice_t> Hq[1];
    builder_do_all_huge_slices(mb, &main_i0, Hq, bigscrapsize);
    builder_push_huge_slices(mm, Hq);

    /* done, at last ! */
    free(mb->data[1]);
    mb->data[1] = NULL;

    return mm;
}

struct matmul_bucket_data_s * matmul_bucket_reload_cache(abobj_ptr xx, const char * filename, param_list pl, int optimized_direction)/*{{{*/
{
    struct matmul_bucket_data_s * mm;
    FILE * f;

    mm = matmul_bucket_init(xx, pl, optimized_direction);
    f = matmul_common_reload_cache_fopen(abbytes(xx,1), mm->public_, filename, MM_EXTENSION, MM_MAGIC);
    if (f == NULL) { free(mm); return NULL; }

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

    return mm;
}/*}}}*/

void matmul_bucket_save_cache(struct matmul_bucket_data_s * mm, const char * filename)/*{{{*/
{
    FILE * f;

    f = matmul_common_save_cache_fopen(abbytes(mm->xab,1), mm->public_, filename, MM_EXTENSION, MM_MAGIC);

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

    ASM_COMMENT("multiplication code -- dense slices"); /* {{{ */
    uint16_t nhstrips = *pos->q16++;
    pos->q16++;        // alignment.
    if (d == !mm->public_->store_transposed) {
        for(uint16_t s = 0 ; s < nhstrips ; s++) {
            uint32_t nrows_packed = read32(pos->q16);
            ASSERT(pos->i + nrows_packed <= pos->nrows_t);
            abt * where = dst + aboffset(x, pos->i);
            abzero(x, where, nrows_packed);
            ASM_COMMENT("critical loop");
            abt const * from = src;
            uint32_t ncoeffs_slice = read32(pos->q16);
            uint32_t j = 0;
            for(uint32_t c = 0 ; c < ncoeffs_slice ; c++) {
                j += *pos->q16++;
                ASSERT(j < pos->ncols_t);
                uint32_t di = *pos->q16++;
                abadd(x, where + aboffset(x, di), from + aboffset(x, j));
            }
            pos->i += nrows_packed;
            ASM_COMMENT("end of critical loop");
        }
    } else {
        /* d == mm->public_->store_transposed */
        for(uint16_t s = 0 ; s < nhstrips ; s++) {
            uint32_t j = 0;
            uint32_t nrows_packed = read32(pos->q16);
            ASM_COMMENT("critical loop");
            uint32_t ncoeffs_slice = read32(pos->q16);
            for(uint32_t c = 0 ; c < ncoeffs_slice ; c++) {
                j += *pos->q16++;
                uint32_t di = *pos->q16++;
                abadd(x, dst + aboffset(x, j), src + aboffset(x, pos->i + di));
            }
            pos->i += nrows_packed;
            ASM_COMMENT("end of critical loop");
        }
    }
    ASM_COMMENT("end of dense slices"); /* }}} */
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

static inline void
fill_buckets_indirect(abobj_ptr x, abt ** sb, const abt * z, const uint8_t * q, unsigned int n)
{
    /* Dispatch data found in z[0]...z[f(n-1)] such that z[f(i)] is in
     * array pointed to by sb[q[2*i+1]]. The function f(i) is given by
     * the sum q[0]+q[2]+...+q[2*(i-1)]. Exactly 2n coefficients are
     * expected in q[] All the sb[] pointers are increased */
    for(unsigned int c = 0 ; c < n ; c++) {
        z += aboffset(x, *q);
        q++;
        abcopy(x, sb[*q], z, 1);
        sb[*q]++;
        q++;
    }
}

static inline void
unfill_buckets_indirect(abobj_ptr x, abt ** sb, abt * z, const uint8_t * q, unsigned int n)
{
    /* Does the converse of the above */
    for(unsigned int c = 0 ; c < n ; c++) {
        z += aboffset(x, *q);
        q++;
        abadd(x, z, sb[*q]);
        sb[*q]++;
        q++;
    }
}

static inline void apply_small_buckets(abobj_ptr x, abt * dst, const abt * z, const uint8_t * q, const unsigned int * ql)
{
    /* This ``applies'' the 256 small buckets whose respective lengths
     * are given by ql[0] to ql[255].
     *
     * For ql[k] <= i < ql[k+1], z[i] is added to dst[k*256+qk[i]], with
     * qk = q + ql[0] + ... + ql[k-1].
     */
    for(int k = 0 ; k < 256 ; k++) {
        unsigned int l = ql[k];
        for( ; l-- ; ) {
            /* For padding coeffs, the assertion can fail if
             * we choose a row not equal to (0,0) -- first in
             * the first bucket.
             */
            abadd(x, dst + aboffset(x, *q), z);
            z++;
            q++;
        }
        dst += aboffset(x, 256);
    }
}

static inline void unapply_small_buckets(abobj_ptr x, const abt * src, abt * z, const uint8_t * q, const unsigned int * ql)
{
    /* converse of the above */
    for(int k = 0 ; k < 256 ; k++) {
        unsigned int l = ql[k];
        for( ; l-- ; ) {
            /* For padding coeffs, the assertion can fail if
             * we choose a row not equal to (0,0) -- first in
             * the first bucket.
             */
            abcopy(x, z, src + aboffset(x, *q), 1);
            z++;
            q++;
        }
        src += aboffset(x, 256);
    }
}

static inline void
fill_buckets_direct(abobj_ptr x, abt ** sb, const abt * z, const uint8_t * q, unsigned int n)
{
    /* Dispatch data found in z[0]...z[n] such that z[i] is in array
     * pointed to by sb[q[i]]. Exactly n coefficients are expected. All
     * the sb[] pointers are increased */
    for(unsigned int c = 0 ; c < n ; c++) {
        abcopy(x, sb[q[c]], z + c, 1);
        sb[q[c]]++;
    }
}

static inline void
unfill_buckets_direct(abobj_ptr x, abt ** sb, abt * z, const uint8_t * q, unsigned int n)
{
    /* Does the converse of the above */
    for(unsigned int c = 0 ; c < n ; c++) {
        abcopy(x, z + c, sb[q[c]], 1);
        sb[q[c]]++;
    }
}

static inline void matmul_bucket_mul_large(struct matmul_bucket_data_s * mm, abt * dst, abt const * src, int d, struct pos_desc * pos)
{
    abobj_ptr x = mm->xab;

    abt * scrap = (abt *) abinit(x, mm->scrapsize);
    ASM_COMMENT("multiplication code -- large (sparse) slices"); /* {{{ */
    unsigned int nlarge = *pos->ql++;
    if (d == !mm->public_->store_transposed) {
        for(unsigned int s = 0 ; s < nlarge ; s++) {
            uint32_t di = *pos->ql++;
            ASSERT(pos->i + di <= pos->nrows_t);
            uint32_t j = 0;
            abzero(x, dst + aboffset(x, pos->i), di);
            for( ; j < pos->ncols_t ; ) {
                uint32_t j1 = j + *pos->ql++;
                uint32_t n = *pos->ql++;
                abt * bucket[256];
                abt const * inp = src + aboffset(x, j);
                abt * outp = dst + aboffset(x, pos->i);
                prepare_buckets(x,bucket,scrap,pos->ql,256);
                fill_buckets_indirect(x, bucket, inp, pos->q8, n);
                apply_small_buckets(x, outp, scrap, pos->q8+2*n, pos->ql);
                pos->q8 += 3*n;
                pos->ql += 256;
                j = j1;
            }
            pos->i += di;
        }
    } else {
        for(unsigned int s = 0 ; s < nlarge ; s++) {
            uint32_t di = *pos->ql++;
            ASSERT(pos->i + di <= pos->nrows_t);
            uint32_t j = 0;
            for( ; j < pos->ncols_t ; ) {
                uint32_t j1 = j + *pos->ql++;
                uint32_t n = *pos->ql++;
                abt * bucket[256];
                abt * outp = dst + aboffset(x, j);
                abt const *  inp = src + aboffset(x, pos->i);
                prepare_buckets(x,bucket,scrap,pos->ql,256);
                unapply_small_buckets(x, inp, scrap, pos->q8+2*n, pos->ql);
                unfill_buckets_indirect(x, bucket, outp, pos->q8, n);
                pos->q8 += 3 * n;
                pos->ql += 256;
                j = j1;
            }
            pos->i += di;
        }
    }
    ASM_COMMENT("end of large (sparse) slices"); /* }}} */
    abclear(x, scrap, mm->scrapsize);
}

static inline void matmul_bucket_mul_huge(struct matmul_bucket_data_s * mm, abt * dst, abt const * src, int d, struct pos_desc * pos)
{
    abobj_ptr x = mm->xab;

    abt * scrap = (abt *) abinit(x, mm->scrapsize);
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
                abt * bigscrap = (abt *) abinit(x, n);
                abt const * inp = src + aboffset(x, j);
                abt * bucket[HUGE_MPLEX_MAX];
                const unsigned int * Lsizes = pos->ql;
                prepare_buckets(x,bucket,bigscrap,pos->ql,nlarge);
                pos->ql += nlarge;
                fill_buckets_indirect(x, bucket, inp, pos->q8, n);
                pos->q8 += 2 * n;
                for(unsigned int k = 0 ; k < nlarge ; k++) {
                    abt * sbucket[256];
                    prepare_buckets(x,sbucket,scrap,pos->ql,256);
                    bucket[k] -= Lsizes[k];
                    fill_buckets_direct(x, sbucket, bucket[k], pos->q8, Lsizes[k]);
                    pos->q8 += Lsizes[k];
                    abt * outp = dst + aboffset(x, pos->i + k * di_sub);
                    apply_small_buckets(x, outp, scrap, pos->q8, pos->ql);
                    pos->q8 += Lsizes[k];
                    pos->ql += 256;
                }
                abclear(x, bigscrap, n);
                j = j1;
            }
            pos->i += di;
        }
    } else {
#if 1
        for(unsigned int s = 0 ; s < nhuge ; s++) {
            uint32_t di = *pos->ql++;
            ASSERT(pos->i + di <= pos->nrows_t);
            uint32_t j = 0;
            unsigned int nlarge = *pos->ql++;
            uint32_t di_sub = iceildiv(di, nlarge);
            for( ; j < pos->ncols_t ; ) {
                uint32_t j1 = j + *pos->ql++;
                unsigned int n = *pos->ql++;
                abt * bigscrap = (abt *) abinit(x, n);
                abt * outp = dst + aboffset(x, j);
                abt * bucket[HUGE_MPLEX_MAX];
                const unsigned int * Lsizes = pos->ql;
                prepare_buckets(x,bucket,bigscrap,pos->ql,nlarge);
                pos->ql += nlarge;
                const uint8_t * q8_saved = pos->q8;
                pos->q8 += 2 * n;
                for(unsigned int k = 0 ; k < nlarge ; k++) {
                    abt * sbucket[256];
                    prepare_buckets(x,sbucket,scrap,pos->ql,256);
                    const uint8_t * fill = pos->q8;
                    const uint8_t * apply = pos->q8 + Lsizes[k];
                    const abt * inp = src + aboffset(x, pos->i + k * di_sub);
                    unapply_small_buckets(x, inp, scrap, apply, pos->ql);

                    unfill_buckets_direct(x, sbucket, bucket[k], fill, Lsizes[k]);
                    pos->q8 += 2 * Lsizes[k];
                    pos->ql += 256;
                }
                swap(pos->q8, q8_saved);
                unfill_buckets_indirect(x, bucket, outp, pos->q8, n);
                pos->q8 += 2 * n;
                swap(pos->q8, q8_saved);
                abclear(x, bigscrap, n);
                j = j1;
            }
            pos->i += di;
        }
#endif
    }
    ASM_COMMENT("end of huge (very sparse) slices"); /* }}} */
    abclear(x, scrap, mm->scrapsize);
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

void matmul_bucket_report(struct matmul_bucket_data_s * mm MAYBE_UNUSED) { }

void matmul_bucket_auxv(struct matmul_bucket_data_s * mm MAYBE_UNUSED, int op MAYBE_UNUSED, va_list ap MAYBE_UNUSED) { }

void matmul_bucket_aux(struct matmul_bucket_data_s * mm, int op, ...)
{
    va_list ap;
    va_start(ap, op);
    matmul_bucket_auxv (mm, op, ap);
    va_end(ap);
}
/* vim: set sw=4: */
