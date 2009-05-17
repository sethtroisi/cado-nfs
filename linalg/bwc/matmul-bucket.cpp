#ifndef __cplusplus
#define _GNU_SOURCE         /* asprintf */
#endif
#define _DARWIN_C_SOURCE    /* for asprintf. _ANSI_SOURCE must be undefined */

#define __STDC_LIMIT_MACROS
#define __STDC_FORMAT_MACROS
/* Manage the in-memory data for the matrix */
/* It's in C++ because the STL is handy, but that's really all there is
 * to it... */

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
/* Same, but for switching from large to very large blocks */
#define DJ_CUTOFF2   1.0e20

/* How many large slices in a very large slice ? */
#define VL_MPLEX        16

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

#define MM_MAGIC_VERSION       0x1007UL
#define MM_MAGIC (MM_MAGIC_FAMILY << 16 | MM_MAGIC_VERSION)

/* see matmul-basic.c */
#define MM_DIR0_PREFERS_TRANSP_MULT   1

struct matmul_bucket_data_s {
    /* repeat the fields from the public_ interface */
    struct matmul_public_s public_[1];
    /* now our private fields */
    abobj_t xab;
    size_t scrapsize;
    vector<uint16_t> t16;       /* For small (dense) slices */
    vector<uint8_t> t8;         /* For large (sparse) slices */
    vector<unsigned int> lsl;   /* Descriptors for large slices */
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

void builder_do_all_small_slices(builder * mb, uint32_t * p_i0, list<small_slice_t> * Sq, unsigned int npack)
{
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

void builder_push_small_slices(struct matmul_bucket_data_s * mm, list<small_slice_t> * Sq)
{
    mm->t16.push_back(Sq->size());
    mm->t16.push_back(0);       // for alignment.

    printf("Flushing %zu small slices\n", Sq->size());

    for(unsigned int s = 0 ; !Sq->empty() ; s++, Sq->pop_front()) {
        small_slice_t * S = &(Sq->front());

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

struct large_slice_vblock_t {
    uint32_t j0;
    uint32_t j1;
    unsigned int n;             // number of coeffs.
    unsigned int pad;           // (half the) number of padding coeffs.
    vector<uint8_t> t8c;        // read: contribution to t8.
    unsigned int ind_sizes[256];// contribution to lsl.
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
    uint8_t * mp = &R->main.front();
    uint8_t * ip[256];
    for(int k = 0 ; k < 256 ; k++) {
        ip[k] = &R->ind[k].front();
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
        /* We have n coeffs, counting padding. This makes 2n bytes in
         * main, and n bytes in the bucket blocks */
        V.n = n;
        V.pad = np;
        size_t nf = n + 2 * np;
        V.t8c.resize(3 * nf);
        uint8_t * q = &V.t8c.front();
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
        memset(V.ind_sizes, 0, sizeof(V.ind_sizes));
        q += 2*nf;
        for(size_t k = 0 ; k < nf ; k++) {
            V.ind_sizes[mp[2 * k + 1]]++;
        }
        mp += 2 * nf;
        for(int k = 0 ; k < 256 ; k++) {
            ASSERT(ip[k]-(&R->ind[k].front())+V.ind_sizes[k] <= (ptrdiff_t) R->ind[k].size());
            memcpy(q, ip[k], V.ind_sizes[k] * sizeof(uint8_t));
            q += V.ind_sizes[k];
            ip[k] += V.ind_sizes[k];
        }
        ASSERT(q == V.t8c.size() + &V.t8c.front());

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
        printf("-> too sparse\n");
        /* We won't keep this slice */
        for(uint32_t j = 0 ; j < mb->ncols_t ; j++) {
            mb->colptrs[j] -= R->col_sizes[j];
        }
        return 0;
    }

    printf("\n");

    split_large_slice_in_vblocks(mb, L, R, scrapsize);

    return 1;
}

void builder_push_one_large_slice(struct matmul_bucket_data_s * mm, large_slice_t * L)
{
    mm->lsl.push_back(L->i1 - L->i0);
    for( ; ! L->vbl.empty() ; L->vbl.pop_front()) {
        large_slice_vblock_t & V(L->vbl.front());
        mm->lsl.push_back(V.j1 - V.j0);
        mm->lsl.push_back(V.n + 2 * V.pad);
        mm->lsl.insert(mm->lsl.end(), &(V.ind_sizes[0]), &(V.ind_sizes[256]));
        mm->t8.insert(mm->t8.end(), V.t8c.begin(), V.t8c.end());
    }
    mm->public_->ncoeffs += L->ncoeffs;
}

void builder_push_large_slices(struct matmul_bucket_data_s * mm, list<large_slice_t> * Lq)
{
    printf("Flushing %zu large slices\n", Lq->size());

    for(unsigned int l = 0 ; !Lq->empty() ; l++, Lq->pop_front()) {
        large_slice_t& L(Lq->front());
        builder_push_one_large_slice(mm, &L);
    }
}

void builder_do_all_large_slices(builder * mb, uint32_t * p_i0, list<large_slice_t> * Lq, unsigned int scrapsize)
{
    unsigned int rem_nrows = mb->nrows_t - *p_i0;
    unsigned int nlarge_slices = iceildiv(rem_nrows,65536);
    unsigned int s;
    for(s = 0 ; s < nlarge_slices ; s++) {
        large_slice_t L[1];
        uint32_t i0 = * p_i0 +  s      * rem_nrows / nlarge_slices;
        uint32_t i1 = * p_i0 + (s + 1) * rem_nrows / nlarge_slices;

        printf("Lsl %u %ss %u..%u\n", s, mb->rowname, i0, i1);
        fflush(stdout);
        builder_do_large_slice(mb, L, i0, i1, scrapsize);
        transfer(Lq, L);
    }
    *p_i0 += s * rem_nrows / nlarge_slices;
}

struct matmul_bucket_data_s * matmul_bucket_build(abobj_ptr xx, const char * filename, param_list pl, int optimized_direction)
{
    struct matmul_bucket_data_s * mm;
    builder mb[1];
    mm = matmul_bucket_init(xx, pl, optimized_direction);
    builder_init(mb, filename, mm);


    /* {{{ Do small slices. */
    /* First take a guess for their expected size, and
     * then arrange for slices of approximately equal size.
     */
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

    free(mb->data[1]);
    mb->data[1] = NULL;

    builder_push_large_slices(mm, Lq);

#if 0
    /* Time to flatten the queue of large slices. First we compute the
     * required storage for the mm->lsl and mm->lsl arrays.
     */
    ASSERT(Lq->size() == nlarge_slices);
    printf("Total %zu large slices\n", Lq->size());
    list<large_slice_t>::iterator ll;
    size_t switch_vlsl = 0;
    for(ll = Lq->begin() ; ll != Lq->end() ; ll++) {
        if (ll->dj_sum / (double) ll->ncoeffs > DJ_CUTOFF2)
            break;
        switch_vlsl++;
    }
    printf("%zu large slices saved for very large slices\n",
            Lq->size() - switch_vlsl);
    printf("Flushing %zu large slices\n", switch_vlsl);
    size_t lsl_size = 1;
    size_t t8_size = 0;
    size_t k = 0;
    for(ll = Lq->begin() ; ll != Lq->end() ; ll++, k++) {
        if (k == switch_vlsl)
            break;
        lsl_size += ll->lslc.size();
        t8_size += ll->t8c.size();
    }
    mm->lsl.reserve(lsl_size);
    mm->t8.reserve(t8_size);
    /* NEW: Now we also give the number of lsl before going vl */
    mm->lsl.push_back(switch_vlsl);
    for(k = 0 ; k < switch_vlsl ; k++) {
        large_slice_t const& L(Lq->front());
        mm->lsl.insert(mm->lsl.end(), L.lslc.begin(), L.lslc.end());
        mm->t8.insert(mm->t8.end(), L.t8c.begin(), L.t8c.end());
        mm->public_->ncoeffs += L.ncoeffs;
        Lq->pop_front();
    }
    printf("Assembling %zu large slices into very large slices\n",
            Lq->size() - switch_vlsl);
    for( ; !Lq->empty() ; ) {
        /* Pick a list of at most VL_MPLEX large slices. Do so with no
         * copy */
        size_t mplexing = 0;
        vector<large_slice_t> fan(VL_MPLEX);
        for( ; !Lq->empty() && mplexing < VL_MPLEX ; mplexing++) {
            swap(fan[mplexing], Lq->front());
            Lq->pop_front();
        }
    }
#endif
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
    size_t n16, n8, nlsl;
    MATMUL_COMMON_READ_ONE32(n16, f);
    MATMUL_COMMON_READ_ONE32(n8, f);
    MATMUL_COMMON_READ_ONE32(nlsl, f);
    mm->t16.resize(n16);
    mm->t8.resize(n8);
    mm->lsl.resize(nlsl);
    MATMUL_COMMON_READ_MANY16(&(mm->t16.front()), n16, f);
    MATMUL_COMMON_READ_MANY8(&(mm->t8.front()), n8, f);
    MATMUL_COMMON_READ_MANY32(&(mm->lsl.front()), nlsl, f);

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
    size_t nlsl = mm->lsl.size();
    MATMUL_COMMON_WRITE_ONE32(n16, f);
    MATMUL_COMMON_WRITE_ONE32(n8, f);
    MATMUL_COMMON_WRITE_ONE32(nlsl, f);
    MATMUL_COMMON_WRITE_MANY16(&(mm->t16.front()), n16, f);
    MATMUL_COMMON_WRITE_MANY8(&(mm->t8.front()), n8, f);
    MATMUL_COMMON_WRITE_MANY32(&(mm->lsl.front()), nlsl, f);

    fclose(f);
}/*}}}*/

static inline uint32_t read32(uint16_t const * & q) /* {{{ */
{
    uint32_t res;
    res = *q++;
    res |= ((uint32_t) *q++) << 16;
    return res;
} /* }}} */

void matmul_bucket_mul(struct matmul_bucket_data_s * mm, abt * dst, abt const * src, int d)
{
    const uint16_t * q16 = &(mm->t16.front());
    const uint8_t * q8 = &(mm->t8.front());
    const unsigned int * ql = &(mm->lsl.front());

    uint32_t i = 0;
    abobj_ptr x = mm->xab;
    uint32_t nrows_t = mm->public_->dim[ mm->public_->store_transposed];
    uint32_t ncols_t = mm->public_->dim[!mm->public_->store_transposed];

    ASM_COMMENT("multiplication code -- dense slices"); /* {{{ */
    uint16_t nhstrips = *q16++;
    q16++;        // alignment.
    if (d == !mm->public_->store_transposed) {
        for(uint16_t s = 0 ; s < nhstrips ; s++) {
            uint32_t nrows_packed = read32(q16);
            ASSERT(i + nrows_packed <= nrows_t);
            abt * where = dst + aboffset(x, i);
            abzero(x, where, nrows_packed);
            ASM_COMMENT("critical loop");
            abt const * from = src;
            uint32_t ncoeffs_slice = read32(q16);
            uint32_t j = 0;
            for(uint32_t c = 0 ; c < ncoeffs_slice ; c++) {
                j += *q16++;
                ASSERT(j < ncols_t);
                uint32_t di = *q16++;
                abadd(x, where + aboffset(x, di), from + aboffset(x, j));
            }
            ASM_COMMENT("end of critical loop");
            i += nrows_packed;
        }
    } else {
        /* d == mm->public_->store_transposed */
        /* BEWARE, it's a priori sub-optimal ! In practice, the
         * difference isn't so striking though. */
        if (mm->public_->iteration[d] == 10) {
            fprintf(stderr, "Warning: Doing many iterations with bad code\n");
        }

        abzero(x, dst, ncols_t);
        for(uint16_t s = 0 ; s < nhstrips ; s++) {
            uint32_t j = 0;
            uint32_t nrows_packed = read32(q16);
            ASM_COMMENT("critical loop");
            uint32_t ncoeffs_slice = read32(q16);
            for(uint32_t c = 0 ; c < ncoeffs_slice ; c++) {
                j += *q16++;
                uint32_t di = *q16++;
                abadd(x, dst + aboffset(x, j), src + aboffset(x, i + di));
            }
            i += nrows_packed;
            ASM_COMMENT("end of critical loop");
        }
    }
    ASM_COMMENT("end of dense slices"); /* }}} */

    abt * scrap = (abt *) abinit(x, mm->scrapsize);

    ASM_COMMENT("multiplication code -- sparse slices"); /* {{{ */
    if (d == !mm->public_->store_transposed) {
        for( ; i < nrows_t ; ) {
            uint32_t di = *ql++;
            uint32_t j = 0;
            abzero(x, dst + aboffset(x, i), di);
            for( ; j < ncols_t ; ) {
                uint32_t dj = *ql++;
                uint32_t nread = *ql++;
                abt const * inp = src + aboffset(x, j);
                abt * bucket[256];
                abt * fence[256];
                abt * z = scrap;
                for(int k = 0 ; k < 256 ; k++) {
                    bucket[k] = z;
                    fence[k] = bucket[k] + ql[k];
                    z += aboffset(x, ql[k]);
                }
                /* fill buckets */
                for( ; nread-- ; ) {
                    ASSERT(inp - src + *q8 < aboffset(x, j+dj));
                    inp += aboffset(x, *q8);
                    q8++;
                    ASSERT(bucket[*q8] < fence[*q8]);
                    ASSERT(bucket[*q8] - scrap < (ptrdiff_t) mm->scrapsize);
                    abcopy(x, bucket[*q8]++, inp, 1);
                    q8++;
                }
                /* apply buckets */
                z = scrap;
                abt * outp = dst + aboffset(x, i);
                for(int k = 0 ; k < 256 ; k++) {
                    unsigned int l = ql[k];
                    for( ; l-- ; ) {
                        /* For padding coeffs, the assertion can fail if
                         * we choose a row not equal to (0,0) -- first in
                         * the first bucket.
                         */
                        ASSERT(outp + *q8 < dst + nrows_t);
                        abadd(x, outp + aboffset(x, *q8), z);
                        z++;
                        q8++;
                    }
                    outp += aboffset(x, 256);
                }
                ql += 256;
                j += dj;
            }
            i += di;
        }
    } else {
        for( ; i < nrows_t ; ) {
            uint32_t di = *ql++;
            uint32_t j = 0;
            for( ; j < ncols_t ; ) {
                uint32_t dj = *ql++;
                uint32_t nread = *ql++;
                abt * bucket[256];
                abt * fence[256];
                abt * z = scrap;
                for(int k = 0 ; k < 256 ; k++) {
                    bucket[k] = z;
                    fence[k] = bucket[k] + ql[k];
                    z += aboffset(x, ql[k]);
                }
                abt * outp = dst + aboffset(x, j);
                /* reverse-fill buckets */
                const uint8_t * q8_saved = q8;
                q8 += 2 * nread;
                z = scrap;
                abt const *  inp = src + aboffset(x, i);
                for(int k = 0 ; k < 256 ; k++) {
                    unsigned int l = ql[k];
                    for( ; l-- ; ) {
                        /* For padding coeffs, the assertion can fail if
                         * we choose a row not equal to (0,0) -- first in
                         * the first bucket.
                         */
                        ASSERT(inp + *q8 < src + nrows_t);
                        abcopy(x, z, inp + aboffset(x, *q8), 1);
                        z++;
                        q8++;
                    }
                    inp += aboffset(x, 256);
                }
                /* reverse-apply buckets */
                std::swap(q8, q8_saved);
                for( ; nread-- ; ) {
                    ASSERT(outp - dst + *q8 < aboffset(x, j+dj));
                    outp += aboffset(x, *q8);
                    q8++;
                    ASSERT(bucket[*q8] < fence[*q8]);
                    ASSERT(bucket[*q8] - scrap < (ptrdiff_t) mm->scrapsize);
                    abadd(x, outp, bucket[*q8]++);
                    q8++;
                }
                q8 = q8_saved;
                ql += 256;
                j += dj;
            }
            i += di;
        }
    }
    ASM_COMMENT("end of sparse slices"); /* }}} */

    abclear(x, scrap, mm->scrapsize);


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
