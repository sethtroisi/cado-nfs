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
#define DJ_LARGE_JUMP   4.0

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
 * buckets'', so as to be able to reuse our scrap space buffer for the
 * following column. Thus a ``large slice'' is typically split into
 * several ``vertical blocks''. Extremely sparse large slices will have
 * only one vertical block.
 *
 * ``Applying the buckets'' means taking a list of column coefficients
 * which have been directed to the bucket because it is known that they
 * affect at least one of the 256 corresponding row indices on output. So
 * in addition the the N stored coefficients in the bucket, we have a
 * list of N 8-bit integers indicating which row index is to be modified.
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

#define UPDATE_DIFFERENCE(i,di,mem) do { \
    di = mem - i; i = mem; mem = di;     \
} while (0)

struct matmul_bucket_data_s * matmul_bucket_build(abobj_ptr xx, const char * filename, param_list pl, int optimized_direction)
{
    struct matmul_bucket_data_s * mm;
    mm = matmul_bucket_init(xx, pl, optimized_direction);

    /* {{{ read matrix data in proper order */
    uint32_t * data[2];

    /* we'll call ``rows'' the streams of data found in data[0] */
    read_easy(filename,
            &data[ mm->public_->store_transposed],
            &data[!mm->public_->store_transposed],
            &mm->public_->dim[0],
            &mm->public_->dim[1]);

    uint32_t nrows_t = mm->public_->dim[ mm->public_->store_transposed];
    uint32_t ncols_t = mm->public_->dim[!mm->public_->store_transposed];
    const char * rowname = rowcol[ mm->public_->store_transposed];
    const char * colname = rowcol[!mm->public_->store_transposed];
    /* }}} */

    /* {{{ Do small slices. */
    /* First take a guess for their expected size, and
     * then arrange for slices of approximately equal size.
     */
    unsigned int npack = L1_CACHE_SIZE;
    if (pl) param_list_parse_uint(pl, "l1_cache_size", &npack);
    npack /= abbytes(mm->xab,1);
    unsigned int nslices = iceildiv(nrows_t, npack);

    unsigned int nslices_index = mm->t16.size();
    mm->t16.push_back(0);
    mm->t16.push_back(0);        // placeholder for alignment

    unsigned int s;

    /* Get some info on the columns */
    vector<uint32_t *> colheads(ncols_t, NULL);
    vector<uint32_t *> colptrs(ncols_t, NULL);
    { /* {{{ initialize column pointers */
        uint32_t * ptr = data[1];
        for(uint32_t j = 0 ; j < ncols_t ; j++) {
            colheads[j] = ptr;
            colptrs[j] = ptr + 1;
            ptr += 1 + *ptr;
        }
    } /* }}} */

    uint32_t * ptr = data[0];
    for(s = 0 ; s < nslices ; s++) {
        uint32_t i0 =  s    * nrows_t / nslices;
        uint32_t i1 = (s+1) * nrows_t / nslices;

        /* We're doing a new slice */
        typedef std::vector<std::pair<uint32_t, uint32_t> > L_t;
        typedef L_t::const_iterator Lci_t;
        typedef L_t::iterator Li_t;
        L_t L;
        uint32_t * ptr0 = ptr;
        for(uint32_t i = i0 ; i < i1 ; i++) {
            for(unsigned int j = 0 ; j < *ptr ; j++) {
                L.push_back(std::make_pair(ptr[1+j],i-i0));
            }
            ptr += 1 + *ptr;
        }
        /* L is the list of (column index, row index) of all
         * coefficients in the current horizontal slice */
        std::sort(L.begin(), L.end());

        /* There are two possible reasons for this slice to be discarded.
         * Either the max dj is too large -- larger than 16 bits -- or the
         * average dj is too large -- larger than some guessed value.
         */
        /* {{{ To start with, convert all j indices to differences */
        uint32_t j = 0;
        uint32_t dj_max = 0;
        double dj_sum=0;
        uint32_t weight=0;
        for(Li_t lp = L.begin() ; lp != L.end() ; ++lp) {
            uint32_t dj;
            UPDATE_DIFFERENCE(j,dj,lp->first);
            dj_sum += dj;
            if (dj > dj_max) { dj_max = dj; }
            weight++;
        }
        /*}}}*/
        double dj_avg = dj_sum / weight;

        if (dj_max >= (1UL << 16) || dj_avg >= DJ_LARGE_JUMP) {
            ptr = ptr0;
            break;
        }

        /* validate this slice. */

        printf("Ssl %u: %ss %u..%u ; w %" PRIu32 " ; avg dj %.1f ; max dj %u\n",
                s, rowname, i0, i1, weight, dj_avg,dj_max);
        
        mm->t16.push_back(i1-i0);
        mm->t16.push_back(0);
        mm->t16.push_back(L.size() & 0xffff);
        mm->t16.push_back(L.size() >> 16);
        for(Lci_t lp = L.begin() ; lp != L.end() ; ++lp) {
            mm->t16.push_back(lp->first);
            mm->t16.push_back(lp->second);
        }
        mm->t16[nslices_index]++;
        mm->public_->ncoeffs += weight;
    }

    /*}}}*/

    unsigned int end_sslices = s * nrows_t / nslices;
    unsigned int rem_nrows = nrows_t - end_sslices;
    unsigned int nlarge_slices = iceildiv(rem_nrows,65536);
    /* {{{ advance all column pointers. Takes time, of course */
    for(uint32_t j = 0 ; j < ncols_t ; j++) {
        uint32_t len = * colheads[j];
        uint32_t * v0 = colheads[j] + 1;
        uint32_t * c = colptrs[j];
        for( ; c - v0 < (ptrdiff_t) len && *c < (uint32_t) end_sslices ; c++) ;
        colptrs[j] = c;
    } /* }}} */

    unsigned int scrapsize = L2_CACHE_SIZE;
    if (pl) param_list_parse_uint(pl, "l2_cache_size", &scrapsize);
    scrapsize /= abbytes(mm->xab,1);

    mm->scrapsize = scrapsize;

    for(s = 0 ; s < nlarge_slices ; s++) {
        uint32_t i0 = end_sslices +  s      * rem_nrows / nlarge_slices;
        uint32_t i1 = end_sslices + (s + 1) * rem_nrows / nlarge_slices;
        mm->lsl.push_back(i1 - i0);
        /* {{{ What is the expected weight for this batch of vstrips ? */
        unsigned int weight = 0;
        uint32_t * q = ptr;
        for(uint32_t i = i0 ; i < i1 ; i++) {
            weight += *q;
            q += 1 + *q;
        }
        ptr = q;
        /* }}} */

        printf("Lsl %u %ss %u..%u : w %u, expected bucket hit rate 1/%.1f\n",
                s, rowname, i0, i1, weight,
                256.0 * ncols_t / (double) weight);


        /* We'll create vblocks one by one, and flush them as soon as the
         * scrap size is exceeded. It's better than trying to guess in
         * advance the number of vblocks since the distribution of the
         * density can't be expected to be uniform across vblocks */
        unsigned int v = 0;

        for(uint32_t j = 0 ; j < ncols_t ; ) {
            /* Start a new block */
            vector<uint8_t> ind[256];
            vector<uint8_t> main;
            uint32_t j0 = j;
            unsigned int locweight = 0;
            uint32_t lastj = j0;
            uint32_t dj_max = 0;
            double dj_sum=0;
            int pad=0;
            // unsigned int used_input = 0;
            for( ; j < ncols_t ; j++) {
                uint32_t len = * colheads[j];
                uint32_t * v0 = colheads[j] + 1;
                uint32_t * c = colptrs[j];
                vector<uint8_t> my_ind[256];
                vector<uint8_t> my_main;
                unsigned int my_weight = 0;
                unsigned int my_pad = 0;
                uint32_t my_lastj = lastj;

                for( ; c - v0 < (ptrdiff_t) len && *c < i1 ; c++) {/*{{{*/
                    uint32_t i = *c - i0;
                    uint8_t w = i / 256;
                    uint32_t diff = j - my_lastj;
                    for( ; diff > 255 ; ) {
                        /* This is equivalent to adding two coefficients
                         * to the matrix, but in the same position. */
                        my_main.push_back(255);
                        my_main.push_back(0);
                        my_ind[0].push_back(0);
                        my_main.push_back(0);
                        my_main.push_back(0);
                        my_ind[0].push_back(0);
                        my_lastj += 255;
                        diff -= 255;
                        my_pad++;
                    }
                    my_main.push_back(diff);
                    my_main.push_back(w);
                    my_ind[w].push_back((uint8_t) i);
                    my_lastj = j;
                    my_weight++;
                }/*}}}*/

                if (locweight + my_weight + 2*(pad + my_pad) > scrapsize) {
                    /* Cannot fit this :-( */
                    break;
                }

                // used_input += (my_weight != 0);

                /* This one fits in the current vblock, copy it */
                main.insert(main.end(), my_main.begin(), my_main.end());
                for(int k = 0 ; k < 256 ; k++) {
                    ind[k].insert(ind[k].end(),
                            my_ind[k].begin(), my_ind[k].end());
                }
                locweight += my_weight;
                pad += my_pad;
                uint32_t diff = my_lastj - lastj;
                dj_sum += diff;
                if (diff > dj_max) dj_max = diff;
                lastj = my_lastj;
                colptrs[j] = c;
            }
            ASSERT_ALWAYS(!main.empty());
            double avg_fill = (double) locweight / (j - j0);
            printf(" vbl%u", v);
            printf(": %ss %u..%u ", colname, j0, j);
            printf("; w %u", locweight);
            printf("; f %.2f", avg_fill);
            printf("; avg dj %.1f", dj_sum / locweight);
            printf("; max dj %u", dj_max);
            if (pad) printf("; pad 6*%d", pad);
            printf("\n");
            mm->lsl.push_back(j - j0);
            mm->lsl.push_back(locweight + 2 * pad);

            /* flush */
            mm->t8.insert(mm->t8.end(), main.begin(), main.end());
            for(int k = 0 ; k < 256 ; k++) {
                mm->lsl.push_back(ind[k].size());
                mm->t8.insert(mm->t8.end(), ind[k].begin(), ind[k].end());
            }
#ifdef  DEBUG_BUCKETS
            /* re-count the number of occurences of each. */
            unsigned int recount[256] = {0, };
            for(unsigned int z = 1 ; z < main.size() ; z+=2) {
                recount[main[z]]++;
            }
            for(int k = 0 ; k < 256 ; k++) {
                ASSERT_ALWAYS(recount[k] == mm->lsl[mm->lsl.size()-256+k]);
            }
#endif
            mm->public_->ncoeffs += locweight;
            v++;
        }
    }
    free(data[0]);
    free(data[1]);
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
        ASSERT_ALWAYS(0);
        /* That's going to be a nightmare. */
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
