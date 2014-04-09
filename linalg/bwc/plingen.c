/* Copyright (C) 1999--2007 Emmanuel Thom'e --- see LICENSE file */
#include "cado.h"

#include <sys/time.h>
#include <sys/types.h>
#include <dirent.h>
#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>

#include <assert.h>

#include "portability.h"
#include "macros.h"
#include "utils.h"
#include "abase.h"
#include "memusage.h"
#include "lingen-matpoly.h"

#include "lingen-bigmatpoly.h"

#ifdef HAVE_MPIR
#include "lingen-bigmatpoly-ft.h"
#endif

#include "bw-common.h"		/* Handy. Allows Using global functions
                                 * for recovering parameters */
#include "bw-common-mpi.h"
#include "filenames.h"
#include "plingen.h"
#include "plingen-tuning.h"

/* Call tree for methods within this program:
 *
 * Two separate entry points.
 *
 * bw_biglingen_collective        when need to go collective
 *     |<--> bw_biglingen_recursive , loops to bw_biglingen_collective
 *     \->bw_lingen_single               when it makes sense to do this locally again
 * bw_lingen_single                      when computation can be done locally
 *    |<->bw_lingen_recursive
 *    |->bw_lingen_basecase
 *
 */
static unsigned int display_threshold = 10;
static int with_timings = 0;

/* This is an indication of the number of bytes we read at a time for A
 * (input) and F (output) */
static unsigned int io_block_size = 1 << 20;

/* If non-zero, then reading from A is actually replaced by reading from
 * a random generator */
static unsigned int random_input_length = 0;

gmp_randstate_t rstate;

#ifdef HAVE_MPIR
static int caching = 1;
#else
static int caching = 0;
#endif

struct bmstatus_s {
    dims d[1];
    unsigned int t;
    int * lucky;

    double t_basecase;
    double t_mp;
    double t_mul;

    unsigned int lingen_threshold;
    unsigned int lingen_mpi_threshold;
    int mpi_dims[2]; /* mpi_dims[0] = mpi[0] * thr[0] */
    MPI_Comm world;     /* reordered, in fact */
};
typedef struct bmstatus_s bmstatus[1];
typedef struct bmstatus_s *bmstatus_ptr;

/* {{{ col sorting */
/* We sort only with respect to the global delta[] parameter. As it turns
 * out, we also access the column index in the same aray and sort with
 * respect to it, but this is only cosmetic.
 *
 * Note that unlike what is asserted in other coiped of the code, sorting
 * w.r.t. the local delta[] value is completely useless. Code which
 * relies on this should be fixed.
 */

typedef int (*sortfunc_t) (const void*, const void*);

static int col_cmp(const int x[2], const int y[2])
{
    for(int i = 0 ; i < 2 ; i++) {
        int d = x[i] - y[i];
        if (d) return d;
    }
    return 0;
}

/* }}} */

static inline unsigned int expected_pi_length(dims * d, unsigned int len)/*{{{*/
{
    /* The idea is that we want something which may account for something
     * exceptional, bounded by probability 2^-64. This corresponds to a
     * column in e (matrix of size m*b) to be spontaneously equal to
     * zero. This happens with probability (#K)^-m.
     * The solution to
     * (#K)^(-m*x) > 2^-64
     * is m*x*log_2(#K) < 64
     *
     * We thus need to get an idea of the value of log_2(#K).
     *
     * (Note that we know that #K^abgroupsize(ab) < 2^64, but that bound
     * might be very gross).
     *
     * The same esitmate can be used to appreciate what we mean by
     * ``luck'' in the end. If a column happens to be zero more than
     * expected_pi_length(d,0) times in a row, then the cause must be
     * more than sheer luck, and we use it to detect generating rows.
     */
    unsigned int m = d->m;
    unsigned int n = d->n;
    unsigned int b = m + n;
    abdst_field ab MAYBE_UNUSED = d->ab;
    unsigned int res = iceildiv(len * m, b);
    mpz_t p;
    mpz_init(p);
    abfield_characteristic(ab, p);
    unsigned int l = mpz_sizeinbase(p, 2);
    l *= abfield_degree(ab);    /* roughtly log_2(#K) */
    mpz_clear(p);
    // unsigned int safety = iceildiv(abgroupsize(ab), m * sizeof(abelt));
    unsigned int safety = iceildiv(64, m * l);
    return res + safety;
}/*}}}*/

/* Forward declaration, it's used by the recursive version */
static int bw_lingen_single(bmstatus_ptr bm, matpoly pi, matpoly E, unsigned int *delta);
static int bw_biglingen_collective(bmstatus_ptr bm, bigmatpoly pi, bigmatpoly E, unsigned int *delta);

/* This destructively cancels the first len coefficients of E, and
 * computes the appropriate matrix pi which achieves this. The
 * elimination is done in accordance with the nominal degrees found in
 * delta.
 *
 * The result is expected to have degree ceil(len*m/b) coefficients, so
 * that E*pi is divisible by X^len.
 */

static int bw_lingen_basecase(bmstatus_ptr bm, matpoly pi, matpoly E, unsigned int *delta) /*{{{*/
{
    dims * d = bm->d;
    unsigned int m = d->m;
    unsigned int n = d->n;
    unsigned int b = m + n;
    abdst_field ab = d->ab;
    ASSERT(E->m == m);
    ASSERT(E->n == b);

    ASSERT(pi->m == 0);
    ASSERT(pi->n == 0);
    ASSERT(pi->alloc == 0);

    /* Allocate something large enough for the result. This will be
     * soon freed anyway. Set it to identity. */
    unsigned int pi_room_base = expected_pi_length(d, E->size);

    matpoly_init(ab, pi, b, b, pi_room_base);
    pi->size = pi_room_base;

    /* Also keep track of the
     * number of coefficients for the columns of pi. Set pi to Id */

    unsigned int *pi_lengths = malloc(b * sizeof(unsigned int));
    for(unsigned int i = 0 ; i < b ; i++) {
        abset_ui(ab, matpoly_coeff(ab, pi, i, i, 0), 1);
        pi_lengths[i] = 1;
    }

    /* This is used below */
    int generator_found = 0;

    /* Keep a list of columns which have been used as pivots at the
     * previous iteration */
    unsigned int * pivots = malloc(m * sizeof(unsigned int));
    int * is_pivot = malloc(b * sizeof(int));
    memset(is_pivot, 0, b * sizeof(int));

    matpoly e;
    matpoly_init(ab, e, m, b, 1);
    e->size = 1;

    for (unsigned int t = 0; t < E->size ; t++, bm->t++) {

        /* {{{ Update the columns of e for degree t. Save computation
         * time by not recomputing those which can easily be derived from
         * previous iteration. Notice that the columns of e are exactly
         * at the physical positions of the corresponding columns of pi.
         */

        abelt_ur tmp_ur;
        abelt_ur_init(ab, &tmp_ur);

        abvec_ur e_ur;
        abvec_ur_init(ab, &e_ur, m*b);
        abvec_ur_set_zero(ab, e_ur, m*b);
        for(unsigned int j = 0 ; j < b ; j++) {
            if (is_pivot[j]) continue;
            /* We should never have to recompute from pi using discarded
             * columns. Discarded columns should always correspond to
             * pivots */
            ASSERT_ALWAYS(bm->lucky[j] >= 0);
            for(unsigned int s = 0 ; s < pi_lengths[j] ; s++) {
                for(unsigned int i = 0 ; i < m ; i++) {
                    for(unsigned int k = 0 ; k < b ; k++) {
                        abmul_ur(ab, tmp_ur,
                                matpoly_coeff(ab, E, i, k, t - s),
                                matpoly_coeff(ab, pi, k, j, s));
                        abelt_ur_add(ab,
                                abvec_ur_coeff_ptr(ab, e_ur, i*b+j),
                                abvec_ur_coeff_ptr_const(ab, e_ur, i*b+j),
                                tmp_ur);
                    }
                }
            }
        }
        unsigned int newluck = 0;
        for(unsigned int j = 0 ; j < b ; j++) {
            if (is_pivot[j]) continue;
            unsigned int nz = 0;
            for(unsigned int i = 0 ; i < m ; i++) {
                abreduce(ab, matpoly_coeff(ab, e, i, j, 0),
                        abvec_ur_coeff_ptr(ab, e_ur, i*b+j));
                nz += abcmp_ui(ab, matpoly_coeff(ab, e, i, j, 0), 0) == 0;
            }
            if (nz == m) {
                newluck++, bm->lucky[j]++;
            } else if (bm->lucky[j] > 0) {
                bm->lucky[j] = 0;
            }
        }
        abelt_ur_clear(ab, &tmp_ur);
        abvec_ur_clear(ab, &e_ur, m*b);
        if (newluck) {
            /* If newluck == n, then we probably have a generator. We add an
             * extra guarantee. newluck==n, for a total of k iterations in a
             * row, means m*n*k coefficients cancelling magically. We would
             * like this to be impossible by mere chance. Thus we want n*k >
             * luck_mini, which can easily be checked */

            int luck_mini = expected_pi_length(d, 0);
            unsigned int luck_sure = 0;

            printf("t=%d, canceled columns:", bm->t);
            for(unsigned int j = 0 ; j < b ; j++) {
                if (bm->lucky[j] > 0) {
                    printf(" %u", j);
                    luck_sure += bm->lucky[j] > luck_mini;
                }
            }

            if (newluck == n && luck_sure == n) {
                if (!generator_found) {
                    printf(", complete generator found, for sure");
                }
                generator_found = 1;
            }
            printf(".\n");
        }
        /* }}} */

        if (generator_found) break;

        int (*ctable)[2] = malloc(b * 2 * sizeof(int));
        /* {{{ Now see in which order I may look at the columns of pi, so
         * as to keep the nominal degrees correct. In contrast with what
         * we used to do before, we no longer apply the permutation to
         * delta. So the delta[] array keeps referring to physical
         * indices, and we'll tune this in the end. */
        for(unsigned int j = 0; j < b; j++) {
            ctable[j][0] = delta[j];
            ctable[j][1] = j;
        }
        qsort(ctable, b, 2 * sizeof(int), (sortfunc_t) & col_cmp);
        /* }}} */

        /* {{{ Now do Gaussian elimination */
        memset(is_pivot, 0, b * sizeof(int));
        unsigned int r = 0;
        /* Loop through logical indices */
        for(unsigned int jl = 0; jl < b; jl++) {
            unsigned int j = ctable[jl][1];
            unsigned int u = 0;
            /* {{{ Find the pivot */
            for( ; u < m ; u++) {
                if (abcmp_ui(ab, matpoly_coeff(ab, e, u, j, 0), 0) != 0)
                    break;
            }
            if (u == m) continue;
            assert(r < m);
            /* }}} */
            pivots[r++] = j;
            is_pivot[j] = 1;
            /* {{{ Cancel this coeff in all other columns. */
            abelt inv;
            abinit(ab, &inv);
            int rc = abinv(ab, inv, matpoly_coeff(ab, e, u, j, 0));
            if (!rc) {
                fprintf(stderr, "Error, found a factor of the modulus: ");
                abfprint(ab, stderr, inv);
                fprintf(stderr, "\n");
                exit(EXIT_FAILURE);
            }
            abneg(ab, inv, inv);
            for (unsigned int kl = jl + 1; kl < b ; kl++) {
                unsigned int k = ctable[kl][1];
                if (abcmp_ui(ab, matpoly_coeff(ab, e, u, k, 0), 0) == 0)
                    continue;
                // add lambda = e[u,k]*-e[u,j]^-1 times col j to col k.
                abelt lambda;
                abinit(ab, &lambda);
                abmul(ab, lambda, inv, matpoly_coeff(ab, e, u, k, 0));

                assert(delta[j] <= delta[k]);
                /* {{{ Apply on both e and pi */
                abelt tmp;
                abinit(ab, &tmp);
                for(unsigned int i = 0 ; i < m ; i++) {
                    /* TODO: Would be better if mpfq had an addmul */
                    abmul(ab, tmp, lambda, matpoly_coeff(ab, e, i, j, 0));
                    abadd(ab,
                            matpoly_coeff(ab, e, i, k, 0),
                            matpoly_coeff(ab, e, i, k, 0),
                            tmp);
                }
                if (bm->lucky[k] < 0) {
                    /* This column is already discarded, don't bother */
                    continue;
                }
                if (bm->lucky[j] < 0) {
                    /* This column is discarded. This is going to
                     * invalidate another column of pi. Not a problem,
                     * unless it's been marked as lucky previously ! */
                    ASSERT_ALWAYS(bm->lucky[k] <= 0);
                    printf("Column %u discarded from now on (through addition from column %u)\n", k, j);
                    bm->lucky[k] = -1;
                    continue;
                }
                for(unsigned int i = 0 ; i < b ; i++) {
                    /* Beware. One may be tempted to think that the code
                     * above is dubious. It is, in fact, not a problem.
                     * As long as the delta[] array undergoes no
                     * disturbing modification, everything is ok.
                     */
                    if (pi_lengths[k] < pi_lengths[j])
                        pi_lengths[k] = pi_lengths[j];
                    for(unsigned int s = 0 ; s < pi_lengths[j] ; s++) {
                        /* TODO: Would be better if mpfq had an addmul */
                        abmul(ab, tmp, lambda,
                                matpoly_coeff(ab, pi, i, j, s));
                        abadd(ab,
                                matpoly_coeff(ab, pi, i, k, s),
                                matpoly_coeff(ab, pi, i, k, s),
                                tmp);
                    }
                }
                abclear(ab, &tmp);
                /* }}} */
                abclear(ab, &lambda);
            }
            abclear(ab, &inv); /* }}} */
        }
        /* }}} */
        free(ctable);

        ASSERT_ALWAYS(r == m);

        /* {{{ Now for all pivots, multiply column in pi by x */
        for (unsigned int j = 0; j < b ; j++) {
            if (!is_pivot[j]) continue;
            if (pi_lengths[j] >= pi->alloc) {
                if (!generator_found) {
                    matpoly_realloc(ab, pi, pi->alloc + MAX(pi->alloc / (m+n), 1));
                    printf("t=%u, expanding allocation for pi (now %zu%%) ; lengths: ",
                            bm->t,
                            100 * pi->alloc / pi_room_base);
                    for(unsigned int j = 0; j < b; j++)
                        printf(" %u", pi_lengths[j]);
                    printf("\n");
                } else {
                    ASSERT_ALWAYS(bm->lucky[j] <= 0);
                    if (bm->lucky[j] == 0)
                        printf("t=%u, column %u discarded from now on\n",
                                bm->t, j);
                    bm->lucky[j] = -1;
                    pi_lengths[j]++;
                    delta[j]++;
                    continue;
                }
            }
            matpoly_multiply_column_by_x(ab, pi, j, pi_lengths[j]);
            pi_lengths[j]++;
            delta[j]++;
        }
        /* }}} */
    }
    pi->size = 0;
    for(unsigned int j = 0; j < b; j++) {
        if (pi_lengths[j] > pi->size)
            pi->size = pi_lengths[j];
    }
    pi->size = MIN(pi->size, pi->alloc);
    matpoly_clear(ab, e);
    free(is_pivot);
    free(pivots);
    free(pi_lengths);   /* What shall we do with this one ??? */

    return generator_found;
}/*}}}*/

double start_time = -1;
void info_init_timer()
{
    start_time = wct_seconds();
}

static const char *size_disp(size_t s, char buf[16])/*{{{*/
{
    char *prefixes = "bkMGT";
    double ds = s;
    const char *px = prefixes;
    for (; px[1] && ds > 500.0;) {
        ds /= 1024.0;
        px++;
    }
    snprintf(buf, 10, "%.1f%c", ds, *px);
    return buf;
}/*}}}*/

void print_info_mp(unsigned int t, matpoly_ptr A, matpoly_ptr B)/*{{{*/
{
    if (with_timings) {
        char buf1[16];
        char buf2[16];
        size_disp(Memusage2(), buf1);
        size_disp(PeakMemusage(), buf2);
        printf("[%.3f %s %s] t=%u, MP(%zu, %zu) --> %zu",
                wct_seconds() - start_time, buf1, buf2,
                t,
                A->size, B->size, MAX(A->size, B->size) - MIN(A->size, B->size) + 1);
    } else {
        printf("t=%u, MP(%zu, %zu) --> %zu\n",
                t,
                A->size, B->size, MAX(A->size, B->size) - MIN(A->size, B->size) + 1);
    }
}/*}}}*/

void print_info_mul(unsigned int t, matpoly_ptr A, matpoly_ptr B)/*{{{*/
{
    if (with_timings) {
        char buf1[16];
        char buf2[16];
        size_disp(Memusage2(), buf1);
        size_disp(PeakMemusage(), buf2);
        printf("[%.3f %s %s] t=%u, MUL(%zu, %zu) --> %zu",
                wct_seconds() - start_time, buf1, buf2,
                t,
                A->size, B->size, A->size + B->size - 1);
    } else {
        printf("t=%u, MUL(%zu, %zu) --> %zu\n",
                t,
                A->size, B->size, A->size + B->size - 1);
    }
}/*}}}*/

void print_info_mpi_mp(unsigned int t, bigmatpoly_ptr A, bigmatpoly_ptr B)/*{{{*/
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank) return;
    if (with_timings) {
        char buf1[16];
        char buf2[16];
        size_disp(Memusage2(), buf1);
        size_disp(PeakMemusage(), buf2);
        printf("[%.3f %s %s] t=%u, MPI-MP%s(%zu, %zu) --> %zu",
                wct_seconds() - start_time, buf1, buf2,
                t,
                caching?"-caching":"",
                A->size, B->size, MAX(A->size, B->size) - MIN(A->size, B->size) + 1);
    } else {
        printf("t=%u, MPI-MP%s(%zu, %zu) --> %zu\n",
                t,
                caching?"-caching":"",
                A->size, B->size, MAX(A->size, B->size) - MIN(A->size, B->size) + 1);
    }
}/*}}}*/

void print_info_mpi_mul(unsigned int t, bigmatpoly_ptr A, bigmatpoly_ptr B)/*{{{*/
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank) return;
    if (with_timings) {
        char buf1[16];
        char buf2[16];
        size_disp(Memusage2(), buf1);
        size_disp(PeakMemusage(), buf2);
        printf("[%.3f %s %s] t=%u, MPI-MUL%s(%zu, %zu) --> %zu",
                wct_seconds() - start_time, buf1, buf2,
                t,
                caching?"-caching":"",
                A->size, B->size, A->size + B->size - 1);
    } else {
        printf("t=%u, MPI-MUL%s(%zu, %zu) --> %zu\n",
                t,
                caching?"-caching":"",
                A->size, B->size, A->size + B->size - 1);
    }
}/*}}}*/

static int bw_lingen_recursive(bmstatus_ptr bm, matpoly pi, matpoly E, unsigned int *delta) /*{{{*/
{
    dims * d = bm->d;
    abdst_field ab = d->ab;
    int done;
    double tt;

    /* XXX I think we have to start with something large enough to get
     * all coefficients of E_right correct */
    size_t half = E->size - (E->size / 2);

    // fprintf(stderr, "Enter %s\n", __func__);
    /* declare an lazy-alloc all matrices */
    matpoly E_left;
    matpoly pi_left;
    matpoly pi_right;
    matpoly E_right;
    matpoly_init(ab, E_left, 0, 0, 0);
    matpoly_init(ab, pi_left, 0, 0, 0);
    matpoly_init(ab, pi_right, 0, 0, 0);
    matpoly_init(ab, E_right, 0, 0, 0);

    matpoly_truncate(ab, E_left, E, half);
    done = bw_lingen_single(bm, pi_left, E_left, delta);
    ASSERT_ALWAYS(pi_left->size);
    matpoly_clear(ab, E_left);

    if (done) {
        matpoly_swap(pi_left, pi);
        matpoly_clear(ab, pi_left);
        // fprintf(stderr, "Leave %s\n", __func__);
        return 1;
    }

    tt = seconds();
    matpoly_rshift(ab, E, E, half - pi_left->size + 1);
    if (E->size > display_threshold) print_info_mp(bm->t, E, pi_left);
    matpoly_mp(ab, E_right, E, pi_left);
    tt = seconds() - tt;
    if (with_timings && E->size > display_threshold) printf("\t[%.2f]\n", tt);
    bm->t_mp += tt;
    done = bw_lingen_single(bm, pi_right, E_right, delta);
    matpoly_clear(ab, E_right);

    tt = seconds();
    if (E->size > display_threshold) print_info_mul(bm->t, pi_left, pi_right);
    matpoly_mul(ab, pi, pi_left, pi_right);
    tt = seconds() - tt;
    if (with_timings && E->size > display_threshold) printf("\t[%.2f]\n", tt);
    bm->t_mul += tt;
    matpoly_clear(ab, pi_left);
    matpoly_clear(ab, pi_right);

    // fprintf(stderr, "Leave %s\n", __func__);
    return done;
}/*}}}*/

/* This version works over MPI */
static int bw_biglingen_recursive(bmstatus_ptr bm, bigmatpoly pi, bigmatpoly E, unsigned int *delta) /*{{{*/
{
    dims * d = bm->d;
    abdst_field ab = d->ab;
    int done;
    double tt;

    int rank;
    MPI_Comm_rank(bm->world, &rank);

    // fprintf(stderr, "Enter %s\n", __func__);
    /* XXX I think we have to start with something large enough to get
     * all coefficients of E_right correct */
    size_t half = E->size - (E->size / 2);

    /* declare an lazy-alloc all matrices */
    bigmatpoly E_left;
    bigmatpoly E_right;
    bigmatpoly pi_left;
    bigmatpoly pi_right;
    bigmatpoly_init(ab, E_left, E, 0, 0, 0);
    bigmatpoly_init(ab, pi_left, pi, 0, 0, 0);
    bigmatpoly_init(ab, pi_right, pi, 0, 0, 0);
    bigmatpoly_init(ab, E_right, E, 0, 0, 0);

    bigmatpoly_truncate_loc(ab, E_left, E, half);
    done = bw_biglingen_collective(bm, pi_left, E_left, delta);
    bigmatpoly_clear(ab, E_left);

    if (done) {
        bigmatpoly_swap(pi_left, pi);
        bigmatpoly_clear(ab, pi_left);
        // fprintf(stderr, "Leave %s\n", __func__);
        return 1;
    }

    bigmatpoly_rshift(ab, E, E, half - pi_left->size + 1);
    tt = wct_seconds();
    if (E->size > display_threshold) print_info_mpi_mp(bm->t, E, pi_left);
#ifdef  HAVE_MPIR
    if (caching) {
        bigmatpoly_mp_caching(ab, E_right, E, pi_left);
    } else {
        bigmatpoly_mp(ab, E_right, E, pi_left);
    }
#else
        bigmatpoly_mp(ab, E_right, E, pi_left);
#endif
    tt = wct_seconds() - tt;
    if (!rank && with_timings && E->size > display_threshold) printf("\t[%.2f]\n", tt);
    bm->t_mp = tt;
    done = bw_biglingen_collective(bm, pi_right, E_right, delta);
    bigmatpoly_clear(ab, E_right);

    tt = wct_seconds();
    if (E->size > display_threshold) print_info_mpi_mul(bm->t, pi_left, pi_right);
#ifdef  HAVE_MPIR
    if (caching) {
        bigmatpoly_mul_caching(ab, pi, pi_left, pi_right);
    } else {
        bigmatpoly_mul(ab, pi, pi_left, pi_right);
    }
#else
    bigmatpoly_mul(ab, pi, pi_left, pi_right);
#endif
    tt = wct_seconds() - tt;
    if (!rank && with_timings && E->size > display_threshold) printf("\t[%.2f]\n", tt);
    bm->t_mul += tt;
    bigmatpoly_clear(ab, pi_left);
    bigmatpoly_clear(ab, pi_right);

    // fprintf(stderr, "Leave %s\n", __func__);
    return done;
}/*}}}*/

static int bw_lingen_single(bmstatus_ptr bm, matpoly pi, matpoly E, unsigned int *delta) /*{{{*/
{
    int rank;
    MPI_Comm_rank(bm->world, &rank);
    ASSERT_ALWAYS(!rank);

    // ASSERT_ALWAYS(E->size < bm->lingen_mpi_threshold);

    // fprintf(stderr, "Enter %s\n", __func__);
    int done;
    if (E->size < bm->lingen_threshold) {
        bm->t_basecase -= seconds();
        done = bw_lingen_basecase(bm, pi, E, delta);
        bm->t_basecase += seconds();
    } else {
        done = bw_lingen_recursive(bm, pi, E, delta);
    }
    // fprintf(stderr, "Leave %s\n", __func__);
    return done;
}/*}}}*/

static int bw_biglingen_collective(bmstatus_ptr bm, bigmatpoly pi, bigmatpoly E, unsigned int *delta)/*{{{*/
{
    dims *d = bm->d;
    abdst_field ab = d->ab;
    unsigned int m = d->m;
    unsigned int n = d->n;
    unsigned int b = m + n;
    int done;
        int rank;
        MPI_Comm_rank(bm->world, &rank);

    // fprintf(stderr, "Enter %s\n", __func__);
    if (E->size >= bm->lingen_mpi_threshold)  {
        done = bw_biglingen_recursive(bm, pi, E, delta);
    } else {
        /* Fall back to local code */
        /* This entails gathering E locally, computing pi locally, and
         * dispathing it back. */

        matpoly sE, spi;
        matpoly_init(ab, sE, m, b, E->size);
        matpoly_init(ab, spi, 0, 0, 0);
        bigmatpoly_gather_mat_alt2(ab, sE, E);
        /* Only the master node does the local computation */
        if (!rank)
            done = bw_lingen_single(bm, spi, sE, delta);
        MPI_Bcast(&done, 1, MPI_INT, 0, bm->world);
        MPI_Bcast(delta, b, MPI_UNSIGNED, 0, bm->world);
        MPI_Bcast(bm->lucky, b, MPI_UNSIGNED, 0, bm->world);
        MPI_Bcast(&(bm->t), 1, MPI_UNSIGNED, 0, bm->world);
        bigmatpoly_scatter_mat_alt2(ab, pi, spi);
        /* Don't forget to broadcast delta from root node to others ! */
        matpoly_clear(ab, spi);
        matpoly_clear(ab, sE);
    }
    // fprintf(stderr, "Leave %s\n", __func__);

    return done;
}/*}}}*/

/**********************************************************************/

void get_minmax_delta_on_solutions(bmstatus_ptr bm, unsigned int * res, unsigned int * delta)/*{{{*/
{
    dims * d = bm->d;
    unsigned int m = d->m;
    unsigned int n = d->n;
    unsigned int maxdelta = 0;
    unsigned int mindelta = UINT_MAX;
    for(unsigned int j = 0 ; j < m + n ; j++) {
        if (bm->lucky[j] <= 0)
            continue;
        if (delta[j] > maxdelta)
            maxdelta = delta[j];
        if (delta[j] < mindelta)
            mindelta = delta[j];
    }
    res[0] = mindelta;
    res[1] = maxdelta;
}/*}}}*/

unsigned int get_max_delta_on_solutions(bmstatus_ptr bm, unsigned int * delta)/*{{{*/
{
    unsigned int res[2];
    get_minmax_delta_on_solutions(bm, res, delta);
    return res[1];
}/*}}}*/

/**********************************************************************/

struct bm_io_s {/*{{{*/
    bmstatus_ptr bm;
    unsigned int t0;
    FILE * f;
    const char * input_file;
    const char * output_file;
    int ascii;
    /* This is only a rolling window ! */
    matpoly A, F;
    unsigned int (*fdesc)[2];
    /* This k is the coefficient in A(X) div X of the next coefficient to
     * be read. This is thus the total number of coefficients of A(X) div
     * X which have been read so far.
     * In writing mode, k is the number of coefficients of F which have
     * been written so far.
     */
    unsigned int k;

    unsigned int guessed_length;
};

typedef struct bm_io_s bm_io[1];
typedef struct bm_io_s * bm_io_ptr;/*}}}*/

/* The reading mode of bm_io is streaming, but with a look-back
 * functionality: we want to be able to access coefficients a few places
 * earlier, so we keep them in memory.
 *
 * The writing mode has a write-ahead feature. Coefficients of the result
 * are written at various times, some earlier than others. The time span
 * is the same as for the reading mode.
 */

double avg_matsize(abdst_field ab, unsigned int m, unsigned int n, int ascii)/*{{{*/
{
    if (!ascii) {
        /* Easy case first. If we have binary input, then we know a priori
         * that the input data must have size a multiple of the element size.
         */
        size_t elemsize = abvec_elt_stride(ab, 1);
        size_t matsize = elemsize * m * n;
        return matsize;
    }

    /* Ascii is more complicated. We're necessarily fragile here.
     * However, assuming that each coefficient comes with only one space,
     * and each matrix with an extra space (this is how the GPU program
     * prints data -- not that this ends up having a considerable impact
     * anyway...), we can guess the number of bytes per matrix. */

    /* Formula for the average number of digits of an integer mod p,
     * written in base b:
     *
     * (k-((b^k-1)/(b-1)-1)/p)  with b = Ceiling(Log(p)/Log(b)).
     */
    double avg;
    mpz_t p, a;
    mpz_init(p);
    mpz_init(a);
    abfield_characteristic(ab, p);
    unsigned long k = ceil(log(mpz_get_d(p))/log(10));
    unsigned long b = 10;
    mpz_ui_pow_ui(a, b, k);
    mpz_sub_ui(a, a, 1);
    mpz_fdiv_q_ui(a, a, b-1);
    avg = k - mpz_get_d(a) / mpz_get_d(p);
    mpz_clear(p);
    mpz_clear(a);
    // printf("Expect roughly %.2f decimal digits for integers mod p.\n", avg);
    double matsize = (avg + 1) * m * n + 1;
    // printf("Expect roughly %.2f bytes for each sequence matrix.\n", matsize);
    return matsize;
}/*}}}*/


/* We write the coefficients of the reversed polynomial \hat{F*\pi}, in
 * increasing degree order. Thus the first coefficients written
 * corresponds to high degree coefficients in \pi.
 * This is mostly for historical reasons, since in fact, we would prefer
 * having the coefficients in the reverse order.
 */

/* let mindelta and maxdelta be (as their name suggests) the minimum and
 * maximum over all deltas corresponding to solution columns.
 *
 * For 0<=i<n, coeff i,j,k of pi becomes coeff i,j,k' of the final f,
 * with k'=k-(maxdelta-delta[j]).
 *
 * For 0<=i<m, coeff n+i,j,k of pi becomes coeff c[i],j,k' of the final f,
 * with k'=k-(maxdelta-delta[j])-(t0-e[j]).
 *
 * Therefore the maximum write-behind distance is (maxdelta-mindelta)+t0.
 * We need one coeff more (because offset goes from 0 to maxoffset).
 */
unsigned int bm_io_set_write_behind_size(bm_io_ptr aa, unsigned int * delta)
{
    bmstatus_ptr bm = aa->bm;
    dims * d = bm->d;
    abdst_field ab = d->ab;
    unsigned int res[2];
    get_minmax_delta_on_solutions(bm, res, delta);
    unsigned int mindelta = res[0];
    unsigned int maxdelta = res[1];
    unsigned int window = maxdelta - mindelta + aa->t0 + 1;
    matpoly_realloc(ab, aa->F, window);
    matpoly_zero(ab, aa->F);
    aa->F->size = window;
    return window;
}

/* Note that coefficients of F are written transposed */
void bm_io_write_one_F_coeff(bm_io_ptr aa, unsigned int deg)
{
    bmstatus_ptr bm = aa->bm;
    dims * d = bm->d;
    abdst_field ab = d->ab;
    unsigned int n = d->n;
    deg = deg % aa->F->alloc;
    // if (random_input_length) return;
    if (aa->ascii && !random_input_length) {
            for(unsigned int j = 0 ; j < n ; j++) {
                for(unsigned int i = 0 ; i < n ; i++) {
                    if (i) fprintf(aa->f, " ");
                    abfprint(ab, aa->f, matpoly_coeff(ab, aa->F, i, j, deg));
                }
                fprintf(aa->f, "\n");
            }
            fprintf(aa->f, "\n");
    } else {
        abelt tmp;
        abinit(ab, &tmp);
            for(unsigned int j = 0 ; j < n ; j++) {
                for(unsigned int i = 0 ; i < n ; i++) {
                    abset(ab, tmp, matpoly_coeff(ab, aa->F, i, j, deg));
                    size_t elemsize = abvec_elt_stride(ab, 1);
                    fwrite(tmp, elemsize, 1, aa->f);
                }
            }
        abclear(ab, &tmp);
    }
}

void bm_io_clear_one_F_coeff(bm_io_ptr aa, unsigned int deg)
{
    bmstatus_ptr bm = aa->bm;
    dims * d = bm->d;
    abdst_field ab = d->ab;
    unsigned int n = d->n;
    deg = deg % aa->F->alloc;
    for(unsigned int j = 0 ; j < n ; j++) {
        for(unsigned int i = 0 ; i < n ; i++) {
            abset_zero(ab, matpoly_coeff(ab, aa->F, i, j, deg));
        }
    }
}

void bm_io_compute_final_F(bm_io_ptr aa, bigmatpoly_ptr xpi, unsigned int * delta)/*{{{*/
{
    bmstatus_ptr bm = aa->bm;
    dims * d = bm->d;
    unsigned int m = d->m;
    unsigned int n = d->n;
    unsigned int b = m + n;
    abdst_field ab = d->ab;
    int rank;
    MPI_Comm_rank(bm->world, &rank);


    /* We are not interested by xpi->size, but really by the number of
     * coefficients for the columns which give solutions. */
    unsigned int maxdelta = get_max_delta_on_solutions(bm, delta);
    unsigned int pilen = maxdelta + 1 - aa->t0;

    if (!rank) printf("Final f(X)=f0(X)pi(X) has degree %u\n", maxdelta);

    /* Decide on the temp storage size */
    double avg = avg_matsize(ab, n, n, aa->ascii);
    unsigned int B = iceildiv(io_block_size, avg);
    if (!rank && !random_input_length)
        printf("Writing F by blocks of %u coefficients (%.1f bytes each)\n",
                B, avg);
    /* This is only temp storage ! */
    matpoly pi;
    matpoly_init(ab, pi, b, b, B);

    unsigned int window = aa->F->alloc;

    /* Which columns of F*pi will make the final generator ? */
    unsigned int * sols = malloc(n * sizeof(unsigned int));
    for(unsigned int j = 0, jj=0 ; j < m + n ; j++) {
        if (bm->lucky[j] <= 0)
            continue;
        sols[jj++]=j;
    }

    double tt0 = wct_seconds();
    double next_report_t = tt0 + 10;
    unsigned next_report_k = pilen / 100;

#if 0 /* {{{ save the code which reads pi forwards. */
    /* To be clear: it's untested, but probably a better base to work on
     * than the version reading backwards... Someday, we might fix the
     * way we write F. This will ease the job for mksol if we consider
     * forcing it to evaluate with Horner. (presently, after much turmoil
     * for writing F*pi reversed, doing so would require mksol to read
     * the coefficients reversed too...)
     */
    unsigned int kpi;
    for(unsigned int kpi0 = 0 ; kpi0 < pilen ; kpi0 += B) {
        matpoly_zero(ab, pi);
        pi->size = B;

        bigmatpoly_gather_mat_partial(ab, pi, xpi, kpi0, B);

        if (rank) continue;

        kpi = kpi0;

        for(unsigned int s = 0 ; s < B ; s++, kpi++) {
            /* coefficient of degree kpi-window is now complete */
            if (kpi >= window) {
                bm_io_write_one_F_coeff(aa, kpi - window);
                bm_io_clear_one_F_coeff(aa, kpi - window);
            }
            /* Now use our coefficient. This might tinker with
             * coefficients up to degree kpi-(window-1) in the file F */

            if (kpi > maxdelta + aa->t0 ) {
                /* this implies that we always have kF > delta[jpi],
                 * whence we expect a zero contribution */
                for(unsigned int i = 0 ; i < m + n ; i++) {
                    for(unsigned int jF = 0 ; jF < n ; jF++) {
                        unsigned int jpi = sols[jF];
                        absrc_elt src = matpoly_coeff_const(ab, pi, i, jpi, s);
                        ASSERT_ALWAYS(abis_zero(ab, src));
                    }
                }
                break;
            }

            for(unsigned int i = 0 ; i < n ; i++) {
                for(unsigned int jF = 0 ; jF < n ; jF++) {
                    unsigned int jpi = sols[jF];
                    unsigned int subtract = maxdelta - delta[jpi];
                    ASSERT(subtract < window);
                    if (kpi < subtract) continue;
                    unsigned int kF = kpi - subtract;
                    absrc_elt src = matpoly_coeff_const(ab, pi, i, jpi, s);
                    abdst_elt dst = matpoly_coeff(ab, aa->F, i, jF, kF % window);
                    ASSERT_ALWAYS(kF <= delta[jpi] || abis_zero(ab, src));
                    abadd(ab, dst, dst, src);
                }
            }
            for(unsigned int ipi = n ; ipi < m + n ; ipi++) {
                for(unsigned int jF = 0 ; jF < n ; jF++) {
                    unsigned int jpi = sols[jF];
                    unsigned int iF = aa->fdesc[ipi-n][1];
                    unsigned int offset = aa->t0 - aa->fdesc[ipi-n][0];
                    unsigned int subtract = maxdelta - delta[jpi] + offset;
                    ASSERT(subtract < window);
                    if (kpi < subtract) continue;
                    unsigned int kF = kpi - subtract;
                    abdst_elt dst = matpoly_coeff(ab, aa->F, iF, jF, kF % window);
                    absrc_elt src = matpoly_coeff_const(ab, pi, ipi, jpi, s);
                    ASSERT_ALWAYS(kF <= delta[jpi] || abis_zero(ab, src));
                    abadd(ab, dst, dst, src);
                }
            }
        }

        if (kpi != kpi0 + B)
            break;

        if (!rank && kpi > next_report_k) {
            double tt = wct_seconds();
            if (tt > next_report_t) {
                fprintf(stderr,
                        "Written %u coefficients (%.1f%%) in %.1f s\n",
                        kpi, 100.0 * kpi / pilen,
                        tt-tt0);
                next_report_t = tt + 10;
                next_report_k = kpi + pilen/ 100;
            }
        }
    }
    /* flush the pipe */
    for(unsigned int s = 0 ; s < window ; s++, kpi++)
        bm_io_write_one_F_coeff(aa, kpi - window);
#endif /* }}} */

    /* we need to read pi backwards. The number of coefficients in pi is
     * pilen = maxdelta + 1 - t0. Hence the first interesting index is
     * maxdelta - t0. However, for notational ease, we'll access
     * coefficients from index maxdelta downwards.
     */

    unsigned int kpi;
    /* index kpi1 is end fence for window */
    for(unsigned int kpi1 = iceildiv(maxdelta + 1, B) * B ; kpi1 > 0 ; kpi1 -= B) {
        unsigned int kpi0 = kpi1 - B;
        /* rare corner case */
        if (kpi0 >= pilen) continue;

        matpoly_zero(ab, pi);
        pi->size = B;

        bigmatpoly_gather_mat_partial(ab, pi, xpi, kpi0, MIN(B, pilen - kpi0));

        if (rank) continue;

        kpi = kpi1 - 1;

        for(unsigned int s = B ; s-- > 0 ; kpi--) {
            if (kpi + window == maxdelta) {
                /* Avoid writing zero coefficients. This can occur !
                 * Example:
                 * tt=(2*4*1200) mod 1009, a = (tt cat tt cat * tt[0..10])
                 */
                for(unsigned int j = 0 ; j < n ; j++) {
                    int z = 1;
                    for(unsigned int i = 0 ; z && i < n ; i++) {
                        absrc_elt src = matpoly_coeff_const(ab, aa->F, i, j, 0);
                        z = abis_zero(ab, src);
                    }

                    if (z) {
                        /* This is a bit ugly. Given that we're going to
                         * shift one column of F, we'll have a
                         * potentially deepedr write-back buffer. Columns
                         * which seemed to be ready still are, but they
                         * will now be said so only at the next step.
                         */
                        printf("Reduced solution column #%u from"
                                " delta=%u to delta=%u\n",
                                sols[j], delta[sols[j]], delta[sols[j]]-1);
                        window++;
                        matpoly_realloc(ab, aa->F, window);
                        delta[sols[j]]--;
                        /* shift this column */
                        for(unsigned int k = 1 ; k < window ; k++) {
                            matpoly_extract_column(ab,
                                    aa->F, j, k-1, aa->F, j, k);
                        }
                        matpoly_zero_column(ab, aa->F, j, window - 1);
                        break;
                    }
                }
            }
            /* coefficient of degree maxdelta-kpi-window is now complete */
            if (kpi + window <= maxdelta) {
                bm_io_write_one_F_coeff(aa, (maxdelta-kpi) - window);
                bm_io_clear_one_F_coeff(aa, (maxdelta-kpi) - window);
            }
            /* Now use our coefficient. This might tinker with
             * coefficients up to degree kpi-(window-1) in the file F */

            if (kpi > maxdelta + aa->t0 ) {
                /* this implies that we always have kF > delta[jpi],
                 * whence we expect a zero contribution */
                for(unsigned int i = 0 ; i < m + n ; i++) {
                    for(unsigned int jF = 0 ; jF < n ; jF++) {
                        unsigned int jpi = sols[jF];
                        absrc_elt src = matpoly_coeff_const(ab, pi, i, jpi, s);
                        ASSERT_ALWAYS(abis_zero(ab, src));
                    }
                }
                continue;
            }

            for(unsigned int i = 0 ; i < n ; i++) {
                for(unsigned int jF = 0 ; jF < n ; jF++) {
                    unsigned int jpi = sols[jF];
                    unsigned int subtract = maxdelta - delta[jpi];
                    ASSERT(subtract < window);
                    if ((maxdelta-kpi) < subtract) continue;
                    unsigned int kF = (maxdelta - kpi) - subtract;
                    absrc_elt src = matpoly_coeff_const(ab, pi, i, jpi, s);
                    abdst_elt dst = matpoly_coeff(ab, aa->F, i, jF, kF % window);
                    ASSERT_ALWAYS(kF <= delta[jpi] || abis_zero(ab, src));
                    abadd(ab, dst, dst, src);
                }
            }
            for(unsigned int ipi = n ; ipi < m + n ; ipi++) {
                for(unsigned int jF = 0 ; jF < n ; jF++) {
                    unsigned int jpi = sols[jF];
                    unsigned int iF = aa->fdesc[ipi-n][1];
                    unsigned int offset = aa->t0 - aa->fdesc[ipi-n][0];
                    unsigned int subtract = maxdelta - delta[jpi] + offset;
                    ASSERT(subtract < window);
                    if ((maxdelta-kpi) < subtract) continue;
                    unsigned int kF = (maxdelta-kpi) - subtract;
                    abdst_elt dst = matpoly_coeff(ab, aa->F, iF, jF, kF % window);
                    absrc_elt src = matpoly_coeff_const(ab, pi, ipi, jpi, s);
                    ASSERT_ALWAYS(kF <= delta[jpi] || abis_zero(ab, src));
                    abadd(ab, dst, dst, src);
                }
            }
        }

        if (!rank && kpi > next_report_k) {
            double tt = wct_seconds();
            if (tt > next_report_t) {
                fprintf(stderr,
                        "Written %u coefficients (%.1f%%) in %.1f s\n",
                        kpi, 100.0 * kpi / pilen,
                        tt-tt0);
                next_report_t = tt + 10;
                next_report_k = kpi + pilen/ 100;
            }
        }
    }
    /* flush the pipe */
    if (!rank) {
        for(unsigned int s = window ; s-- > 0 ; kpi--)
            bm_io_write_one_F_coeff(aa, maxdelta - kpi - window);
    }

    free(sols);
}/*}}}*/

/* read 1 coefficient into the sliding window of input coefficients of
 * the input series A. The io_window parameter controls the size of the
 * sliding window. There are in fact two behaviours:
 *  - io_window == 0: there is no sliding window, really, and the new
 *    coefficient is appended as the last coefficient of aa->A.
 *  - io_window > 0: there, we really have a sliding window. Coeffs
 *    occupy places in a circular fashion within the buffer.
 */
int bm_io_read1(bm_io_ptr aa, unsigned int io_window)/*{{{*/
{
    bmstatus_ptr bm = aa->bm;
    dims *d = bm->d;
    abdst_field ab = d->ab;
    unsigned int m = d->m;
    unsigned int n = d->n;
    matpoly_ptr A = aa->A;
    static unsigned int generated_random_coefficients = 0;

    unsigned int pos = aa->k;
    if (io_window) {
        pos = pos % io_window;
        ASSERT_ALWAYS(A->size == io_window);
    } else {
        ASSERT_ALWAYS(A->size == aa->k);
        if (aa->k >= A->alloc) {
            matpoly_realloc(aa->bm->d->ab, A, A->alloc + 1);
        }
        A->size++;
    }
    if (random_input_length) {
        if (generated_random_coefficients >= random_input_length)
            return 0;

        for (unsigned int i = 0; i < m ; i++) {
            for (unsigned int j = 0; j < n ; j++) {
                abdst_elt x = matpoly_coeff(ab, A, i, j, pos);
                abrandom(ab, x, rstate);
            }
        }
        generated_random_coefficients++;
        aa->k++;
        return 1;
    }

    for (unsigned int i = 0; i < m ; i++) {
        for (unsigned int j = 0; j < n ; j++) {
            abdst_elt x = matpoly_coeff(ab, A, i, j, pos);
            int rc;
            if (aa->ascii) {
                rc = abfscan(ab, aa->f, x);
                rc = rc == 1;
            } else {
                size_t elemsize = abvec_elt_stride(ab, 1);
                rc = fread(x, elemsize, 1, aa->f);
                rc = rc == 1;
                abnormalize(ab, x);
            }
            if (!rc) {
                if (i == 0 && j == 0) {
                    return 0;
                }
                fprintf(stderr,
                        "Parse error while reading coefficient (%d,%d,%d)%s\n",
                        i, j, 1 + aa->k,
                        aa->ascii ? "" : " (forgot --ascii?)");
                exit(EXIT_FAILURE);
            }
        }
    }
    aa->k++;
    return 1;
}/*}}}*/

void bm_io_init(bm_io_ptr aa, bmstatus_ptr bm, const char * input_file, const char * output_file, int ascii)/*{{{*/
{
    memset(aa, 0, sizeof(*aa));
    aa->bm = bm;
    aa->ascii = ascii;
    aa->input_file = input_file;
    aa->output_file = output_file;
}/*}}}*/

void bm_io_begin_read(bm_io_ptr aa)/*{{{*/
{
    bmstatus_ptr bm = aa->bm;
    dims *d = bm->d;
    abdst_field ab = d->ab;
    unsigned int m = d->m;
    unsigned int n = d->n;
    int rank;
    MPI_Comm_rank(aa->bm->world, &rank);
    if (rank) return;

    matpoly_init(ab, aa->A, m, n, 1);

    if (random_input_length)
        return;

    aa->f = fopen(aa->input_file, aa->ascii ? "r" : "rb");
    DIE_ERRNO_DIAG(aa->f == NULL, "fopen", aa->input_file);

    /* read and discard the first coefficient */
    if (!bm_io_read1(aa, 0)) {
        fprintf(stderr, "Read error from %s\n", aa->input_file);
        exit(EXIT_FAILURE);
    }
    /* discard ! */
    aa->A->size = aa->k = 0;
}/*}}}*/

void bm_io_end_read(bm_io_ptr aa)/*{{{*/
{
    int rank;
    MPI_Comm_rank(aa->bm->world, &rank);
    if (rank) return;
    matpoly_clear(aa->bm->d->ab, aa->A);
    if (random_input_length) return;
    fclose(aa->f);
    aa->f = NULL;
}/*}}}*/

void bm_io_begin_write(bm_io_ptr aa)/*{{{*/
{
    bmstatus_ptr bm = aa->bm;
    dims *d = bm->d;
    abdst_field ab = d->ab;
    unsigned int n = d->n;
    int rank;
    MPI_Comm_rank(aa->bm->world, &rank);
    if (rank) return;

    matpoly_init(ab, aa->F, n, n, aa->t0 + 1);

    if (random_input_length) {
        /* It's horrible. I really want "checksum" to be printed at the
         * exact right moment, hence the following hack...  */
        aa->f = popen("(echo \"checksum:\" ; sha1sum) | xargs echo", "w");
    } else {
        aa->f = fopen(aa->output_file, aa->ascii ? "w" : "wb");
    }
    DIE_ERRNO_DIAG(aa->f == NULL, "fopen", aa->output_file);
}/*}}}*/

void bm_io_end_write(bm_io_ptr aa)/*{{{*/
{
    int rank;
    MPI_Comm_rank(aa->bm->world, &rank);
    if (rank) return;
    matpoly_clear(aa->bm->d->ab, aa->F);
    if (random_input_length) {
        pclose(aa->f);
    } else {
        fclose(aa->f);
    }
    aa->f = NULL;
}/*}}}*/

void bm_io_clear(bm_io_ptr aa)/*{{{*/
{
    if (aa->fdesc) free(aa->fdesc);
    memset(aa, 0, sizeof(*aa));
}/*}}}*/

void bm_io_guess_length(bm_io_ptr aa)/*{{{*/
{
    bmstatus_ptr bm = aa->bm;
    dims * d = bm->d;
    unsigned int m = d->m;
    unsigned int n = d->n;
    abdst_field ab = d->ab;
    int rank;
    MPI_Comm_rank(aa->bm->world, &rank);

    if (random_input_length) {
        aa->guessed_length = random_input_length;
        return;
    }

    if (!rank) {
        struct stat sbuf[1];
        int rc = fstat(fileno(aa->f), sbuf);
        if (rc < 0) {
            perror(aa->input_file);
            exit(EXIT_FAILURE);
        }

        size_t filesize = sbuf->st_size;

        if (!aa->ascii) {
            size_t avg = avg_matsize(ab, m, n, aa->ascii);
            if (filesize % avg) {
                fprintf(stderr, "File %s has %zu bytes, while its size should be amultiple of %zu bytes (assuming binary input; perhaps --ascii is missing ?).\n", aa->input_file, filesize, avg);
                exit(EXIT_FAILURE);
            }
            aa->guessed_length = filesize / avg;
        } else {
            double avg = avg_matsize(ab, m, n, aa->ascii);
            double expected_length = filesize / avg;
            if (!rank)
                printf("Expect roughly %.2f items in the sequence.\n", expected_length);

            /* First coefficient is always lighter, so we add a +1. The
             * 5% are here really only to take into account the
             * deviations, but we don't expect much */
            size_t guessed_length = 1 + ceil(1.05 * expected_length);
            if (!rank)
                printf("With safety margin: expect length %zu at most\n", guessed_length);

            aa->guessed_length = guessed_length;
        }
#if 0
        /* we don't have the struct bw at hand here... */
        if (bw->end || bw->start) {
            printf("(Note: from bw parameters, we expect %u).\n",
                    bw->end - bw->start);
        }
        printf(".\n");
#endif
    }
    MPI_Bcast(&(aa->guessed_length), 1, MPI_UNSIGNED, 0, aa->bm->world);
}/*}}}*/

void bm_io_compute_initial_F(bm_io_ptr aa) /*{{{ */
{
    bmstatus_ptr bm = aa->bm;
    dims *d = bm->d;
    abdst_field ab = d->ab;
    unsigned int m = d->m;
    unsigned int n = d->n;
    matpoly_ptr A = aa->A;

    unsigned int (*fdesc)[2] = malloc(2 * m * sizeof(unsigned int));

    int rank;
    MPI_Comm_rank(aa->bm->world, &rank);
    if (!rank) {
        /* read the first few coefficients. Expand A accordingly as we are
         * doing the read */

        ASSERT(A->m == m);
        ASSERT(A->n == n);

        abelt tmp;
        abinit(ab, &tmp);

        /* First try to create the initial F matrix */
        printf("Computing t0\n");

        /* We want to create a full rank m*m matrix M, by extracting columns
         * from the first coefficients of A */

        matpoly M;
        matpoly_init(ab, M, m, m, 1);
        M->size = 1;

        /* For each integer i between 0 and m-1, we have a column, picked
         * from column cnum[i] of coeff exponent[i] of A which, once reduced modulo
         * the other ones, has coefficient at row pivots[i] unequal to zero.
         */
        unsigned int *pivots = malloc(m * sizeof(unsigned int));
        unsigned int *exponents = malloc(m * sizeof(unsigned int));
        unsigned int *cnum = malloc(m * sizeof(unsigned int));
        unsigned int r = 0;

        for (unsigned int k = 0; r < m ; k++) {
            /* read a new coefficient */
            bm_io_read1(aa, 0);

            for (unsigned int j = 0; r < m && j < n; j++) {
                /* Extract a full column into M */
                matpoly_extract_column(ab, M, r, 0, A, j, k);

                /* Now reduce it modulo all other columns */
                for (unsigned int v = 0; v < r; v++) {
                    unsigned int u = pivots[v];
                    /* the v-th column in the M is known to
                     * kill coefficient u (more exactly, to have a -1 as u-th
                     * coefficient, and zeroes for the other coefficients
                     * referenced in the pivots[0] to pivots[v-1] indices).
                     */
                    /* add M[u,r]*column v of M to column r of M */
                    for(unsigned int i = 0 ; i < m ; i++) {
                        if (i == u) continue;
                        abmul(ab, tmp,
                                  matpoly_coeff(ab, M, i, v, 0),
                                  matpoly_coeff(ab, M, u, r, 0));
                        abadd(ab, matpoly_coeff(ab, M, i, r, 0),
                                  matpoly_coeff(ab, M, i, r, 0),
                                  tmp);
                    }
                    abset_zero(ab,
                            matpoly_coeff(ab, M, u, r, 0));
                }
                unsigned int u = 0;
                for( ; u < m ; u++) {
                    if (abcmp_ui(ab, matpoly_coeff(ab, M, u, r, 0), 0) != 0)
                        break;
                }
                if (u == m) {
                    printf("[X^%d] A, col %d does not increase rank (still %d)\n",
                           k, j, r);

                    /* we need at least m columns to get as starting matrix
                     * with full rank. Given that we have n columns per
                     * coefficient, this means at least m/n matrices.
                     */

                    if (k * n > m + 40) {
                        printf("The choice of starting vectors was bad. "
                               "Cannot find %u independent cols within A\n", m);
                        exit(EXIT_FAILURE);
                    }
                    continue;
                }

                /* Bingo, it's a new independent col. */
                pivots[r] = u;
                cnum[r] = j;
                exponents[r] = k;

                /* Multiply the column so that the pivot becomes -1 */
                int rc = abinv(ab, tmp, matpoly_coeff(ab, M, u, r, 0));
                if (!rc) {
                    fprintf(stderr, "Error, found a factor of the modulus: ");
                    abfprint(ab, stderr, tmp);
                    fprintf(stderr, "\n");
                    exit(EXIT_FAILURE);
                }
                abneg(ab, tmp, tmp);
                for(unsigned int i = 0 ; i < m ; i++) {
                    abmul(ab, matpoly_coeff(ab, M, i, r, 0),
                              matpoly_coeff(ab, M, i, r, 0),
                              tmp);
                }

                r++;

                // if (r == m)
                    printf
                        ("[X^%d] A, col %d increases rank to %d (head row %d)\n",
                         aa->k, j, r, u);
            }
        }

        if (r != m) {
            printf("This amount of data is insufficient. "
                   "Cannot find %u independent cols within A\n", m);
            exit(EXIT_FAILURE);
        }

        aa->t0 = exponents[r - 1] + 1;
        ASSERT_ALWAYS(aa->t0 == aa->k);

        printf("Found satisfying init data for t0=%d\n", aa->t0);

        /* We've also got some adjustments to make: room for one extra
         * coefficient is needed in A. Reading of further coefficients will
         * pickup where they stopped, and will always leave the last t0+1
         * coefficients readable. */
        matpoly_realloc(ab, A, aa->t0 + 1);
        A->size++;


        for(unsigned int j = 0 ; j < m ; j++) {
            fdesc[j][0] = exponents[j];
            fdesc[j][1] = cnum[j];
            ASSERT_ALWAYS(exponents[j] < aa->t0);
        }
        free(pivots);
        free(exponents);
        free(cnum);
        matpoly_clear(ab, M);
        abclear(ab, &tmp);
    }
    aa->fdesc = fdesc;
    MPI_Bcast(aa->fdesc, 2*m, MPI_UNSIGNED, 0, aa->bm->world);
    MPI_Bcast(&(aa->t0), 1, MPI_UNSIGNED, 0, aa->bm->world);
    bm->t = aa->t0;
}				/*}}} */

void bm_io_compute_E(bm_io_ptr aa, bigmatpoly_ptr xE)/*{{{*/
{
    // F0 is exactly the n x n identity matrix, plus the
    // X^(s-exponent)e_{cnum} vectors. fdesc has the (exponent, cnum)
    // pairs
    bmstatus_ptr bm = aa->bm;
    dims * d = bm->d;
    unsigned int m = d->m;
    unsigned int n = d->n;
    unsigned int b = m + n;
    abdst_field ab = d->ab;
    int rank;
    MPI_Comm_rank(aa->bm->world, &rank);

    unsigned int guess = aa->guessed_length;
    ASSERT(bigmatpoly_check_pre_init(xE));
    bigmatpoly_finish_init(ab, xE, m, b, guess);

    /* Decide on the temp storage size */
    double avg = avg_matsize(ab, m, n, aa->ascii);
    unsigned int B = iceildiv(io_block_size, avg);
    if (!rank)
        printf("Reading A by blocks of %u coefficients (%.1f bytes each)\n",
                B, avg);
    /* This is only temp storage ! */
    matpoly E;
    matpoly_init(ab, E, m, b, B);

    unsigned int window = aa->t0 + 1;

    double tt0 = wct_seconds();
    double next_report_t = tt0 + 10;
    unsigned next_report_k = guess / 100;

    unsigned int kE;

    for(unsigned int kE0 = 0 ; kE0 + aa->t0 < aa->guessed_length ; kE0 += B) {
        /* Setting E->size is rather artificial, since we use E essentially
         * as an area where we may write coefficients freely. The only aim is
         * to escape some safety checks involving ->size in matpoly_part */
        matpoly_zero(ab, E);
        E->size = B;

        kE = kE0;

        if (!rank) {
            for(unsigned int s = 0 ; s < B ; s++, kE++) {
                if (!bm_io_read1(aa, window)) {
                    /* the +1 is because we're playing a nasty trick by
                     * forcibly discarding the first coefficient in the
                     * sequence. Just don't convey bad information... */
                    fprintf(stderr, "EOF met after reading %u coefficients\n",
                            aa->k + 1);
                    break;
                }

                if (kE + aa->t0 > aa->guessed_length) {
                    fprintf(stderr, "Going past guessed length ???\n");
                }

                for(unsigned int j = 0 ; j < n ; j++) {
                    unsigned int kA = kE + aa->t0;
                    matpoly_extract_column(ab, E, j, s, aa->A, j, kA % window);
                }
                for(unsigned int jE = n ; jE < m + n ; jE++) {
                    unsigned int jA = aa->fdesc[jE-n][1];
                    unsigned int offset = aa->fdesc[jE-n][0];
                    unsigned int kA = kE + offset;
                    matpoly_extract_column(ab, E, jE, s, aa->A, jA, kA % window);
                }
            }
        }
        MPI_Bcast(&(aa->k), 1, MPI_UNSIGNED, 0, aa->bm->world);
        MPI_Bcast(&(kE), 1, MPI_UNSIGNED, 0, aa->bm->world);
        E->size = kE - kE0;

        bigmatpoly_scatter_mat_partial(ab, xE, E, kE0, kE - kE0);

        if (!rank) {
            ASSERT_ALWAYS(kE + aa->t0 == aa->k);
            if (aa->k > next_report_k) {
                double tt = wct_seconds();
                if (tt > next_report_t) {
                    fprintf(stderr,
                            "Read %u coefficients (%.1f%%)"
                            " in %.1f s (%.1f MB/s)\n",
                            aa->k, 100.0 * aa->k / guess,
                            tt-tt0, aa->k * avg / (tt-tt0)/1.0e6);
                    next_report_t = tt + 10;
                    next_report_k = aa->k + guess/ 100;
                }
            }
        }
        if (kE < kE0 + B)
            break;
    }
    bigmatpoly_set_size(xE, kE);
    matpoly_clear(ab, E);
}/*}}}*/

void usage()
{
    fprintf(stderr, "Usage: ./plingen [options, to be documented]\n");
    fprintf(stderr,
	    "General information about bwc options is in the README file\n");
    exit(EXIT_FAILURE);
}

void bmstatus_init(bmstatus_ptr bm, unsigned int m, unsigned int n)/*{{{*/
{
    /* Fields will not all be significant right now. The initialization
     * is a long process */
    memset(bm, 0, sizeof(bmstatus));
    bm->d->m = m;
    bm->d->n = n;
    bm->lucky = malloc((m + n) * sizeof(int));
    memset(bm->lucky, 0, (m + n) * sizeof(int));
}/*}}}*/

void bmstatus_clear(bmstatus_ptr bm)/*{{{*/
{
    free(bm->lucky);
    memset(bm, 0, sizeof(bmstatus));
}/*}}}*/


unsigned int count_lucky_columns(bmstatus_ptr bm)/*{{{*/
{
    dims * d = bm->d;
    unsigned int m = d->m;
    unsigned int n = d->n;
    unsigned int b = m + n;
    int luck_mini = expected_pi_length(d, 0);
    MPI_Bcast(bm->lucky, b, MPI_UNSIGNED, 0, bm->world);
    unsigned int nlucky = 0;
    for(unsigned int j = 0 ; j < b ; nlucky += bm->lucky[j++] >= luck_mini) ;
    return nlucky;
}/*}}}*/

void check_luck_condition(bmstatus_ptr bm)/*{{{*/
{
    dims * d = bm->d;
    unsigned int m = d->m;
    unsigned int n = d->n;
    unsigned int nlucky = count_lucky_columns(bm);

    int rank;
    MPI_Comm_rank(bm->world, &rank);

    if (nlucky == n)
        return;

    if (!rank) {
        fprintf(stderr, "Could not find the required set of solutions (nlucky=%u)\n", nlucky);
    }
    if (random_input_length) {
        static int once=0;
        if (once++) {
            if (!rank) {
                fprintf(stderr, "Solution-faking loop crashed\n");
            }
            MPI_Abort(bm->world, EXIT_FAILURE);
        }
        if (!rank) {
            printf("Random input: faking successful computation\n");
        }
        for(unsigned int j = 0 ; j < n ; j++) {
            bm->lucky[(j * 1009) % (m+n)] = expected_pi_length(d, 0);
        }
        check_luck_condition(bm);
        return;
    }

    MPI_Abort(bm->world, EXIT_FAILURE);
}/*}}}*/

void display_deltas(bmstatus_ptr bm, unsigned int * delta)/*{{{*/
{
    dims * d = bm->d;
    unsigned int m = d->m;
    unsigned int n = d->n;

    int rank;
    MPI_Comm_rank(bm->world, &rank);

    if (!rank) {
        printf("Final, t=%u: delta =", bm->t);
        for(unsigned int j = 0; j < m + n; j++) {
            printf(" %u", delta[j]);
            if (bm->lucky[j] < 0) {
                printf("(*)");
            }
        }
        printf("\n");
    }
}/*}}}*/


int main(int argc, char *argv[])
{
    bmstatus bm;
    dims * d = bm->d;
    int tune = 0;
    int ascii = 0;

    setvbuf(stdout, NULL, _IONBF, 0);
    setvbuf(stderr, NULL, _IONBF, 0);

    gmp_randinit_default(rstate);

    /* {{{ Parameter list. Feed the common bw[] struct */
    param_list pl;
    param_list_init(pl);

    param_list_configure_switch(pl, "--tune", &tune);
    param_list_configure_switch(pl, "--ascii", &ascii);
    param_list_configure_switch(pl, "--timings", &with_timings);
    bw_common_init_mpi(bw, pl, &argc, &argv);

    int rank;
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    info_init_timer();

    param_list_parse_uint(pl, "random-input-with-length", &random_input_length);
    param_list_parse_int(pl, "caching", &caching);

    const char * afile = param_list_lookup_string(pl, "afile");

    if (bw->m == -1) {
	fprintf(stderr, "no m value set\n");
	exit(EXIT_FAILURE);
    }
    if (bw->n == -1) {
	fprintf(stderr, "no n value set\n");
	exit(EXIT_FAILURE);
    }
    if (!tune && !(afile || random_input_length)) {
        fprintf(stderr, "No afile provided\n");
        exit(EXIT_FAILURE);
    }
    if (!param_list_lookup_string(pl, "prime")) {
	fprintf(stderr, "no prime set\n");
	exit(EXIT_FAILURE);
    }

    /* we allow ffile and ffile to be both NULL */
    const char * tmp = param_list_lookup_string(pl, "ffile");
    char * ffile = NULL;
    if (tmp) {
        ffile = strdup(tmp);
    } else if (afile) {
        int rc = asprintf(&ffile, "%s.gen", afile);
        ASSERT_ALWAYS(rc >= 0);
    }
    ASSERT_ALWAYS((afile==NULL) == (ffile == NULL));

#ifndef HAVE_MPIR
    if (caching) {
        fprintf(stderr, "--caching=1 only supported with MPIR\n");
        exit(EXIT_FAILURE);
    }
#endif

    bmstatus_init(bm, bw->m, bw->n);

    unsigned int m = d->m;
    unsigned int n = d->n;
    abdst_field ab = d->ab;

    abfield_init(ab);
    {
	mpz_t p;
	mpz_init_set_ui(p, 2);
	param_list_parse_mpz(pl, "prime", p);
	abfield_specify(ab, MPFQ_PRIME_MPZ, p);
	mpz_clear(p);
    }
    abmpi_ops_init(ab);

    bm->lingen_threshold = 10;
    bm->lingen_mpi_threshold = 1000;
    param_list_parse_uint(pl, "lingen-threshold", &(bm->lingen_threshold));
    param_list_parse_uint(pl, "display-threshold", &(display_threshold));
    param_list_parse_uint(pl, "lingen-mpi-threshold", &(bm->lingen_mpi_threshold));
    param_list_parse_uint(pl, "io-block-size", &(io_block_size));
    {
        unsigned long random_seed = 0;
        param_list_parse_ulong(pl, "random_seed", &random_seed);
        gmp_randseed_ui(rstate, random_seed);
    }


#if defined(FAKEMPI_H_)
    bm->lingen_mpi_threshold = UINT_MAX;
#endif

    /* }}} */

    /* {{{ Parse MPI args. Make bm->world a better mpi communicator */
    bm->mpi_dims[0] = 1;
    bm->mpi_dims[1] = 1;
    param_list_parse_intxint(pl, "mpi", bm->mpi_dims);
    {
        /* Reorder all mpi nodes so that each node gets the given number
         * of jobs, but close together.
         */
        int mpi[2] = { bm->mpi_dims[0], bm->mpi_dims[1], };
#ifdef  FAKEMPI_H_
        if (mpi[0]*mpi[1] > 1) {
            fprintf(stderr, "non-trivial option mpi= can't be used with fakempi. Please do an MPI-enabled build (MPI=1)\n");
            exit(EXIT_FAILURE);
        }
#endif
        int thr[2] = {1,1};
        param_list_parse_intxint(pl, "thr", thr);
        /* Asssume numbering in MPI_COMM_WORLD is all mpi jobs from the
         * same node together. So we pick them by bunches of size
         * thr[0]*thr[1].
         */
        if (!rank)
            printf("size=%d mpi=%dx%d thr=%dx%d\n", size, mpi[0], mpi[1], thr[0], thr[1]);
        ASSERT_ALWAYS(size == mpi[0] * mpi[1] * thr[0] * thr[1]);
        /* The mpi and thr command line argument lead to the same number
         * of working processes than with krylov. Except that here, we're
         * not really doing threads, but running real mpi jobs with their
         * own address space (hence mpirun will need a bigger -n
         * argument).  In a sense, here we strive to follow the semantics
         * used in bwc otherwise.  In fact, we should probably use the
         * pi_wiring structures as well, so that we can say that we
         * consistently use the same infrastructure.  For simplicity of
         * the code, we don't.
         */
        bm->mpi_dims[0] *= thr[0];
        bm->mpi_dims[1] *= thr[1];
        if (bm->mpi_dims[0] != bm->mpi_dims[1]) {
            if (!rank)
                fprintf(stderr, "The current plingen code is limited to square splits ; here, we received a %d x %d split, which will not work\n",
                    bm->mpi_dims[0], bm->mpi_dims[1]);
            abort();
        } else if ((m % bm->mpi_dims[0] != 0) || (n % bm->mpi_dims[0] != 0)) {
            if (!rank)
                fprintf(stderr, "The process grid dimensions must divide gcd(m,n)\n");
            abort();
        }
        int tl_grank = rank % (thr[0] * thr[1]); // thread-level global rank
        int tl_irank = tl_grank / thr[1];
        int tl_jrank = tl_grank % thr[1];
        int ml_grank = rank / (thr[0] * thr[1]); // thread-level global rank
        int ml_irank = ml_grank / mpi[1];
        int ml_jrank = ml_grank % mpi[1];
        int irank = ml_irank * thr[0] + tl_irank;
        int jrank = ml_jrank * thr[1] + tl_jrank;
        int newrank = irank * mpi[1] * thr[1] + jrank;
        MPI_Comm_split(MPI_COMM_WORLD, 0, newrank, &(bm->world));
    }
    /* }}} */


    /* plingen tuning accepts some arguments. We look them up so as to
     * avoid failures down the line */
    param_list_lookup_string(pl, "B");
    param_list_lookup_string(pl, "catchsig");
    

    if (param_list_warn_unused(pl))
	usage();

    if (tune) {
        plingen_tuning(bm->d->ab, bm->d->m, bm->d->n, bm->world, pl);
        MPI_Finalize();
        return 0;
    }


    /* We now have a protected structure for a_reading task which does
     * the right thing concerning parallelism among MPI nodes (meaning
     * that non-root nodes essentially do nothing while the master job
     * does the I/O stuff) */
    bm_io aa;
    bm_io_init(aa, bm, afile, ffile, ascii);
    bm_io_begin_read(aa);
    bm_io_guess_length(aa);

    bm_io_compute_initial_F(aa);

    /* this is somewhat ugly, too */
    unsigned int t0 = bm->t;
    unsigned int * delta = malloc((m + n) * sizeof(unsigned int));
    for(unsigned int j = 0 ; j < m + n ; delta[j++]=t0);

    int go_mpi = aa->guessed_length >= bm->lingen_mpi_threshold;

    if (go_mpi && !rank) {
        if (size) {
            printf("Expected length %u exceeds MPI threshold %u, going MPI now.\n", aa->guessed_length, bm->lingen_mpi_threshold);
        } else {
            printf("Expected length %u exceeds MPI threshold %u, but the process is not running  in an MPI context.\n", aa->guessed_length, bm->lingen_mpi_threshold);
        }
    }

    bigmatpoly model;
    bigmatpoly xpi, xE;
    bigmatpoly_init_model(model, bm->world, bm->mpi_dims[0], bm->mpi_dims[1]);
    bigmatpoly_init(ab, xE, model, 0, 0, 0);
    bigmatpoly_init(ab, xpi, model, 0, 0, 0);   /* pre-init for now */

    /* This is a dispatching read */
    bm_io_compute_E(aa, xE);
    bm_io_end_read(aa);

    /* If we're running over MPI, we have done a dispatching read of E.
     * Probably there's a good reason for it, since we expect our
     * mpi_threshold to be below the size of E. In the degenerate case
     * where we're running over MPI while the mpi threshold is too large,
     * then: 1) it's silly 2) we'll run the MPI code still, since this is
     * the one which is capable of running gather & scatter. There will
     * be copies, we somehow we asked for it.
     */
    if (size) {
        bw_biglingen_collective(bm, xpi, xE, delta);
    } else {
        /* When communicator has size 1 (hence no mpi), we may save a few
         * copies by taking a shortcut */
        bw_lingen_single(bm, bigmatpoly_my_cell(xpi), bigmatpoly_my_cell(xE), delta);
        bigmatpoly_set_size(xpi, bigmatpoly_my_cell(xpi)->size);
    }

    /* This possibly aborts */
    check_luck_condition(bm);
    display_deltas(bm, delta);
    if (!rank) printf("(pi->alloc = %zu)\n", bigmatpoly_my_cell(xpi)->alloc);

    /* and this is a gathering write */
    bm_io_begin_write(aa);
    /* this reallocates F */

    bm_io_set_write_behind_size(aa, delta);
    bm_io_compute_final_F(aa, xpi, delta);
    bm_io_end_write(aa);

    if (!rank) {
        printf("t_basecase = %.2f\n", bm->t_basecase);
        printf("t_mp = %.2f\n", bm->t_mp);
        printf("t_mul = %.2f\n", bm->t_mul);
    }

    /* clear everything */
    bigmatpoly_clear(ab, xE);
    bigmatpoly_clear(ab, xpi);
    bigmatpoly_clear_model(model);
    bm_io_clear(aa);
    free(delta);
    abmpi_ops_clear(ab);
    abfield_clear(ab);
    bmstatus_clear(bm);
    bw_common_clear(bw);
    param_list_clear(pl);
    if (ffile) free(ffile);

    gmp_randclear(rstate);

    MPI_Finalize();
    return 0;
}

/* vim:set sw=4 sta et: */
