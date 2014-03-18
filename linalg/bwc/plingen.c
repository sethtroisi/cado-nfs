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
 * bw_biglingen_single  [complete data on master node]
 * |->bw_biglingen_collective  when need to go collective
 *     |<--> bw_lingen_bigrecursive , loops to bw_biglingen_collective
 *     \->bw_biglingen_single     when it makes sense to do this locally again
 * \->bw_lingen                   when computation can be done locally
 *    |<->bw_lingen_recursive
 *    |->bw_lingen_basecase
 *
 */
static unsigned int display_threshold = 10;
static int with_timings = 0;

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
static int bw_lingen(bmstatus_ptr bm, matpoly pi, matpoly E, unsigned int *delta);
static int bw_biglingen_collective(bmstatus_ptr bm, bigmatpoly pi, bigmatpoly E, unsigned int *delta);

static int bw_biglingen_single(bmstatus_ptr bm, matpoly pi, matpoly E, unsigned int *delta);


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

static const char *size_disp(size_t s, char buf[16])
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
}

void print_info_mp(unsigned int t, matpoly_ptr A, matpoly_ptr B)
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
}

void print_info_mul(unsigned int t, matpoly_ptr A, matpoly_ptr B)
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
}

void print_info_mpi_mp(unsigned int t, bigmatpoly_ptr A, bigmatpoly_ptr B)
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
}

void print_info_mpi_mul(unsigned int t, bigmatpoly_ptr A, bigmatpoly_ptr B)
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
}

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
    done = bw_lingen(bm, pi_left, E_left, delta);
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
    done = bw_lingen(bm, pi_right, E_right, delta);
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
static int bw_lingen_bigrecursive(bmstatus_ptr bm, bigmatpoly pi, bigmatpoly E, unsigned int *delta) /*{{{*/
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

static int bw_lingen(bmstatus_ptr bm, matpoly pi, matpoly E, unsigned int *delta) /*{{{*/
{
    int rank;
    MPI_Comm_rank(bm->world, &rank);
    ASSERT_ALWAYS(rank == 0);
    ASSERT_ALWAYS(E->size < bm->lingen_mpi_threshold);

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

    // fprintf(stderr, "Enter %s\n", __func__);
    if (E->size >= bm->lingen_mpi_threshold)  {
        done = bw_lingen_bigrecursive(bm, pi, E, delta);
    } else {
        /* Fall back to local code */
        /* This entails gathering E locally, computing pi locally, and
         * dispathing it back. */

        matpoly sE, spi;
        matpoly_init(ab, sE, m, b, E->size);
        matpoly_init(ab, spi, 0, 0, 0);
        bigmatpoly_gather_mat(ab, sE, E);
        /* Only the master node does the local computation */
        done = bw_biglingen_single(bm, spi, sE, delta);
        bigmatpoly_scatter_mat(ab, pi, spi);
        /* Don't forget to broadcast delta from root node to others ! */
        matpoly_clear(ab, spi);
        matpoly_clear(ab, sE);
    }
    // fprintf(stderr, "Leave %s\n", __func__);

    return done;
}/*}}}*/

static int bw_biglingen_single(bmstatus_ptr bm, matpoly pi, matpoly E, unsigned int *delta)/*{{{*/
{
    /* This version of the code is called collectively from all nodes,
     * but with the input E on rank 0. This normally happens _only_ at
     * the start. */
    dims *d = bm->d;
    unsigned int m = d->m;
    unsigned int n = d->n;
    unsigned int b = m + n;
    int done;
    int go_mpi = 0;
    go_mpi = E->size >= bm->lingen_mpi_threshold;

    // fprintf(stderr, "Enter %s\n", __func__);
    if (!go_mpi) {
        int rank;
        MPI_Comm_rank(bm->world, &rank);
        if (rank == 0)
            done = bw_lingen(bm, pi, E, delta);
        MPI_Bcast(&done, 1, MPI_INT, 0, bm->world);
        MPI_Bcast(delta, b, MPI_UNSIGNED, 0, bm->world);
        MPI_Bcast(bm->lucky, b, MPI_UNSIGNED, 0, bm->world);
        MPI_Bcast(&(bm->t), 1, MPI_UNSIGNED, 0, bm->world);
    } else {
        /* We are going to do this collectively */
        bigmatpoly model;
        abdst_field ab = d->ab;
        bigmatpoly xpi, xE;

        bigmatpoly_init_model(model, bm->world, bm->mpi_dims[0], bm->mpi_dims[1]);
        // ???
        /* We prefer to allocate soon. The interface doesn't really like
         * lazy allocation at the moment */
        // bigmatpoly_init(ab, xE, model, m, b, E->size);
        bigmatpoly_init(ab, xE, model, 0, 0, 0);
        bigmatpoly_init(ab, xpi, model, 0, 0, 0);   /* pre-init for now */
        bigmatpoly_scatter_mat_alt(ab, xE, E);
        ASSERT_ALWAYS(xE->size);
        done = bw_biglingen_collective(bm, xpi, xE, delta);
        bigmatpoly_gather_mat(ab, pi, xpi);
        bigmatpoly_clear(ab, xE);
        bigmatpoly_clear(ab, xpi);
        bigmatpoly_clear_model(model);
    }
    // fprintf(stderr, "Leave %s\n", __func__);
    return done;
}/*}}}*/

/**********************************************************************/

unsigned int get_max_delta_on_solutions(bmstatus_ptr bm, unsigned int * delta)/*{{{*/
{
    dims * d = bm->d;
    unsigned int m = d->m;
    unsigned int n = d->n;
    unsigned int maxdelta = 0;
    for(unsigned int j = 0 ; j < m + n ; j++) {
        if (bm->lucky[j] <= 0)
            continue;
        if (delta[j] > maxdelta)
            maxdelta = delta[j];
    }
    return maxdelta;
}/*}}}*/

void compute_final_F_red(bmstatus_ptr bm, matpoly f, unsigned int (*fdesc)[2], unsigned int t0, matpoly pi, unsigned int * delta)/*{{{*/
{
    dims * d = bm->d;
    unsigned int m = d->m;
    unsigned int n = d->n;
    unsigned int b = m + n;
    abdst_field ab = d->ab;

    unsigned int maxdelta = get_max_delta_on_solutions(bm, delta);

    // F0 is exactly the n x n identity matrix, plus the
    // X^(t0-exponent)e_{cnum} vectors. fdesc has the (exponent, cnum)
    // pairs.
    // Note that the degree of F0 is t0, whenceforth its length is t0+1.
    // Our argument really is t0.
    // Note also that the degree of the result is maxdelta, hence the
    // length is maxdelta + 1. Coefficients are read from pi as long as
    // the condition k+t0 <= maxdelta is satisfied. Which means length
    // pilen as below.
    unsigned int flen = maxdelta + 1;
    unsigned int pilen = maxdelta + 1 - t0;

    /* pi->size is not necessarily correct. I think it's always at least
     * the value pilen computed above, but I can conceive that it be less
     * (it's not clear, as the deltas give nominal degrees). In any case,
     * having pi->size somewhat larger than pilen is entirely possible,
     * because once the solution columns are fixed, the rest of pi has
     * its columns gaining degrees at each step.
     */

    ASSERT(f->m == 0);
    ASSERT(f->n == 0);
    ASSERT(f->alloc == 0);

    matpoly_init(ab, f, n, n, flen);
    f->size = flen;

    printf("Computing value of f(X)=f0(X)pi(X) (degree %u)\n", maxdelta);
    printf("Final, t=%u: delta =", bm->t);
    for(unsigned int j = 0; j < b; j++) {
        printf(" %u", delta[j]);
        if (bm->lucky[j] < 0) {
            printf("(*)");
        }
    }
    printf(", (alloc=%zu)", pi->alloc);
    printf("\n");

    unsigned int * pi_colidx = malloc(n * sizeof(unsigned int));
    unsigned int nlucky = 0;
    for(unsigned int j = 0 ; j < b ; j++) {
        if (bm->lucky[j] > 0) {
            pi_colidx[nlucky++] = j;
        }
    }
    ASSERT_ALWAYS(nlucky == n);

    /* First treat the part of F0 which is the identity matrix. */
    for(unsigned int i = 0 ; i < n ; i++) {
        /* We can't do a long copy because we have holes */
        for(unsigned int k = 0 ; k < pilen ; k++) {
            for(unsigned int j = 0 ; j < n ; j++) {
                unsigned int j1 = pi_colidx[j];
                abset(ab,
                        matpoly_coeff(ab, f, i, j, k),
                        matpoly_coeff(ab, pi, i, j1, k));
            }
        }
    }
    
    /* Then add some parts of pi */
    for(unsigned int i = 0 ; i < m ; i++) {
        unsigned int c = fdesc[i][1];
        unsigned int e = fdesc[i][0];
        /* Add X^(t0-e) times row i+n of pi to row c of f */
        for(unsigned int k = 0 ; k < pilen ; k++) {
            for(unsigned int j = 0 ; j < n ; j++) {
                unsigned int j1 = pi_colidx[j];
                abadd(ab,
                        matpoly_coeff(ab, f, c, j, k+t0-e),
                        matpoly_coeff(ab, f, c, j, k+t0-e),
                        matpoly_coeff(ab, pi, i+n, j1, k));
            }
        }
    }

    free(pi_colidx);
}/*}}}*/

void write_f(bmstatus_ptr bm, const char * filename, matpoly f_red, unsigned int * delta, int ascii)/*{{{*/
{
    dims * d = bm->d;
    unsigned int m = d->m;
    unsigned int n = d->n;
    abdst_field ab = d->ab;
    FILE * f = fopen(filename, ascii ? "w" : "wb");
    DIE_ERRNO_DIAG(f == NULL, "fopen", filename);
    unsigned int maxdelta = get_max_delta_on_solutions(bm, delta);
    unsigned int flen = maxdelta + 1;
    unsigned int * sols = malloc(n * sizeof(unsigned int));
    for(unsigned int j = 0, jj=0 ; j < m + n ; j++) {
        if (bm->lucky[j] <= 0)
            continue;
        sols[jj++]=j;
    }
    /* Apparently it occurs that F can have a zero coefficient in
     * delta[j], or at least some columns of F. This seems to have a
     * slight potential of missing part of the solution space. Hence we
     * check whether, for each column, the coefficients of degree
     * delta[j] in F are zero or not. If they are all zero, we reduce
     * delta[j] by one unit for the purpose of this function, as this
     * will lead to the same sequence of vectors being computed in the
     * end. Whether or not this is the _real_ delta[j], corresponding to
     * the right-hand side having zero coefficients for degrees
     * [delta[j]..t[, is out of our concern.
     */
    for(unsigned int jj = 0 ; jj < n ; jj++) {
        unsigned int j = sols[jj];
        unsigned int delta_orig = delta[j];
        for(int z = 1; z && delta[j] ; delta[j]-=z) {
            for(unsigned int i = 0 ; z && i < n ; i++) {
                z = abis_zero(ab, matpoly_coeff(ab, f_red, i, jj, delta[j]));
            }
        }
        if (delta_orig > delta[j]) {
            printf("Reduced solution column #%u from delta=%u to delta=%u\n",
                    j, delta_orig, delta[j]);
        }
    }

    /* Do as we do in lingen-binary.cpp: transpose the solutions */
    if (ascii) {
        for(unsigned int k = 0 ; k < flen ; k++) {
            for(unsigned int jj = 0 ; jj < n ; jj++) {
                for(unsigned int i = 0 ; i < n ; i++) {
                    if (i) fprintf(f, " ");
                    unsigned int j = sols[jj];
                    if (k <= delta[j]) {
                        abfprint(ab, f, matpoly_coeff(ab, f_red, i, jj, delta[j]-k));
                    } else {
                        fprintf(f, "0");
                    }
                }
                fprintf(f, "\n");
            }
            fprintf(f, "\n");
        }
    } else {
        abelt tmp;
        abinit(ab, &tmp);
        for(unsigned int k = 0 ; k < flen ; k++) {
            for(unsigned int jj = 0 ; jj < n ; jj++) {
                unsigned int j = sols[jj];
                for(unsigned int i = 0 ; i < n ; i++) {
                    abset_zero(ab, tmp);
                    if (k <= delta[j])
                        abset(ab, tmp, matpoly_coeff(ab, f_red, i, jj, delta[j]-k));
                    size_t elemsize = abvec_elt_stride(ab, 1);
                    fwrite(tmp, elemsize, 1, f);
                }
            }
        }
        abclear(ab, &tmp);
    }
    fclose(f);
}/*}}}*/

struct a_reading_task_s {/*{{{*/
    bmstatus_ptr bm;
    unsigned int t0;
    FILE * f;
    const char * input_file;
    int ascii;
    /* This is only a rolling window ! */
    matpoly A;
    unsigned int (*fdesc)[2];
    /* This k is the coefficient in A(X) div X of the next coefficient to
     * be read. This is thus the total number of coefficients of A(X) div
     * X which have been read so far */
    unsigned int k;

    unsigned int guessed_length;
    double avg_matsize;
};

typedef struct a_reading_task_s a_reading_task[1];
typedef struct a_reading_task_s * a_reading_task_ptr;/*}}}*/

int a_reading_task_read1(a_reading_task_ptr aa, unsigned int modulus)/*{{{*/
{
    bmstatus_ptr bm = aa->bm;
    dims *d = bm->d;
    abdst_field ab = d->ab;
    unsigned int m = d->m;
    unsigned int n = d->n;
    matpoly_ptr A = aa->A;

    unsigned int pos = aa->k;
    if (modulus) {
        pos = pos % modulus;
        ASSERT_ALWAYS(A->size == modulus);
    } else {
        ASSERT_ALWAYS(A->size == aa->k);
        if (aa->k >= A->alloc) {
            matpoly_realloc(aa->bm->d->ab, A, A->alloc + 1);
        }
        A->size++;
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

void a_reading_task_init(a_reading_task_ptr aa, bmstatus_ptr bm, const char * input_file, int ascii)/*{{{*/
{
    memset(aa, 0, sizeof(*aa));
    aa->bm = bm;
    aa->ascii = ascii;
    aa->input_file = input_file;
    dims * d = bm->d;
    unsigned int m = d->m;
    unsigned int n = d->n;
    abdst_field ab = d->ab;

    aa->f = fopen(input_file, ascii ? "r" : "rb");
    DIE_ERRNO_DIAG(aa->f == NULL, "fopen", input_file);

    matpoly_init(ab, aa->A, m, n, 1);
    /* read and discard the first coefficient */
    if (!a_reading_task_read1(aa, 0)) {
        fprintf(stderr, "Read error from %s\n", input_file);
        exit(EXIT_FAILURE);
    }
    /* discard ! */
    aa->A->size = aa->k = 0;
}/*}}}*/

void a_reading_task_clear(a_reading_task_ptr aa)/*{{{*/
{
    fclose(aa->f);
    /* caller code actually steals the pointer from us (and replaces it
     * by a NULL */
    if (aa->fdesc) free(aa->fdesc);
    matpoly_clear(aa->bm->d->ab, aa->A);
}/*}}}*/

void a_reading_task_guess_length(a_reading_task_ptr aa)/*{{{*/
{
    bmstatus_ptr bm = aa->bm;
    dims * d = bm->d;
    unsigned int m = d->m;
    unsigned int n = d->n;
    abdst_field ab = d->ab;

    struct stat sbuf[1];
    int rc = fstat(fileno(aa->f), sbuf);
    if (rc < 0) {
        perror(aa->input_file);
        exit(EXIT_FAILURE);
    }

    size_t filesize = sbuf->st_size;

    /* Easy case first. If we have binary input, then we know a priori
     * that the input data must have size a multiple of the element size.
     */
    if (!aa->ascii) {
        size_t elemsize = abvec_elt_stride(ab, 1);
        size_t matsize = elemsize * m * n;
        if (filesize % matsize) {
            fprintf(stderr, "File %s has %zu bytes, while its size should be amultiple of %zu bytes (assuming binary input; perhaps --ascii is missing ?).\n", aa->input_file, filesize, matsize);
            exit(EXIT_FAILURE);
        }
        aa->avg_matsize = matsize;
        aa->guessed_length = filesize / matsize;
        return;
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
    {
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
        printf("Expect roughly %.2f decimal digits for integers mod p.\n", avg);
    }

    double matsize = (avg + 1) * m * n + 1;
    printf("Expect roughly %.2f bytes for each sequence matrix.\n", matsize);

    double expected_length = filesize / matsize;
    printf("Expect roughly %.2f items in the sequence.\n", expected_length);

    /* First coefficient is always lighter, so we add a +1. The 5% are
     * here really only to take into account the deviations, but we don't
     * expect much */
    size_t guessed_length = 1 + ceil(1.05 * expected_length);
    printf("With safety margin: expect length %zu at most\n", guessed_length);

    aa->avg_matsize = matsize;
    aa->guessed_length = guessed_length;

#if 0
    /* we don't have the struct bw at hand here... */
    if (bw->end || bw->start) {
        printf("(Note: from bw parameters, we expect %u).\n",
                bw->end - bw->start);
    }
    printf(".\n");
#endif
}/*}}}*/

void a_reading_task_compute_F(a_reading_task_ptr aa) /*{{{ */
{
    bmstatus_ptr bm = aa->bm;
    dims *d = bm->d;
    abdst_field ab = d->ab;
    unsigned int m = d->m;
    unsigned int n = d->n;
    matpoly_ptr A = aa->A;

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
        a_reading_task_read1(aa, 0);

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
 
    bm->t = aa->t0;

    /* We've also got some adjustments to make: room for one extra
     * coefficient is needed in A. Reading of further coefficients will
     * pickup where they stopped, and will always leave the last t0+1
     * coefficients readable. */
    matpoly_realloc(ab, A, aa->t0 + 1);
    A->size++;
    

    unsigned int (*fdesc)[2] = malloc(2 * m * sizeof(unsigned int));
    for(unsigned int j = 0 ; j < m ; j++) {
        fdesc[j][0] = exponents[j];
        fdesc[j][1] = cnum[j];
        ASSERT_ALWAYS(exponents[j] < aa->t0);
    }
    aa->fdesc = fdesc;
    free(pivots);
    free(exponents);
    free(cnum);
    matpoly_clear(ab, M);
    abclear(ab, &tmp);
}				/*}}} */

int a_reading_task_compute_E_oneblock(a_reading_task_ptr aa, matpoly E, unsigned int pos, unsigned int block_size)/*{{{*/
{
    // F0 is exactly the n x n identity matrix, plus the
    // X^(s-exponent)e_{cnum} vectors. fdesc has the (exponent, cnum)
    // pairs
    bmstatus_ptr bm = aa->bm;
    dims * d = bm->d;
    unsigned int m = d->m;
    unsigned int n = d->n;
    unsigned int t0 = bm->t;
    abdst_field ab = d->ab;

    unsigned int s;

    for(s = 0 ; s < block_size ; s++) {
        /* keep this value of k, since this is the one of interest to us */
        unsigned int k = aa->k;

        if (!a_reading_task_read1(aa, t0 + 1)) {
            fprintf(stderr, "EOF met after reading %u coefficients\n", k);
            E->size = k - t0;
            break;
        }

        unsigned int kr = k % (t0 + 1);

        if (k > aa->guessed_length) {
            fprintf(stderr, "Going past guessed length ???\n");
        }

        for(unsigned int j = 0 ; j < n ; j++) {
            /* Column j of E simply comes from multiplication by the
             * identity matrix. Only coefficients of degree t0 and
             * above in A contribute, here. In effect, coefficient of
             * degree k in e is exactly coefficient of degree k + t0
             * in E.
             */
            matpoly_extract_column(ab, E, j, pos + s, aa->A, j, kr);
        }
        for(unsigned int j = n ; j < m + n ; j++) {
            /* Because the initial matrix F has powers of x in its
             * right part, we use columns of A from earlier reads
             * here.
             *
             * Let c, e be as below. Column j of E is x^(t0-e) times
             * column c of A. Coefficient of degree k in E is thus
             * coefficient of degree k + t0 in the product A *
             * x^(t0-e), which is thus coefficient of degree k + e in
             * A. Since e is in [0..t0-1]. Therefore, our cyclic
             * buffer of depth t0 + 1 still has the required
             * coefficient.
             */
            unsigned int c = aa->fdesc[j-n][1];
            unsigned int e = aa->fdesc[j-n][0];
            matpoly_extract_column(ab, E, j, pos + s, aa->A, c, (k - t0 + e) % (t0 + 1));
        }
    }
    return s;
}/*}}}*/

void a_reading_task_compute_E_local(a_reading_task_ptr aa, matpoly E)/*{{{*/
{
    // F0 is exactly the n x n identity matrix, plus the
    // X^(s-exponent)e_{cnum} vectors. fdesc has the (exponent, cnum)
    // pairs
    bmstatus_ptr bm = aa->bm;
    dims * d = bm->d;
    unsigned int m = d->m;
    unsigned int n = d->n;
    unsigned int b = m + n;
    unsigned int t0 = bm->t;
    abdst_field ab = d->ab;

    unsigned int guess = aa->guessed_length;
    ASSERT(!E->m && !E->n && !E->alloc);
    matpoly_init(ab, E, m, b, guess);
    E->size = guess;

    double tt0 = wct_seconds();
    double next_t = tt0 + 10;
    unsigned int k0 = aa->k;
    ASSERT_ALWAYS(k0 == t0);

    for(unsigned next_k = guess / 100 ; aa->k < guess ; ) {
        /* read E by blocks of B coefficients */
        unsigned int B = 1000;
        unsigned int nr = a_reading_task_compute_E_oneblock(aa, E, aa->k - t0, B);
        if (aa->k > next_k) {
            double tt = wct_seconds();
            if (tt > next_t) {
                fprintf(stderr,
                        "Read %u coefficients (%.1f%%) in %.1f s (%.1f MB/s)\n",
                        aa->k, 100.0 * aa->k / guess,
                        tt-tt0, (aa->k - k0)*aa->avg_matsize/(tt-tt0)/1.0e6);
                next_t = tt + 10;
                next_k = aa->k + guess/ 100;
            }
        }
        if (nr < B)
            break;
    }
}/*}}}*/

void set_random_input(bmstatus_ptr bm, matpoly A, unsigned int length) /* {{{ */
{
    dims * d = bm->d;
    unsigned int m = d->m;
    unsigned int n = d->n;
    abdst_field ab = d->ab;

    matpoly_init(ab, A, m, n, length);
    A->size = length;

    gmp_randstate_t rstate;
    gmp_randinit_default(rstate);
    for(unsigned int k = 0 ; k < length ; k++) {
	for (unsigned int i = 0; i < m ; i++) {
	    for (unsigned int j = 0; j < n ; j++) {
                abdst_elt x = matpoly_coeff(ab, A, i, j, k);
                abrandom(ab, x, rstate);
            }
        }
    }

    /* we force an existing generator. This will be a fairly trivial one,
     * but that's not really an issue as far as I can tell */
    if (length > 20) {
        for(unsigned int k = length-10 ; k < length ; k++) {
            for (unsigned int i = 0; i < m ; i++) {
                for (unsigned int j = 0; j < n ; j++) {
                    absrc_elt x = matpoly_coeff(ab, A, i, (j*1009)%n, k - (length - 10));
                    abdst_elt y = matpoly_coeff(ab, A, i, j, k);
                    abset(ab, y, x);
                }
            }
        }
    }

    gmp_randclear(rstate);
} /* }}} */


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



int main(int argc, char *argv[])
{
    bmstatus bm;
    dims * d = bm->d;
    int tune = 0;
    int ascii = 0;
    unsigned int random_input_length = 0;


    setvbuf(stdout, NULL, _IONBF, 0);
    setvbuf(stderr, NULL, _IONBF, 0);

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

#ifndef HAVE_MPIR
    if (caching) {
        fprintf(stderr, "--caching=1 only supported with MPIR\n");
        exit(EXIT_FAILURE);
    }
#endif

    bmstatus_init(bm, bw->m, bw->n);

    unsigned int m = d->m;
    unsigned int n = d->n;
    unsigned int b = m + n;
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
        if (rank == 0)
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
            if (rank == 0)
                fprintf(stderr, "The current plingen code is limited to square splits ; here, we received a %d x %x split, which will not work\n",
                    bm->mpi_dims[0], bm->mpi_dims[1]);
            abort();
        } else if ((m % bm->mpi_dims[0] != 0) || (n % bm->mpi_dims[0] != 0)) {
            if (rank == 0)
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

    unsigned int (*fdesc)[2] = NULL;    /* gcc is stupid */

    /* For the moment we read the complete thing in local memory before
     * dispatching, which is admittedly somewhat wasteful... Of course we
     * should rather do an mpi-level read. Either with only one
     * jobreading, and handing over data immediately to its peers, or
     * with several accesses to the file(s). */
    matpoly E;
    matpoly_init(ab, E, 0, 0, 0);


    MPI_Bcast(&(bm->t), 1, MPI_UNSIGNED, 0, bm->world);

    if (rank == 0) { /* {{{ compute E by reading A on the fly */
        // compute_initial_E(bm, E, A, fdesc);
        if (!random_input_length) {
            a_reading_task aa;
            a_reading_task_init(aa, bm, afile, ascii);
            a_reading_task_guess_length(aa);
            a_reading_task_compute_F(aa);
            fdesc = aa->fdesc;

            /* This is one possible way to go */
            a_reading_task_compute_E_local(aa, E);

            // compute_initial_E_from_A_ondisk(bm, E, afile, fdesc, ascii);
            aa->fdesc = NULL;
            a_reading_task_clear(aa);
        } else {
            matpoly_init(ab, E, m, b, random_input_length);
            gmp_randstate_t rstate;
            gmp_randinit_default(rstate);
            abvec_random(ab, E->x, m * b * random_input_length, rstate);
            gmp_randclear(rstate);
        }

    } /* }}} */


    /* This will quite probably be changed. We are playing nasty games
     * here, with E->size being used to draw decisions even though there
     * is no corresponding allocation.
     */
    MPI_Bcast(&(E->size), 1, MPI_UNSIGNED, 0, bm->world);

    // bw_bbpoly piL, piR, piP;

    unsigned int t0 = bm->t;

    unsigned int * delta = malloc((m + n) * sizeof(unsigned int));
    for(unsigned int j = 0 ; j < m + n ; j++) {
        delta[j] = t0;
    }

    matpoly pi;
    matpoly_init(ab, pi, 0, 0, 0);
    bw_biglingen_single(bm, pi, E, delta);
    matpoly_clear(ab, E);

    unsigned int nlucky = 0;
    int luck_mini = expected_pi_length(d, 0);
    for(unsigned int j = 0 ; j < b ; nlucky += bm->lucky[j++] >= luck_mini) ;

    if (rank == 0) {
        /* TODO: consider luck only below probability 2^-64 */
        if (nlucky == n) {
            matpoly f_red;
            matpoly_init(ab, f_red, 0, 0, 0);
            compute_final_F_red(bm, f_red, fdesc, t0, pi, delta);
            if (random_input_length) {
                fprintf(stderr, "Not writing result for random data\n");
            } else {
                char * f_filename;
                int rc = asprintf(&f_filename, "%s.gen", afile);
                ASSERT_ALWAYS(rc >= 0);
                write_f(bm, f_filename, f_red, delta, ascii);
                free(f_filename);
            }
            matpoly_clear(ab, f_red);
        } else {
            fprintf(stderr, "Could not find the required set of solutions (nlucky=%u)\n", nlucky);
        }

        printf("t_basecase = %.2f\n", bm->t_basecase);
        printf("t_mp = %.2f\n", bm->t_mp);
        printf("t_mul = %.2f\n", bm->t_mul);
        free(delta);
        matpoly_clear(ab, pi);
        free(fdesc);
    }

    abmpi_ops_clear(ab);
    abfield_clear(ab);
    bmstatus_clear(bm);
    bw_common_clear(bw);
    param_list_clear(pl);

    MPI_Finalize();
    return 0;
}

/* vim:set sw=4 sta et: */
