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

#include <assert.h>

#include "portability.h"
#include "macros.h"
#include "utils.h"
#include "abase.h"
#include "polymat.h"
#include "bigpolymat.h"
#include "bw-common.h"		/* Handy. Allows Using global functions
                                 * for recovering parameters */
#include "filenames.h"
#include "plingen.h"

// #include "types.h"
// #include "macros.h"
// #include "auxfuncs.h"
// #include "bw_scalar.h"
// #include "timer.h"
// #include "variables.h"
// #include "structure.h"
// #include "gmp-hacks.h"
// #include "modulus_hacks.h"
// #include "e_polynomial.h"
// #include "twisting_polynomials.h"
// #include "fft_on_matrices.hpp"
// #include "field_def.h"
// #include "field_prime.h"
// #include "field_quad.h"
// #include "field_usage.h"
// /* #include "version.h" */
// #include "master_common.h"
// #include "ops_poly.hpp"

static unsigned int display_threshold = 100;

struct bmstatus_s {
    dims d[1];
    unsigned int t;
    int * lucky;

    double t_basecase;
    double t_mp;

    unsigned int lingen_threshold;
    unsigned int lingen_mpi_threshold;
    int mpi_dims[2];
    MPI_Comm world;     /* reordered, in fact */
};
typedef struct bmstatus_s bmstatus[1];
typedef struct bmstatus_s *bmstatus_ptr;

#if 0/*{{{*/
bw_nbpoly f_poly;
int t_counter;
int *global_delta;
int global_sum_delta;
int dontread_pi = 0;
unsigned int *chance_list;
static int recursion_level;	/* static, hence 0 on init */
static int rec_threshold = 1;
static int print_min = 10;
static int preferred_quadratic_algorithm;	/* Defaults to 0 (old) */
static int t_init;
static int check_input = 1;

const char f_base_filename[] = "F_INIT";

#define SAVE_LEVEL_THRESHOLD	4

void reclevel_prolog(void)
{
    int i;
    printf("%2d [ t=%6d ] ", recursion_level, t_counter);
    for (i = 0; i < recursion_level; i++)
	printf("  ");
}

int sum_delta(dims * d, int * delta)
{
    int i, res = 0;
    for (i = 0; i < d->b; i++)
	res += delta[i];
    return res;
}

int max_delta(dims * d, int *delta)
{
    int i, res = -1;
    for (i = 0; i < d->b; i++)
	if (delta[i] > res)
	    res = delta[i];
    return res;
}


/* {{{ low priority */
#if 0
static int save_pi(struct t_poly *pi, int t_start, int t_middle, int t_end)
{
    char filename[FILENAME_LENGTH];
    FILE *f;
    int res;

    sprintf(filename, pi_meta_filename, t_start, t_end);
    f = fopen(filename, "wb");
    if (f == NULL)
	return -1;
    res = tp_write(f, pi);
    if (fclose(f) < 0)
	res = -1;

    if (res == 0) {
	printf("Saved file %s\n", filename);
    } else {
	fprintf(stderr, "Failure to save %s: %s\n",
		filename, strerror(errno));
	return -1;
    }

    if (t_middle == -1)
	return res;

#if 0
    /* Unlinking not done, for safety */
    sprintf(filename, pi_meta_filename, t_start, t_middle);

    if (unlink(filename) < 0) {
	fprintf(stderr, "Cannot unlink %s: %s\n", filename, strerror(errno));
    } else {
	printf("Unlinked %s\n", filename);
    }

    sprintf(filename, pi_meta_filename, t_middle, t_end);

    if (unlink(filename) < 0) {
	fprintf(stderr, "Cannot unlink %s: %s\n", filename, strerror(errno));
    } else {
	printf("Unlinked %s\n", filename);
    }
#endif
    return res;
}


static int retrieve_pi_files(struct t_poly **p_pi, int t_start)
{
    DIR *pi_dir;
    struct dirent *curr;
    int n_pi_files;
    struct couple {
	int s;
	int e;
    } *pi_files;
    const char *pattern;
    int i;
    struct t_poly *left = NULL, *right = NULL;
    ft_order_t order;
    int o_i;
    struct dft_bb *dft_left, *dft_right, *dft_prod;
    double tt;

    *p_pi = NULL;

    if (dontread_pi)
	return t_start;

    if ((pi_dir = opendir(".")) == NULL) {
	perror(".");
	return t_start;
    }

    printf("Scanning directory %s for pi files\n", ".");

    pattern = strrchr(pi_meta_filename, '/');
    if (pattern == NULL) {
	pattern = pi_meta_filename;
    } else {
	pattern++;
    }

    for (n_pi_files = 0; (curr = readdir(pi_dir)) != NULL;) {
	int s, e;
	if (sscanf(curr->d_name, pattern, &s, &e) == 2) {
	    printf("Found %s\n", curr->d_name);
	    if (s > e) {
		printf("but that's a stupid one\n");
		continue;
	    }
	    n_pi_files++;
	}
    }

    if (n_pi_files == 0) {
	printf("Found no pi files\n");
	return t_start;
    }

    pi_files = (struct couple *) malloc(n_pi_files * sizeof(struct couple));

    rewinddir(pi_dir);

    for (i = 0; i < n_pi_files && (curr = readdir(pi_dir)) != NULL;) {
	if (sscanf(curr->d_name,
		   pattern, &(pi_files[i].s), &(pi_files[i].e)) == 2) {
	    if (pi_files[i].s > pi_files[i].e)
		continue;
	    i++;
	}
    }
    n_pi_files = i;
    closedir(pi_dir);

    /* The rule is: only look at the best candidate. It's not worth
     * bothering about more subtle cases */
    for (;;) {
	int t_max;
	int best = -1;
	FILE *f;
	char filename[FILENAME_LENGTH];

	printf("Scanning for data starting at t=%d\n", t_start);
	t_max = -1;
	for (i = 0; i < n_pi_files; i++) {
	    if (pi_files[i].s == t_start && pi_files[i].e != -1) {
		printf("candidate : ");
		printf(pattern, pi_files[i].s, pi_files[i].e);
		printf("\n");
		if (pi_files[i].e > t_max) {
		    t_max = pi_files[i].e;
		    best = i;
		}
	    }
	}
	if (t_max == -1) {
	    printf("Could not find such data\n");
	    break;
	}

	sprintf(filename, pi_meta_filename, t_start, t_max);
	printf("trying %s\n", filename);
	f = fopen(filename, "rb");
	if (f == NULL) {
	    perror(filename);
	    pi_files[best].e = -1;
	    continue;
	}
	/* Which degree can we expect for t_start..t_max ?
	 */

	unsigned int pideg;
	pideg = iceildiv(m_param * (t_max - t_start), bigdim);
	pideg += 10;
	if (t_max > bm->len) {
	    pideg += t_max - bm->len;
	}

	right = tp_read(f, pideg);
	fclose(f);

	if (right == NULL) {
	    printf("%s : bad or nonexistent data\n", filename);
	    pi_files[best].e = -1;
	    continue;
	}

	if (left == NULL) {
	    left = right;
	    right = NULL;
	    t_start = t_max;
	    continue;
	}

	printf("Beginning multiplication\n");
	*p_pi = tp_comp_alloc(left, right);
	core_if_null(*p_pi, "*p_pi");

	order.set((*p_pi)->degree + 1);
	o_i = order;

	dft_left = fft_tp_dft(left, order, &tt);
	printf("DFT(pi_left,%d) : %.2fs\n", o_i, tt);
	core_if_null(dft_left, "dft_left");

	dft_right = fft_tp_dft(right, order, &tt);
	printf("DFT(pi_right,%d) : %.2fs\n", o_i, tt);
	core_if_null(dft_right, "dft_right");

	dft_prod = fft_bbb_conv(dft_left, dft_right, &tt);
	printf("CONV(pi_left,pi_right,%d) : %.2fs\n", o_i, tt);
	core_if_null(dft_prod, "dft_prod");

	fft_tp_invdft(*p_pi, dft_prod, &tt);
	printf("IDFT(pi,%d) : %.2fs\n", o_i, tt);

	tp_free(left);
	tp_free(right);
	dft_bb_free(dft_left);
	dft_bb_free(dft_right);
	dft_bb_free(dft_prod);
	left = *p_pi;
	right = NULL;
	t_start = t_max;
    }
    free(pi_files);

    *p_pi = left;
    return t_start;
}
#endif /* }}} */

#endif/*}}}*/

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
static int bw_lingen(bmstatus_ptr bm, polymat pi, polymat E, unsigned int *delta);


/* This destructively cancels the first len coefficients of E, and
 * computes the appropriate matrix pi which achieves this. The
 * elimination is done in accordance with the nominal degrees found in
 * delta.
 *
 * The result is expected to have degree ceil(len*m/b) coefficients, so
 * that E*pi is divisible by X^len.
 */

static int bw_lingen_basecase(bmstatus_ptr bm, polymat pi, polymat E, unsigned int *delta) /*{{{*/
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

    polymat_init(pi, b, b, pi_room_base);

    /* Also keep track of the
     * number of coefficients for the columns of pi. Set pi to Id */

    unsigned int *pi_lengths = malloc(b * sizeof(unsigned int));
    for(unsigned int i = 0 ; i < b ; i++) {
        abset_ui(ab, polymat_coeff(pi, i, i, 0), 1);
        pi_lengths[i] = 1;
    }

    /* This is used below */
    int generator_found = 0;

    /* Keep a list of columns which have been used as pivots at the
     * previous iteration */
    unsigned int * pivots = malloc(m * sizeof(unsigned int));
    int * is_pivot = malloc(b * sizeof(int));
    memset(is_pivot, 0, b * sizeof(int));

    polymat e;
    polymat_init(e, m, b, 1);

    for (unsigned int t = 0; t < E->size ; t++, bm->t++) {

        /* {{{ Update the columns of e for degree t. Save computation
         * time by not recomputing those which can easily be derived from
         * previous iteration. Notice that the columns of e are exactly
         * at the physical positions of the corresponding columns of pi.
         */

        abelt_ur tmp_ur;
        abelt_ur_init(ab, &tmp_ur);
        polymat_ur e_ur;
        polymat_ur_init(e_ur, m, b, 1);
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
                                polymat_coeff(E, i, k, t - s),
                                polymat_coeff(pi, k, j, s));
                        abelt_ur_add(ab,
                                polymat_ur_coeff(e_ur, i, j, 0),
                                polymat_ur_coeff(e_ur, i, j, 0),
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
                abreduce(ab,
                        polymat_coeff(e, i, j, 0),
                        polymat_ur_coeff(e_ur, i, j, 0)
                        );
                nz += abcmp_ui(ab, polymat_coeff(e, i, j, 0), 0) == 0;
            }
            if (nz == m) {
                newluck++, bm->lucky[j]++;
            } else if (bm->lucky[j] > 0) {
                bm->lucky[j] = 0;
            }
        }
        abelt_ur_clear(ab, &tmp_ur);
        polymat_ur_clear(e_ur);
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
                if (abcmp_ui(ab, polymat_coeff(e, u, j, 0), 0) != 0)
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
            int rc = abinv(ab, inv, polymat_coeff(e, u, j, 0));
            if (!rc) {
                fprintf(stderr, "Error, found a factor of the modulus: ");
                abfprint(ab, stderr, inv);
                fprintf(stderr, "\n");
                exit(1);
            }
            abneg(ab, inv, inv);
            for (unsigned int kl = jl + 1; kl < b ; kl++) {
                unsigned int k = ctable[kl][1];
                if (abcmp_ui(ab, polymat_coeff(e, u, k, 0), 0) == 0)
                    continue;
                // add lambda = e[u,k]*-e[u,j]^-1 times col j to col k.
                abelt lambda;
                abinit(ab, &lambda);
                abmul(ab, lambda, inv, polymat_coeff(e, u, k, 0));

                assert(delta[j] <= delta[k]);
                /* {{{ Apply on both e and pi */
                abelt tmp;
                abinit(ab, &tmp);
                for(unsigned int i = 0 ; i < m ; i++) {
                    /* TODO: Would be better if mpfq had an addmul */
                    abmul(ab, tmp, lambda, polymat_coeff(e, i, j, 0));
                    abadd(ab,
                            polymat_coeff(e, i, k, 0),
                            polymat_coeff(e, i, k, 0),
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
                                polymat_coeff(pi, i, j, s));
                        abadd(ab,
                                polymat_coeff(pi, i, k, s),
                                polymat_coeff(pi, i, k, s),
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
                    polymat_realloc(pi, pi->alloc + pi->alloc / (m+n));
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
            bwmat_move_coeffs(ab,
                    polymat_part(pi, 0, j, 1), b,
                    polymat_part(pi, 0, j, 0), b,
                    b * pi_lengths[j]);
            for(unsigned int i = 0 ; i < b ; i++) {
                abset_ui(ab, polymat_coeff(pi, i, j, 0), 0);
            }
            pi_lengths[j]++;
            delta[j]++;
        }
        /* }}} */

        /*
        printf("t=%u:", bm->t);
        for(unsigned int j = 0; j < b; j++) {
            printf(" %u", delta[j]);
        }
        printf("\n");
        */
    }
    /*
        printf("t=%u: delta =", bm->t);
        for(unsigned int j = 0; j < b; j++) {
            printf(" %u", delta[j]);
        }
        printf("\n");
        printf("t=%u: pi_length =", bm->t + len);
        for(unsigned int j = 0; j < b; j++) {
            printf(" %u", pi_lengths[j]);
        }
        printf("\n");
        */
    for(unsigned int j = 0; j < b; j++) {
        if (pi_lengths[j] > pi->size)
            pi->size = pi_lengths[j];
    }
    pi->size = MIN(pi->size, pi->alloc);
    polymat_clear(e);
    free(is_pivot);
    free(pivots);
    free(pi_lengths);   /* What shall we do with this one ??? */

    return generator_found;
}/*}}}*/

static int bw_lingen_recursive(bmstatus_ptr bm, polymat pi, polymat E, unsigned int *delta) /*{{{*/
{
    dims * d = bm->d;
    unsigned int m = d->m;
    unsigned int n = d->n;
    unsigned int b = m + n;
    abdst_field ab = d->ab;
    int done;

    /* XXX I think we have to start with something large enough to get
     * all coefficients of E_right correct */
    size_t half = E->size - (E->size / 2);
    polymat E_left;
    polymat_init(E_left, m, b, half);

    bwmat_copy_coeffs(ab,
            polymat_part(E_left,0,0,0),1,
            polymat_part(E,0,0,0),1,
            m*b*half);
    E_left->size = half;

    polymat pi_left;
    polymat_init(pi_left, 0, 0, 0);

    done = bw_lingen(bm, pi_left, E_left, delta);

    polymat_clear(E_left);

    if (done) {
        polymat_swap(pi_left, pi);
        polymat_clear(pi_left);
        return 1;
    }

    /* Do a naive middle product for the moment, just to make sure I'm
     * not speaking nonsense */
    /* First get coefficient bounds in E */
    /* E_i0 + (max degree in pi) = half */
    unsigned int E_i0 = half - (pi_left->size - 1);
    unsigned int E_i1 = E->size;

    /* length of the middle product is the difference of lengths + 1 */
    unsigned mp_len = E_i1 - E_i0 - (pi_left->size - 1);
    if (E->size > display_threshold)
        printf("t=%u, MP(%zu, %u) --> %u\n", bm->t, pi_left->size, E_i1 - E_i0, mp_len);

    polymat E_right;
    polymat_init(E_right, m, b, mp_len);
    E_right->size = mp_len;

    bm->t_mp -= seconds();
    polymat_mp_raw(ab,
            E_right, 0,
            E, E_i0, E_i1 - E_i0,
            pi_left, 0, pi_left->size, 0, 0);
    bm->t_mp += seconds();

    polymat pi_right;
    polymat_init(pi_right, 0, 0, 0);

    done = bw_lingen(bm, pi_right, E_right, delta);

    polymat_clear(E_right);

    if (E->size > display_threshold)
        printf("t=%u, MUL(%zu, %zu) --> %zu\n", bm->t, pi_left->size, pi_right->size, pi_left->size + pi_right ->size - 1);
    bm->t_mp -= seconds();
    polymat_mul(ab, pi, pi_left, pi_right);
    bm->t_mp += seconds();
    
    polymat_clear(pi_left);
    polymat_clear(pi_right);

    return done;
}/*}}}*/

void debug_matrix_print(polymat_srcptr M, const char * fmt, ...)
{
    return;
    va_list ap;
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    va_start(ap, fmt);
    /* We have three buffers: initial message, with matrix, and with rank
     * info prepended.
     */
    char * buf[3];
    int len[3];
    int alloc[3];
    int pos[3];

    len[0] = vasprintf(&(buf[0]), fmt, ap);
    ASSERT_ALWAYS(len[0] >= 0);

    alloc[1] = len[0] + 16 * M->m * M->n + 4*M->m;
    buf[1] = malloc(alloc[1]);
    pos[1] = 0;

    memcpy(buf[1], buf[0], len[0] + 1);
    pos[1] = len[0];

    /* For the moment we're assuming only one unsigned long */
    for(unsigned int i = 0 ; i < M->m ; i++) {
        int v;
        for(unsigned int j = 0 ; j < M->n ; j++) {
            ASSERT_ALWAYS(pos[1] < alloc[1] - 16);
            v = snprintf(buf[1] + pos[1], alloc[1]-pos[1],  " %4lu", ((unsigned long*)(M->x))[i*M->n + j]);
            ASSERT_ALWAYS(v >= 0);
            pos[1] += v;
        }
        ASSERT_ALWAYS(pos[1] < alloc[1] - 10);
        v = snprintf(buf[1] + pos[1], alloc[1]-pos[1], "\n");
        ASSERT_ALWAYS(v >= 0);
        pos[1] += v;
    }
    len[1] = pos[1];

    /* Prepend rank info, now */
    alloc[2] = len[1] + 4 * M->m + len[0];
    buf[2] = malloc(alloc[2]);
    pos[2] = 0;

    for(pos[1] = 0 ; pos[1] < len[1] ; ) {
        int v;
        ASSERT_ALWAYS(pos[2] < alloc[2] - 4);
        v = snprintf(buf[2] + pos[2], alloc[2] - pos[2], "%2d: ", rank);
        ASSERT_ALWAYS(v >= 0);
        pos[2] += v;
        for( ; pos[1] < len[1] ; ) {
            ASSERT_ALWAYS(pos[2] < alloc[2]);
            if ((buf[2][pos[2]++] = buf[1][pos[1]++]) == '\n')
                break;
        }
    }

    fputs(buf[2], stdout);
    free(buf[2]);
    free(buf[1]);
    free(buf[0]);
}

/* This version works over MPI */
static int bw_lingen_bigrecursive(bmstatus_ptr bm, bigpolymat pi, bigpolymat E, unsigned int *delta) /*{{{*/
{
    dims * d = bm->d;
    unsigned int m = d->m;
    unsigned int n = d->n;
    unsigned int b = m + n;
    unsigned int m0 = E->m0;
    unsigned int b0 = E->n0;    /* E is m * b */
    // unsigned int n0 = b0 - m0;

    abdst_field ab = d->ab;
    int done;

        int irank;
        int jrank;
        MPI_Comm_rank(E->col, &irank);
        MPI_Comm_rank(E->row, &jrank);

    debug_matrix_print(bigpolymat_my_cell(E), 
            "Entering bw_lingen_bigrecursive, size %zu/%zu\n"
            "local matrix E at degree 0:\n", E->size, bigpolymat_my_cell(E)->size);

    if (E->size < bm->lingen_mpi_threshold) {
        /* Fall back to local code */
        /* This entails gathering E locally, computing pi locally, and
         * dispathing it back. */

        int done;
        polymat sE, spi;
        polymat_init(sE, m, b, E->size);
        polymat_init(spi, 0, 0, 0);
        bigpolymat_gather_mat(ab, sE, E);
        /* Only the master node does the local computation */
        if (!irank && !jrank) {
            debug_matrix_print(sE, "Global matrix [X^0]E at root (size %zu):\n", sE->size);
            done = bw_lingen_recursive(bm, spi, sE, delta);
            debug_matrix_print(spi, "Output: global matrix [X^0]pi at root (size %zu):\n", spi->size);
        }
        bigpolymat_scatter_mat(ab, pi, spi);
        debug_matrix_print(bigpolymat_my_cell(pi), "Output: local matrix [X^0]pi (complete size %zu/%zu):\n", pi->size, bigpolymat_my_cell(pi)->size);
        /* Don't forget to broadcast delta from root node to others ! */
        MPI_Bcast(&done, 1, MPI_INT, 0, bm->world);
        MPI_Bcast(delta, b, MPI_UNSIGNED, 0, bm->world);
        MPI_Bcast(bm->lucky, b, MPI_UNSIGNED, 0, bm->world);
        MPI_Bcast(&(bm->t), 1, MPI_UNSIGNED, 0, bm->world);
        polymat_clear(spi);
        polymat_clear(sE);
        return done;
    }


    /* XXX I think we have to start with something large enough to get
     * all coefficients of E_right correct */
    size_t half = E->size - (E->size / 2);
    bigpolymat E_left;
    bigpolymat_init(E_left, E, m, b, half);

    /* should be part of the interface, I guess. OTOH it remains simple
     * enough */
    bwmat_copy_coeffs(ab,
            polymat_part(bigpolymat_my_cell(E_left),0,0,0),1,
            polymat_part(bigpolymat_my_cell(E),0,0,0),1,
            m0*b0*half);

    bigpolymat_set_size(E_left, half);

    bigpolymat pi_left;
    bigpolymat_init(pi_left, pi, 0, 0, 0);      /* pre-init */

    done = bw_lingen_bigrecursive(bm, pi_left, E_left, delta);

    bigpolymat_clear(E_left);

    if (done) {
        bigpolymat_swap(pi_left, pi);
        bigpolymat_clear(pi_left);
        return 1;
    }

    /* Do a naive middle product for the moment, just to make sure I'm
     * not speaking nonsense */
    /* First get coefficient bounds in E */
    /* E_i0 + (max degree in pi) = half */
    unsigned int E_i0 = half - (pi_left->size - 1);
    unsigned int E_i1 = E->size;

    /* length of the middle product is the difference of lengths + 1 */
    unsigned mp_len = E_i1 - E_i0 - (pi_left->size - 1);
    if (!irank && !jrank && E->size > display_threshold)
        printf("t=%u, MPI-MP(%zu, %u) --> %u\n", bm->t, pi_left->size, E_i1 - E_i0, mp_len);

    bigpolymat E_right;
    bigpolymat_init(E_right, E, m, b, mp_len);
    E_right->size = mp_len;

    bm->t_mp -= seconds();
    bigpolymat_mp_raw(ab,
            E_right, 0,
            E, E_i0, E_i1 - E_i0,
            pi_left, 0, pi_left->size, 0, 0);
    bm->t_mp += seconds();

    bigpolymat pi_right;
    bigpolymat_init(pi_right, pi, 0, 0, 0);     /* pre-init */

    done = bw_lingen_bigrecursive(bm, pi_right, E_right, delta);

    bigpolymat_clear(E_right);

    if (!irank && !jrank && E->size > display_threshold)
        printf("t=%u, MPI-MUL(%zu, %zu) --> %zu\n", bm->t, pi_left->size, pi_right->size, pi_left->size + pi_right ->size - 1);
    bm->t_mp -= seconds();
    bigpolymat_mul(ab, pi, pi_left, pi_right);
    bm->t_mp += seconds();
    
    bigpolymat_clear(pi_left);
    bigpolymat_clear(pi_right);

    return done;
}/*}}}*/


static int/*{{{*/
bw_lingen(bmstatus_ptr bm, polymat pi, polymat E, unsigned int *delta)
{
    dims *d = bm->d;
    abdst_field ab = d->ab;
    unsigned int m = d->m;
    unsigned int n = d->n;
    unsigned int b = m + n;

    if (E->size >= bm->lingen_mpi_threshold) {
        /* We are going to delegate to the MPI code */
        bigpolymat model;
        bigpolymat xpi, xE;

        bigpolymat_init_model(model, bm->world, bm->mpi_dims[0], bm->mpi_dims[1]);
        /* We prefer to allocate soon. The interface doesn't really like
         * lazy allocation at the moment */
        bigpolymat_init(xE, model, m, b, E->size);
        bigpolymat_init(xpi, model, 0, 0, 0);   /* pre-init for now */
        bigpolymat_scatter_mat(ab, xE, E);
        ASSERT_ALWAYS(xE->size);
        int res = bw_lingen_bigrecursive(bm, xpi, xE, delta);
        bigpolymat_gather_mat(ab, pi, xpi);
        bigpolymat_clear(xE);
        bigpolymat_clear(xpi);
        bigpolymat_clear_model(model);
        return res;
    } else if (E->size < bm->lingen_threshold) {
        bm->t_basecase -= seconds();
        int res = bw_lingen_basecase(bm, pi, E, delta);
        bm->t_basecase += seconds();
        return res;
    } else {
        return bw_lingen_recursive(bm, pi, E, delta);
    }
}/*}}}*/

#if 0/*{{{*/
int check_zero_and_advance(struct e_coeff *ec, unsigned int kill)
{
    unsigned int i;
    for (i = 0; i < kill; i++) {
	int res = 1;
	int j;
	for (j = 0; res && j < bigdim; j++) {
	    if (!mcol_is_zero(mbmat_col(mbpoly_coeff(ec->p, 0), j))) {
		die("argh, not zero !\n", 1);
		return 1;
	    }
	}
	if (!res) {
	    return 0;
	}
	ec_advance(ec, 1);
    }
    return 1;
}

static double
bw_recursive_algorithm(struct e_coeff *ec, int *delta, struct t_poly **p_pi)
{
    struct t_poly *pi_left, *pi_right;
    int deg, ldeg, rdeg;
    struct dft_mb *dft_e_left, *dft_e_middle;
    struct dft_bb *dft_pi_left, *dft_pi_right;
    struct dft_bb *dft_pi;
    ft_order_t sub_order;
    int so_i;
    int expected_pi_deg;
    double t_dft_e_l, t_dft_pi_l, t_conv_e, t_idft_e,
	t_dft_pi_r, t_conv_pi, t_idft_pi, t_ft, t_cv, t_sub;
    int kill;
    int t0 = t_counter;
    double tt = 0;

    tt -= seconds();

    deg = ec->degree;

    /* Repartition of the job:
     *
     *          left            right
     * deg==0   1               0       (never recursive)
     * deg==1   1               1
     * deg==2   2               1
     * deg==n   n/2 + 1         (n+1)/2
     * 
     * The figures are for the number of steps, each one corres-
     * ponding to a m/(m+n) increase of the average degree of pi.
     */

    ldeg = (deg / 2) + 1;
    rdeg = (deg + 1) / 2;

    assert(ldeg && rdeg && ldeg + rdeg == deg + 1);

    /* We aim at computing ec * pi / X^ldeg. The degree of this
     * product will be
     *
     * ec->degree + pi->degree - ldeg
     *
     * (We are actually only interested in the low (ec->degree-ldeg)
     * degree part of the product, but the whole thing is required)
     *
     * The expected value of pi->degree is 
     *  ceil(ldeg*m/(m+n))
     *
     * The probability that pi exceeds this expected degree
     * depends on the base field, but is actually low.
     * However, by the end of the computations, this does
     * happen because the degrees increase unevenly.
     *
     * The DFTs of e and pi can be computed using only the
     * number of points given above, *even if their actual
     * degree is higher*. The FFT routines need to have
     * provision for this.
     *
     * The number of points will then be the smallest power
     * of 2 above deg+ceil(ldeg*m/(m+n))-ldeg+1
     */

    expected_pi_deg = 10 + iceildiv(ldeg * m_param, bigdim);
#ifdef	HAS_CONVOLUTION_SPECIAL
    kill = ldeg;
#else
    kill = 0;
#endif
    sub_order.set(deg + expected_pi_deg - kill + 1);

    so_i = sub_order;

    dft_e_left = fft_ec_dft(ec, sub_order, &t_dft_e_l);
    reclevel_prolog();
    printf("DFT(e,%d) : %.2fs\n", so_i, t_dft_e_l);
    core_if_null(dft_e_left, "dft_e_left");

    ec->degree = ldeg - 1;
    t_sub = bw_lingen(ec, delta, &pi_left);

    if (t_counter < t0 + ldeg) {
	printf("Exceptional situation, small generator ; escaping\n");
	*p_pi = pi_left;
	dft_mb_free(dft_e_left);
	return tt + seconds();
    }

    printf("deg(pi_l)=%d, bound is %d\n", pi_left->degree, expected_pi_deg);
    if (!sub_order.fits(deg + pi_left->degree - kill + 1)) {
	printf("Warning : pi grows above its expected degree...\n");
	printf("order %d , while :\n"
	       "deg=%d\n"
	       "deg(pi_left)=%d\n"
	       "ldeg-1=%d\n"
	       "hence, %d is too big\n",
	       so_i,
	       deg, pi_left->degree, ldeg - 1,
	       deg + pi_left->degree - kill + 1);

	int n_exceptional = 0;
	for (int i = 0; i < bigdim; i++) {
	    n_exceptional += chance_list[i] * m_param;
	}
	if (!n_exceptional) {
	    die("This should only happen at the end of the computation\n", 1);
	}
	*p_pi = pi_left;
	dft_mb_free(dft_e_left);
	return tt + seconds();
    }

    dft_pi_left = fft_tp_dft(pi_left, sub_order, &t_dft_pi_l);
    reclevel_prolog();
    printf("DFT(pi_l,%d) : %.2fs\n", so_i, t_dft_pi_l);
    core_if_null(dft_pi_left, "dft_pi_left");

#ifdef  HAS_CONVOLUTION_SPECIAL
    dft_e_middle = fft_mbb_conv_sp(dft_e_left, dft_pi_left, ldeg, &t_conv_e);
#else
    dft_e_middle = fft_mbb_conv(dft_e_left, dft_pi_left, &t_conv_e);
#endif
    reclevel_prolog();
    printf("CONV(e*pi_l,%d) : %.2fs\n", so_i, t_conv_e);
    core_if_null(dft_e_middle, "dft_e_middle");

    /* This is a special convolution in the sense that we
     * compute f(w)*g(w) / w^k for k=ldeg, since we are
     * interested in fg div X^k (we know fg mod X^k==0)
     */

    ec_park(ec);
    ec_untwist(ec);

    fft_mb_invdft(ec->p, dft_e_middle, deg - kill, &t_idft_e);
    ec->degree = deg - kill;

    check_zero_and_advance(ec, ldeg - kill);
    reclevel_prolog();
    printf("IDFT(e,%d) : %.2fs\n", (int) (dft_e_middle->order), t_idft_e);

    dft_mb_free(dft_e_middle);
    dft_mb_free(dft_e_left);

    assert(ec->degree == rdeg - 1);

    t_sub += bw_lingen(ec, delta, &pi_right);
    printf("deg(pi_r)=%d, bound is %d\n", pi_right->degree, expected_pi_deg);

    *p_pi = tp_comp_alloc(pi_left, pi_right);
    core_if_null(*p_pi, "*p_pi");

    printf("deg(pi_prod)=%d (order %d)\n", (*p_pi)->degree, so_i);
    assert(sub_order.fits((*p_pi)->degree + 1));

    dft_pi_right = fft_tp_dft(pi_right, sub_order, &t_dft_pi_r);
    reclevel_prolog();
    printf("DFT(pi_r,%d) : %.2fs\n", so_i, t_dft_pi_r);
    core_if_null(dft_pi_right, "dft_pi_right");

    dft_pi = fft_bbb_conv(dft_pi_left, dft_pi_right, &t_conv_pi);
    reclevel_prolog();
    printf("CONV(pi_l*pi_r,%d) : %.2fs\n", so_i, t_conv_pi);
    core_if_null(dft_pi, "dft_pi");


    fft_tp_invdft(*p_pi, dft_pi, &t_idft_pi);
    reclevel_prolog();
    printf("IDFT(pi,%d) : %.2fs\n", (int) (dft_pi->order), t_idft_pi);

    dft_bb_free(dft_pi);
    dft_bb_free(dft_pi_right);
    dft_bb_free(dft_pi_left);
    tp_free(pi_left);
    tp_free(pi_right);

    reclevel_prolog();

    t_ft = t_dft_e_l + t_dft_pi_l + t_idft_e + t_dft_pi_r + t_idft_pi;
    t_cv = t_conv_e + t_conv_pi;

    printf("proper : %.2fs (%.2fs FT + %.2fs CV), sub : %.2fs\n",
	   t_ft + t_cv, t_ft, t_cv, t_sub);
    printf("constants : c_ft=%.4e c_cv=%.4e		# %d,%d\n",
	   t_ft / (double) (so_i << so_i),
	   t_cv / (double) (1 << so_i), deg, so_i);
    printf("Different values for M1:");
    printf("   e_left: M1=%.3e\n",
	   t_dft_e_l / (so_i << so_i) / (m_param * bigdim));
    printf("  pi_left: M1=%.3e\n",
	   t_dft_pi_l / (so_i << so_i) / (bigdim * bigdim));
    printf("    e_inv: M1=%.3e\n",
	   t_idft_e / (so_i << so_i) / (m_param * bigdim));
    printf(" pi_right: M1=%.3e\n",
	   t_dft_pi_r / (so_i << so_i) / (bigdim * bigdim));
    printf("   pi_inv: M1=%.3e\n",
	   t_idft_pi / (so_i << so_i) / (bigdim * bigdim));
    printf("   e_conv: M1=%.3e\n",
	   t_conv_e / (1 << so_i) / (m_param * bigdim * bigdim));
    printf("  pi_conv: M1=%.3e\n",
	   t_conv_pi / (1 << so_i) / (bigdim * bigdim * bigdim));
    return tt + seconds();
}

static double bw_lingen(struct e_coeff *ec, int *delta, struct t_poly **p_pi)
{
    /* int check_chance; */
    double inner;
    int deg;
    int t_before;
    int did_rec = 0;

    t_before = t_counter;
    deg = ec->degree;
    reclevel_prolog();
    printf("Degree %d\n", deg);

    recursion_level++;
    /* check_chance=(t_counter > bm->len - 5); */
    if (ec->degree < rec_threshold) {
	/* check_chance=(t_counter > bm->len - 5 - ec->degree); */
	inner = bw_traditional_algorithm(ec, delta, p_pi, 1);
    } else {
	did_rec = 1;
	inner = bw_recursive_algorithm(ec, delta, p_pi);
    }
    recursion_level--;

    reclevel_prolog();
    printf("Degree %d : took %.2fs\n", deg, inner);

    if (recursion_level <= SAVE_LEVEL_THRESHOLD) {
	if (!did_rec || recursion_level == SAVE_LEVEL_THRESHOLD) {
	    save_pi(*p_pi, t_before, -1, t_counter);
	} else {
	    int t_middle;
	    t_middle = t_before + (deg / 2) + 1;
	    save_pi(*p_pi, t_before, t_middle, t_counter);
	}
    }

    return inner;
}

void showuse(void)
{
    die("Usage : bw-master <bank#>\n", 1);
}

#endif/*}}}*/

/**********************************************************************/

unsigned int (*compute_initial_F(bmstatus_ptr bm, polymat A))[2] /*{{{ */
{				
    dims *d = bm->d;
    abdst_field ab = d->ab;
    unsigned int m = d->m;
    unsigned int n = d->n;

    ASSERT(A->m == m);
    ASSERT(A->n == n);

    abelt tmp;
    abinit(ab, &tmp);

    /* First try to create the initial F matrix */
    printf("Computing t0\n");

    /* We want to create a full rank m*m matrix M, by extracting columns
     * from the first coefficients of A */

    polymat M;
    polymat_init(M, m, m, 1);

    /* For each integer i between 0 and m-1, we have a column, picked
     * from column cnum[i] of coeff exponent[i] of A which, once reduced modulo
     * the other ones, has coefficient at row pivots[i] unequal to zero.
     */
    unsigned int *pivots = malloc(m * sizeof(unsigned int));
    unsigned int *exponents = malloc(m * sizeof(unsigned int));
    unsigned int *cnum = malloc(m * sizeof(unsigned int));
    unsigned int r = 0;

    for (unsigned int k = 0; r < m && k < A->size; k++) {
	for (unsigned int j = 0; r < m && j < n; j++) {
	    /* Extract a full column into M */
	    bwmat_copy_coeffs(ab,
			      polymat_part(M, 0, r, 0),
			      m, polymat_part(A, 0, j, k), n, m);

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
                              polymat_coeff(M, i, v, 0),
                              polymat_coeff(M, u, r, 0));
                    abadd(ab, polymat_coeff(M, i, r, 0),
                              polymat_coeff(M, i, r, 0),
                              tmp);
                }
                abset_zero(ab,
                        polymat_coeff(M, u, r, 0));
	    }
            unsigned int u = 0;
            for( ; u < m ; u++) {
                if (abcmp_ui(ab, polymat_coeff(M, u, r, 0), 0) != 0)
                    break;
            }
            if (u == m) {
		printf("[X^%d] A, col %d does not increase rank (still %d)\n",
		       k, j, r);
		if (k * n > m + 40) {
		    printf("The choice of starting vectors was bad. "
			   "Cannot find %u independent cols within A\n", m);
		    exit(1);
		}
		continue;
	    }

	    /* Bingo, it's a new independent col. */
	    pivots[r] = u;
	    cnum[r] = j;
	    exponents[r] = k;

	    /* Multiply the column so that the pivot becomes -1 */
            int rc = abinv(ab, tmp, polymat_coeff(M, u, r, 0));
            if (!rc) {
                fprintf(stderr, "Error, found a factor of the modulus: ");
                abfprint(ab, stderr, tmp);
                fprintf(stderr, "\n");
                exit(1);
            }
            abneg(ab, tmp, tmp);
            for(unsigned int i = 0 ; i < m ; i++) {
                abmul(ab, polymat_coeff(M, i, r, 0),
                          polymat_coeff(M, i, r, 0),
                          tmp);
            }

	    r++;

	    // if (r == m)
		printf
		    ("[X^%d] A, col %d increases rank to %d (head row %d)\n",
		     k, j, r, u);
	}
    }

    if (r != m) {
	printf("This amount of data is insufficient. "
	       "Cannot find %u independent cols within A\n", m);
	exit(1);
    }

    unsigned int t0 = exponents[r - 1] + 1;
    printf("Found satisfying init data for t0=%d\n", t0);
 
    bm->t = t0;

    unsigned int (*fdesc)[2] = malloc(2 * m * sizeof(unsigned int));
    for(unsigned int j = 0 ; j < m ; j++) {
        fdesc[j][0] = exponents[j];
        fdesc[j][1] = cnum[j];
    }
    free(pivots);
    free(exponents);
    free(cnum);
    polymat_clear(M);
    abclear(ab, &tmp);

    return fdesc;
}				/*}}} */

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


void compute_final_F_red(bmstatus_ptr bm, polymat f, unsigned int (*fdesc)[2], unsigned int t0, polymat pi, unsigned int * delta)/*{{{*/
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

    polymat_init(f, n, n, flen);

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
                        polymat_coeff(f, i, j, k),
                        polymat_coeff(pi, i, j1, k));
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
                        polymat_coeff(f, c, j, k+t0-e),
                        polymat_coeff(f, c, j, k+t0-e),
                        polymat_coeff(pi, i+n, j1, k));
            }
        }
    }

    free(pi_colidx);
}/*}}}*/


void write_f(bmstatus_ptr bm, const char * filename, polymat f_red, unsigned int * delta, int ascii)/*{{{*/
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
    if (ascii) {
        for(unsigned int k = 0 ; k < flen ; k++) {
            for(unsigned int i = 0 ; i < n ; i++) {
                for(unsigned int jj = 0 ; jj < n ; jj++) {
                    if (jj) fprintf(f, " ");
                    unsigned int j = sols[jj];
                    if (k <= delta[j]) {
                        abfprint(ab, f, polymat_coeff(f_red, i, jj, delta[j]-k));
                    } else {
                        printf("0");
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
            for(unsigned int i = 0 ; i < n ; i++) {
                for(unsigned int jj = 0 ; jj < n ; jj++) {
                    unsigned int j = sols[jj];
                    abset_zero(ab, tmp);
                    if (k <= delta[j])
                        abset(ab, tmp, polymat_coeff(f_red, i, jj, delta[j]-k));
                    fwrite(tmp, sizeof(abelt), 1, f);
                }
            }
        }
        abclear(ab, &tmp);
    }
    fclose(f);
}/*}}}*/


void compute_initial_E(bmstatus_ptr bm, polymat E, polymat A, unsigned int (*fdesc)[2])/*{{{*/
{
    // F0 is exactly the n x n identity matrix, plus the
    // X^(s-exponent)e_{cnum} vectors. fdesc has the (exponent, cnum)
    // pairs
    dims * d = bm->d;
    unsigned int m = d->m;
    unsigned int n = d->n;
    unsigned int b = m + n;
    unsigned int t0 = bm->t;
    ASSERT(A->m == m);
    ASSERT(A->n == n);

    ASSERT(!E->m && !E->n && !E->alloc);
    abdst_field ab = d->ab;
    /* Now we're ready to compute E, which is nothing more than a rewrite
     * of A, of course. */
    polymat_init(E, m, b, A->size - t0);
    for(unsigned int k = t0 ; k < A->size ; k++) {
        for(unsigned int j = 0 ; j < n ; j++) {
            /* Take column j of A, shifted by t0 positions */
            bwmat_copy_coeffs(ab,
                    polymat_part(E, 0, j, k-t0), b,
                    polymat_part(A, 0, j, k), n,
                    m);
        }
        for(unsigned int j = n ; j < m + n ; j++) {
            /* Take column cnum[j-n] of coeff exponents[j-n] of A, to
             * reach position t0. Which means that since we're
             * effectively computing E from coefficient t0 onwards, we
             * shift by exponents[j-n] coefficients */
            unsigned int c = fdesc[j-n][1];
            unsigned int e = fdesc[j-n][0];
            bwmat_copy_coeffs(ab,
                    polymat_part(E, 0, j, k-t0), b,
                    polymat_part(A, 0, c, k-t0+e), n,
                    m);
        }
    }
    E->size = A->size - t0;
}/*}}}*/


void read_data_for_series(bmstatus_ptr bm, polymat A, /* {{{ */
			  const char *input_file, int ascii_input)
{
    dims * d = bm->d;
    unsigned int m = d->m;
    unsigned int n = d->n;
    abdst_field ab = d->ab;

    unsigned int guess_len = 1000;

    ASSERT(!A->m && !A->n && !A->alloc);
    polymat_init(A, m, n, guess_len);

    FILE *f = fopen(input_file, ascii_input ? "r" : "rb");
    DIE_ERRNO_DIAG(f == NULL, "fopen", input_file);

    unsigned int k = 0;
    int eof_met = 0;
    for( ; !eof_met ; k++) {
        if (k == A->alloc) {
            polymat_realloc(A, A->alloc + A->alloc / 10);
        }

	/* coefficient 0 will be read to position 0 (k-!!k=0), but all
	 * other coefficients from position 1 onwards will be stored to
	 * position k-1 (!!k=1), effectively chopping the first
	 * coefficient.
	 */
	int k1 = k - ! !k;
	for (unsigned int i = 0; i < m && !eof_met ; i++) {
	    for (unsigned int j = 0; j < n && !eof_met ; j++) {
                abdst_elt x = polymat_coeff(A, i, j, k1);
		int rc;
                if (ascii_input) {
                    rc = abfscan(ab, f, x);
                    rc = rc == 1;
                } else {
                    rc = fread(x, sizeof(abelt), 1, f);
                    rc = rc == 1;
                    abreduce(ab, x, x);
                }
		if (!rc) {
                    if (i == 0 && j == 0) {
                        eof_met = 1;
                        break;
                    }
		    fprintf(stderr,
			    "Parse error in %s while reading coefficient (%d,%d,%d)\n",
			    input_file, i, j, k);
		    exit(1);
		}
	    }
	}
    }
    /* We can bail out only via eof_met, so this incurs an extra k++ */
    k--;
    fclose(f);
    printf("Using A(X) div X in order to consider Y as starting point\n");
    A->size = k - 1;
} /* }}} */


void usage()
{
    fprintf(stderr, "Usage: ./plingen [options, to be documented]\n");
    fprintf(stderr,
	    "General information about bwc optiosn is in the README file\n");
    exit(1);
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
    gmp_randstate_t rstate;

    MPI_Init(&argc, &argv);
    int rank;
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    gmp_randinit_default(rstate);
    gmp_randseed_ui(rstate, 1);

    setvbuf(stdout, NULL, _IONBF, 0);
    setvbuf(stderr, NULL, _IONBF, 0);

    /* {{{ Parameter list. Feed the common bw[] struct */
    param_list pl;
    param_list_init(pl);

    param_list_configure_switch(pl, "--tune", &tune);
    param_list_configure_switch(pl, "--ascii", &ascii);
    bw_common_init(bw, pl, &argc, &argv);

    const char * afile = param_list_lookup_string(pl, "afile");

    if (bw->m == -1) {
	fprintf(stderr, "no m value set\n");
	exit(1);
    }
    if (bw->n == -1) {
	fprintf(stderr, "no n value set\n");
	exit(1);
    }
    if (!tune && !afile) {
        fprintf(stderr, "No afile provided\n");
        exit(1);
    }
    if (!param_list_lookup_string(pl, "prime")) {
	fprintf(stderr, "no prime set\n");
	exit(1);
    }

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

    bm->lingen_threshold = 10;
    bm->lingen_mpi_threshold = 1000;
    param_list_parse_uint(pl, "lingen-threshold", &(bm->lingen_threshold));
    param_list_parse_uint(pl, "lingen-mpi-threshold", &(bm->lingen_mpi_threshold));

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
            exit(1);
        }
#endif
        int thr[2] = {1,1};
        param_list_parse_intxint(pl, "thr", thr);
        /* Asssume numbering in MPI_COMM_WORLD is all mpi jobs from the
         * same node together. So we pick them by bunches of size
         * thr[0]*thr[1].
         */
        printf("size=%d mpi=%dx%d thr=%dx%d\n", size, mpi[0], mpi[1], thr[0], thr[1]);
        ASSERT_ALWAYS(size == mpi[0] * mpi[1] * thr[0] * thr[1]);
        /* Keep the semantics which exist for krylov and so on --
         * even though we are not really using threads here. */
        bm->mpi_dims[0] *= thr[0];
        bm->mpi_dims[1] *= thr[1];
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

    if (param_list_warn_unused(pl))
	usage();
    /* }}} */

    /* TODO: delegate them to other functions, elsewhere... */
    int tune_bm_basecase = tune;
    int tune_mp = tune;

    if (tune_bm_basecase) {
        unsigned int maxtune = 10000 / (m * n);
        for(unsigned int k = 10 ; k < maxtune ; k += k/10) {
            unsigned int * delta = malloc((m + n) * sizeof(unsigned int));
            polymat E, pi;
            polymat_init(pi, 0, 0, 0);
            polymat_init(E, m, m+n, maxtune);
            E->size = k;
            for(unsigned int v = 0 ; v < E->m * E->n * E->size ; v++) {
                abrandom(ab, E->x[v], rstate);
            }
            double tt = seconds();
            for(unsigned int j = 0 ; j < m + n ; delta[j++]=1);
            bm->t = 1;
            bw_lingen_basecase(bm, pi, E, delta);
            printf("%zu %.2e\n", E->size, (seconds()-tt) / (k * k));
            polymat_clear(pi);
            polymat_clear(E);
            free(delta);
        }
    }

    if (tune_mp) {
        /* Now for benching mp and plain mul */
        unsigned int maxtune = 10000 / (m*n);
        /* Bench the products which would come together with a k-steps
         * basecase algorithm. IOW a 2k, one-level recursive call incurs
         * twice the k-steps basecase, plus once the timings counted here
         * (presently, this has Karatsuba complexity)
         */
        for(unsigned int k = 10 ; k < maxtune ; k+= k/10) {
            polymat E, piL, piR, pi, Er;
            unsigned int sE = k*(m+2*n)/(m+n);
            unsigned int spi = k*m/(m+n);
            polymat_init(E, m, m+n, sE);
            polymat_init(piL, m+n, m+n, spi);
            polymat_init(piR, m+n, m+n, spi);
            polymat_init(pi, m+n, m+n, spi*2);
            polymat_init(Er, m, m+n, sE-spi+1);
            E->size = sE;
            for(unsigned int v = 0 ; v < E->m * E->n * E->size ; v++) {
                abrandom(ab, E->x[v], rstate);
            }
            piL->size = spi;
            piR->size = spi;
            for(unsigned int v = 0 ; v < piL->m * piL->n * piL->size ; v++) {
                abrandom(ab, piL->x[v], rstate);
                abrandom(ab, piR->x[v], rstate);
            }
            double ttmp = 0, ttmul = 0;
            ttmp -= seconds();
            polymat_mp(ab, Er, E, piL);
            ttmp += seconds();
            ttmul -= seconds();
            polymat_mul(ab, pi, piL, piR);
            ttmul += seconds();
            double ttmpq = ttmp / (k*k);
            double ttmulq = ttmul / (k*k);
            double ttmpk = ttmp / pow(k, 1.58);
            double ttmulk = ttmul / pow(k, 1.58);
            printf("%u [%.2e+%.2e = %.2e] [%.2e+%.2e = %.2e]\n",
                    k,
                    ttmpq, ttmulq, ttmpq + ttmulq,
                    ttmpk, ttmulk, ttmpk + ttmulk
                    );
            // (seconds()-tt) / (k*k)); // ((sE-spi) * spi) / (m*(m+n)*(m+n)));
            // printf("%zu %.2e\n", E->size, (seconds()-tt) / (k*k)); // (spi * spi) / ((m+n)*(m+n)*(m+n)));
            polymat_clear(E);
            polymat_clear(piL);
            polymat_clear(piR);
            polymat_clear(pi);
            polymat_clear(Er);
        }
    }

    if (tune)
        return 0;

    unsigned int (*fdesc)[2] = NULL;    /* gcc is stupid */

    /* For the moment we read the complete thing in local memory before
     * dispatching, which is admittedly somewhat wasteful... Of course we
     * should rather do an mpi-level read. Either with only one
     * jobreading, and handing over data immediately to its peers, or
     * with several accesses to the file(s). */
    polymat E;
    polymat_init(E, 0, 0, 0);

    if (rank == 0) { /* {{{ Read A, compute F0 and E, and keep only E */
        polymat A;
        polymat_init(A, 0, 0, 0);
        printf("Reading scalar data in polynomial ``a'' from %s\n", afile);
        read_data_for_series(bm, A, afile, ascii);

        printf("Read %zu+1=%zu iterations",
                A->size, A->size+ 1);
        if (bw->end || bw->start) {
            printf(" (bw parameters: expect %u)",
                    bw->end - bw->start);
        }
        printf(".\n");
        /* Data read stage completed. */

        fdesc = compute_initial_F(bm, A);

        compute_initial_E(bm, E, A, fdesc);

        printf("Throwing out a(X)\n");
        polymat_clear(A);
    } /* }}} */
    MPI_Bcast(&(bm->t), 1, MPI_UNSIGNED, 0, bm->world);
    /* This will quite probably be changed. We are playing nasty games
     * here, with E->size being used to draw decisions even though the is
     * no corresponding allocation.
     */
    MPI_Bcast(&(E->size), 1, MPI_UNSIGNED, 0, bm->world);

    // bw_bbpoly piL, piR, piP;

    unsigned int t0 = bm->t;

    unsigned int * delta = malloc((m + n) * sizeof(unsigned int));
    for(unsigned int j = 0 ; j < m + n ; j++) {
        delta[j] = t0;
    }

    polymat pi;
    polymat_init(pi, 0, 0, 0);
    /* At this point, we're not propagating the mpi info here */
    bw_lingen(bm, pi, E, delta);
    polymat_clear(E);

    unsigned int nlucky = 0;
    int luck_mini = expected_pi_length(d, 0);
    for(unsigned int j = 0 ; j < b ; nlucky += bm->lucky[j++] >= luck_mini) ;

    if (rank == 0) {
        /* TODO: consider luck only below probability 2^-64 */
        if (nlucky == n) {
            polymat f_red;
            polymat_init(f_red, 0, 0, 0);
            compute_final_F_red(bm, f_red, fdesc, t0, pi, delta);
            char * f_filename;
            int rc = asprintf(&f_filename, "%s.gen", afile);
            ASSERT_ALWAYS(rc >= 0);
            write_f(bm, f_filename, f_red, delta, ascii);
            free(f_filename);
            polymat_clear(f_red);
        } else {
            fprintf(stderr, "Could not find the required set of solutions (nlucky=%u)\n", nlucky);
        }

        printf("t_basecase = %.2f\n", bm->t_basecase);
        printf("t_mp = %.2f\n", bm->t_mp);
        free(delta);
        polymat_clear(pi);
        free(fdesc);
    }

#if 0
    struct e_coeff *ec;
    struct t_poly *pi_left, *pi_right, *pi_prod;

    printf("ec->degree=%d, new_t=%d, t_counter=%d\n",
	   ec->degree, new_t, t_counter);

    global_sum_delta = sum_delta(d, global_delta);

    if (ec->degree >= 0) {
	bw_lingen(ec, global_delta, &pi_right);
    } else {
	pi_right = pi_left;
	pi_left = NULL;
    }

    compute_f_final(pi_prod);
    bw_commit_f(f_poly, global_delta);
    print_chance_list(bm->len, chance_list);
    /* Now I can clean up everything if I want... if I want... */
    ft_order_t::cleanup();
#endif

    abfield_clear(ab);
    bmstatus_clear(bm);
    bw_common_clear(bw);
    param_list_clear(pl);
    gmp_randclear(rstate);

    MPI_Finalize();
    return 0;
}

/* vim:set sw=4 sta et: */
