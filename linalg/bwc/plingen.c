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

#include "macros.h"
#include "utils.h"
#include "abase.h"
#include "bw-common.h"		/* Handy. Allows Using global functions for recovering
				   parameters */
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

struct bmstatus_s {
    dims d[1];
    int * delta;
    int * lucky;
    // bw_nbpoly E;
    // bw_nbmat e0;        /* Useful ??? */
    int t;
    unsigned int len;
    /* Not clear we need the following: t0, F0. */
};
typedef struct bmstatus_s bmstatus[1];
typedef struct bmstatus_s *bmstatus_ptr;

/* External debug code I use is sensible to zero matrices appearing in
 * polynomial expansions, and uses this for detecting termination of the
 * expansion. Thus for gdb'ing the code, I set 1 below.
 */
#define DEBUG_EXTRA_ROOM        0

#if 0
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
    f = fopen(filename, "w");
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
	f = fopen(filename, "r");
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

#endif

/* {{{ col sorting */
/* XXX FIXME.
 *
 * The text below is plain wrong we should not be relying on the columns
 * of pi be in any degree order.
 */
/* It might seem merely cosmetic and useless to sort w.r.t both the
 * global and local nominal degrees. In fact, it is crucial for the
 * corectness of the computations. (Imagine a 2-step increase, starting
 * with uneven global deltas, and hitting an even situation in the
 * middle. One has to sort out the local deltas to prevent trashing the
 * whole picture).
 *
 * The positional sort, however, *is* cosmetic (makes debugging easier).
 */

typedef int (*sortfunc_t) (const void*, const void*);

static int col_cmp(const int x[3], const int y[3])
{
    for(int i = 0 ; i < 3 ; i++) {
        int d = x[i] - y[i];
        if (d) return d;
    }
    return 0;
}

/* }}} */

#if 0
static int bw_check_chance(bw_mbmat e, unsigned int *clist)
{
    int i, j;
    unsigned int maxchance;

    maxchance = 0;

    for (j = 0; j < bigdim; j++) {
	for (i = 0; i < m_param; i++) {
	    bw_reduce_short_scalar(mbmat_scal(e, i, j));
	}
    }

    for (j = 0; j < bigdim; j++) {
	if (mcol_is_zero(mbmat_col(e, clist ? clist[j] : j))) {
	    if (++chance_list[j] > maxchance)
		maxchance = chance_list[j];
	    printf("Column %d happens to be zero ! (%d)\n",
		   j, chance_list[j]);
	    if (t_counter < 30) {
		fprintf(stderr, "surely a degenerate case :-(\n");
		BUG();
	    }
	} else
	    chance_list[j] = 0;
    }

    return maxchance;
}
#endif

static inline unsigned int expected_pi_length(dims * d, unsigned int len)
{
    /* The idea is that we want something which may account for something
     * exceptional, bounded by probability 2^-64. This corresponds to a
     * column in e (matrix of size m*b) to be spontaneously equal to
     * zero. This happens with probability (#K)^-m. Now (2^64)^sizeof(abelt) is
     * roughly #K^abgroupsize(ab). Thus the solution to
     * (#K)^(-m*x) > 2^-64
     * is (#K)^(-m*x*sizeof(abelt)) > #K^(-abgroupsize(ab))
     * i.e. m*x*sizeof(abelt) < abgroupsize(ab).
     */
    unsigned int m = d->m;
    unsigned int n = d->n;
    unsigned int b = m + n;
    abdst_field ab MAYBE_UNUSED = d->ab;
    return iceildiv(len * m, b) + iceildiv(abgroupsize(ab), m * sizeof(abelt));
}

/* Forward declaration, it's used by the recursive version */
// static double bw_lingen(struct e_coeff *, int *, struct t_poly **);

/* This destructively cancels the first len coefficients of E, and
 * computes the appropriate matrix pi which achieves this. The
 * elimination is done in accordance with the nominal degrees found in
 * delta.
 *
 * The result is expected to have degree ceil(len*m/b) coefficients, so
 * that E*pi is divisible by X^len.
 */

static bw_bbpoly
bw_lingen_basecase(bmstatus_ptr bm, bw_mbpoly E, unsigned int len, int *delta)
{
    dims * d = bm->d;
    unsigned int m = d->m;
    unsigned int n = d->n;
    unsigned int b = m + n;
    abdst_field ab = d->ab;

    /* Allocate something large enough for the result. This will be
     * soon freed anyway. Set it to identity. */
    unsigned int pi_room = expected_pi_length(d, len);
    bw_bbpoly pi = polymat_alloc(b, b, pi_room + DEBUG_EXTRA_ROOM);

    /* Keep a permutation on [0..b[ indicating which is the physical
     * position of column j in the matrix pi. Also keep track of the
     * number of coefficients for the columns of pi. Set pi to Id */
    unsigned int *pi_perm = malloc(b * sizeof(unsigned int));
    unsigned int *pi_lengths = malloc(b * sizeof(unsigned int));
    polymat_zero(pi, b, b, pi_room + DEBUG_EXTRA_ROOM);
    for(unsigned int i = 0 ; i < b ; i++) {
        abset_ui(ab, polymat_coeff(pi, b, b, i, i, 0), 1);
        pi_perm[i] = i;
        pi_lengths[i] = 1;
    }

    /* Keep a list of columns which have been used as pivots at the
     * previous iteration */
    unsigned int * pivots = malloc(m * sizeof(unsigned int));
    int * is_pivot = malloc(b * sizeof(int));
    memset(is_pivot, 0, b * sizeof(int));

    bw_mbmat e = polymat_alloc(m, b, 1);
    polymat_zero(e, m, b, 1);

    for (unsigned int t = 0; t < len ; t++) {

        /* {{{ Update the columns of e for degree t. Save computation
         * time by not recomputing those which can easily be derived from
         * previous iteration. Notice that the columns of e are exactly
         * at the physical positions of the corresponding columns of pi.
         */

        abelt_ur tmp_ur;
        abelt_ur_init(ab, &tmp_ur);
        bw_mbmat_ur e_ur = polymat_alloc_ur(m, b, 1);
        polymat_zero_ur(e_ur, m, b, 1);
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
                                polymat_coeff(E, m, b, i, k, t - s),
                                polymat_coeff(pi, b, b, k, j, s));
                        abelt_ur_add(ab,
                                polymat_coeff_ur(e_ur, m, b, i, j, 0),
                                polymat_coeff_ur(e_ur, m, b, i, j, 0),
                                tmp_ur);
                    }
                }
            }
        }
        int newluck = 0;
        for(unsigned int j = 0 ; j < b ; j++) {
            if (is_pivot[j]) continue;
            unsigned int nz = 0;
            for(unsigned int i = 0 ; i < m ; i++) {
                abreduce(ab,
                        polymat_coeff(e, m, b, i, j, 0),
                        polymat_coeff_ur(e_ur, m, b, i, j, 0)
                        );
                nz += abcmp_ui(ab, polymat_coeff(e, m, b, i, j, 0), 0) == 0;
            }
            if (nz == m) {
                newluck++, bm->lucky[j]++;
            } else if (bm->lucky[j] > 0) {
                bm->lucky[j] = 0;
            }
        }
        abelt_ur_clear(ab, &tmp_ur);
        polymat_free_ur(e_ur /* , m, b, 1 */);
        if (newluck) {
            printf("t=%d, canceled columns:", t + bm->t);
            for(unsigned int j = 0 ; j < b ; j++) {
                if (bm->lucky[j] > 0)
                    printf(" %u", j);
            }
            printf("\n");
        }
        /* }}} */

        /* {{{ Now see in which order I may look at the columns of pi, so
         * as to keep the nominal degrees correct. In contrast with what
         * we used to do before, we no longer apply the permutation to
         * delta. So the delta[] array keeps referring to physical
         * indices, and we'll tune this in the end. */
        int (*ctable)[3] = malloc(b * 3 * sizeof(int));
        for(unsigned int j = 0; j < b; j++) {
            ctable[j][0] = delta[j];
            ctable[j][1] = pi_lengths[j];
            ctable[j][2] = j;
        }
        qsort(ctable, b, 3 * sizeof(int), (sortfunc_t) & col_cmp);
        for(unsigned int j = 0; j < b; j++) {
            pi_perm[j] = ctable[j][2];
        }
        free(ctable);
        /* }}} */

        /* {{{ Now do Gaussian elimination */
        memset(is_pivot, 0, b * sizeof(int));
        unsigned int r = 0;
        /* Loop through logical indices */
        for(unsigned int jl = 0; jl < b; jl++) {
            unsigned int j = pi_perm[jl];
            unsigned int u = 0;
            /* {{{ Find the pivot */
            for( ; u < m ; u++) {
                if (abcmp_ui(ab, polymat_coeff(e, m, b, u, j, 0), 0) != 0)
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
            abinv(ab, inv, polymat_coeff(e, m, b, u, j, 0));
            abneg(ab, inv, inv);
            for (unsigned int kl = jl + 1; kl < b ; kl++) {
                unsigned int k = pi_perm[kl];
                if (abcmp_ui(ab, polymat_coeff(e, m, b, u, k, 0), 0) == 0)
                    continue;
                // add lambda = e[u,k]*-e[u,j]^-1 times col j to col k.
                abelt lambda;
                abinit(ab, &lambda);
                abmul(ab, lambda, inv, polymat_coeff(e, m, b, u, k, 0));

                assert(delta[j] <= delta[k]);
                /* {{{ Apply on both e and pi */
                abelt tmp;
                abinit(ab, &tmp);
                for(unsigned int i = 0 ; i < m ; i++) {
                    /* TODO: Would be better if mpfq had an addmul */
                    abmul(ab, tmp, lambda, polymat_coeff(e, m, b, i, j, 0));
                    abadd(ab,
                            polymat_coeff(e, m, b, i, k, 0),
                            polymat_coeff(e, m, b, i, k, 0),
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
                    /* XXX FIXME. The assert() below is wrong. */
                    /* Not clear whether this assert() is correct or not.
                     * If it isn't, then why do we sort w.r.t this in
                     * ctable ? */
                    ASSERT_ALWAYS(pi_lengths[j] <= pi_lengths[k]);
                    /* If we are to allow length[j] > length[k], then we
                     * should update length[k], of course */
                    for(unsigned int s = 0 ; s < pi_lengths[j] ; s++) {
                        /* TODO: Would be better if mpfq had an addmul */
                        abmul(ab, tmp, lambda,
                                polymat_coeff(pi, b, b, i, j, s));
                        abadd(ab,
                                polymat_coeff(pi, b, b, i, k, s),
                                polymat_coeff(pi, b, b, i, k, s),
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

        ASSERT_ALWAYS(r == m);

        /* {{{ Now for all pivots, multiply column in pi by x */
        for (unsigned int j = 0; j < b ; j++) {
            if (!is_pivot[j]) continue;
            if (pi_lengths[j] >= pi_room) {
                ASSERT_ALWAYS(bm->lucky[j] <= 0);
                if (bm->lucky[j] == 0)
                    printf("Column %u discarded from now on\n", j);
                bm->lucky[j] = -1;
                pi_lengths[j]++;
                delta[j]++;
                continue;
            }
            bwmat_move_coeffs(ab,
                    polymat_part(pi, b, b, 0, j, 1), b,
                    polymat_part(pi, b, b, 0, j, 0), b,
                    b * pi_lengths[j]);
            for(unsigned int i = 0 ; i < b ; i++) {
                abset_ui(ab, polymat_coeff(pi, b, b, i, j, 0), 0);
            }
            pi_lengths[j]++;
            delta[j]++;
        }
        /* }}} */

        /*
        printf("t=%u:", bm->t + t);
        for(unsigned int j = 0; j < b; j++) {
            printf(" %u", delta[j]);
        }
        printf("\n");
        */
    }
        printf("t=%u: delta =", bm->t + len);
        for(unsigned int j = 0; j < b; j++) {
            printf(" %u", delta[j]);
        }
        printf("\n");
        printf("t=%u: pi_length =", bm->t + len);
        for(unsigned int j = 0; j < b; j++) {
            printf(" %u", pi_lengths[j]);
        }
        printf("\n");
    polymat_free(e /*, m, b, 1 */);
    free(is_pivot);
    free(pivots);
    free(pi_perm);
    free(pi_lengths);   /* What shall we do with this one ??? */

    return pi;
}


#if 0
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

#endif

/**********************************************************************/

unsigned int (*compute_initial_F(bmstatus_ptr bm, bw_mnpoly A))[2] /*{{{ */
{				
    dims *d = bm->d;
    abdst_field ab = d->ab;
    unsigned int m = d->m;
    unsigned int n = d->n;
    unsigned int len = bm->len;
    abelt tmp;
    abinit(ab, &tmp);

    /* First try to create the initial F matrix */
    printf("Computing t0\n");

    /* We want to create a full rank m*m matrix M, by extracting columns
     * from the first coefficients of A */

    bw_mmmat M = polymat_alloc(m, m, 1);
    /* For each integer i between 0 and m-1, we have a column, picked
     * from column cnum[i] of coeff exponent[i] of A which, once reduced modulo
     * the other ones, has coefficient at row pivots[i] unequal to zero.
     */
    unsigned int *pivots = malloc(m * sizeof(unsigned int));
    unsigned int *exponents = malloc(m * sizeof(unsigned int));
    unsigned int *cnum = malloc(m * sizeof(unsigned int));
    unsigned int r = 0;

    for (unsigned int k = 0; r < m && k < len; k++) {
	for (unsigned int j = 0; r < m && j < n; j++) {
	    /* Extract a full column into M */
	    bwmat_copy_coeffs(ab,
			      polymat_part(M, m, m, 0, r, 0),
			      m, polymat_part(A, m, n, 0, j, k), n, m);

            /* Now reduce it modulo all other columns */
	    for (unsigned int v = 0; v < r; v++) {
		unsigned int u = pivots[v];
		/* the v-th column in the M is known to
		 * kill coefficient u (more exactly, to have a -1 as u-th
		 * coefficient, and zeroes for the other coefficients
		 * referenced in the pivots[0] to pivots[v-1] indices).
		 */
                for(unsigned int i = 0 ; i < m ; i++) {
                    if (i == u) continue;
                    abmul(ab, tmp,
                              polymat_coeff(M, m, m, i, v, 0),
                              polymat_coeff(M, m, m, u, r, 0));
                    abadd(ab, polymat_coeff(M, m, m, i, v, 0),
                              polymat_coeff(M, m, m, i, v, 0),
                              tmp);
                }
                abset_zero(ab,
                        polymat_coeff(M, m, m, u, r, 0));
	    }
            unsigned int u = 0;
            for( ; u < m ; u++) {
                if (abcmp_ui(ab, polymat_coeff(M, m, m, u, r, 0), 0) != 0)
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
            abinv(ab, tmp, polymat_coeff(M, m, m, u, r, 0));
            abneg(ab, tmp, tmp);
            for(unsigned int i = 0 ; i < m ; i++) {
                abmul(ab, polymat_coeff(M, m, m, i, r, 0),
                          polymat_coeff(M, m, m, i, r, 0),
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
 
    for(unsigned int j = 0 ; j < m + n ; j++) {
        bm->delta[j] = t0;
    }
    bm->t = t0;

    unsigned int (*fdesc)[2] = malloc(2 * m * sizeof(unsigned int));
    for(unsigned int j = 0 ; j < m ; j++) {
        fdesc[j][0] = exponents[j];
        fdesc[j][1] = cnum[j];
    }
    free(pivots);
    free(exponents);
    free(cnum);
    polymat_free(M /* , m, m, 1 */);
    abclear(ab, &tmp);

    return fdesc;
}				/*}}} */

unsigned int get_max_delta_on_solutions(bmstatus_ptr bm)/*{{{*/
{
    dims * d = bm->d;
    unsigned int m = d->m;
    unsigned int n = d->n;
    unsigned int maxdelta = 0;
    for(unsigned int j = 0 ; j < m + n ; j++) {
        if (bm->lucky[j] <= 0)
            continue;
        unsigned int delta = bm->delta[j];
        if (delta > maxdelta)
            maxdelta = delta;
    }
    return maxdelta;
}/*}}}*/


bw_nnpoly compute_final_F_red(bmstatus_ptr bm, unsigned int (*fdesc)[2], unsigned int t0, bw_bbpoly pi)/*{{{*/
{
    dims * d = bm->d;
    unsigned int m = d->m;
    unsigned int n = d->n;
    unsigned int b = m + n;
    abdst_field ab = d->ab;

    unsigned int maxdelta = get_max_delta_on_solutions(bm);

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
    bw_nbpoly f = polymat_alloc(n, n, flen);
    polymat_zero(f, n, n, flen);

    printf("Computing value of f(X)=f0(X)pi(X) (degree %u)\n", maxdelta);
    printf("t0=%u pilen=%u\n", t0, pilen);

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
                        polymat_coeff(f, n, n, i, j, k),
                        polymat_coeff(pi, b, b, i, j1, k));
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
                        polymat_coeff(f, n, n, c, j, k+t0-e),
                        polymat_coeff(f, n, n, c, j, k+t0-e),
                        polymat_coeff(pi, b, b, i+n, j1, k));
            }
        }
    }

    free(pi_colidx);

    return f;
}/*}}}*/


void write_f(bmstatus_ptr bm, bw_nnpoly f_red)
{
    dims * d = bm->d;
    unsigned int n = d->n;
    abdst_field ab = d->ab;
    FILE * f = fopen(LINGEN_F_FILE, "w");
    DIE_ERRNO_DIAG(f == NULL, "fopen", LINGEN_F_FILE);
    unsigned int maxdelta = get_max_delta_on_solutions(bm);
    unsigned int flen = maxdelta + 1;
    for(unsigned int k = 0 ; k < flen ; k++) {
        for(unsigned int i = 0 ; i < n ; i++) {
            for(unsigned int j = 0 ; j < n ; j++) {
                if (j) fprintf(f, " ");
                abfprint(ab, f, polymat_coeff(f_red, n, n, i, j, k));
            }
            fprintf(f, "\n");
        }
        fprintf(f, "\n");
    }
    fclose(f);
}


bw_mbpoly compute_initial_E(bmstatus_ptr bm, bw_mnpoly A, unsigned int (*fdesc)[2])/*{{{*/
{
    // F0 is exactly the n x n identity matrix, plus the
    // X^(s-exponent)e_{cnum} vectors. fdesc has the (exponent, cnum)
    // pairs
    dims * d = bm->d;
    unsigned int m = d->m;
    unsigned int n = d->n;
    unsigned int b = m + n;
    unsigned int len = bm->len;
    unsigned int t0 = bm->t;
    abdst_field ab = d->ab;
    /* Now we're ready to compute E, which is nothing more than a rewrite
     * of A, of course. */
    bw_mbpoly E = polymat_alloc(m, b, bm->len - t0 + DEBUG_EXTRA_ROOM);
    polymat_zero(E, m, b, bm->len - t0 + DEBUG_EXTRA_ROOM);
    for(unsigned int k = t0 ; k < len ; k++) {
        for(unsigned int j = 0 ; j < n ; j++) {
            /* Take column j of A, shifted by t0 positions */
            bwmat_copy_coeffs(ab,
                    polymat_part(E, m, b, 0, j, k-t0), b,
                    polymat_part(A, m, n, 0, j, k), n,
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
                    polymat_part(E, m, b, 0, j, k-t0), b,
                    polymat_part(A, m, n, 0, c, k-t0+e), n,
                    m);
        }
    }
    return E;
}/*}}}*/


void read_data_for_series(bmstatus_ptr bm, bw_mnpoly a_poly, /* {{{ */
			  const char *input_file)
{
    dims * d = bm->d;
    unsigned int m = d->m;
    unsigned int n = d->n;
    abdst_field ab = d->ab;
    FILE *f = fopen(input_file, "r");
    DIE_ERRNO_DIAG(f == NULL, "fopen", input_file);
    for (unsigned int k = 0; k < bm->len; k++) {
	/* coefficient 0 will be read to position 0 (k-!!k=0), but all
	 * other coefficients from position 1 onwards will be stored to
	 * position k-1 (!!k=1), effectively chopping the first
	 * coefficient.
	 */
	int k1 = k - ! !k;
	for (unsigned int i = 0; i < m; i++) {
	    for (unsigned int j = 0; j < n; j++) {
		int rc =
		    abfscan(ab, f, polymat_coeff(a_poly, m, n, i, j, k1));
		if (rc != 1) {
		    fprintf(stderr,
			    "Parse error in %s while reading coefficient (%d,%d,%d)\n",
			    input_file, i, j, k);
		    exit(1);
		}
	    }
	}
    }
    fclose(f);
    printf("Using A(X) div X in order to consider Y as starting point\n");
    bm->len--;
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
    bm->delta = malloc((m + n) * sizeof(int));
    bm->lucky = malloc((m + n) * sizeof(int));
    memset(bm->delta, 0, (m + n) * sizeof(int));
    memset(bm->lucky, 0, (m + n) * sizeof(int));
}/*}}}*/

void bmstatus_clear(bmstatus_ptr bm)/*{{{*/
{
    free(bm->delta);
    free(bm->lucky);
    memset(bm, 0, sizeof(bmstatus));
}/*}}}*/

int main(int argc, char *argv[])
{
    bmstatus bm;
    dims * d = bm->d;

    setvbuf(stdout, NULL, _IONBF, 0);
    setvbuf(stderr, NULL, _IONBF, 0);

    /* {{{ Parameter list. Feed the common bw[] struct */
    param_list pl;
    param_list_init(pl);

    bw_common_init(bw, pl, &argc, &argv);

    if (bw->m == -1) {
	fprintf(stderr, "no m value set\n");
	exit(1);
    }

    bmstatus_init(bm, bw->m, bw->n);

    unsigned int m = d->m;
    unsigned int n = d->n;
    unsigned int b = m + n;
    abdst_field ab = d->ab;

    if (!param_list_lookup_string(pl, "prime")) {
	fprintf(stderr, "no prime set\n");
	exit(1);
    }

    abfield_init(ab);
    {
	mpz_t p;
	mpz_init_set_ui(p, 2);
	param_list_parse_mpz(pl, "prime", p);
	abfield_specify(ab, MPFQ_PRIME_MPZ, p);
	mpz_clear(p);
    }

    if (param_list_warn_unused(pl))
	usage();
    param_list_clear(pl);
    /* }}} */

    char input_file[FILENAME_MAX];

    /* {{{ detect the input file -- there must be only one file. */
    unsigned int n0, n1, j0, j1;
    {
	DIR * dir = opendir(".");
	struct dirent * de;
	input_file[0] = '\0';
	for (; (de = readdir(dir)) != NULL;) {
	    int len;
	    int rc = sscanf(de->d_name, A_FILE_PATTERN "%n",
			    &n0, &n1, &j0, &j1, &len);
	    /* rc is expected to be 4 or 5 depending on our reading of the
	     * standard */
	    if (rc < 4 || len != (int) strlen(de->d_name)) {
		continue;
	    }
	    if (input_file[0] != '\0') {
		fprintf(stderr, "Found two possible file names %s and %s\n",
			input_file, de->d_name);
		exit(1);
	    }
	    size_t clen = MAX((size_t) len, sizeof(input_file));
	    memcpy(input_file, de->d_name, clen);
	}
	closedir(dir);
        if (!input_file[0]) {
            fprintf(stderr, "No input file found for pattern %s\n", A_FILE_PATTERN);
            exit(1);
        }
    }				/* }}} */

    /* {{{ More argument checking. Learn input size */
    if (bw->n == 0) {
	bw->n = n1 - n0;
    } else if (bw->n != (int) (n1 - n0)) {
	fprintf(stderr, "n value mismatch (config says %d, A file says %u)\n",
		bw->n, n1 - n0);
	exit(1);
    }

    if (bw->end == 0) {
	ASSERT_ALWAYS(bw->start == 0);
	bw->end = j1 - j0;
    } else if (bw->end - bw->start > (int) (j1 - j0)) {
	fprintf(stderr, "sequence file %s is too short\n", input_file);
	exit(1);
    }

    bm->len = bw->end - bw->start;
    /* }}} */

    unsigned int (*fdesc)[2];

    bw_mbpoly E;

    { /* {{{ Read A, compute F0 and E, and keep only E */
        bw_mnpoly a_poly = polymat_alloc(m, n, bm->len);
        polymat_zero(a_poly, m, n, bm->len);

        printf("Reading scalar data in polynomial ``a'' from %s\n", input_file);
        read_data_for_series(bm, a_poly, input_file);

        /* Data read stage completed. */

        fdesc = compute_initial_F(bm, a_poly);

        E = compute_initial_E(bm, a_poly, fdesc);

        printf("Throwing out a(X)\n");
        polymat_free(a_poly /* , m, n, bm->len */);
    } /* }}} */

    // bw_bbpoly piL, piR, piP;

    unsigned int len = bm->len - bm->t;
    unsigned int t0 = bm->t;

    bw_bbpoly piL = bw_lingen_basecase(bm, E, len, bm->delta);
    polymat_free(E /* , m, b, len */);
    bm->t += len;

    unsigned int nlucky = 0;
    for(unsigned int j = 0 ; j < b ; nlucky += bm->lucky[j++] > 0) ;

    if (nlucky == n) {
        bw_nnpoly f_red = compute_final_F_red(bm, fdesc, t0, piL);
        write_f(bm, f_red);
        polymat_free(f_red);
    } else {
        fprintf(stderr, "Could not find the required set of solutions (nlucky=%u)\n", nlucky);
    }

    polymat_free(piL /* , b, b, pilen */);
    free(fdesc);


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
    return 0;
}

/* vim:set sw=4 sta et: */
