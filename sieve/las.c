#include "cado.h"
#include <stdint.h>     /* AIX wants it first (it's a bug) */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h> /* for PRIx64 macro and strtoumax */
#include <math.h>   // for ceiling, floor in cfrac
#include <ctype.h>
#include <float.h>
#include <pthread.h>
#include <fcntl.h>   /* for _O_BINARY */
#include <stdarg.h> /* Required so that GMP defines gmp_vfprintf() */
#include "fb.h"
#include "portability.h"
#include "utils.h"           /* lots of stuff */
#include "ecm/facul.h"
#include "bucket.h"
#include "trialdiv.h"
#include "las-config.h"
#include "las-types.h"
#include "las-coordinates.h"
#include "las-debug.h"
#include "las-duplicate.h"
#include "las-report-stats.h"
#include "las-norms.h"
#include "las-unsieve.h"
#include "las-arith.h"
#include "las-qlattice.h"
#include "las-smallsieve.h"
#include "las-descent-helpers.h"
#include "las-cofactor.h"
#include "las-fill-in-buckets.h"

#ifdef HAVE_SSE41
/* #define SSE_SURVIVOR_SEARCH 1 */
#include <smmintrin.h>
#endif

// #define HILIGHT_START   "\e[01;31m"
// #define HILIGHT_END   "\e[00;30m"

#define HILIGHT_START   ""
#define HILIGHT_END   ""

/* This global mutex should be locked in multithreaded parts when a
 * thread does a read / write, especially on stdout, stderr...
 */
pthread_mutex_t io_mutex = PTHREAD_MUTEX_INITIALIZER;

int create_descent_hints = 0;
double tt_qstart;

#define BUCKET_REGION (1 << LOG_BUCKET_REGION)

/* {{{ for cofactorization statistics */
int stats = 0; /* 0: nothing, 1: write stats file, 2: read stats file,
                  the stats file can be used with gnuplot, for example:
                  splot "stats.dat" u 1:2:3, "stats.dat" u 1:2:4 */
double stats_prob = 2e-4;
FILE *stats_file;
FILE *sievestats_file;
uint32_t **cof_call; /* cof_call[r][a] is the number of calls of the
                        cofactorization routine with a cofactor of r bits on
                        the rational side, and a bits on the algebraic side */
uint32_t **cof_succ; /* cof_succ[r][a] is the corresponding number of
                        successes, i.e., of call that lead to a relation */
/* }}} */

/* }}} */

/*****************************/

/* siever_config stuff */

void siever_config_display(FILE * o, siever_config_srcptr sc)/*{{{*/
{
    fprintf(o, "# Sieving parameters for q~2^%d on the %s side\n",
            sc->bitsize, sidenames[sc->side]);
    /* Strive to keep these output lines untouched */
    fprintf(o,
	    "# Sieving parameters: rlim=%lu alim=%lu lpbr=%d lpba=%d\n",
	    sc->sides[RATIONAL_SIDE]->lim,
            sc->sides[ALGEBRAIC_SIDE]->lim,
            sc->sides[RATIONAL_SIDE]->lpb,
	    sc->sides[ALGEBRAIC_SIDE]->lpb);
    fprintf(o,
	    "#                     mfbr=%d mfba=%d rlambda=%1.1f alambda=%1.1f\n",
	    sc->sides[RATIONAL_SIDE]->mfb,
            sc->sides[ALGEBRAIC_SIDE]->mfb,
            sc->sides[RATIONAL_SIDE]->lambda,
	    sc->sides[ALGEBRAIC_SIDE]->lambda);
    fprintf(o, "#                     skewness=%1.1f\n",
	    sc->skewness);
}/*}}}*/

int siever_config_cmp(siever_config_srcptr a, siever_config_srcptr b)/*{{{*/
{
    return memcmp(a, b, sizeof(siever_config));
}/*}}}*/

/* siever_config_match only tells whether this config entry is one we
 * want to use for a given problem size/side */
int siever_config_match(siever_config_srcptr a, siever_config_srcptr b)/*{{{*/
{
    return a->side == b->side && a->bitsize == b->bitsize;
}/*}}}*/



/* sieve_info stuff */
static void sieve_info_init_trialdiv(sieve_info_ptr si, int side)/*{{{*/
{
    /* Our trial division needs odd divisors, 2 is handled by mpz_even_p().
       If the FB primes to trial divide contain 2, we skip over it.
       We assume that if 2 is in the list, it is the first list entry,
       and that it appears at most once. */

    /* XXX This function consider the full contents of the list
     * si->sides[side]->fb to be trialdiv primes, therefore that
     * sieve_info_split_bucket_fb_for_threads must have already been
     * called, so that the only primes which are still in s->fb are the
     * trial-divided ones */
    sieve_side_info_ptr s = si->sides[side];
    unsigned long pmax = MIN((unsigned long) si->conf->bucket_thresh,
                             trialdiv_get_max_p());
    s->trialdiv_primes = fb_extract_bycost (s->fb, pmax, si->conf->td_thresh);
    int n;
    for (n = 0; s->trialdiv_primes[n] != FB_END; n++);
    int skip2 = n > 0 && s->trialdiv_primes[0] == 2;
    s->trialdiv_data = trialdiv_init (s->trialdiv_primes + skip2, n - skip2);
}/*}}}*/

static void sieve_info_clear_trialdiv(sieve_info_ptr si, int side)/*{{{*/
{
        trialdiv_clear (si->sides[side]->trialdiv_data);
        free (si->sides[side]->trialdiv_primes);
}/*}}}*/

/* {{{ Factor base handling */
/* {{{ Initialize the factor bases */
void sieve_info_init_factor_bases(las_info_ptr las, sieve_info_ptr si, param_list pl)
{
    double tfb;
    const char * fbcfilename = param_list_lookup_string(pl, "fbc");

    if (fbcfilename != NULL) {
        /* Try to read the factor base cache file. If that fails, because
           the file does not exist or is not compatible with our parameters,
           it will be written after we generate the factor bases. */
        verbose_output_print(0, 1, "# Mapping memory image of factor base from file %s\n",
               fbcfilename);
        if (fb_mmap_fbc(&si->sides[0]->fb, &si->sides[0]->fb_bucket_threads,
                        &si->sides[1]->fb, &si->sides[1]->fb_bucket_threads,
                        fbcfilename, las->nb_threads)) {
            si->sides[0]->fb_is_mmapped = 1;
            si->sides[1]->fb_is_mmapped = 1;
            verbose_output_print(0, 1, "# Finished mapping memory image of factor base\n");
            return;
        } else {
            verbose_output_print(0, 1, "# Could not map memory image of factor base\n");
        }
    }

    for(int side = 0 ; side < 2 ; side++) {
        mpz_poly_ptr pol = las->cpoly->pols[side];
        sieve_side_info_ptr sis = si->sides[side];
        unsigned long lim = si->conf->sides[side]->lim;
        if (pol->deg > 1) {
            tfb = seconds ();
            char fbparamname[4];
            snprintf(fbparamname, sizeof(fbparamname), "fb%d", side);
            const char * fbfilename = param_list_lookup_string(pl, fbparamname);
            verbose_output_print(0, 1, "# Reading %s factor base from %s\n", sidenames[side], fbfilename);
            int ok = fb_read (&sis->fb, &sis->fb_bucket_threads, fbfilename,
                              si->conf->bucket_thresh, las->nb_threads,
                              lim, si->conf->sides[side]->powlim);
            FATAL_ERROR_CHECK(!ok, "Error reading factor base file");
            ASSERT_ALWAYS(sis->fb != NULL);
            sis->fb_is_mmapped = 0;
            tfb = seconds () - tfb;
            verbose_output_print(0, 1,
                    "# Reading %s factor base of %zuMb took %1.1fs\n",
                    sidenames[side],
                    fb_size (sis->fb) >> 20, tfb);
        } else {
            tfb = seconds ();
            int ok = fb_make_linear (&sis->fb, &sis->fb_bucket_threads,
                                     (const mpz_t *) pol->coeff, (fbprime_t) lim,
                                     si->conf->bucket_thresh, las->nb_threads,
                                     si->conf->sides[side]->powlim, 1);
            FATAL_ERROR_CHECK(!ok, "Error creating rational factor base");
            sis->fb_is_mmapped = 0;
            tfb = seconds () - tfb;
            verbose_output_print(0, 1, "# Creating rational factor base of %zuMb took %1.1fs\n",
                     fb_size (sis->fb) >> 20, tfb);
        }
    }
    if (fbcfilename != NULL) {
        verbose_output_print(0, 1, "# Writing memory image of factor base to file %s\n", fbcfilename);
        fb_dump_fbc(si->sides[0]->fb,
                    (const factorbase_degn_t **) (si->sides[0]->fb_bucket_threads),
                    si->sides[1]->fb,
                    (const factorbase_degn_t **) (si->sides[1]->fb_bucket_threads),
                    fbcfilename, las->nb_threads);
        verbose_output_print(0, 1, "# Finished writing memory image of factor base\n");
    }
}
/*}}}*/

/* {{{ reordering of the small factor base
 *
 * We split the small factor base in several non-overlapping, contiguous
 * zones:
 *
 *      - powers of 2 (up until the pattern sieve limit)
 *      - powers of 3 (up until the pattern sieve limit)
 *      - trialdiv primes (not powers)
 *      - resieved primes
 *      (- powers of trialdiv primes)
 *      - rest.
 *
 * Problem: bad primes may in fact be pattern sieved, and we might want
 * to pattern-sieve more than just the ``do it always'' cases where p is
 * below the pattern sieve limit.
 *
 * The answer to this is that such primes are expected to be very very
 * rare, so we don't really bother. If we were to do something, we could
 * imagine setting up a schedule list for projective primes -- e.g. a
 * priority queue. But it feels way overkill.
 *
 * Note that the pre-treatment (splitting the factor base in chunks) can
 * be done once and for all.
 */

void reorder_fb(sieve_info_ptr si, int side)
{
    factorbase_degn_t * fb_pow2, * fb_pow2_base;
    factorbase_degn_t * fb_pow3, * fb_pow3_base;
    factorbase_degn_t * fb_td, * fb_td_base;
    // factorbase_degn_t * fb_pow_td, * fb_pow_td_base;
    factorbase_degn_t * fb_rs, * fb_rs_base;
    factorbase_degn_t * fb_rest, * fb_rest_base;

    factorbase_degn_t * fb_base = si->sides[side]->fb;
    factorbase_degn_t * fb = fb_base;

    size_t sz = fb_size(fb);

    fb_pow2 = fb_pow2_base = (factorbase_degn_t *) malloc(sz);
    fb_pow3 = fb_pow3_base = (factorbase_degn_t *) malloc(sz);
    fb_td = fb_td_base = (factorbase_degn_t *) malloc(sz);
    // fb_pow_td = fb_pow_td_base = (factorbase_degn_t *) malloc(sz);
    fb_rs = fb_rs_base = (factorbase_degn_t *) malloc(sz);
    fb_rest = fb_rest_base = (factorbase_degn_t *) malloc(sz);

    fbprime_t plim = si->conf->bucket_thresh;
    fbprime_t costlim = si->conf->td_thresh;

#define PUSH_LIST(x) do {						\
            memcpy(fb_## x, fb, fb_entrysize(fb));			\
            fb_## x = fb_next(fb_## x);					\
} while (0)

    size_t pattern2_size = sizeof(unsigned long) * 2;
    for( ; fb->p != FB_END ; fb = fb_next(fb)) {
        /* The extra conditions on powers of 2 and 3 are related to how
         * pattern-sieving is done.
         */
        if ((fb->p%2)==0 && fb->p <= pattern2_size) {
            PUSH_LIST(pow2);
        } else if (fb->p == 3) {
            PUSH_LIST(pow3);
        } else if (fb->p <= plim && fb->p <= costlim * fb->nr_roots) {
            if (!is_prime_power(fb->p)) {
                PUSH_LIST(td);
            } else {
                // PUSH_LIST(pow_td);
                PUSH_LIST(rest);
            }
        } else {
            if (!is_prime_power(fb->p)) {
                PUSH_LIST(rs);
            } else {
                PUSH_LIST(rest);
            }
        }
    }
#undef PUSH_LIST

#define APPEND_LIST(x) do {						\
    char * pb = (char*) (void*) fb_ ## x ## _base;			\
    char * p  = (char*) (void*) fb_ ## x;				\
    si->sides[side]->fb_parts->x[0] = fb;                               \
    si->sides[side]->fb_parts_x->x[0] = n;                              \
    memcpy(fb, pb, p - pb);						\
    fb = fb_skip(fb, p - pb);						\
    n += fb_diff(fb_ ## x, fb_ ## x ## _base);                          \
    si->sides[side]->fb_parts->x[1] = fb;                               \
    si->sides[side]->fb_parts_x->x[1] = n;                              \
} while (0)
    unsigned int n = 0;
    fb = fb_base;

    APPEND_LIST(pow2);
    APPEND_LIST(pow3);
    APPEND_LIST(td);
    APPEND_LIST(rs);
    APPEND_LIST(rest);
    fb->p = FB_END;

    free(fb_pow2_base);
    free(fb_pow3_base);
    free(fb_td_base);
    free(fb_rs_base);
    free(fb_rest_base);

#undef  APPEND_LIST

}

/* }}} */

/* {{{ Print some statistics about the factor bases
 * This also fills the field si->sides[*]->max_bucket_fill_ratio, which
 * is used to verify that per-thread allocation for buckets is
 * sufficient.
 * */
void sieve_info_print_fb_statistics(las_info_ptr las, sieve_info_ptr si, int side)
{
    const int n = las->nb_threads;
    sieve_side_info_ptr s = si->sides[side];
    double bucket_fill_ratio[n];

    verbose_output_print(0, 1, "# Number of small-sieved primes in %s factor base = %zu\n", sidenames[side], fb_nroots_total(s->fb));

    /* Counting the bucket-sieved primes per thread.  */
    unsigned long * nn = (unsigned long *) malloc(n * sizeof(unsigned long));
    ASSERT_ALWAYS(nn);
    memset(nn, 0, n * sizeof(unsigned long));
    for (int i = 0; i < n; ++i) {
        bucket_fill_ratio[i] = 0;
        factorbase_degn_t *fb = s->fb_bucket_threads[i];
        for (; fb->p != FB_END; fb = fb_next(fb)) {
            nn[i] += fb->nr_roots;
            bucket_fill_ratio[i] += fb->nr_roots / (double) fb->p;
        }
    }
    verbose_output_print(0, 1, "# Number of bucket-sieved primes in %s factor base per thread =", sidenames[side]);
    for(int i = 0 ; i < n ; i++)
        verbose_output_print(0, 1, " %lu", nn[i]);
    verbose_output_print(0, 1, "\n");
    verbose_output_print(0, 1, "# Inverse sum of bucket-sieved primes in %s factor base per thread =", sidenames[side]);
    for(int i = 0 ; i < n ; i++)
        verbose_output_print(0, 1, " %.5f", bucket_fill_ratio[i]);

    double min_bucket_fill_ratio = bucket_fill_ratio[n-1];
    double max_bucket_fill_ratio = bucket_fill_ratio[0];
    for(int i = 0 ; i < n ; i++) {
        double r = bucket_fill_ratio[i];
        if (r < min_bucket_fill_ratio) min_bucket_fill_ratio = r;
        if (r > max_bucket_fill_ratio) max_bucket_fill_ratio = r;
    }
    verbose_output_print(0, 1, " [hit jitter %g]\n",
            (max_bucket_fill_ratio / min_bucket_fill_ratio - 1));
    /* enable some margin in the bucket size */
    s->max_bucket_fill_ratio = max_bucket_fill_ratio * 1.07;
    free(nn);
}
/*}}}*/

/*}}}*/

/* {{{ sieve_info_init_... */
static void
sieve_info_init_from_siever_config(las_info_ptr las, sieve_info_ptr si, siever_config_srcptr sc, param_list pl)
{
    memset(si, 0, sizeof(sieve_info));
    si->cpoly = las->cpoly;
    memcpy(si->conf, sc, sizeof(siever_config));
    ASSERT_ALWAYS(sc->logI > 0);
    si->I = 1 << sc->logI;
    si->J = 1 << (sc->logI - 1);

    /* Initialize the number of buckets */

    /* If LOG_BUCKET_REGION == (sc->logI-1), then one bucket (whose size is the
     * L1 cache size) is actually one line. This changes some assumptions
     * in sieve_small_bucket_region and resieve_small_bucket_region, where
     * we want to differentiate on the parity on j.
     */
    ASSERT_ALWAYS(LOG_BUCKET_REGION >= (sc->logI - 1));

#ifndef SUPPORT_I17
    if (sc->logI >= 17) {
        fprintf(stderr,
                "Error: -I 17 requires setting the SUPPORT_I17 flag at compile time\n");
        abort();
    }
#endif

    /* this is the maximal value of the number of buckets (might be less
       for a given special-q if J is smaller) */
    si->nb_buckets_max = 1 + ((si->J << si->conf->logI) - 1) / BUCKET_REGION;
    verbose_output_print(0, 1, "# bucket_region = %u\n", BUCKET_REGION);
    if (si->nb_buckets_max < THRESHOLD_K_BUCKETS)
      verbose_output_print(0, 1, "# nb_buckets_max = %u, one pass for the buckets sort\n", si->nb_buckets_max);
    else if (si->nb_buckets_max < THRESHOLD_M_BUCKETS)
      verbose_output_print(0, 1, "# nb_buckets_max = %u, two passes for the buckets sort\n", si->nb_buckets_max);
    else
      verbose_output_print(0, 1, "# nb_buckets_max = %u, three passes for the buckets sort\n", si->nb_buckets_max);

    si->j_div = init_j_div(si->J);
    si->us = init_unsieve_data(si->I);
    mpz_init(si->doing->p);
    mpz_init(si->doing->r);

    /* Allocate memory for transformed polynomials */
    sieve_info_init_norm_data(si);

    /* This function in itself is too expensive if called often.
     * In the descent, where we initialize a sieve_info for each pair
     * (bitsize_of_q, side), we would call it dozens of time.
     * For the moment, we will assume that all along the hint file, the
     * values of I, alim, rlim are the constant, so that the factor bases
     * are exactly the same and can be shared.
     */
    if (las->sievers == si) {
        sieve_info_init_factor_bases(las, si, pl);
    } else {
        // We are in descent mode, it seems, so let's not duplicate the
        // factor base data.
        // A few sanity checks, first.
        ASSERT_ALWAYS(las->sievers->conf->logI == si->conf->logI);
        ASSERT_ALWAYS(las->sievers->conf->bucket_thresh == si->conf->bucket_thresh);
        ASSERT_ALWAYS(las->sievers->conf->sides[0]->lim == si->conf->sides[0]->lim);
        ASSERT_ALWAYS(las->sievers->conf->sides[0]->powlim == si->conf->sides[0]->powlim);
        ASSERT_ALWAYS(las->sievers->conf->sides[1]->lim == si->conf->sides[1]->lim);
        ASSERT_ALWAYS(las->sievers->conf->sides[1]->powlim == si->conf->sides[1]->powlim);
        // Then, copy relevant data from the first sieve_info
        verbose_output_print(0, 1, "# Do not regenerate factor base data: copy it from first siever\n");
        for (int side = 0; side < 2; side++) {
            sieve_side_info_ptr sis = si->sides[side];
            sieve_side_info_ptr sis0 = las->sievers->sides[side];
            sis->fb = sis0->fb;
            sis->fb_bucket_threads = sis0->fb_bucket_threads;
            sis->fb_is_mmapped = sis0->fb_is_mmapped;
        }
    }


    /* TODO: We may also build a strategy book, given that several
     * strategies will be similar. Presently we spend some time creating
     * each of them for the descent case.
     */

    for(int s = 0 ; s < 2 ; s++) {
        sieve_info_print_fb_statistics(las, si, s);
        /* init_norms (si, s); */ /* only depends on scale, logmax, lognorm_table */
        sieve_info_init_trialdiv(si, s); /* Init refactoring stuff */

        mpz_init (si->BB[s]);
        mpz_init (si->BBB[s]);
        mpz_init (si->BBBB[s]);
        unsigned long lim = si->conf->sides[s]->lim;
        mpz_ui_pow_ui (si->BB[s], lim, 2);
        mpz_mul_ui (si->BBB[s], si->BB[s], lim);
        mpz_mul_ui (si->BBBB[s], si->BBB[s], lim);

        /* The strategies also depend on the special-q used within the
         * descent, assuming lim / lpb depend on the sq bitsize */
        verbose_output_print(0, 1, "# Creating strategy for %d%c/%s [lim=%lu lpb=%u]\n",
                sc->bitsize, sidenames[sc->side][0], sidenames[s],
                sc->sides[s]->lim, sc->sides[s]->lpb);
        verbose_output_print(0, 1, "# Using %d+3 P-1/P+1/ECM curves\n",
                nb_curves (sc->sides[s]->lpb));
        si->sides[s]->strategy = facul_make_strategy(
                sc->sides[s]->lim, sc->sides[s]->lpb, 0);
        reorder_fb(si, s);

        verbose_output_print(0, 2, "# small %s factor base", sidenames[s]);
        factorbase_degn_t ** q;
        q = si->sides[s]->fb_parts->pow2;
        verbose_output_print(0, 2, ": %d pow2", fb_diff(q[1], q[0]));
        q = si->sides[s]->fb_parts->pow3;
        verbose_output_print(0, 2, ", %d pow3", fb_diff(q[1], q[0]));
        q = si->sides[s]->fb_parts->td;
        verbose_output_print(0, 2, ", %d td", fb_diff(q[1], q[0]));
        q = si->sides[s]->fb_parts->rs;
        verbose_output_print(0, 2, ", %d rs", fb_diff(q[1], q[0]));
        q = si->sides[s]->fb_parts->rest;
        verbose_output_print(0, 2, ", %d rest", fb_diff(q[1], q[0]));
        verbose_output_print(0, 2, " (total %zu)\n", fb_nroots_total(si->sides[s]->fb));
    }
}
/* }}} */

/* This is an ownership transfer from the current head of the todo list
 * into si->doing.  The field todo->next is not accessed, since it does
 * not make sense for si->doing, which is only a single item. The head of
 * the todo list is pruned */
void sieve_info_pick_todo_item(sieve_info_ptr si, las_todo_ptr * todo)
{
    ASSERT_ALWAYS(mpz_poly_is_root(si->cpoly->pols[(*todo)->side], (*todo)->r, (*todo)->p));
    mpz_clear(si->doing->p);
    mpz_clear(si->doing->r);
    memcpy(si->doing, *todo, sizeof(las_todo));
    free(*todo);
    *todo = si->doing->next;
    si->doing->next = 0;
    /* sanity check */
    if (!mpz_probab_prime_p(si->doing->p, 1)) {
        verbose_output_vfprint(1, 0, gmp_vfprintf, "Error, %Zd is not prime\n",
                               si->doing->p);
        exit(1);
    }
    ASSERT_ALWAYS(si->conf->side == si->doing->side);
}

static void sieve_info_update (FILE *output, sieve_info_ptr si, int nb_threads)/*{{{*/
{
  /* essentially update the fij polynomials and J value */
  sieve_info_update_norm_data(output, si, nb_threads);

  /* update number of buckets */
  si->nb_buckets = 1 + ((si->J << si->conf->logI) - 1) / BUCKET_REGION;

  /* Update the steps of the log of factor base primes */
  for(int side = 0 ; side < 2 ; side++) {
      sieve_side_info_ptr sis = si->sides[side];
      unsigned long lim = si->conf->sides[side]->lim;
      sis->log_steps_max = fb_make_steps(sis->log_steps, lim, sis->scale * LOG_SCALE);
  }

}/*}}}*/

static void sieve_info_clear (las_info_ptr las, sieve_info_ptr si)/*{{{*/
{
    clear_unsieve_data(si->us);
    si->us = NULL;
    clear_j_div(si->j_div);
    si->j_div = NULL;

    for(int s = 0 ; s < 2 ; s++) {
        facul_clear_strategy (si->sides[s]->strategy);
        si->sides[s]->strategy = NULL;
        sieve_info_clear_trialdiv(si, s);
        if (si == las->sievers) {
            free(si->sides[s]->fb);
            if (! si->sides[s]->fb_is_mmapped) {
                for(int i = 0 ; i < las->nb_threads ; i++) {
                    free(si->sides[s]->fb_bucket_threads[i]);
                }
            }
            free(si->sides[s]->fb_bucket_threads);
        }

        mpz_clear (si->BB[s]);
        mpz_clear (si->BBB[s]);
        mpz_clear (si->BBBB[s]);
    }
    sieve_info_clear_norm_data(si);
    mpz_clear(si->doing->p);
    mpz_clear(si->doing->r);
}/*}}}*/

/* las_info stuff */

static void las_info_init_hint_table(las_info_ptr las, param_list pl)/*{{{*/
{
    const char * filename = param_list_lookup_string(pl, "descent-hint");
    if (filename == NULL) return;
    char line[1024];
    unsigned int hint_alloc = 0;
    unsigned int hint_size = 0;
    FILE * f;
    f = fopen(filename, "r");
    if (f == NULL) {
        fprintf(stderr, "%s: %s\n", filename, strerror(errno));
        /* There's no point in proceeding, since it would really change
         * the behaviour of the program to do so */
        exit(1);
    }
    las->max_hint_bitsize[0] = 0;
    las->max_hint_bitsize[1] = 0;
    for(;;) {
        char * x = fgets(line, sizeof(line), f);
        double t;
        unsigned long z;
        /* Tolerate comments and blank lines */
        if (x == NULL) break;
        if (*x == '#') continue;
        for( ; *x && isspace(*x) ; x++) ;
        if (!*x) continue;

        /* We have a new entry to parse */
        if (hint_size >= hint_alloc) {
            hint_alloc = 2 * hint_alloc + 8;
            las->hint_table = realloc(las->hint_table, hint_alloc * sizeof(descent_hint));
        }

        descent_hint_ptr h = las->hint_table[hint_size++];
        siever_config_ptr sc = h->conf;

        z = strtoul(x, &x, 10);
        ASSERT_ALWAYS(z > 0);
        sc->bitsize = z;
        switch(*x++) {
            case 'a' : sc->side = ALGEBRAIC_SIDE; break;
            case 'r' : sc->side = RATIONAL_SIDE; break;
            default:
                       fprintf(stderr, "%s: parse error at %s\n", filename, line);
                       exit(1);
        }
        for( ; *x && isspace(*x) ; x++) ;
        t = strtod(x, &x); ASSERT_ALWAYS(t > 0);
        h->expected_time = t;
        for( ; *x && isspace(*x) ; x++) ;
        t = strtod(x, &x); ASSERT_ALWAYS(t > 0);
        h->expected_success = t;
        for( ; *x && isspace(*x) ; x++) ;
        for( ; *x && !isdigit(*x) ; x++) ;
        z = strtoul(x, &x, 10); ASSERT_ALWAYS(z > 0);
        sc->logI = z;
        
        for(int s = 0 ; s < 2 ; s++) {
            for( ; *x && isspace(*x) ; x++) ;
            z = strtoul(x, &x, 10); ASSERT_ALWAYS(z > 0);
            sc->sides[s]->lim = z;
            for( ; *x && !isdigit(*x) ; x++) ;
            z = strtoul(x, &x, 10); ASSERT_ALWAYS(z > 0);
            sc->sides[s]->lpb = z;
            for( ; *x && !isdigit(*x) ; x++) ;
            z = strtoul(x, &x, 10); ASSERT_ALWAYS(z > 0);
            sc->sides[s]->mfb = z;
            for( ; *x && !isdigit(*x) ; x++) ;
            t = strtod(x, &x); ASSERT_ALWAYS(t > 0);
            sc->sides[s]->lambda = t;
        }
        for( ; *x ; x++) ASSERT_ALWAYS(isspace(*x));

        if (sc->bitsize > las->max_hint_bitsize[sc->side])
            las->max_hint_bitsize[sc->side] = sc->bitsize;
        // Copy default value for non-given parameters
        sc->skewness = las->default_config->skewness;
        sc->bucket_thresh = las->default_config->bucket_thresh;
        sc->td_thresh = las->default_config->td_thresh;
        sc->unsieve_thresh = las->default_config->unsieve_thresh;
        sc->sides[0]->powlim = las->default_config->sides[0]->powlim;
        sc->sides[1]->powlim = las->default_config->sides[1]->powlim;
    }
    if (las->hint_table == NULL) {
        fprintf(stderr, "%s: no data ??\n", filename);
        exit(1);
    }
    if (hint_size >= hint_alloc) {
        hint_alloc = 2 * hint_alloc + 8;
        las->hint_table = realloc(las->hint_table, hint_alloc * sizeof(descent_hint));
    }
    las->hint_table[hint_size++]->conf->bitsize = 0;
    las->hint_table = realloc(las->hint_table, hint_size * sizeof(descent_hint));

    /* Allocate the quick lookup tables */

    for(int s = 0 ; s < 2 ; s++) {
        unsigned int n = las->max_hint_bitsize[s] + 1;
        las->hint_lookups[s] = malloc(n * sizeof(unsigned int));
        for(unsigned int i = 0 ; i < n ; i++) {
            las->hint_lookups[s][i] = -1;
        }
    }
    for(descent_hint_ptr h = las->hint_table[0] ; h->conf->bitsize ; h++) {
        siever_config_ptr sc = h->conf;
        int d = las->hint_lookups[sc->side][sc->bitsize];
        if (d >= 0) {
            fprintf(stderr, "Error: two hints found for %d%c\n",
                    sc->bitsize, sidenames[sc->side][0]);
            exit(1);
        }
        las->hint_lookups[sc->side][sc->bitsize] = h - las->hint_table[0];

        /* We must make sure that the default factor base bounds are
         * larger than what we have for all the hint cases, so that it
         * remains reasonable to base our work on the larger factor base
         * (thus doing incomplete sieving).
         */
        /* But now, the .poly file does not contain lim data anymore.
         * This is up to the user to do this check, anyway, because
         * it is not sure that makefb was run with the alim given in the
         * poly file.
         */
        /*
        for(int s = 0 ; s < 2 ; s++) {
            ASSERT_ALWAYS(sc->sides[s]->lim <= las->cpoly->pols[s]->lim);
        }
        */
    }

    fclose(f);
}/*}}}*/

static void las_info_clear_hint_table(las_info_ptr las)/*{{{*/
{
    if (las->hint_table == NULL) return;
    for(int s = 0 ; s < 2 ; s++) {
        free(las->hint_lookups[s]);
        las->hint_lookups[s] = NULL;
        las->max_hint_bitsize[s] = 0;
    }
    free(las->hint_table);
    las->hint_table = NULL;
}/*}}}*/

static void las_info_init(las_info_ptr las, param_list pl)/*{{{*/
{
    memset(las, 0, sizeof(las_info));
    /* {{{ Parse and install general operational flags */
    las->outputname = param_list_lookup_string(pl, "out");
    /* Init output file */
    las->output = stdout;
    if (las->outputname) {
	if (!(las->output = fopen_maybe_compressed(las->outputname, "w"))) {
	    fprintf(stderr, "Could not open %s for writing\n", las->outputname);
	    exit(EXIT_FAILURE);
	}
    }

    las->verbose = param_list_parse_switch(pl, "-v");

    verbose_output_init(3);
    verbose_output_add(0, las->output, las->verbose + 1);
    verbose_output_add(1, stderr, 1);

    /* Channel 2 is for statistics. We always print them to las' normal output */
    verbose_output_add(2, las->output, 1);
    if (param_list_parse_switch(pl, "-stats-stderr")) {
        /* If we should also print stats to stderr, add stderr to channel 2 */
        verbose_output_add(2, stderr, 1);
    }

    verbose_set_enabled_flags(pl);
    param_list_print_command_line(las->output, pl);
    las_display_config_flags();

    las->suppress_duplicates = param_list_parse_switch(pl, "-dup");
    las->nb_threads = 1;		/* default value */
    param_list_parse_int(pl, "t", &las->nb_threads);
    if (las->nb_threads <= 0) {
	fprintf(stderr,
		"Error, please provide a positive number of threads\n");
	param_list_clear(pl);
	exit(EXIT_FAILURE);
    }
    /* }}} */
    /* {{{ Parse polynomial */
    cado_poly_init(las->cpoly);
    const char *tmp;
    if ((tmp = param_list_lookup_string(pl, "poly")) == NULL) {
        fprintf(stderr, "Error: -poly is missing\n");
        param_list_print_usage(pl, NULL, stderr);
	cado_poly_clear(las->cpoly);
	param_list_clear(pl);
        exit(EXIT_FAILURE);
    }

    if (!cado_poly_read(las->cpoly, tmp)) {
	fprintf(stderr, "Error reading polynomial file %s\n", tmp);
	cado_poly_clear(las->cpoly);
	param_list_clear(pl);
	exit(EXIT_FAILURE);
    }

    /* {{{ Parse default siever config (fill all possible fields) */

    siever_config_ptr sc = las->default_config;
    memset(sc, 0, sizeof(siever_config));

    sc->skewness = las->cpoly->skew;
    /* -skew (or -S) may override (or set) the skewness given in the
     * polynomial file */
    param_list_parse_double(pl, "skew", &(sc->skewness));

    if (sc->skewness <= 0.0) {
	fprintf(stderr, "Error, please provide a positive skewness\n");
	cado_poly_clear(las->cpoly);
	param_list_clear(pl);
	exit(EXIT_FAILURE);
    }
    /* }}} */

    /* The default config is not necessarily a complete bit of
     * information.
     *
     * The field with the bitsize of q is here filled with q0 as found in
     * the command line. Note though that this is only one choice among
     * several possible (q1, or just no default).
     * For the descent, we do not intend to have a default config, thus
     * specifying q0 makes no sense. Likewise, for file-based todo lists,
     * we have no default either, and no siever is configured to provide
     * this ``default'' behaviour (yet, the default bits here are used to
     * pre-fill the config data later on).
     */
    mpz_t q0;
    mpz_init(q0);
    if (param_list_parse_mpz(pl, "q0", q0)) {
        sc->bitsize = mpz_sizeinbase(q0, 2);
    }
    mpz_clear(q0);
    sc->side = param_list_parse_switch(pl, "-ratq") ? RATIONAL_SIDE : ALGEBRAIC_SIDE;
    param_list_parse_double(pl, "lambda0", &(sc->sides[RATIONAL_SIDE]->lambda));
    param_list_parse_double(pl, "lambda1", &(sc->sides[ALGEBRAIC_SIDE]->lambda));
    int seen = 1;
    seen  = param_list_parse_int   (pl, "I",       &(sc->logI));
    seen &= param_list_parse_ulong (pl, "lim0",    &(sc->sides[RATIONAL_SIDE]->lim));
    seen &= param_list_parse_int   (pl, "lpb0",    &(sc->sides[RATIONAL_SIDE]->lpb));
    seen &= param_list_parse_int   (pl, "mfb0",    &(sc->sides[RATIONAL_SIDE]->mfb));
    seen &= param_list_parse_ulong (pl, "lim1",    &(sc->sides[ALGEBRAIC_SIDE]->lim));
    seen &= param_list_parse_int   (pl, "lpb1",    &(sc->sides[ALGEBRAIC_SIDE]->lpb));
    seen &= param_list_parse_int   (pl, "mfb1",    &(sc->sides[ALGEBRAIC_SIDE]->mfb));
    if (!seen) {
        fprintf(stderr, "Error: options -I, -lim0, -lpb0, -mfb0, "
                " -lim1, -lpb1, -mfb1 are mandatory.\n");
	cado_poly_clear(las->cpoly);
	param_list_clear(pl);
        exit(EXIT_FAILURE);
    }

    /* compute lambda0 from mfb0 and lpb0 if not given */
    if (sc->sides[RATIONAL_SIDE]->lambda == 0.0)
      sc->sides[RATIONAL_SIDE]->lambda = (double) sc->sides[RATIONAL_SIDE]->mfb
        / (double) sc->sides[RATIONAL_SIDE]->lpb;

    /* compute lambda1 from mfb1 and lpb1 if not given */
    if (sc->sides[ALGEBRAIC_SIDE]->lambda == 0.0)
      sc->sides[ALGEBRAIC_SIDE]->lambda = (double) sc->sides[ALGEBRAIC_SIDE]->mfb
        / (double) sc->sides[ALGEBRAIC_SIDE]->lpb;

    /* Parse optional siever configuration parameters */
    sc->td_thresh = 1024;	/* default value */
    param_list_parse_uint(pl, "tdthresh", &(sc->td_thresh));

    sc->unsieve_thresh = 100;
    if (param_list_parse_uint(pl, "unsievethresh", &(sc->unsieve_thresh))) {
        verbose_output_print(0, 1, "# Un-sieving primes > %u\n",
                sc->unsieve_thresh);
    }

    sc->bucket_thresh = 1 << sc->logI;	/* default value */
    /* overrides default only if parameter is given */
    param_list_parse_ulong(pl, "bkthresh", &(sc->bucket_thresh));

    const char *powlim_params[2] = {"powlim0", "powlim1"};
    for (int side = 0; side < 2; side++) {
        if (!param_list_parse_ulong(pl, powlim_params[side], &sc->sides[side]->powlim)) {
            sc->sides[side]->powlim = sc->bucket_thresh - 1;
            verbose_output_print(0, 1, "# Using default value of %lu for -%s\n",
                     sc->sides[side]->powlim, powlim_params[side]);
        } else if (sc->sides[side]->powlim >= sc->bucket_thresh) {
            sc->sides[side]->powlim = sc->bucket_thresh - 1;
            verbose_output_print(0, 1, "# -%s reduced to %lu\n",
                    powlim_params[side], sc->sides[side]->powlim);
        }
    }

    /* }}} */

    /* {{{ Init and parse info regarding work to be done by the siever */
    /* Actual parsing of the command-line fragments is done within
     * las_todo_feed, but this is an admittedly contrived way to work */
    mpz_init_set_ui(las->todo_q0, 0);
    mpz_init_set_ui(las->todo_q1, 0);
    const char * filename = param_list_lookup_string(pl, "todo");
    if (filename) {
        las->todo_list_fd = fopen(filename, "r");
        if (las->todo_list_fd == NULL) {
            fprintf(stderr, "%s: %s\n", filename, strerror(errno));
            /* There's no point in proceeding, since it would really change
             * the behaviour of the program to do so */
            cado_poly_clear(las->cpoly);
            param_list_clear(pl);
            exit(EXIT_FAILURE);
        }
    }
    /* }}} */
    las_info_init_hint_table(las, pl);
    /* Allocate room for only one sieve_info */
    las->sievers = malloc(sizeof(sieve_info));
    memset(las->sievers, 0, sizeof(sieve_info));
}/*}}}*/

void las_info_clear(las_info_ptr las)/*{{{*/
{
    las_info_clear_hint_table(las);

    for(sieve_info_ptr si = las->sievers ; si->conf->bitsize ; si++) {
        sieve_info_clear(las, si);
    }
    free(las->sievers);
    if (las->outputname) {
        fclose_maybe_compressed(las->output, las->outputname);
    }
    verbose_output_clear();
    mpz_clear(las->todo_q0);
    mpz_clear(las->todo_q1);
    if (las->todo_list_fd)
        fclose(las->todo_list_fd);
    cado_poly_clear(las->cpoly);
}/*}}}*/

/* Look for an existing sieve_info in las->sievers with configuration matching
   that in sc; if none exists, create one. */
sieve_info_ptr get_sieve_info_from_config(las_info_ptr las, siever_config_srcptr sc, param_list pl)/*{{{*/
{
    int n = 0;
    sieve_info_ptr si;
    for(si = las->sievers ; si->conf->bitsize ; si++, n++) {
        if (siever_config_cmp(si->conf, sc) == 0)
            break;
    }
    if (si->conf->bitsize)
        return si;
    /* We've hit the end marker. Need to add a new config. */
    las->sievers = realloc(las->sievers, (n+2) * sizeof(sieve_info));
    si = las->sievers + n;
    verbose_output_print(0, 1, "# Creating new sieve configuration for q~2^%d on the %s side\n",
            sc->bitsize, sidenames[sc->side]);
    sieve_info_init_from_siever_config(las, si, sc, pl);
    memset(si + 1, 0, sizeof(sieve_info));
    siever_config_display(las->output, sc);
    return si;
}/*}}}*/

void las_todo_push_withdepth(las_todo_ptr * d, mpz_srcptr p, mpz_srcptr r, int side, int depth)/*{{{*/
{
    las_todo_ptr nd = malloc(sizeof(las_todo));
    memset(nd, 0, sizeof(las_todo));
    mpz_init_set(nd->p, p);
    mpz_init_set(nd->r, r);
    nd->side = side;
    nd->depth = depth;
    nd->next = *d;
    *d = nd;
}
void las_todo_push(las_todo_ptr * d, mpz_srcptr p, mpz_srcptr r, int side)
{
    las_todo_push_withdepth(d, p, r, side, 0);
}
void las_todo_push_closing_brace(las_todo_ptr * d, int depth)
{
    las_todo_ptr nd = malloc(sizeof(las_todo));
    memset(nd, 0, sizeof(las_todo));
    mpz_init(nd->p);
    mpz_init(nd->r);
    nd->side = -1;
    nd->depth = depth;
    nd->next = *d;
    *d = nd;
}

/*}}}*/

void las_todo_pop(las_todo_ptr * d)/*{{{*/
{
    las_todo_ptr nd = *d;
    *d = nd->next;
    mpz_clear(nd->p);
    mpz_clear(nd->r);
    free(nd);
}

int las_todo_pop_closing_brace(las_todo_ptr * d)
{
    if ((*d)->side >= 0)
        return 0;
    las_todo_pop(d);
    return 1;
}
/*}}}*/

/* {{{ Populating the todo list */
/* See below in main() for documentation about the q-range and q-list
 * modes */
/* These functions return non-zero if the todo list is not empty */
int las_todo_feed_qrange(las_info_ptr las, param_list pl)
{
    if (las->todo) return 1; /* keep going */
    /* handy aliases */
    mpz_ptr q = las->todo_q0;
    mpz_ptr q1 = las->todo_q1;

    int pushed = 0;

    int qside = las->default_config->side;

    if (mpz_cmp_ui(q, 0) == 0) {
        ASSERT_ALWAYS(param_list_parse_mpz(pl, "q0", q));
        if (!param_list_parse_mpz(pl, "q1", q1)) {
            mpz_t r;
            mpz_init(r);
            if (!param_list_parse_mpz(pl, "rho", r)) {
                if (qside != RATIONAL_SIDE) {
                    fprintf(stderr, "Error: single special-q requires -rho on algebraic side\n");
                    exit(EXIT_FAILURE);
                }
                if (!mpz_probab_prime_p(q, 2)) {
                    mpz_t q2;
                    mpz_init(q2);
                    mpz_nextprime(q2, q);
                    verbose_output_vfprint(1, 0, gmp_vfprintf, "Warning: fixing q=%Zd to next prime q=%Zd\n", q, q2);
                    mpz_set(q, q2);
                    mpz_clear(q2);
                }
                /* Arrange so that the normal code handles this q and
                 * computes rho. */
                mpz_add_ui(q1, q, 1);
                mpz_sub_ui(q, q, 1);
            } else {
                ASSERT_ALWAYS(mpz_probab_prime_p(q, 2));
                las_todo_push(& las->todo, q, r, qside);
                pushed++;
                /* Anyway we're not subtracting anything from q, so that
                 * the tail code won't do anything */
                mpz_add_ui(q1, q, 1);
            }
            mpz_clear(r);
        } else {
            /* Nextprime should return q0 if q0 is prime */
            mpz_sub_ui(q, q, 1);
        }
    }
    /* Otherwise we're going to process the next few sq's and put them
     * into the list */
    mpz_t * roots;
    int deg = MAX(las->cpoly->rat->deg, las->cpoly->alg->deg);
    roots = (mpz_t *) malloc (deg * sizeof (mpz_t));
    for(int i = 0 ; i < deg  ; i++) {
        mpz_init(roots[i]);
    }

    las_todo_ptr * pnext = &(las->todo);
    for( ; pushed < 10 && mpz_cmp(q, q1) < 0 ; ) {
        mpz_nextprime(q, q);
        if (mpz_cmp(q, q1) >= 0)
            break;
        mpz_poly_ptr f = las->cpoly->pols[qside];
        int nroots = mpz_poly_roots (roots, f, q);
        if (param_list_parse_switch(pl, "-galois")) {
            if (nroots % 2) {
                fprintf(stderr, "Number of roots modulo q is odd. Don't know how to interpret -galois.\n");
                ASSERT_ALWAYS(0);
            }
            // Keep only one root among {r, 1/r} orbits.
            modulusul_t mm;
            unsigned long qq = mpz_get_ui(q);
            modul_initmod_ul(mm, qq);
            residueul_t r1, r2;
            modul_init(r1, mm);
            modul_init(r2, mm);
            for (int k = 0; k < nroots; k++) {
                unsigned long rr = mpz_get_ui(roots[k]);
                modul_set_ul(r1, rr, mm);
                int kk = 0;
                for (int l = k+1; l < nroots; ++l) {
                    unsigned long ss = mpz_get_ui(roots[l]);
                    modul_set_ul(r2, ss, mm);
                    modul_mul(r2, r2, r1, mm);
                    if (modul_is1(r2, mm)) {
                        kk = l;
                        break;
                    }
                }
                ASSERT_ALWAYS(kk != 0); // Should always find an inverse.
                // Remove it from the list
                for (int l = kk; l < nroots-1; ++l) {
                    mpz_set(roots[l], roots[l+1]);
                }
                nroots--;
            }
            modul_clear(r1, mm);
            modul_clear(r2, mm);
            modul_clearmod(mm);
        }
        /* {{{ This print is now dead (or we're going to print it at
         * weird times)
         *
                if (nroots > 0) {
                    gmp_fprintf (las->output, "### q=%Zd: root%s", q0,
                            (nroots == 1) ? "" : "s");
                    for (i = 1; i <= (int) nroots; i++)
                        gmp_fprintf (las->output, " %Zd", roots[nroots-i]);
                    fprintf (las->output, "\n");
                }
         * }}} */

        for(int i = 0 ; i < nroots ; i++) {
            las_todo_push(pnext, q, roots[i], qside);
            pnext = &((*pnext)->next);
            pushed++;
        }
    }

    for(int i = 0 ; i < deg  ; i++) {
        mpz_clear(roots[i]);
    }
    free(roots);
    return pushed;
}

/* Format of a file with a list of special-q (-todo option):
 *   - Comments are allowed (start line with #)
 *   - Blank lines are ignored
 *   - Each valid line must have the form
 *       s q r
 *     where s is "a" or "r", giving the "algebraic" or "rational" side
 *     of the special q, and q and r are as usual.
 */
int las_todo_feed_qlist(las_info_ptr las, param_list pl)
{
    if (las->todo) return 1; /* keep going */
    char line[1024];
    FILE * f = las->todo_list_fd;
    /* The fgets call below is blocking, so flush las->output here just to
     * be sure. */
    fflush(las->output);
    char * x;
    for( ; ; ) {
        x = fgets(line, sizeof(line), f);
        /* Tolerate comments and blank lines */
        if (x == NULL) return 0;
        if (*x == '#') continue;
        for( ; *x && isspace(*x) ; x++) ;
        if (!*x) continue;
        break;
    }

    /* We have a new entry to parse */
    mpz_t p, r;
    int side = -1;
    mpz_init(p);
    mpz_init(r);
    int nread;
    int rc;
    switch(*x++) {
        case 'a' : side = ALGEBRAIC_SIDE;
                   for( ; *x && !isdigit(*x) ; x++) ;
                   rc = gmp_sscanf(x, "%Zi %Zi%n", p, r, &nread);
                   x+=nread;
                   ASSERT_ALWAYS(rc == 2);  /* %n does not count */
                   ASSERT_ALWAYS(mpz_probab_prime_p(p, 2));
                   break;
        case 'r' : side = RATIONAL_SIDE; 
                   mpz_set_ui(r, 0);
                   for( ; *x && !isdigit(*x) ; x++) ;
                   rc = gmp_sscanf(x, "%Zi %Zi%n", p, r, &nread);
                   x+=nread;
                   ASSERT_ALWAYS(rc >= 1);  /* %n does not count */
                   ASSERT_ALWAYS(mpz_probab_prime_p(p, 2));
                   if (rc == 2)
                       break;
                   // If the root is not specified, then we assume that
                   // the side is really rational, and then we compute
                   // the root.
                   mpz_poly_ptr f = las->cpoly->pols[RATIONAL_SIDE];
                   ASSERT_ALWAYS(f->deg == 1);
                   int nroots = mpz_poly_roots (&r, f, p);
                   ASSERT_ALWAYS(nroots == 1);
                   break;
        default:
                   /* We may as well default on the command-line switch
                    */
                   fprintf(stderr, "%s: parse error at %s\n",
                           param_list_lookup_string(pl, "todo"), line);
                   exit(1);
    }
    for( ; *x ; x++) ASSERT_ALWAYS(isspace(*x));
    las_todo_push(&(las->todo), p, r, side);
    mpz_clear(p);
    mpz_clear(r);
    return 1;
}


int las_todo_feed(las_info_ptr las, param_list pl)
{
    if (las->todo) return 1; /* keep going */
    if (las->todo_list_fd)
        return las_todo_feed_qlist(las, pl);
    else
        return las_todo_feed_qrange(las, pl);
}
/* }}} */

/* Only when tracing. This function gets called once per
 * special-q only. Here we compute the two norms corresponding to
 * the traced (a,b) pair, and start by dividing out the special-q
 * from the one where it should divide */
void init_trace_k(sieve_info_srcptr si, param_list pl)
{
    struct trace_ab_t ab;
    struct trace_ij_t ij;
    struct trace_Nx_t Nx;
    int have_trace_ab = 0, have_trace_ij = 0, have_trace_Nx = 0;

    const char *abstr = param_list_lookup_string(pl, "traceab");
    if (abstr != NULL) {
        if (sscanf(abstr, "%"SCNd64",%"SCNu64, &ab.a, &ab.b) == 2)
            have_trace_ab = 1;
        else {
            fprintf (stderr, "Invalid value for parameter: -traceab %s\n",
                     abstr);
            exit (EXIT_FAILURE);
        }
    }

    const char *ijstr = param_list_lookup_string(pl, "traceij");
    if (ijstr != NULL) {
        if (sscanf(ijstr, "%d,%u", &ij.i, &ij.j) == 2) {
            have_trace_ij = 1;
        } else {
            fprintf (stderr, "Invalid value for parameter: -traceij %s\n",
                     ijstr);
            exit (EXIT_FAILURE);
        }
    }

    const char *Nxstr = param_list_lookup_string(pl, "traceNx");
    if (Nxstr != NULL) {
        if (sscanf(Nxstr, "%u,%u", &Nx.N, &Nx.x) == 2)
            have_trace_Nx = 1;
        else {
            fprintf (stderr, "Invalid value for parameter: -traceNx %s\n",
                     Nxstr);
            exit (EXIT_FAILURE);
        }
    }

    trace_per_sq_init(si, have_trace_Nx ? &Nx : NULL,
        have_trace_ab ? &ab : NULL, 
        have_trace_ij ? &ij : NULL);
}

/* utility. Can go elsewhere */
void thread_do(thread_data * thrs, void * (*f) (thread_data_ptr), int n)/*{{{*/
{
    if (n == 1) {
        /* Then don't bother with pthread calls */
        (*f)(thrs[0]);
        return;
    }
    pthread_t * th = malloc(n * sizeof(pthread_t)); 
    ASSERT_ALWAYS(th);

#if 0
    /* As a debug measure, it's possible to activate this branch instead
     * of the latter. In effect, this causes las to run in a
     * non-multithreaded way, albeit strictly following the code path of
     * the multithreaded case.
     */
    for (int i = 0; i < n ; ++i) {
        (*f)(thrs[i]);
    }
#else
    for (int i = 0; i < n ; ++i) {
        int ret = pthread_create(&(th[i]), NULL, 
		(void * (*)(void *)) f,
                (void *)(thrs[i]));
        ASSERT_ALWAYS(ret == 0);
    }
    for (int i = 0; i < n ; ++i) {
        int ret = pthread_join(th[i], NULL);
        ASSERT_ALWAYS(ret == 0);
    }
#endif

    free(th);
}/*}}}*/

/* {{{ apply_buckets */
#ifndef TRACE_K
/* backtrace display can't work for static symbols (see backtrace_symbols) */
NOPROFILE_STATIC
#endif

void
apply_one_bucket (unsigned char *S, bucket_array_t BA, const int i,
        where_am_I_ptr w)
{
  unsigned int next_logp_j;
  unsigned char logp;
  bucket_update_t *next_logp_change, *next_align, **pnlc, *read_ptr;

  read_ptr = BA.bucket_read[i];
  pnlc = BA.logp_idx + i;
  next_logp_j = 0;
  WHERE_AM_I_UPDATE(w, p, 0);
  while (read_ptr < BA.bucket_write[i]) {
    logp = BA.logp_val[next_logp_j++];
    if (LIKELY(next_logp_j < (unsigned int) BA.nr_logp)) {
      pnlc = (bucket_update_t **)((size_t) pnlc + BA.size_b_align);
      next_logp_change = *pnlc;
    }
    else
      next_logp_change = BA.bucket_write[i];
    next_align = (bucket_update_t *) (((size_t) read_ptr + 0x3F) & ~((size_t) 0x3F));
    if (UNLIKELY(next_align > next_logp_change)) next_align = next_logp_change;
    while (read_ptr < next_align) {
      uint16_t x;
      x = (read_ptr++)->x;
      WHERE_AM_I_UPDATE(w, x, x);
      sieve_increase(S + x, logp, w);
    }
    while (read_ptr + 16 <= next_logp_change) {
      uint64_t x0, x1, x2, x3, x4, x5, x6, x7;
      uint16_t x;
#ifdef HAVE_SSE2
      _mm_prefetch(((void *) read_ptr)+256, _MM_HINT_NTA);
#endif
      x0 = ((uint64_t *) read_ptr)[0];
      x1 = ((uint64_t *) read_ptr)[1];
      x2 = ((uint64_t *) read_ptr)[2];
      x3 = ((uint64_t *) read_ptr)[3];
      x4 = ((uint64_t *) read_ptr)[4];
      x5 = ((uint64_t *) read_ptr)[5];
      x6 = ((uint64_t *) read_ptr)[6];
      x7 = ((uint64_t *) read_ptr)[7];
      read_ptr += 16;
      __asm__ __volatile__ (""); /* To be sure all x? are read together in one operation */
#ifdef CADO_LITTLE_ENDIAN
#define INSERT_2_VALUES(X) do {						\
	(X) >>= 16; x = (uint16_t) (X);					\
	WHERE_AM_I_UPDATE(w, x, x); sieve_increase(S + x, logp, w);	\
	(X) >>= 32;							\
	WHERE_AM_I_UPDATE(w, x, X); sieve_increase(S + (X), logp, w);	\
      } while (0);
#else
#define INSERT_2_VALUES(X) do {						\
	x = (uint16_t) (X);						\
	WHERE_AM_I_UPDATE(w, x, x); sieve_increase(S + x, logp, w);	\
	(X) >>= 32; x = (uint16_t) (X);					\
	WHERE_AM_I_UPDATE(w, x, x); sieve_increase(S + x, logp, w);	\
      } while (0);

#endif
      INSERT_2_VALUES(x0); INSERT_2_VALUES(x1); INSERT_2_VALUES(x2); INSERT_2_VALUES(x3);
      INSERT_2_VALUES(x4); INSERT_2_VALUES(x5); INSERT_2_VALUES(x6); INSERT_2_VALUES(x7);
    }
    while (read_ptr < next_logp_change) {
      uint16_t x;
      x = (read_ptr++)->x;
      WHERE_AM_I_UPDATE(w, x, x);
      sieve_increase(S + x, logp, w);
    }
  }
}

/* {{{ Trial division */
typedef struct {
    uint64_t *fac;
    int n;
} factor_list_t;

#define FL_MAX_SIZE 200

void factor_list_init(factor_list_t *fl) {
    fl->fac = (uint64_t *) malloc (FL_MAX_SIZE * sizeof(uint64_t));
    ASSERT_ALWAYS(fl->fac != NULL);
    fl->n = 0;
}

void factor_list_clear(factor_list_t *fl) {
    free(fl->fac);
}

static void 
factor_list_add(factor_list_t *fl, const uint64_t p)
{
  ASSERT_ALWAYS(fl->n < FL_MAX_SIZE);
  fl->fac[fl->n++] = p;
}

// print a comma-separated list of factors.
// returns the number of factor printed (in particular, a comma is needed
// after this output only if the return value is non zero)
int factor_list_fprint(FILE *f, factor_list_t fl) {
    int i;
    for (i = 0; i < fl.n; ++i) {
        if (i) fprintf(f, ",");
        fprintf(f, "%" PRIx64, fl.fac[i]);
    }
    return i;
}


static const int bucket_prime_stats = 0;
static long nr_bucket_primes = 0;
static long nr_div_tests = 0;
static long nr_composite_tests = 0;
static long nr_wrap_was_composite = 0;
/* The entries in BP must be sorted in order of increasing x */
static void
divide_primes_from_bucket (factor_list_t *fl, mpz_t norm, const unsigned int N MAYBE_UNUSED, const int x,
                           bucket_primes_t *BP, const unsigned long fbb)
{
  bucket_prime_t prime;
  while (!bucket_primes_is_end (BP)) {
      prime = get_next_bucket_prime (BP);
      if (prime.x > x)
        {
          rewind_primes_by_1 (BP);
          break;
        }
      if (prime.x == x) {
          if (bucket_prime_stats) nr_bucket_primes++;
          unsigned long p = prime.p;
          while (p <= fbb) {
              if (bucket_prime_stats) nr_div_tests++;
              if (LIKELY(mpz_divisible_ui_p (norm, p))) {
                  int isprime;
                  modulusul_t m; 
                  modul_initmod_ul (m, p);
                  if (bucket_prime_stats) nr_composite_tests++;
                  isprime = modul_isprime (m);
                  modul_clearmod (m);
                  if (LIKELY(isprime)) {
                      break;
                  } else {
                    if (bucket_prime_stats) nr_wrap_was_composite++;
                  }
              }

              /* It may have been a case of incorrectly reconstructing p
                 from bits 1...16, so let's try if a bigger prime works.

                 Warning: this strategy may fail, since we might find a
                 composite p+k1*BUCKET_P_WRAP dividing the norm, while we
                 really want a larger prime p+k2*BUCKET_P_WRAP. In that case,
                 if a prime dividing p+k1*BUCKET_P_WRAP also divides the norm,
                 it might lead to a bucket error (p = ... does not divide),
                 moreover the wanted prime p+k2*BUCKET_P_WRAP will not be found
                 and we might miss some relations. */
              p += BUCKET_P_WRAP;
          }
          if (UNLIKELY(p > fbb)) {
              verbose_output_print(1, 0,
                       "# Error, p = %lu does not divide at (N,x) = (%u,%d)\n",
                       (unsigned long) prime.p, N, x);
              abort();
          }
          do {
              factor_list_add(fl, p);
              mpz_divexact_ui (norm, norm, p);
          } while (mpz_divisible_ui_p (norm, p));
      }
  }
}


NOPROFILE_STATIC void
trial_div (factor_list_t *fl, mpz_t norm, const unsigned int N, int x,
           factorbase_degn_t *fb, bucket_primes_t *primes,
	   trialdiv_divisor_t *trialdiv_data, const unsigned long fbb,
           int64_t a MAYBE_UNUSED, uint64_t b MAYBE_UNUSED)
{
    const int trial_div_very_verbose = 0;
    int nr_factors;
    fl->n = 0; /* reset factor list */

    if (trial_div_very_verbose)
        verbose_output_vfprint(1, 0, gmp_vfprintf, "# trial_div() entry, x = %d, norm = %Zd\n", x, norm);

    // handle 2 separately, if it is in fb
    if (fb->p == 2) {
        int bit = mpz_scan1(norm, 0);
        int i;
        for (i = 0; i < bit; ++i) {
            fl->fac[fl->n] = 2;
            fl->n++;
        }
        if (trial_div_very_verbose)
            verbose_output_vfprint(1, 0, gmp_vfprintf, "# x = %d, dividing out 2^%d, norm = %Zd\n", x, bit, norm);
        mpz_tdiv_q_2exp(norm, norm, bit);
        fb = fb_next (fb); // cannot do fb++, due to variable size !
    }

    // remove primes in "primes" that map to x
    divide_primes_from_bucket (fl, norm, N, x, primes, fbb);
#ifdef TRACE_K /* {{{ */
    if (trace_on_spot_ab(a,b) && fl->n) {
        verbose_output_print(1, 0, "# divided by 2 + primes from bucket that map to %u: ", x);
        if (!factor_list_fprint(stderr, *fl))
            verbose_output_print(1, 0, "(none)");
        verbose_output_vfprint(1, 0, gmp_vfprintf, ", remaining norm is %Zd\n", norm);
    }
#endif /* }}} */
    if (trial_div_very_verbose)
        verbose_output_vfprint(1, 0, gmp_vfprintf, "# x = %d, after dividing out bucket/resieved norm = %Zd\n", x, norm);

    do {
      /* Trial divide primes with precomputed tables */
#define TRIALDIV_MAX_FACTORS 32
      int i;
      unsigned long factors[TRIALDIV_MAX_FACTORS];
      if (trial_div_very_verbose) {
          /* FIXME: Multi-threading can garble this */
          verbose_output_print(1, 0, "# Trial division by ");
          for (i = 0; trialdiv_data[i].p != 1; i++)
              verbose_output_print(1, 0, " %lu", trialdiv_data[i].p);
          verbose_output_print(1, 0, "\n");
      }

      nr_factors = trialdiv (factors, norm, trialdiv_data, TRIALDIV_MAX_FACTORS);

      for (i = 0; i < MIN(nr_factors, TRIALDIV_MAX_FACTORS); i++)
      {
          if (trial_div_very_verbose)
              verbose_output_print (1, 0, " %lu", factors[i]);
          factor_list_add (fl, factors[i]);
      }
      if (trial_div_very_verbose)
          verbose_output_vfprint(1, 0, gmp_vfprintf, "\n# After trialdiv(): norm = %Zd\n", norm);
    } while (nr_factors == TRIALDIV_MAX_FACTORS + 1);
}
/* }}} */

/* Compute a checksum over the bucket region.

   We import the bucket region into an mpz_t and take it modulo
   checksum_prime. The checksums for different bucket regions are added up,
   modulo checksum_prime. This makes the combined checksum independent of
   the order in which buckets are processed, but it is dependent on size of
   the bucket region. Note that the selection of the sieve region, i.e., of J
   depends somewhat on the number of threads, as we want an equal number of
   bucket regions per thread. Thus the checksums are not necessarily
   comparable between runs with different numbers of threads. */

static const unsigned int checksum_prime = 4294967291; /* < 2^32 */

/* Combine two checksums. Simply (checksum+checksum2) % checksum_prime,
   but using modul_*() to handle sums >= 2^32 correctly. */
static unsigned int
combine_checksum(const unsigned int checksum1, const unsigned int checksum2)
{
    modulusul_t m;
    residueul_t r1, r2;
    unsigned int checksum;

    modul_initmod_ul(m, checksum_prime);
    modul_init(r1, m);
    modul_set_ul(r1, checksum1, m);
    modul_init(r2, m);
    modul_set_ul(r2, checksum2, m);
    modul_add(r1, r1, r2, m);
    checksum = modul_get_ul(r1, m);
    modul_clear(r1, m);
    modul_clear(r2, m);
    modul_clearmod(m);
    return checksum;
}

static unsigned int
bucket_checksum(const unsigned char *bucket, const unsigned int prev_checksum)
{
    mpz_t mb;
    unsigned long checksum;

    mpz_init(mb);
    mpz_import(mb, BUCKET_REGION, -1, sizeof(unsigned char), -1, 0, bucket);
    checksum = mpz_tdiv_ui(mb, checksum_prime);
    mpz_clear(mb);
    checksum = combine_checksum (checksum, prev_checksum);

    return checksum;
}

/* Adds the number of sieve reports to *survivors,
   number of survivors with coprime a, b to *coprimes */

NOPROFILE_STATIC int
factor_survivors (thread_data_ptr th, int N, unsigned char * S[2], where_am_I_ptr w MAYBE_UNUSED)
{
    las_info_ptr las = th->las;
    sieve_info_ptr si = th->si;
    int cpt = 0;
    int surv = 0, copr = 0;
    mpz_t norm[2];
    factor_list_t factors[2];
    mpz_array_t *f[2] = { NULL, };
    uint32_array_t *m[2] = { NULL, }; /* corresponding multiplicities */
    bucket_primes_t primes[2];
    mpz_t BLPrat;       /* alone ? */
#ifdef  DLP_DESCENT
    double deadline = DBL_MAX;
    relation_t winner[1];
    relation_init(winner);
#endif  /* DLP_DESCENT */
    uint32_t cof_rat_bitsize = 0; /* placate gcc */
    uint32_t cof_alg_bitsize = 0; /* placate gcc */
    const unsigned int first_j = N << (LOG_BUCKET_REGION - si->conf->logI);
    const unsigned long nr_lines = 1U << (LOG_BUCKET_REGION - si->conf->logI);

    for(int side = 0 ; side < 2 ; side++) {
        f[side] = alloc_mpz_array (1);
        m[side] = alloc_uint32_array (1);

        factor_list_init(&factors[side]);
        mpz_init (norm[side]);
    }

    mpz_init (BLPrat);
    mpz_set_ui (BLPrat, si->conf->sides[RATIONAL_SIDE]->lim);
    mpz_mul_2exp (BLPrat, BLPrat, si->conf->sides[RATIONAL_SIDE]->lpb); /* fb bound * lp bound */

    unsigned char * alg_S = S[ALGEBRAIC_SIDE];
    unsigned char * rat_S = S[RATIONAL_SIDE];
#ifdef TRACE_K /* {{{ */
    if (trace_on_spot_Nx(N, trace_Nx.x)) {
        fprintf(stderr, "# When entering factor_survivors for bucket %u, alg_S[%u]=%u, rat_S[%u]=%u\n",
                trace_Nx.N, trace_Nx.x, alg_S[trace_Nx.x], trace_Nx.x, rat_S[trace_Nx.x]);
        verbose_output_vfprint(1, 0, gmp_vfprintf, "# Remaining norms which have not been accounted for in sieving: (%Zd, %Zd)\n", traced_norms[0], traced_norms[1]);
    }
#endif  /* }}} */

    if (las->verbose >= 2) {
        /* Update the checksums over the bucket regions */
        for (int side = 0; side < 2; side++)
            th->checksum_post_sieve[side] = bucket_checksum(S[side], th->checksum_post_sieve[side]);
    }

    /* This is the one which gets the merged information in the end */
    unsigned char * SS = S[0];

#ifdef TRACE_K /* {{{ */
    sieve_side_info_ptr rat = si->sides[RATIONAL_SIDE];
    sieve_side_info_ptr alg = si->sides[ALGEBRAIC_SIDE];
    for (int x = 0; x < 1 << LOG_BUCKET_REGION; x++) {
        if (trace_on_spot_Nx(N, x)) {
            fprintf(stderr, "# alg->Bound[%u]=%u, rat->Bound[%u]=%u\n",
                    alg_S[trace_Nx.x], alg_S[x] <= alg->bound ? 0 : alg->bound,
                    rat_S[trace_Nx.x], rat_S[x] <= rat->bound ? 0 : rat->bound);
        }
    }
#endif /* }}} */

    for (unsigned int j = 0; j < nr_lines; j++)
    {
        unsigned char * const both_S[2] = {
            S[0] + (j << si->conf->logI), 
            S[1] + (j << si->conf->logI)
        };
        const unsigned char both_bounds[2] = {
            si->sides[0]->bound,
            si->sides[1]->bound,
        };
        surv += search_survivors_in_line(both_S, both_bounds,
                                         si->conf->logI, j + first_j, N,
                                         si->j_div, si->conf->unsieve_thresh,
                                         si->us);
        /* Make survivor search create a list of x-coordinates that survived
           instead of changing sieve array? More localized accesses in
           purge_bucket() that way */
    }

    /* Copy those bucket entries that belong to sieving survivors and
       store them with the complete prime */
    /* FIXME: choose a sensible size here */

    for(int side = 0 ; side < 2 ; side++) {
        WHERE_AM_I_UPDATE(w, side, side);
        primes[side] = init_bucket_primes (BUCKET_REGION);

        for (int i = 0; i < las->nb_threads; ++i) {
            thread_data_ptr other = th + i - th->id;
            purge_bucket (&primes[side], other->sides[side]->BA, N, SS);
        }

        /* Resieve small primes for this bucket region and store them 
           together with the primes recovered from the bucket updates */
        resieve_small_bucket_region (&primes[side], N, SS, th->si->sides[side]->rsd, th->sides[side]->rsdpos, si, w);

        /* Sort the entries to avoid O(n^2) complexity when looking for
           primes during trial division */
        bucket_sortbucket (&primes[side]);
    }

    /* Scan array one long word at a time. If any byte is <255, i.e. if
       the long word is != 0xFFFF...FF, examine the bytes 
       FIXME: We can use SSE to scan 16 bytes at a time, but have to make 
       sure that SS is 16-aligned first, thus currently disabled. */
#if defined(HAVE_SSE41) && defined(SSE_SURVIVOR_SEARCH)
    const size_t together = sizeof(__m128i);
    __m128i ones128 = (__m128i) {-1,-1};
    const __m128i * restrict SS_lw = (const __m128i *)SS;
#else
    const int together = sizeof(unsigned long);
    const unsigned long * restrict SS_lw = (const unsigned long *)SS;
#endif

    for ( ; (unsigned char *) SS_lw < SS + BUCKET_REGION; SS_lw++) {
#ifdef TRACE_K
        size_t trace_offset = (const unsigned char *) SS_lw - SS;
        if ((unsigned int) N == trace_Nx.N && (unsigned int) trace_offset <= trace_Nx.x && 
            (unsigned int) trace_offset + together > trace_Nx.x) {
            fprintf(stderr, "# Slot [%u] in bucket %u has value %u\n",
                    trace_Nx.x, trace_Nx.N, SS[trace_Nx.x]);
        }
#endif
#if defined(HAVE_SSE41) && defined(SSE_SURVIVOR_SEARCH)
        if (LIKELY(_mm_testc_si128(*SS_lw, ones128)))
            continue;
#else
        if (LIKELY(*SS_lw == (unsigned long)(-1L)))
            continue;
#endif
        size_t offset = (const unsigned char *) SS_lw - SS;
        for (size_t x = offset; x < offset + together; ++x) {
            if (SS[x] == 255) continue;

            th->rep->survivor_sizes[rat_S[x]][alg_S[x]]++;
            
            /* For factor_leftover_norm, we need to pass the information of the
             * sieve bound. If a cofactor is less than the square of the sieve
             * bound, it is necessarily prime. we implement this by keeping the
             * log to base 2 of the sieve limits on each side, and compare the
             * bitsize of the cofactor with their double.
             */
            int64_t a;
            uint64_t b;

#ifdef  DLP_DESCENT
            double t = seconds();

            if (t >= deadline) {
                /* Also break the outer loop (ugly!) */
#if defined(HAVE_SSE41) && defined(SSE_SURVIVOR_SEARCH)
                SS_lw = (const __m128i *)(SS + BUCKET_REGION);
#else
                SS_lw = (const unsigned long *)(SS + BUCKET_REGION);
#endif
                verbose_output_print(0, 1, "# [descent] Aborting, deadline passed\n");
                break;
            }
#endif  /* DLP_DESCENT */

            // Compute algebraic and rational norms.
            NxToAB (&a, &b, N, x, si);
#ifdef TRACE_K
            if (trace_on_spot_ab(a, b))
              fprintf (stderr, "# about to start cofactorization for (%"
                       PRId64 ",%" PRIu64 ")  %zu %u\n", a, b, x, SS[x]);
#endif
            /* since a,b both even were not sieved, either a or b should
             * be odd. However, exceptionally small norms, even without
             * sieving, may fall below the report bound (see tracker
             * issue #15437). Therefore it is safe to continue here. */
            // ASSERT((a | b) & 1);
            if (UNLIKELY(((a | b) & 1) == 0))
            {
                th->rep->both_even++;
                continue;
                /*
                pthread_mutex_lock(&io_mutex);
                fprintf (stderr, "# Error: a and b both even for N = %d, x = %d,\n"
                        "i = %d, j = %d, a = %ld, b = %lu\n",
                        N, x, ((x + N*BUCKET_REGION) & (si->I - 1))
                        - (si->I >> 1),
                        (x + N*BUCKET_REGION) >> si->conf->logI,
                        (long) a, (unsigned long) b);
                abort();
                pthread_mutex_unlock(&io_mutex);
                */
            }

            /* Since the q-lattice is exactly those (a, b) with
               a == rho*b (mod q), q|b  ==>  q|a  ==>  q | gcd(a,b) */
            /* FIXME: fast divisibility test here! */
            if (b == 0 || (mpz_cmp_ui(si->doing->p, b) <= 0 && b % mpz_get_ui(si->doing->p) == 0))
                continue;

            copr++;

            int pass = 1;

            int i;
            unsigned int j;
            for(int z = 0 ; pass && z < 2 ; z++) {
                int side = RATIONAL_SIDE ^ z;   /* start with rational */
                int lim = si->conf->sides[side]->lim;

                // Trial divide rational norm
                /* Compute the norms using the polynomials transformed to 
                   i,j-coordinates. The transformed polynomial on the 
                   special-q side is already divided by q */
                NxToIJ (&i, &j, N, x, si);
		mpz_poly_homogeneous_eval_siui (norm[side], si->sides[side]->fij, i, j);


#ifdef TRACE_K
                if (trace_on_spot_ab(a, b)) {
                    verbose_output_vfprint(1, 0, gmp_vfprintf, "# start trial division for norm=%Zd on %s side for (%" PRId64 ",%" PRIu64 ")\n",norm[side],sidenames[side],a,b);
                }
#endif
                trial_div (&factors[side], norm[side], N, x,
                        si->sides[side]->fb,
                        &primes[side], si->sides[side]->trialdiv_data,
                        lim, a, b);

                pass = check_leftover_norm (norm[side], si, side);
#ifdef TRACE_K
                if (trace_on_spot_ab(a, b)) {
                    verbose_output_vfprint(1, 0, gmp_vfprintf, "# checked leftover norm=%Zd on %s side for (%" PRId64 ",%" PRIu64 "): %d\n",norm[side],sidenames[side],a,b,pass);
                }
#endif
            }
            if (!pass) continue;

            if (stats != 0)
            {
                cof_rat_bitsize = mpz_sizeinbase (norm[RATIONAL_SIDE], 2);
                cof_alg_bitsize = mpz_sizeinbase (norm[ALGEBRAIC_SIDE], 2);
                if (stats == 1) /* learning phase */
                    /* no need to use a mutex here: either we use one thread only
                       to compute the cofactorization data and if several threads
                       the order is irrelevant. The only problem that can happen
                       is when two threads increase the value at the same time,
                       and it is increased by 1 instead of 2, but this should
                       happen rarely. */
                    cof_call[cof_rat_bitsize][cof_alg_bitsize] ++;
                else /* stats == 2: we use the learning data */
                {
                    /* we store the initial number of cofactorization calls in
                       cof_call[0][0] and the remaining nb in cof_succ[0][0] */
                    cof_call[0][0] ++;
                    /* Warning: the <= also catches cases when succ=call=0 */
                    if ((double) cof_succ[cof_rat_bitsize][cof_alg_bitsize] <
                            (double) cof_call[cof_rat_bitsize][cof_alg_bitsize] *
                            stats_prob)
                        continue;
                    cof_succ[0][0] ++;
                }
            }

            pass = factor_both_leftover_norms(norm, BLPrat, f, m, si);
#ifdef TRACE_K
            if (trace_on_spot_ab(a, b) && pass == 0)
              verbose_output_vfprint(1, 0, gmp_vfprintf, "# factor_leftover_norm failed for (%" PRId64 ",%" PRIu64 "), remains %Zd, %Zd unfactored\n",
                           a, b, norm[0], norm[1]);
#endif

            if (pass <= 0) continue; /* a factor was > 2^lpb, or some
                                        factorization was incomplete */

            /* yippee: we found a relation! */

            if (stats == 1) /* learning phase */
                cof_succ[cof_rat_bitsize][cof_alg_bitsize] ++;

            // ASSERT (bin_gcd_int64_safe (a, b) == 1);

            relation_t rel[1];
            memset(rel, 0, sizeof(rel));
            rel->a = a;
            rel->b = b; 
            for (int side = 0; side < 2; side++) {
                for(int i = 0 ; i < factors[side].n ; i++)
                    relation_add_prime(rel, side, factors[side].fac[i]);
                for (unsigned int i = 0; i < f[side]->length; ++i) {
                    if (!mpz_fits_ulong_p(f[side]->data[i]))
                        fprintf(stderr, "Warning: misprinted relation because of large prime of %zu bits at (%" PRId64 ",%" PRIu64 ")\n",
                                mpz_sizeinbase(f[side]->data[i], 2), a, b);
                    for (unsigned int j = 0; j < m[side]->data[i]; j++) {
                        relation_add_prime(rel, side, mpz_get_ui(f[side]->data[i]));
                    }
                }
            }

            relation_add_prime(rel, si->conf->side, mpz_get_ui(si->doing->p));

            relation_compress_rat_primes(rel);
            relation_compress_alg_primes(rel);

#ifdef TRACE_K
            if (trace_on_spot_ab(a, b)) {
                fprintf(stderr, "# Relation for (%" PRId64 ",%" PRIu64 ") printed\n", a, b);
            }
#endif

            {
                int do_check = th->las->suppress_duplicates;
                int is_dup = do_check && relation_is_duplicate(las->output, rel,
                        las->nb_threads, si);
                const char *comment = is_dup ? "# DUPE " : "";
                if (!is_dup)
                    cpt++;
                pthread_mutex_lock(&io_mutex);
                if (create_descent_hints) {
                    verbose_output_print(0, 1, "(%1.4f) ", seconds() - tt_qstart);
                }
                verbose_output_print(0, 3, "# i=%d, j=%u, lognorms = %hhu, %hhu\n",
                        i, j, S[0][x], S[1][x]);
                fprint_relation(las->output, rel, comment);
                pthread_mutex_unlock(&io_mutex);
            }

            /* Build histogram of lucky S[x] values */
            th->rep->report_sizes[S[RATIONAL_SIDE][x]][S[ALGEBRAIC_SIDE][x]]++;

#ifdef  DLP_DESCENT
            /* For descent mode: compute the expected time to finish given
             * the factor sizes, and deduce a deadline.  Assuming that not
             * all encountered factors are below the factor base bound, if we
             * expect an additional time T to finish the decomposition, we
             * keep looking for a better decomposition for a grace time which
             * is computed as x*T, for some configurable ratio x (one might
             * think of x=0.2 for instance. x is DESCENT_GRACE_TIME_RATIO),
             * which defines a ``deadline'' for next step.  [If all factors
             * happen to be smooth, the deadline is _now_.] If within the
             * grace period, a new relation is found, with an earlier implied
             * deadline, the deadline is updated.
             */
            
            if (las->hint_table) {
                double time_left = 0;
                for(int side = 0 ; side < 2 ; side++) {
                    for (unsigned int i = 0; i < f[side]->length; ++i) {
                        /* We assume that the base of the precomputed log
                         * goes as far as the lpb given on command line.
                         */
                        if (mpz_cmp_ui(f[side]->data[i],
                            (1UL<<las->default_config->sides[side]->lpb)) <= 0)
                            continue;
                        /*
                        if (mpz_cmp_ui(f[side]->data[i], si->conf->sides[side]->lim) <= 0)
                            continue;
                            */
                        unsigned int n = mpz_sizeinbase(f[side]->data[i], 2);
                        int k = las->hint_lookups[side][n];
                        if (k < 0) {
                            /* Having no data is normal. The factor base
                             * bound does not have to have any connection
                             * with the database of prime ideals for
                             * which the virtual log has been saved
                             * already (factor base extension can enter
                             * into play here). */
                            // verbose_output_print(0, 1, "# [descent] Warning: cannot estimate refactoring time for relation involving %d%c\n", n, sidenames[side][0]);
                            continue;
                        }
                        descent_hint_ptr h = las->hint_table[k];
                        time_left += h->expected_time;
                    }
                }
                // verbose_output_print(0, 1, "# [descent] This relation entails an additional time of %.2f for the smoothing process\n", time_left);

                double new_deadline = seconds() + DESCENT_GRACE_TIME_RATIO * time_left;
                if (new_deadline < deadline) {
                    double delta = DBL_MAX;
                    if (deadline != DBL_MAX)
                        delta = (deadline-new_deadline)/DESCENT_GRACE_TIME_RATIO;
                    deadline = new_deadline;
                    relation_copy(winner, rel);
                    if (time_left == 0) {
                        verbose_output_print(0, 1, "# [descent] Yiippee, splitting done\n");
                        /* break both loops */
#if defined(HAVE_SSE41) && defined(SSE_SURVIVOR_SEARCH)
                SS_lw = (const __m128i *)(SS + BUCKET_REGION);
#else
                SS_lw = (const unsigned long *)(SS + BUCKET_REGION);
#endif
                        relation_clear(rel);
                        break;
                    } else if (delta != DBL_MAX) {
                        verbose_output_print(0, 1, "# [descent] Improved ETA by %.2f\n", delta);
                    }
                }
            }
#endif  /* DLP_DESCENT */
            relation_clear(rel);
        }
    }

    th->rep->survivors1 += surv;
    th->rep->survivors2 += copr;

#ifdef  DLP_DESCENT
    if (las->hint_table && cpt) {
        /* If we've been doing the descent, and assuming we succeeded, we
         * have to reschedule the possibly still missing large primes in the
         * todo list */
        mpz_t rho, q;
        mpz_init(rho);
        mpz_init(q);
        for(int i = 0 ; i < winner->nb_rp ; i++) {
            int side = RATIONAL_SIDE;
            unsigned long p = winner->rp[i].p;
            /* See comment above */
            if (p <= (1UL<<las->default_config->sides[side]->lpb))
                continue;
            if (mpz_cmp_ui(si->doing->p, p)==0 && side == si->doing->side)
                continue;
            // If q > 64 bits, then the preivous comparison does not
            // work. Let's do a dirty hack, here, because I don't know
            // what would be the best patch. FIXME
            // So we check only the less significant bits to decide
            // whether we are seeing the current special-q.
            if (side == si->doing->side &&
                mpz_sizeinbase(si->doing->p, 2) >= 8*sizeof(unsigned long)) {
                unsigned long pp = mpz_get_ui(si->doing->p);
                if (pp == p)
                    continue;
            }
            unsigned int n = ULONG_BITS - clzl(p);
            int k = las->hint_lookups[side][n];
            if (k < 0) continue;
            mpz_set_ui(q, p);
            /* Beware, cpoly->m mod p would be wrong ! */
            /* This can't work on 32-bits */
            mpz_set_ui(rho, relation_compute_r (winner->a, winner->b, p));
            verbose_output_vfprint(0, 1, gmp_vfprintf, "# [descent] "HILIGHT_START"pushing %s (%Zd,%Zd) [%d%c]"HILIGHT_END" to todo list\n", sidenames[side], q, rho, mpz_sizeinbase(q, 2), sidenames[side][0]);
            las_todo_push_withdepth(&(las->todo), q, rho, side, si->doing->depth + 1);
        }
        relation_compute_all_r(winner);
        for(int i = 0 ; i < winner->nb_ap ; i++) {
            int side = ALGEBRAIC_SIDE;
            unsigned long p = winner->ap[i].p;
            /* See comment above */
            if (p <= (1UL<<las->default_config->sides[side]->lpb))
                continue;
            if (mpz_cmp_ui(si->doing->p, p)==0 && side == si->doing->side)
                continue;
            // If q > 64 bits, then the preivous comparison does not
            // work. Let's do a dirty hack, here, because I don't know
            // what would be the best patch. FIXME
            // So we check only the less significant bits to decide
            // whether we are seeing the current special-q.
            if (side == si->doing->side &&
                mpz_sizeinbase(si->doing->p, 2) >= 8*sizeof(unsigned long)) {
                unsigned long pp = mpz_get_ui(si->doing->p);
                if (pp == p)
                    continue;
            }
            unsigned int n = ULONG_BITS - clzl(p);
            int k = las->hint_lookups[side][n];
            if (k < 0) continue;
            mpz_set_ui(q, p);
            mpz_set_ui(rho, winner->ap[i].r);
            verbose_output_vfprint(0, 1, gmp_vfprintf, "# [descent] "HILIGHT_START"pushing %s (%Zd,%Zd) [%d%c]"HILIGHT_END" to todo list\n", sidenames[side], q, rho, mpz_sizeinbase(q, 2), sidenames[side][0]);
            las_todo_push_withdepth(&(las->todo), q, rho, side, si->doing->depth + 1);
        }
        mpz_clear(q);
        mpz_clear(rho);

        las_descent_helper_found_relation(las->descent_helper, winner);
    }
#endif  /* DLP_DESCENT */

    mpz_clear (BLPrat);

    for(int side = 0 ; side < 2 ; side++) {
        clear_bucket_primes (&primes[side]);
        mpz_clear(norm[side]);
        factor_list_clear(&factors[side]);
        clear_uint32_array (m[side]);
        clear_mpz_array (f[side]);
    }

    return cpt;
}

/* }}} */

/****************************************************************************/

MAYBE_UNUSED static inline void subusb(unsigned char *S1, unsigned char *S2, ssize_t offset)
{
  int ex = (unsigned int) S1[offset] - (unsigned int) S2[offset];
  if (UNLIKELY(ex < 0)) S1[offset] = 0; else S1[offset] = ex;	     
}

/* S1 = S1 - S2, with "-" in saturated arithmetical,
 * and memset(S2, 0, EndS1-S1).
 */
void SminusS (unsigned char *S1, unsigned char *EndS1, unsigned char *S2) {
#ifndef HAVE_SSE2
  ssize_t mysize = EndS1 - S1;
  unsigned char *cS2 = S2;
  while (S1 < EndS1) {
    subusb(S1,S2,0);
    subusb(S1,S2,1);
    subusb(S1,S2,2);
    subusb(S1,S2,3);
    subusb(S1,S2,4);
    subusb(S1,S2,5);
    subusb(S1,S2,6);
    subusb(S1,S2,7);
    S1 += 8; S2 += 8;
  }
  memset(cS2, 0, mysize);
#else
  __m128i *S1i = (__m128i *) S1, *EndS1i = (__m128i *) EndS1, *S2i = (__m128i *) S2,
    z = _mm_setzero_si128();
  while (S1i < EndS1i) {
    __m128i x0, x1, x2, x3;
    __asm__ __volatile__
      ("prefetcht0 0x1000(%0)\n"
       "prefetcht0 0x1000(%1)\n"
       "movdqa (%0),%2\n"
       "movdqa 0x10(%0),%3\n"
       "movdqa 0x20(%0),%4\n"
       "movdqa 0x30(%0),%5\n"
       "psubusb (%1),%2\n"
       "psubusb 0x10(%1),%3\n"
       "psubusb 0x20(%1),%4\n"
       "psubusb 0x30(%1),%5\n"
       "movdqa %6,(%1)\n"
       "movdqa %6,0x10(%1)\n"
       "movdqa %6,0x20(%1)\n"
       "movdqa %6,0x30(%1)\n"
       "movdqa %2,(%0)\n"
       "movdqa %3,0x10(%0)\n"
       "movdqa %4,0x20(%0)\n"
       "movdqa %5,0x30(%0)\n"
       "add $0x40,%0\n"
       "add $0x40,%1\n"
       : "+&r"(S1i), "+&r"(S2i), "=&x"(x0), "=&x"(x1), "=&x"(x2), "=&x"(x3) : "x"(z));
    /* I prefer use ASM than intrinsics to be sure each 4 instructions which
     * use exactly a cache line are together. I'm 99% sure it's not useful...
     * but it's more beautiful :-)
     */
    /*
    __m128i x0, x1, x2, x3;
    _mm_prefetch(S1i + 16, _MM_HINT_T0); _mm_prefetch(S2i + 16, _MM_HINT_T0);
    x0 = _mm_load_si128(S1i + 0);         x1 = _mm_load_si128(S1i + 1);
    x2 = _mm_load_si128(S1i + 2);         x3 = _mm_load_si128(S1i + 3);
    x0 = _mm_subs_epu8(S2i[0], x0);       x1 = _mm_subs_epu8(S2i[1], x1);
    x2 = _mm_subs_epu8(S2i[2], x2);       x3 = _mm_subs_epu8(S2i[3], x3);
    _mm_store_si128(S2i + 0, z);          _mm_store_si128(S1i + 1, z);
    _mm_store_si128(S2i + 2, z);          _mm_store_si128(S1i + 3, z);
    _mm_store_si128(S1i + 0, x0);         _mm_store_si128(S1i + 1, x1);
    _mm_store_si128(S1i + 2, x2);         _mm_store_si128(S1i + 3, x3);
    S1i += 4; S2i += 4;
    */
  }
#endif 
}

/* Move above ? */
/* {{{ process_bucket_region
 * th->id gives the number of the thread: it is supposed to deal with the set
 * of bucket_regions corresponding to that number, ie those that are
 * congruent to id mod nb_thread.
 *
 * The other threads are accessed by combining the thread pointer th and
 * the thread id: the i-th thread is at th - id + i
 */
void * process_bucket_region(thread_data_ptr th)
{
    where_am_I w MAYBE_UNUSED;
    las_info_ptr las = th->las;
    sieve_info_ptr si = th->si;

    WHERE_AM_I_UPDATE(w, si, si);

    las_report_ptr rep = th->rep;

    WHERE_AM_I_UPDATE(w, N, th->id);

    unsigned char * S[2];

    unsigned int my_row0 = (BUCKET_REGION >> si->conf->logI) * th->id;
    unsigned int skiprows = (BUCKET_REGION >> si->conf->logI)*(las->nb_threads-1);


    /* This is local to this thread */
    for(int side = 0 ; side < 2 ; side++) {
        sieve_side_info_ptr s = si->sides[side];
        thread_side_data_ptr ts = th->sides[side];

        ts->ssdpos = small_sieve_start(s->ssd, my_row0, si);
        ts->rsdpos = small_sieve_copy_start(ts->ssdpos, s->fb_parts_x->rs);

        /* local sieve region */
        S[side] = ts->bucket_region;
    }
    unsigned char *SS = th->SS;
    memset(SS, 0, BUCKET_REGION);

    /* loop over appropriate set of sieve regions */
    for (unsigned int i = th->id; i < si->nb_buckets; i += las->nb_threads) 
      {
        WHERE_AM_I_UPDATE(w, N, i);

        {
            const int side = RATIONAL_SIDE;
            WHERE_AM_I_UPDATE(w, side, side);

            sieve_side_info_ptr s = si->sides[side];
            thread_side_data_ptr ts = th->sides[side];
        
            /* Init rational norms */
            rep->tn[side] -= seconds_thread ();
#ifdef SMART_NORM
	    init_norms_bucket_region(S[side], i, si, side, 1);
#else
	    init_norms_bucket_region(S[side], i, si, side, 0);
#endif
            // Invalidate the first row except (1,0)
            if (!i) {
                int pos10 = 1+((si->I)>>1);
                unsigned char n10 = S[side][pos10];
                memset(S[side], 255, si->I);
                S[side][pos10] = n10;
            }
            rep->tn[side] += seconds_thread ();
#if defined(TRACE_K) 
            if (trace_on_spot_N(w->N))
              fprintf (stderr, "# After rationals init_norms_bucket_region, N=%u S[%u]=%u\n", w->N, trace_Nx.x, S[side][trace_Nx.x]);
#endif

            /* Apply rational buckets */
            rep->ttbuckets_apply -= seconds_thread();
            for (int j = 0; j < las->nb_threads; ++j)  {
                thread_data_ptr ot = th + j - th->id;
                apply_one_bucket(SS, ot->sides[side]->BA, i, w);
            }
	    SminusS(S[side], S[side] + BUCKET_REGION, SS);
            rep->ttbuckets_apply += seconds_thread();

            /* Sieve small rational primes */
            sieve_small_bucket_region(SS, i, s->ssd, ts->ssdpos, si, side, w);
	    SminusS(S[side], S[side] + BUCKET_REGION, SS);
#if defined(TRACE_K) 
            if (trace_on_spot_N(w->N))
              fprintf (stderr, "# Final value on rational side, N=%u rat_S[%u]=%u\n", w->N, trace_Nx.x, S[side][trace_Nx.x]);
#endif
        }

        {
            const int side = ALGEBRAIC_SIDE;
            WHERE_AM_I_UPDATE(w, side, side);

            sieve_side_info_ptr s = si->sides[side];
            thread_side_data_ptr ts = th->sides[side];

            /* Init algebraic norms */
            rep->tn[side] -= seconds_thread ();

#ifdef SMART_NORM
	    init_norms_bucket_region(S[side], i, si, side, 1);
#else
            init_norms_bucket_region(S[side], i, si, side, 0);
#endif
            rep->tn[side] += seconds_thread ();
#if defined(TRACE_K) 
            if (trace_on_spot_N(w->N))
              fprintf (stderr, "# After algebraics init_norms_bucket_region, N=%u S[%u]=%u\n", w->N, trace_Nx.x, S[side][trace_Nx.x]);
#endif

            /* Apply algebraic buckets */
            rep->ttbuckets_apply -= seconds_thread();
            for (int j = 0; j < las->nb_threads; ++j) {
                thread_data_ptr ot = th + j - th->id;
                apply_one_bucket(SS, ot->sides[side]->BA, i, w);
            }
	    SminusS(S[side], S[side] + BUCKET_REGION, SS);
            rep->ttbuckets_apply += seconds_thread();

            /* Sieve small algebraic primes */
            sieve_small_bucket_region(SS, i, s->ssd, ts->ssdpos, si, side, w);
	    SminusS(S[side], S[side] + BUCKET_REGION, SS);
#if defined(TRACE_K) 
            if (trace_on_spot_N(w->N))
              fprintf (stderr, "# Final value on algebraic side, N=%u alg_S[%u]=%u\n", w->N, trace_Nx.x, S[side][trace_Nx.x]);
#endif
        }

        /* Factor survivors */
        rep->ttf -= seconds_thread ();
        rep->reports += factor_survivors (th, i, S, w);
        rep->ttf += seconds_thread ();

        /* For the descent mode, we bail out as early as possible */
        if (las->hint_table && rep->reports) break;

        for(int side = 0 ; side < 2 ; side++) {
            sieve_side_info_ptr s = si->sides[side];
            thread_side_data_ptr ts = th->sides[side];
            small_sieve_skip_stride(s->ssd, ts->ssdpos, skiprows, si);
            int * b = s->fb_parts_x->rs;
            memcpy(ts->rsdpos, ts->ssdpos + b[0], (b[1]-b[0]) * sizeof(int));
        }
      }

    for(int side = 0 ; side < 2 ; side++) {
        thread_side_data_ptr ts = th->sides[side];
        free(ts->ssdpos);
        free(ts->rsdpos);
    }

    return NULL;
}/*}}}*/

/* thread handling */
static thread_data * thread_data_alloc(las_info_ptr las, int n)/*{{{*/
{
    thread_data * thrs = (thread_data *) malloc(n * sizeof(thread_data));
    ASSERT_ALWAYS(thrs);
    memset(thrs, 0, n * sizeof(thread_data));

    for(int i = 0 ; i < n ; i++) {
        thrs[i]->id = i;
        thrs[i]->las = las;
        las_report_init(thrs[i]->rep);
        thrs[i]->checksum_post_sieve[0] = thrs[i]->checksum_post_sieve[1] = 0;

        /* Allocate memory for each thread's two bucket regions (one for each
           side) and for the intermediate sum (only one for both sides) */
        for (int side = 0; side < 2; side++) {
          thread_side_data_ptr ts = thrs[i]->sides[side];
          // printf ("# Allocating thrs[%d]->sides[%d]->bucket_region\n", i, side);
          ts->bucket_region = (unsigned char *) contiguous_malloc(BUCKET_REGION + MEMSET_MIN);
        }
        thrs[i]->SS = (unsigned char *) contiguous_malloc(BUCKET_REGION);

        
    }
    return thrs;
}/*}}}*/

static void thread_data_free(thread_data * thrs, int n)/*{{{*/
{
    for (int i = 0; i < n ; ++i) {
        las_report_clear(thrs[i]->rep);
        ASSERT_ALWAYS(thrs[i]->SS != NULL);

        /* Free the memory for the bucket regions */
        for (int side = 0; side < 2; side++) {
          thread_side_data_ptr ts = thrs[i]->sides[side];
          // printf ("# Freeing thrs[%d]->sides[%d]->bucket_region\n", i, side);
          contiguous_free(ts->bucket_region);
          ts->bucket_region = NULL;
        }
        contiguous_free(thrs[i]->SS);
        thrs[i]->SS = NULL;
    }
    free(thrs);
}/*}}}*/

void thread_pickup_si(thread_data * thrs, sieve_info_ptr si, int n)/*{{{*/
{
    for (int i = 0; i < n ; ++i) {
        thrs[i]->si = si;
        for(int s = 0 ; s < 2 ; s++) {
            thrs[i]->sides[s]->fb_bucket = si->sides[s]->fb_bucket_threads[i];
            thrs[i]->sides[s]->log_steps = si->sides[s]->log_steps;
            thrs[i]->sides[s]->log_steps_max = si->sides[s]->log_steps_max;
        }
    }
}/*}}}*/

/* {{{ thread_buckets_alloc
 * TODO: Allow allocating larger buckets if bucket_fill_ratio ever grows
 * above the initial guess. This can easily be made a permanent choice.
 *
 * Note also that we could consider having bucket_fill_ratio global.
 */
static void thread_buckets_alloc(thread_data *thrs, unsigned int n)/*{{{*/
{
  for (unsigned int i = 0; i < n ; ++i) {
    for(unsigned int side = 0 ; side < 2 ; side++) {
      thread_side_data_ptr ts = thrs[i]->sides[side];
      sieve_info_srcptr si = thrs[i]->si;
      /* We used to re-allocate whenever the number of buckets changed. Now we
         always allocate memory for the max. number of buckets, so that we
         never have to re-allocate */
      uint32_t nb_buckets;
      /* If shell environment variable LAS_REALLOC_BUCKETS is *not* set,
         always allocate memory for the max. number of buckets, so that we
         never have to re-allocate. If it is set, allocate just enough for
         for the current number of buckets, which will re-allocate memory
         if number of buckets changes. */
      if (getenv("LAS_REALLOC_BUCKETS") == NULL) {
        nb_buckets = si->nb_buckets_max;
      } else {
        nb_buckets = si->nb_buckets;
      }
      init_buckets(&(ts->BA), &(ts->kBA), &(ts->mBA),
                   si->sides[side]->max_bucket_fill_ratio, nb_buckets);
    }
  }
}/*}}}*/

static void thread_buckets_free(thread_data * thrs, unsigned int n)/*{{{*/
{
  for (unsigned int i = 0; i < n ; ++i) {
    thread_side_data_ptr ts;
    for(unsigned int side = 0 ; side < 2 ; side++) {
      // fprintf ("# Freeing buckets, thread->id=%d, side=%d\n", thrs[i]->id, side);
      ts = thrs[i]->sides[side];
      /* if there is no special-q in the interval, the arrays are not malloced */
      if (ts->BA.bucket_write != NULL)
        clear_buckets(&(ts->BA), &(ts->kBA), &(ts->mBA));
    }
  }
}/*}}}*/

static double thread_buckets_max_full(thread_data * thrs, int n)/*{{{*/
{
    double mf, mf0 = 0;
    for (int i = 0; i < n ; ++i) {
        mf = buckets_max_full (thrs[i]->sides[RATIONAL_SIDE]->BA);
        if (mf > mf0) mf0 = mf;
        mf = buckets_max_full (thrs[i]->sides[ALGEBRAIC_SIDE]->BA);
        if (mf > mf0) mf0 = mf;
    }
    return mf0;
}/*}}}*/

/* {{{ las_report_accumulate_threads_and_display
 * This function does three distinct things.
 *  - accumulates the timing reports for all threads into a collated report
 *  - display the per-sq timing relative to this report, and the given
 *    timing argument (in seconds).
 *  - merge the per-sq report into a global report
 */
void las_report_accumulate_threads_and_display(las_info_ptr las, sieve_info_ptr si, las_report_ptr report, thread_data * thrs, double qt0)
{
    /* Display results for this special q */
    las_report rep;
    las_report_init(rep);
    unsigned int checksum_post_sieve[2] = {0, 0};
    
    for (int i = 0; i < las->nb_threads; ++i) {
        las_report_accumulate(rep, thrs[i]->rep);
        for (int side = 0; side < 2; side++)
            checksum_post_sieve[side] = combine_checksum(checksum_post_sieve[side], thrs[i]->checksum_post_sieve[side]);
    }
    verbose_output_print(0, 2, "# ");
    /* verbose_output_print(0, 2, "%lu survivors after rational sieve,", rep->survivors0); */
    verbose_output_print(0, 2, "%lu survivors after algebraic sieve, ", rep->survivors1);
    verbose_output_print(0, 2, "coprime: %lu\n", rep->survivors2);
    verbose_output_print(0, 2, "# Checksums over sieve region: after all sieving: %u, %u\n", checksum_post_sieve[0], checksum_post_sieve[1]);
    verbose_output_vfprint(0, 1, gmp_vfprintf, "# %lu relation(s) for %s (%Zd,%Zd)\n", rep->reports, sidenames[si->conf->side], si->doing->p, si->doing->r);
    double qtts = qt0 - rep->tn[0] - rep->tn[1] - rep->ttf;
    if (rep->both_even) {
        verbose_output_print(0, 1, "# Warning: found %lu hits with i,j both even (not a bug, but should be very rare)\n", rep->both_even);
    }
#ifdef HAVE_RUSAGE_THREAD
    int dont_print_tally = 0;
#else
    int dont_print_tally = 1;
#endif
    if (dont_print_tally && las->nb_threads > 1) {
        verbose_output_print(0, 1, "# Time for this special-q: %1.4fs [tally available only in mono-thread]\n", qt0);
    } else {
        verbose_output_print(0, 1, "# Time for this special-q: %1.4fs [norm %1.4f+%1.4f, sieving %1.4f"
            " (%1.4f + %1.4f + %1.4f),"
            " factor %1.4f]\n", qt0,
            rep->tn[RATIONAL_SIDE],
            rep->tn[ALGEBRAIC_SIDE],
            qtts,
            rep->ttbuckets_fill,
            rep->ttbuckets_apply,
            qtts-rep->ttbuckets_fill-rep->ttbuckets_apply,
            rep->ttf);
    }
    las_report_accumulate(report, rep);
    las_report_clear(rep);
}/*}}}*/


/*************************** main program ************************************/


static void declare_usage(param_list pl)
{
  param_list_usage_header(pl,
  "In the names and in the descriptions of the parameters, below there are often\n"
  "aliases corresponding to the convention that 0 is the rational side and 1\n"
  "is the algebraic side. If the two sides are algebraic, then the word\n"
  "'rational' just means the side number 0. Note also that for a rational\n"
  "side, the factor base is recomputed on the fly (or cached), and there is\n"
  "no need to provide a fb0 parameter.\n"
  );

  param_list_decl_usage(pl, "poly", "polynomial file");
  param_list_decl_usage(pl, "fb0",   "factor base file on the rational side");
  param_list_decl_usage(pl, "fb1",   "(alias fb) factor base file on the algebraic side");
  param_list_decl_usage(pl, "fbc",  "factor base cache file");
  param_list_decl_usage(pl, "q0",   "left bound of special-q range");
  param_list_decl_usage(pl, "q1",   "right bound of special-q range");
  param_list_decl_usage(pl, "rho",  "sieve only root r mod q0");
  param_list_decl_usage(pl, "v",    "(switch) verbose mode, also prints sieve-area checksums");
  param_list_decl_usage(pl, "out",  "filename where relations are written, instead of stdout");
  param_list_decl_usage(pl, "t",   "number of threads to use");
  param_list_decl_usage(pl, "ratq", "(switch) use rational special-q");

  param_list_decl_usage(pl, "I",    "set sieving region to 2^I");
  param_list_decl_usage(pl, "skew", "(alias S) skewness");
  param_list_decl_usage(pl, "lim0", "(alias rlim) rational factor base bound");
  param_list_decl_usage(pl, "lim1", "(alias alim) algebraic factor base bound");
  param_list_decl_usage(pl, "lpb0", "(alias lpbr) set rational large prime bound to 2^lpb0");
  param_list_decl_usage(pl, "lpb1", "(alias lpba) set algebraic large prime bound to 2^lpb1");
  param_list_decl_usage(pl, "mfb0", "(alias mfbr) set rational cofactor bound 2^mfb0");
  param_list_decl_usage(pl, "mfb1", "(alias mfba) set algebraic cofactor bound 2^mfb1");
  param_list_decl_usage(pl, "lambda0", "(alias rlambda) rational lambda value");
  param_list_decl_usage(pl, "lambda1", "(alias alambda) algebraic lambda value");
  param_list_decl_usage(pl, "powlim0", "(alias rpowlim) limit on powers on rat side");
  param_list_decl_usage(pl, "powlim1", "(alias apowlim) limit on powers on alg side");
  param_list_decl_usage(pl, "tdthresh", "trial-divide primes p/r <= ththresh (r=number of roots)");
  param_list_decl_usage(pl, "bkthresh", "bucket-sieve primes p >= bkthresh");
  param_list_decl_usage(pl, "unsievethresh", "Unsieve all p > unsievethresh where p|gcd(a,b)");

  param_list_decl_usage(pl, "allow-largesq", "(switch) allows large special-q, e.g. for a DL descent");
  param_list_decl_usage(pl, "stats-stderr", "(switch) print stats to stderr in addition to stdout/out file");
  param_list_decl_usage(pl, "todo", "provide file with a list of special-q to sieve instead of qrange");
  param_list_decl_usage(pl, "descent-hint", "filename with tuned data for the descent, for each special-q bitsize");
  param_list_decl_usage(pl, "mkhint", "(switch) _create_ a descent file, instead of reading one");
  param_list_decl_usage(pl, "no-prepare-hints", "(switch) defer initialization of siever precomputed structures (one per special-q side) to time of first actual use");
  param_list_decl_usage(pl, "dup", "(switch) suppress duplicate relations");
  param_list_decl_usage(pl, "galois", "(switch) for reciprocal polynomials, sieve only half of the q's");
#ifdef TRACE_K
  param_list_decl_usage(pl, "traceab", "Relation to trace, in a,b format");
  param_list_decl_usage(pl, "traceij", "Relation to trace, in i,j format");
  param_list_decl_usage(pl, "traceNx", "Relation to trace, in N,x format");
#endif
  verbose_decl_usage(pl);
}

int main (int argc0, char *argv0[])/*{{{*/
{
    las_info las;
    double t0, tts, wct;
    unsigned long nr_sq_processed = 0;
    int allow_largesq = 0;
    double totJ = 0.0;
    int argc = argc0;
    char **argv = argv0;
    double max_full = 0.;

#ifdef HAVE_MINGW
    _fmode = _O_BINARY;     /* Binary open for all files */
#endif

    param_list pl;
    param_list_init(pl);

    declare_usage(pl);

    /* Passing NULL is allowed here. Find value with
     * param_list_parse_switch later on */
    param_list_configure_switch(pl, "-v", NULL);
    param_list_configure_switch(pl, "-ratq", NULL);
    param_list_configure_switch(pl, "-no-prepare-hints", NULL);
    param_list_configure_switch(pl, "-allow-largesq", &allow_largesq);
    param_list_configure_switch(pl, "-stats-stderr", NULL);
    param_list_configure_switch(pl, "-mkhint", &create_descent_hints);
    param_list_configure_switch(pl, "-dup", NULL);
    param_list_configure_switch(pl, "-galois", NULL);
    param_list_configure_alias(pl, "-skew", "-S");
    param_list_configure_alias(pl, "-fb1", "-fb");
    param_list_configure_alias(pl, "-lim0", "-rlim");
    param_list_configure_alias(pl, "-lim1", "-alim");
    param_list_configure_alias(pl, "-lpb0", "-lpbr");
    param_list_configure_alias(pl, "-lpb1", "-lpba");
    param_list_configure_alias(pl, "-mfb0", "-mfbr");
    param_list_configure_alias(pl, "-mfb1", "-mfba");
    param_list_configure_alias(pl, "-lambda0", "-rlambda");
    param_list_configure_alias(pl, "-lambda1", "-alambda");
    param_list_configure_alias(pl, "-powlim0", "-rpowlim");
    param_list_configure_alias(pl, "-powlim1", "-apowlim");

    argv++, argc--;
    for( ; argc ; ) {
        if (param_list_update_cmdline(pl, &argc, &argv)) { continue; }
        /* Could also be a file */
        FILE * f;
        if ((f = fopen(argv[0], "r")) != NULL) {
            param_list_read_stream(pl, f);
            fclose(f);
            argv++,argc--;
            continue;
        }
        fprintf(stderr, "Unhandled parameter %s\n", argv[0]);
        param_list_print_usage(pl, argv0[0], stderr);
        exit(EXIT_FAILURE);
    }

    las_info_init(las, pl);    /* side effects: prints cmdline and flags */

    /*
    int descent_bootstrap = param_list_lookup_string(pl, "target") != NULL;
    */
    int descent_lower = param_list_lookup_string(pl, "descent-hint") != NULL;

    /* We have the following dependency chain (not sure the account below
     * is exhaustive).
     *
     * q0 -> q0d (double) -> si->sides[*]->{scale,logmax}
     * q0 -> (I, lpb, lambda) for the descent
     * 
     * scale -> logp's in factor base.
     *
     * I -> splittings of the factor base among threads.
     *
     * This is probably enough to justify having separate sieve_info's
     * for the given sizes.
     */

    if (!param_list_parse_switch(pl, "-no-prepare-hints")) {
        /* Create a default siever instance among las->sievers if needed */
        if (las->default_config->bitsize)
            get_sieve_info_from_config(las, las->default_config, pl);

        /* Create all siever configurations from the preconfigured hints */
        /* This can also be done dynamically if needed */
        if (las->hint_table) {
            for(descent_hint_ptr h = las->hint_table[0] ; h->conf->bitsize ; h++) {
                siever_config_ptr sc = h->conf;
                get_sieve_info_from_config(las, sc, pl);
            }
        }
    }

#if defined(HAVE_GCC_STYLE_AMD64_INLINE_ASM) && defined(LAS_MEMSET)
    extern size_t max_cache, min_stos;
    max_cache = direct_write_vs_stos ();
    min_stos = stos_vs_write128 ();
    verbose_output_print(0, 1, "# Las normalisation memset: movaps from 33(0x21) to %zu(0x%zx); rep stosq until %zu(0x%zx); movntps after\n", min_stos, min_stos, max_cache, max_cache);
#endif
    
    thread_data * thrs = thread_data_alloc(las, las->nb_threads);

    las_report report;
    las_report_init(report);

    t0 = seconds ();
    wct = wct_seconds();
    verbose_output_print(0, 1, "#\n");

    where_am_I w MAYBE_UNUSED;
    WHERE_AM_I_UPDATE(w, las, las);

    if (descent_lower)
        las->descent_helper = las_descent_helper_alloc();

    /* {{{ Doc on todo list handling
     * The function las_todo_feed behaves in different
     * ways depending on whether we're in q-range mode or in q-list mode.
     *
     * q-range mode: the special-q's to be handled are specified as a
     * range. Then, whenever the las->todo list almost runs out, it is
     * refilled if possible, up to the limit q1 (-q0 & -rho just gives a
     * special case of this).
     *
     * q-list mode: the special-q's to be handled are always read from a
     * file. Therefore each new special-q to be handled incurs a
     * _blocking_ read on the file, until EOF. This mode is also used for
     * the descent, which has the implication that the read occurs if and
     * only if the todo list is empty. }}} */

    /* pop() is achieved by sieve_info_pick_todo_item */
    for( ; las_todo_feed(las, pl) ; ) {
        if (descent_lower && las_todo_pop_closing_brace(&(las->todo))) {
            las_descent_helper_done_node(las->descent_helper);
            if (las_descent_helper_current_depth(las->descent_helper) == 0)
                las_descent_helper_display_last_tree(las->descent_helper, las->output);
            continue;
        }

        siever_config current_config;
        memcpy(current_config, las->default_config, sizeof(siever_config));
        current_config->bitsize = mpz_sizeinbase(las->todo->p, 2);
        current_config->side = las->todo->side;

        /* Do we have a hint table with specifically tuned parameters,
         * well suited to this problem size ? */
        if (las->hint_table) {
            for(descent_hint_ptr h = las->hint_table[0] ; h->conf->bitsize ; h++) {
                siever_config_ptr sc = h->conf;
                if (!siever_config_match(sc, current_config)) continue;
                verbose_output_print(0, 1, "# Using existing sieving parameters from hint list for q~2^%d on the %s side [%d%c]\n", sc->bitsize, sidenames[sc->side], sc->bitsize, sidenames[sc->side][0]);
                memcpy(current_config, sc, sizeof(siever_config));
            }
        }

        /* Maybe create a new siever ? */
        sieve_info_ptr si = get_sieve_info_from_config(las, current_config, pl);
        WHERE_AM_I_UPDATE(w, si, si);

        sieve_info_pick_todo_item(si, &(las->todo));

        las_descent_helper_new_node(las->descent_helper, si->conf, si->doing->depth);
        if (descent_lower)
            las_todo_push_closing_brace(&(las->todo), si->doing->depth);
#if 0
        /* I think that there is sufficient provision in las_todo_feed now */
        ASSERT_ALWAYS(mpz_cmp_ui(las->todo->r, 0) != 0);
        if (si->cpoly->pols[si->doing->pside]->deg == 1) {
            /* compute the good rho */
            int n;
            n = mpz_poly_roots (&si->doing->r, si->cpoly->pols[si->doing->pside], si->doing->p);
            ASSERT_ALWAYS(n);
        } else {
            mpz_set(si->doing->r, las->todo->r);
        }
#endif

        /* Check whether q is larger than the large prime bound.
         * This can create some problems, for instance in characters.
         * By default, this is not allowed, but the parameter
         * -allow-largesq is a by-pass to this test.
         */
        if (!allow_largesq) {
            if ((int)mpz_sizeinbase(si->doing->p, 2) >
                    si->conf->sides[si->conf->side]->lpb) {
                fprintf(stderr, "ERROR: The special q is larger than the "
                        "large prime bound.\n");
                fprintf(stderr, "       You can disable this check with "
                        "the -allow-largesq argument,\n");
                fprintf(stderr, "       It is for instance useful for the "
                        "descent.\n");
                exit(EXIT_FAILURE);
            }
        }

        double qt0 = seconds();
        tt_qstart = seconds();

        if (SkewGauss (si) != 0)
            continue;

        /* check |a0|, |a1| < 2^31 if we use fb_root_in_qlattice_31bits */
#ifndef SUPPORT_LARGE_Q
        if (si->a0 <= -2147483648L || 2147483648L <= si->a0 ||
            si->a1 <= -2147483648L || 2147483648L <= si->a1)
          {
            fprintf (stderr, "Error, too large special-q, define SUPPORT_LARGE_Q. Skipping this special-q.\n");
            continue;
          }
#endif

        /* FIXME: maybe we can discard some special q's if a1/a0 is too large,
         * see http://www.mersenneforum.org/showthread.php?p=130478 
         *
         * Just for correctness at the moment we're doing it in very
         * extreme cases, see bug 15617
         */
        if (sieve_info_adjust_IJ(si, las->nb_threads) == 0) {
            verbose_output_vfprint(0, 1, gmp_vfprintf, "# "HILIGHT_START"Discarding %s q=%Zd; rho=%Zd;"HILIGHT_END" a0=%" PRId64 "; b0=%" PRId64 "; a1=%" PRId64 "; b1=%" PRId64 "; raw_J=%u;\n",
                    sidenames[si->conf->side],
                    si->doing->p, si->doing->r, si->a0, si->b0, si->a1, si->b1,
                    si->J);
            continue;
        }


        verbose_output_vfprint(0, 1, gmp_vfprintf, "# "HILIGHT_START"Sieving %s q=%Zd; rho=%Zd;"HILIGHT_END" a0=%" PRId64 "; b0=%" PRId64 "; a1=%" PRId64 "; b1=%" PRId64 ";",
                sidenames[si->conf->side],
                si->doing->p, si->doing->r, si->a0, si->b0, si->a1, si->b1);
        if (si->doing->depth) {
            verbose_output_print(0, 1, " # within descent, currently at depth %d", si->doing->depth);
        }
        verbose_output_print(0, 1, "\n");
        nr_sq_processed ++;

        /* checks the value of J,
         * precompute the skewed polynomials of f(x) and g(x), and also
         * their floating-point versions */
        sieve_info_update (las->output, si, las->nb_threads);
        totJ += (double) si->J;
        verbose_output_print(0, 2, "# I=%u; J=%u\n", si->I, si->J);
        if (las->verbose >= 2) {
            fprintf (las->output, "# f_0'(x) = ");
            mpz_poly_fprintf(las->output, si->sides[0]->fij);
            fprintf (las->output, "# f_1'(x) = ");
            mpz_poly_fprintf(las->output, si->sides[1]->fij);
        }

#ifdef TRACE_K
        init_trace_k(si, pl);
#endif

        /* WARNING. We're filling the report info for thread 0 for
         * ttbuckets_fill, while in fact the cost is over all threads.
         * The problem is always the fact that we don't have a proper
         * per-thread timer. So short of a better solution, this is what
         * we do. True, it's not clean.
         * (we need this data to be caught by
         * las_report_accumulate_threads_and_display further down, hence
         * this hack).
         */
        thrs[0]->rep->ttbuckets_fill -= seconds();

        thread_pickup_si(thrs, si, las->nb_threads);

        /* Allocate buckets */
        thread_buckets_alloc(thrs, las->nb_threads);

        /* Fill in rat and alg buckets */
        thread_do(thrs, &fill_in_buckets_both, las->nb_threads);

        max_full = MAX(max_full, thread_buckets_max_full(thrs, las->nb_threads));
        ASSERT_ALWAYS(max_full <= 1.0 || /* see commented code below */
                 fprintf (stderr, "max_full=%f, see #14987\n", max_full) == 0);
#if 0   /* {{{ I no longer believe we can save something if this happens */
        /* See bug #14987 on the tracker */
        if (max_full >= 1.0) {
            for (i = 0; i < las->nb_threads; ++i) {
                fprintf(stderr, "intend to free [%d] max_full=%f %f\n",
                        i,
                        buckets_max_full (thrs[i]->sides[RATIONAL_SIDE]->BA),
                        buckets_max_full (thrs[i]->sides[ALGEBRAIC_SIDE]->BA));
            }
            thread_buckets_free(thrs); /* may crash. See below */

            si->bucket_limit_multiplier *= 1.1 * max_full;
            max_full = 1.0/1.1;
            nroots++;   // ugly: redo the same class
            // when doing one big malloc, there's some chance that the
            // bucket overrun actually stepped over the next bucket. In
            // this case, the freeing of buckets in the code above might
            // have succeeded, so we can hope to resume with this special
            // q. On the other hand, if we have one malloc per bucket,
            // the free() calls above are guaranteed to crash.
            // Thus it's okay to proceed, if we're lucky enough to reach
            // here. Note that increasing bucket_limit will have a
            // permanent effect on the rest of this run.
            // abort();
            continue;
        }
#endif /* }}} */

        thrs[0]->rep->ttbuckets_fill += seconds();

        /* This can now be factored out ! */
        for(int side = 0 ; side < 2 ; side++) {
            sieve_side_info_ptr s = si->sides[side];

            small_sieve_init(s->ssd, las, s->fb, si, side);
            small_sieve_info("small sieve", side, s->ssd);

            small_sieve_extract_interval(s->rsd, s->ssd, s->fb_parts_x->rs);
            small_sieve_info("resieve", side, s->rsd);
        }

        /* FIXME: For the descent, the current logic is not
         * multithread-capable, as the two threads each try to produce a
         * winning relation. There should be a global counter with a
         * read-write lock, or something like that.
         */
        /* Process bucket regions in parallel */
        thread_do(thrs, &process_bucket_region, las->nb_threads);

        /* clear */
        for(int side = 0 ; side < 2 ; side++) {
            small_sieve_clear(si->sides[side]->ssd);
            small_sieve_clear(si->sides[side]->rsd);
        }
        qt0 = seconds() - qt0;
        las_report_accumulate_threads_and_display(las, si, report, thrs, qt0);

#ifdef TRACE_K
        trace_per_sq_clear(si);
#endif

      } // end of loop over special q ideals.

    thread_buckets_free(thrs, las->nb_threads);

    if (descent_lower) {
        verbose_output_print(0, 1, "# Now displaying again the results of all descents\n");
        las_descent_helper_display_all_trees(las->descent_helper, las->output);
        las_descent_helper_free(las->descent_helper);
    }

    t0 = seconds () - t0;
    wct = wct_seconds() - wct;
    verbose_output_print (2, 1, "# Average J=%1.0f for %lu special-q's, max bucket fill %f\n",
            totJ / (double) nr_sq_processed, nr_sq_processed, max_full);
    tts = t0;
    tts -= report->tn[0];
    tts -= report->tn[1];
    tts -= report->ttf;
    if (las->verbose)
        facul_print_stats (las->output);

    /*{{{ Display tally */
#ifdef HAVE_RUSAGE_THREAD
    int dont_print_tally = 0;
#else
    int dont_print_tally = 1;
#endif
    if (bucket_prime_stats) {
        verbose_output_print(2, 1, "# nr_bucket_primes = %lu, nr_div_tests = %lu, nr_composite_tests = %lu, nr_wrap_was_composite = %lu\n",
                 nr_bucket_primes, nr_div_tests, nr_composite_tests, nr_wrap_was_composite);
    }

    if (dont_print_tally && las->nb_threads > 1) 
        verbose_output_print (2, 1, "# Total cpu time %1.2fs [tally available only in mono-thread]\n", t0);
    else
        verbose_output_print (2, 1, "# Total cpu time %1.2fs [norm %1.2f+%1.1f, sieving %1.1f"
                " (%1.1f + %1.1f + %1.1f),"
                " factor %1.1f]\n", t0,
                report->tn[RATIONAL_SIDE],
                report->tn[ALGEBRAIC_SIDE],
                tts,
                report->ttbuckets_fill,
                report->ttbuckets_apply,
                tts-report->ttbuckets_fill-report->ttbuckets_apply,
                report->ttf);

    verbose_output_print (2, 1, "# Total elapsed time %1.2fs, per special-q %gs, per relation %gs\n",
                 wct, wct / (double) nr_sq_processed, wct / (double) report->reports);
    verbose_output_print (2, 1, "# Total %lu reports [%1.3gs/r, %1.1fr/sq]\n",
            report->reports, t0 / (double) report->reports,
            (double) report->reports / (double) nr_sq_processed);
    
    /*}}}*/

    thread_data_free(thrs, las->nb_threads);

    las_report_clear(report);

    las_info_clear(las);

    param_list_clear(pl);

    return 0;
}/*}}}*/
