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
#include "fb.h"
#include "portability.h"
#include "utils.h"           /* lots of stuff */
#include "ecm/facul.h"
#include "bucket.h"
#include "trialdiv.h"
#include "mpz_poly.h"
#include "las-config.h"
#include "las-types.h"
#include "las-coordinates.h"
#include "las-debug.h"
#include "las-report-stats.h"
#include "las-norms.h"
#include "las-unsieve.h"
#include "las-arith.h"
#include "las-qlattice.h"
#include "las-smallsieve.h"
#include "las-descent-helpers.h"
#ifdef HAVE_SSE41
/* #define SSE_SURVIVOR_SEARCH 1 */
#include <smmintrin.h>
#endif

static inline uint64_t cputicks()
{
        uint64_t r;
        __asm__ __volatile__(
                "rdtsc\n\t"
                "shlq $32, %%rdx\n\t"
                "orq %%rdx, %%rax\n\t"
                : "=a"(r)
                :
                : "%rdx", "cc");
        return r;
}

#define LOG_SCALE 1.4426950408889634 /* 1/log(2) to 17 digits, rounded to
                                        nearest. This is enough to uniquely
                                        identify the corresponding IEEE 754
                                        double precision number */
// #define HILIGHT_START   "\e[01;31m"
// #define HILIGHT_END   "\e[00;30m"

#define HILIGHT_START   ""
#define HILIGHT_END   ""

#ifndef HAVE_LOG2
static double MAYBE_UNUSED log2 (double x)
{
  return log (x) / log (2.0);
}
#endif

#ifndef HAVE_EXP2
static double MAYBE_UNUSED exp2 (double x)
{
  return exp (x * log (2.0));
}
#endif

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
double cof_calls[2][256] = {{0},{0}};
double cof_fails[2][256] = {{0},{0}};
/* }}} */

/* For memcpy in fill_in_k_buckets & fill_in_m_buckets.
   When you have to move data which the lenght is a
   static const uint8_t N <= 16 bytes,
   it's faster to move optimal_move[N] with a memcpy
   (if you could do this, of course) :
   it's done with only one or two instructions */
static const uint8_t optimal_move[] = { 0, 1, 2, 4, 4, 8, 8, 8, 8, 16, 16, 16, 16, 16, 16, 16, 16 };

/* Test if entry x in bucket region n is divisible by p */
void test_divisible_x (const fbprime_t p, const unsigned long x, const int n,
		       sieve_info_srcptr si, int side);
int factor_leftover_norm (mpz_t n,
                          mpz_array_t* const factors,
			  uint32_array_t* const multis,
                          sieve_info_ptr si, int side);
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
#if 0   /* disabled for the moment, as we haven't made the skewness part
         * of the siever_config. Should we ? */
    fprintf(o, "#                     skewness=%1.1f\n",
	    las->cpoly->skew);
#endif
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
    s->trialdiv_primes = fb_extract_bycost (s->fb, si->bucket_thresh,
            si->td_thresh);
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
/* {{{ silly code to factor small numbers. Should go to utils/ */
/* Finds prime factors p < lim of n and returns a pointer to a zero-terminated
   list of those factors. Repeated factors are stored only once. */
static fbprime_t *
factor_small (mpz_t n, fbprime_t lim)
{
  unsigned long p;
  unsigned long l; /* number of prime factors */
  fbprime_t *f;

  l = 0;
  f = (fbprime_t*) malloc (sizeof (fbprime_t));
  FATAL_ERROR_CHECK(f == NULL, "malloc failed");
  for (p = 2; p <= lim; p = getprime (p))
    {
      if (mpz_divisible_ui_p (n, p))
        {
          l ++;
          f = (fbprime_t*) realloc (f, (l + 1) * sizeof (fbprime_t));
          FATAL_ERROR_CHECK(f == NULL, "realloc failed");
          f[l - 1] = p;
        }
    }
  f[l] = 0; /* end of list marker */
  getprime (0);
  return f;
}
/*}}}*/

/* {{{ Initialize the factor bases */
void sieve_info_init_factor_bases(las_info_ptr las, sieve_info_ptr si, param_list pl)
{
    double tfb;
    /* TODO should these go into siever_config or not ? */
    int rpow_lim = 0, apow_lim = 0;
    param_list_parse_int(pl, "rpowlim", &rpow_lim);
    param_list_parse_int(pl, "apowlim", &apow_lim);
 
    for(int side = 0 ; side < 2 ; side++) {
        cado_poly_side_ptr pol = las->cpoly->pols[side];
        sieve_side_info_ptr sis = si->sides[side];
        unsigned long lim = si->conf->sides[side]->lim;
        if (pol->degree > 1) {
            fbprime_t *leading_div;
            tfb = seconds ();
            /* Do we actually use leading_div anywhere? */
            leading_div = factor_small (pol->f[pol->degree], lim);
            /* FIXME: fbfilename should allow *distinct* file names, of
             * course, for each side (think about the bi-algebraic case)
             */
            const char * fbfilename = param_list_lookup_string(pl, "fb");
            /* If apowlim is not given, or if it is too large, set it to
             * its maximum allowed value */
            if (apow_lim >= si->bucket_thresh) {
                apow_lim = si->bucket_thresh - 1;
                printf ("# apow_lim reduced to %d\n", apow_lim);
            }
            if (apow_lim == 0) 
                apow_lim = si->bucket_thresh - 1;
            fprintf(las->output, "# Reading %s factor base from %s\n", sidenames[side], fbfilename);
            int ok = fb_read_split (&sis->fb, &sis->fb_bucket_threads, fbfilename,
                                    sis->scale * LOG_SCALE, si->bucket_thresh,
                                    las->nb_threads, las->verbose, lim, apow_lim);
            FATAL_ERROR_CHECK(!ok, "Error reading factor base file");
            ASSERT_ALWAYS(sis->fb != NULL);
            tfb = seconds () - tfb;
            fprintf (las->output, 
                    "# Reading %s factor base of %zuMb took %1.1fs\n",
                    sidenames[side],
                    fb_size (sis->fb) >> 20, tfb);
            free (leading_div);
        } else {
            tfb = seconds ();
            if (rpow_lim >= si->bucket_thresh)
              {
                rpow_lim = si->bucket_thresh - 1;
                printf ("# rpow_lim reduced to %d\n", rpow_lim);
              }
            int ok = fb_make_linear (&sis->fb, &sis->fb_bucket_threads,
                                     (const mpz_t *) pol->f, (fbprime_t) lim,
                                     si->bucket_thresh, las->nb_threads,
                                     rpow_lim, sis->scale * LOG_SCALE,
                                     las->verbose, 1, las->output);
            FATAL_ERROR_CHECK(!ok, "Error creating rational factor base");
            tfb = seconds () - tfb;
            fprintf (las->output, "# Creating rational factor base of %zuMb took %1.1fs\n",
                     fb_size (sis->fb) >> 20, tfb);
        }
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

    fbprime_t plim = si->bucket_thresh;
    fbprime_t costlim = si->td_thresh;

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

    fprintf (las->output, "# Number of small-sieved primes in %s factor base = %zu\n", sidenames[side], fb_nroots_total(s->fb));

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
    fprintf (las->output, "# Number of bucket-sieved primes in %s factor base per thread =", sidenames[side]);
    for(int i = 0 ; i < n ; i++)
        fprintf (las->output, " %lu", nn[i]);
    fprintf(las->output, "\n");
    fprintf (las->output, "# Inverse sum of bucket-sieved primes in %s factor base per thread =", sidenames[side]);
    for(int i = 0 ; i < n ; i++)
        fprintf (las->output, " %.5f", bucket_fill_ratio[i]);

    double min_bucket_fill_ratio = bucket_fill_ratio[n-1];
    double max_bucket_fill_ratio = bucket_fill_ratio[0];
    for(int i = 0 ; i < n ; i++) {
        double r = bucket_fill_ratio[i];
        if (r < min_bucket_fill_ratio) min_bucket_fill_ratio = r;
        if (r > max_bucket_fill_ratio) max_bucket_fill_ratio = r;
    }
    fprintf(las->output, " [hit jitter %.2f%%]\n",
            100 * (max_bucket_fill_ratio / min_bucket_fill_ratio - 1));
    /* enable some margin in the bucket size */
    s->max_bucket_fill_ratio = max_bucket_fill_ratio * 1.05;
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

    si->bucket_thresh = si->I;	/* default value */
    /* overrides default only if parameter is given */
    param_list_parse_int(pl, "bkthresh", &(si->bucket_thresh));

    si->td_thresh = 1024;	/* default value */
    param_list_parse_uint(pl, "tdthresh", &(si->td_thresh));

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
    si->nb_buckets = 1 + ((si->J << si->conf->logI) - 1) / BUCKET_REGION;
    fprintf(las->output, "# bucket_region = %u\n", BUCKET_REGION);
    if (si->nb_buckets < THRESHOLD_K_BUCKETS)
      fprintf(las->output, "# nb_buckets = %u, one pass for the buckets sort\n", si->nb_buckets);
    else if (si->nb_buckets < THRESHOLD_M_BUCKETS)
      fprintf(las->output, "# nb_buckets = %u, two passes for the buckets sort\n", si->nb_buckets);
    else
      fprintf(las->output, "# nb_buckets = %u, three passes for the buckets sort\n", si->nb_buckets);

    sieve_info_init_unsieve_data(si);
    mpz_init(si->doing->p);
    mpz_init(si->doing->r);

    /* This (mildly) depends on the special q, but more importantly also
     * on the lpb/... bounds */
    sieve_info_init_norm_data(las->output, si, exp2(sc->bitsize), sc->side);

    /* TODO: This function in itself is too expensive if called often */
    sieve_info_init_factor_bases(las, si, pl);

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
        fprintf(las->output, "# Creating strategy for %d%c/%s [lim=%lu lpb=%u]\n",
                sc->bitsize, sidenames[sc->side][0], sidenames[s],
                sc->sides[s]->lim, sc->sides[s]->lpb);
        fprintf(las->output, "# Using %d+3 P-1/P+1/ECM curves\n",
                nb_curves (sc->sides[s]->lpb));
        si->sides[s]->strategy = facul_make_strategy(
                sc->sides[s]->lim, sc->sides[s]->lpb);
        reorder_fb(si, s);
        if (las->verbose) {
            fprintf(las->output, "# small %s factor base", sidenames[s]);
            factorbase_degn_t ** q;
            q = si->sides[s]->fb_parts->pow2;
            fprintf(las->output, ": %d pow2", fb_diff(q[1], q[0]));
            q = si->sides[s]->fb_parts->pow3;
            fprintf(las->output, ", %d pow3", fb_diff(q[1], q[0]));
            q = si->sides[s]->fb_parts->td;
            fprintf(las->output, ", %d td", fb_diff(q[1], q[0]));
            q = si->sides[s]->fb_parts->rs;
            fprintf(las->output, ", %d rs", fb_diff(q[1], q[0]));
            q = si->sides[s]->fb_parts->rest;
            fprintf(las->output, ", %d rest", fb_diff(q[1], q[0]));
            fprintf(las->output, " (total %zu)\n", fb_nroots_total(si->sides[s]->fb));
        }
    }
}
/* }}} */

/* This is an ownership transfer from the current head of the todo list
 * into si->doing.  The field todo->next is not accessed, since it does
 * not make sense for si->doing, which is only a single item. The head of
 * the todo list is pruned */
void sieve_info_pick_todo_item(sieve_info_ptr si, las_todo_ptr * todo)
{
    poly_t f;
    poly_alloc(f, si->cpoly->pols[(*todo)->side]->degree);
    poly_set(f,  si->cpoly->pols[(*todo)->side]->f, si->cpoly->pols[(*todo)->side]->degree);
    mpz_t x;
    mpz_init(x);
    poly_eval_mod_mpz(x, f, (*todo)->r, (*todo)->p);
    ASSERT_ALWAYS(mpz_cmp_ui(x, 0) == 0);
    mpz_clear(x);
    poly_free(f);
    mpz_clear(si->doing->p);
    mpz_clear(si->doing->r);
    memcpy(si->doing, *todo, sizeof(las_todo));
    free(*todo);
    *todo = si->doing->next;
    si->doing->next = 0;
    /* sanity check */
    if (!mpz_probab_prime_p(si->doing->p, 1)) {
        gmp_fprintf(stderr, "Error, %Zd is not prime\n", si->doing->p);
        exit(1);
    }
    ASSERT_ALWAYS(si->conf->side == si->doing->side);
}

/* return 0 if we should discard that special-q, in which case we intend
 * to discard this special-q. For this reason, si->J is then set to an
 * unrounded value, for diagnostic.
 *
 * The current check for discarding is whether we do fill one bucket or
 * not. If we don't even achieve that, we should of course discard.
 *
 * Now for efficiency reasons, the ``minimum reasonable'' number of
 * buckets should be more than that.
 */
int sieve_info_adjust_IJ(sieve_info_ptr si, double skewness, int nb_threads)/*{{{*/
{
    /* compare skewed max-norms: let u = [a0, b0] and v = [a1, b1],
       and u' = [a0, s*b0], v' = [a1, s*b1] where s is the skewness.
       Assume |u'| <= |v'|.
       We know from Gauss reduction that u' and v' form an angle of at
       least pi/3, thus since their determinant is q*s, we have
       |u'|^2 <= |u'| * |v'| <= q*s/sin(pi/3) = 2*q*s/sqrt(3).
       Define B := sqrt(2/sqrt(3)*q/s), then |a0| <= s*B and |b0| <= B.

       If we take J <= I/2*min(s*B/|a1|, B/|b1|), then for any j <= J:
       |a1|*J <= I/2*s*B and |b1|*J <= I/2*B, thus
       |a| = |a0*i+a1*j| <= s*B*I and |b| <= |b0*i+b1*j| <= B*I.
    */
    double maxab1, maxab0;
    maxab1 = si->b1 * skewness;
    maxab1 = maxab1 * maxab1 + si->a1 * si->a1;
    maxab0 = si->b0 * skewness;
    maxab0 = maxab0 * maxab0 + si->a0 * si->a0;
    if (maxab0 > maxab1) { /* exchange u and v, thus I and J */
        int64_t oa[2] = { si->a0, si->a1 };
        int64_t ob[2] = { si->b0, si->b1 };
        si->a0 = oa[1]; si->a1 = oa[0];
        si->b0 = ob[1]; si->b1 = ob[0];
    }
    maxab1 = MAX(fabs(si->a1), fabs(si->b1) * skewness);
    /* make sure J does not exceed I/2 */
    /* FIXME: We should not have to compute this B a second time. It
     * appears in sieve_info_init_norm_data already */
    double B = sqrt (2.0 * mpz_get_d(si->doing->p) / (skewness * sqrt (3.0)));
    if (maxab1 >= B * skewness)
        si->J = (uint32_t) (B * skewness / maxab1 * (double) (si->I >> 1));
    else
        si->J = si->I >> 1;

    /* Make sure the bucket region size divides the sieve region size, 
       partly covered bucket regions may lead to problems when 
       reconstructing p from half-empty buckets. */
    /* Compute number of i-lines per bucket region, must be integer */
    ASSERT_ALWAYS(LOG_BUCKET_REGION >= si->conf->logI);
    uint32_t i = 1U << (LOG_BUCKET_REGION - si->conf->logI);
    i *= nb_threads;  /* ensures nb of bucket regions divisible by nb_threads */

    /* Bug 15617: if we round up, we are not true to our promises */
    uint32_t nJ = (si->J / i) * i; /* Round down to multiple of i */

    /* XXX No rounding if we intend to abort */
    if (nJ > 0) si->J = nJ;
    return nJ > 0;
}/*}}}*/

static void sieve_info_update (sieve_info_ptr si, int nb_threads)/*{{{*/
{
  /* essentially update the fij polynomials */
  sieve_info_update_norm_data(si, nb_threads);

  /* update number of buckets */
  si->nb_buckets = 1 + ((si->J << si->conf->logI) - 1) / BUCKET_REGION;
}/*}}}*/

static void sieve_info_clear (las_info_ptr las, sieve_info_ptr si)/*{{{*/
{
    sieve_info_clear_unsieve_data(si);

    for(int s = 0 ; s < 2 ; s++) {
        facul_clear_strategy (si->sides[s]->strategy);
        si->sides[s]->strategy = NULL;
        sieve_info_clear_trialdiv(si, s);
        free(si->sides[s]->fb);
        for(int i = 0 ; i < las->nb_threads ; i++) {
            free(si->sides[s]->fb_bucket_threads[i]);
        }
        free(si->sides[s]->fb_bucket_threads);

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

    param_list_print_command_line(las->output, pl);
    las_display_config_flags(las->output);

    las->verbose = param_list_parse_switch(pl, "-v");
    las->nb_threads = 1;		/* default value */
    param_list_parse_int(pl, "mt", &las->nb_threads);
    if (las->nb_threads <= 0) {
	fprintf(stderr,
		"Error, please provide a positive number of threads\n");
	exit(EXIT_FAILURE);
    }
    /* }}} */
    /* {{{ Parse polynomial */
    cado_poly_init(las->cpoly);
    const char *tmp;
    if ((tmp = param_list_lookup_string(pl, "poly")) == NULL) {
        fprintf(stderr, "Error: -poly is missing\n");
        param_list_print_usage(pl, NULL, stderr);
        exit(EXIT_FAILURE);
    }

    if (!cado_poly_read(las->cpoly, tmp)) {
	fprintf(stderr, "Error reading polynomial file %s\n", tmp);
	exit(EXIT_FAILURE);
    }

    /* -skew (or -S) may override (or set) the skewness given in the
     * polynomial file */
    /* TODO should the skewness go into siever_config or not ? */
    param_list_parse_double(pl, "skew", &(las->cpoly->skew));

    if (las->cpoly->skew <= 0.0) {
	fprintf(stderr, "Error, please provide a positive skewness\n");
	exit(EXIT_FAILURE);
    }
    /* }}} */

    /* {{{ Parse default siever config (fill all possible fields) */

    siever_config_ptr sc = las->default_config;
    memset(sc, 0, sizeof(siever_config));

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
    int seen = 1;
    seen  = param_list_parse_int   (pl, "I",       &(sc->logI));
    seen &= param_list_parse_ulong (pl, "rlim",    &(sc->sides[RATIONAL_SIDE]->lim));
    seen &= param_list_parse_int   (pl, "lpbr",    &(sc->sides[RATIONAL_SIDE]->lpb));
    seen &= param_list_parse_int   (pl, "mfbr",    &(sc->sides[RATIONAL_SIDE]->mfb));
    seen &= param_list_parse_double(pl, "rlambda", &(sc->sides[RATIONAL_SIDE]->lambda));
    seen &= param_list_parse_ulong (pl, "alim",    &(sc->sides[ALGEBRAIC_SIDE]->lim));
    seen &= param_list_parse_int   (pl, "lpba",    &(sc->sides[ALGEBRAIC_SIDE]->lpb));
    seen &= param_list_parse_int   (pl, "mfba",    &(sc->sides[ALGEBRAIC_SIDE]->mfb));
    seen &= param_list_parse_double(pl, "alambda", &(sc->sides[ALGEBRAIC_SIDE]->lambda));
    if (!seen) {
        fprintf(stderr, "Error: options -I, -rlim, -lpbr, -mfbr, -rlambda,"
                " -alim, -lpba, -mfba, -alambda are mandatory.\n");
        exit(EXIT_FAILURE);
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
            exit(1);
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
    if (las->outputname)
        fclose_maybe_compressed(las->output, las->outputname);
    mpz_clear(las->todo_q0);
    mpz_clear(las->todo_q1);
    if (las->todo_list_fd)
        fclose(las->todo_list_fd);
    cado_poly_clear(las->cpoly);
}/*}}}*/

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
    fprintf(las->output, "# Creating new sieve configuration for q~2^%d on the %s side\n",
            sc->bitsize, sidenames[sc->side]);
    sieve_info_init_from_siever_config(las, si, sc, pl);
    memset(si + 1, 0, sizeof(sieve_info));
    siever_config_display(las->output, sc);
    fprintf(las->output, "#                     skewness=%1.1f\n",
            las->cpoly->skew);
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
                    gmp_fprintf(stderr, "Warning: fixing q=%Zd to next prime", q);
                    mpz_nextprime(q, q);
                    gmp_fprintf(stderr, " q=%Zd\n", q);
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
    int deg = MAX(las->cpoly->rat->degree, las->cpoly->alg->degree);
    roots = (mpz_t *) malloc (deg * sizeof (mpz_t));
    for(int i = 0 ; i < deg  ; i++) {
        mpz_init(roots[i]);
    }

    las_todo_ptr * pnext = &(las->todo);
    for( ; pushed < 10 && mpz_cmp(q, q1) < 0 ; ) {
        mpz_nextprime(q, q);
        if (mpz_cmp(q, q1) >= 0)
            break;
        cado_poly_side_ptr f = las->cpoly->pols[qside];
        int nroots = poly_roots(roots, f->f, f->degree, q);

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
                   rc = gmp_sscanf(x, "%Zi%n", p, &nread);
                   x+=nread;
                   ASSERT_ALWAYS(rc == 1);  /* %n does not count */
                   ASSERT_ALWAYS(mpz_probab_prime_p(p, 2));
                   cado_poly_side_ptr f = las->cpoly->pols[RATIONAL_SIDE];
                   ASSERT_ALWAYS(f->degree == 1);
                   int nroots = poly_roots(&r, f->f, f->degree, p);
                   ASSERT_ALWAYS(nroots == 1);
                   /* We may in fact also have the root specified. We
                    * ignore what is in the file, the root is computed
                    * unconditionally above */
                   gmp_sscanf(x, "%*Zi%n", &nread);
                   x+=nread;
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

/*******************************************************************/
/********        Walking in p-lattices                    **********/

/* {{{ p-lattice stuff */


// Compute the root r describing the lattice inside the q-lattice
// corresponding to the factor base prime (p,R).
// Formula: r = - (a1-R*b1)/(a0-R*b0) mod p
// Assumes p < 2^32

/* General version of the lattice transform function. Allows projective
   roots in input and output, and handles prime powers.
   In input, if the root is projective, say s/t (mod p) and t is
   non-invertible (mod p), then we expect R = p + (t/s mod p).
   On output, if the root is projective, say u/v (mod p) and v is
   non-invertible (mod p), then return value r = p + (v/u mod p).
   So projective roots are stored as their reciprocal, and have p added
   to signal the fact that it's a reciprocal value.
*/


/*
 * Algorithm by Franke and Kleinjung for lattice sieving of largish
 * primes.
 */

/* DONT modify this: asm code writes in this with hardcoded deplacement */
typedef struct {
  int32_t a0; uint32_t a1;
  int32_t b0; uint32_t b1;
} plattice_info_t;

// Proposition 1 of [FrKl05]:
// Compute a basis <(alpha, beta), (gamma, delta)> of the p-lattice
// inside the q-lattice, such that
//    beta, delta > 0
//    -I < alpha <= 0 <= gamma < I
//    gamma-alpha >= I
//
// Sizes:
//    p is less than 32 bits and I fits easily in 32 bits.
//    So, alpha and beta fit easily in 32 bits, since they are less than I
//    Now, gamma and delta are also bounded by p, so 32 bits is enough
//    However: a and c can be as large as p*I (not both ?).
//    We still store them in 32 bits, since if they are larger, it means
//    that as soon as they are added to the offset for S, the index will
//    be out of range for S and the loop stops. Hence, this is safe to
//    replace a and c by a large value within 32 bits, when they are
//    larger than 32 bits.
//    Except that the choice of this ``large value'' requires some
//    caution. We need a value which can be used either for a or c, or
//    both, so that adding a, or c, or both, to a value within [0,IJ[ is
//    guaranteed to exceed IJ, but can't wrap around. Up to I=15, it's
//    rather easy. With the rescaling of J, at worst we may have IJ
//    within the range [2^29,2^30[. Thus if a and c are set to 2^30-1,
//    a.k.a INT32_MAX/2, then adding either, or the sum of both, to a
//    valid value x is guaranteed to be at most 3*2^30, which fits within
//    32 bits.
//    For I=16, it's much, much harder. Unless somebody comes up with a
//    nice idea, I see no way to avoid 64-bit arithmetic (which has some
//    cost, sure, but not terribly expensive). For consistency, we make
//    all data types for x, a, and c 64-bit in this case, and choose the
//    overflow constants as UINT32_MAX.
//

// OBSOLETE: the large value is now UINT64MAX>>2. So, plattice_x_t must
// be 64 bits. We could deal until logI <= 31 when logJ = logI -1.

typedef uint64_t plattice_x_t;

// Return value:
// * non-zero if everything worked ok
// * zero when the algorithm failed. This can happen when p is a prime power,
//   and g, gcd(p,r) >= I, since then the subtractive Euclidean algorithm will
//   yield (a0=g, b0=0) at some point --- or the converse --- and the loop
//   while (|a0| >= I) a0 += b0 will loop forever.
//
// Note that on a c166 example, this code alone accounts for almost 20%
// of the computation time.


// 4 reduce_plattice optimized versions.
// - First is the optimized original & use the subtractive Euclidean algorithm
//   to avoid the 32 bits integer division. It's the fastest for old processors
//   and AMD processors.
// - Second is the optimized version with the 32 bits integer division.
//   It's the fastest for Intel beginning from the Nehalem processor.
// - Third is the version with double. No interest for X86 (use the next); maybe
//   for others processors ?
// - Last is the SSE(2) version with double. In Nehalem, it's the slowest version
//   except the previous, of course. Maybe in the Intel lastest procesors ?

#define reduce_plattice reduce_plattice1

NOPROFILE_INLINE int
reduce_plattice0 (plattice_info_t *pli, const fbprime_t p, const fbprime_t r, sieve_info_srcptr si)
{
  int32_t a0 = - (int32_t) p, a1 = 0, b0 = (int32_t) r, b1 = 1, k;
#if MOD2_CLASSES_BS
  const int32_t hI = (int32_t) ((si->I) >> 1);
#else
  const int32_t hI = (int32_t) (si->I);
#endif
  const int32_t mhI = -hI;
    /* subtractive variant of Euclid's algorithm */
    for(;;) {
        /* a0 < 0 <= b0 < -a0 */
        if (b0 < hI) break;
        /* a0 < 0 < b0 < -a0 */
        /* for( ; a0 += b0, a1 += b1, a0 + b0 <= 0 ; ); */
	for (;;) {
	  a0 += b0; a1 += b1; if (a0 + b0 > 0) break;
	  a0 += b0; a1 += b1; if (a0 + b0 > 0) break;
	  a0 += b0; a1 += b1; if (a0 + b0 > 0) break;
	  a0 += b0; a1 += b1; if (a0 + b0 > 0) break;
	  a0 += b0; a1 += b1; if (a0 + b0 > 0) break;
	  a0 += b0; a1 += b1; if (a0 + b0 > 0) break;
	  a0 += b0; a1 += b1; if (a0 + b0 > 0) break;
	  a0 += b0; a1 += b1; if (a0 + b0 > 0) break;
	  a0 += b0; a1 += b1; if (a0 + b0 > 0) break;
	  a0 += b0; a1 += b1; if (a0 + b0 > 0) break;
	  a0 += b0; a1 += b1; if (a0 + b0 > 0) break;
	  a0 += b0; a1 += b1; if (a0 + b0 > 0) break;
	  a0 += b0; a1 += b1; if (a0 + b0 > 0) break;
	  a0 += b0; a1 += b1; if (a0 + b0 > 0) break;
	  a0 += b0; a1 += b1; if (a0 + b0 > 0) break;
	  a0 += b0; a1 += b1; if (a0 + b0 > 0) break;
	}
        /* -b0 < a0 <= 0 < b0 */
        if (a0 > mhI) break;
        /* -b0 < a0 < 0 < b0 */
        /* for( ; b0 += a0, b1 += a1, b0 + a0 >= 0 ; ); */
	for (;;) {
	  b0 += a0; b1 += a1; if (b0 + a0 < 0) break;
	  b0 += a0; b1 += a1; if (b0 + a0 < 0) break;
	  b0 += a0; b1 += a1; if (b0 + a0 < 0) break;
	  b0 += a0; b1 += a1; if (b0 + a0 < 0) break;
	  b0 += a0; b1 += a1; if (b0 + a0 < 0) break;
	  b0 += a0; b1 += a1; if (b0 + a0 < 0) break;
	  b0 += a0; b1 += a1; if (b0 + a0 < 0) break;
	  b0 += a0; b1 += a1; if (b0 + a0 < 0) break;
	  b0 += a0; b1 += a1; if (b0 + a0 < 0) break;
	  b0 += a0; b1 += a1; if (b0 + a0 < 0) break;
	  b0 += a0; b1 += a1; if (b0 + a0 < 0) break;
	  b0 += a0; b1 += a1; if (b0 + a0 < 0) break;
	  b0 += a0; b1 += a1; if (b0 + a0 < 0) break;
	  b0 += a0; b1 += a1; if (b0 + a0 < 0) break;
	  b0 += a0; b1 += a1; if (b0 + a0 < 0) break;
	  b0 += a0; b1 += a1; if (b0 + a0 < 0) break;
	}
        /* a0 < 0 <= b0 < -a0 */
    }
    k = b0 - hI - a0;
    if (b0 > -a0) {
        if (UNLIKELY(!a0)) return 0;
        /* Now that |a0| < hI, we switch to classical division, since
           if say |a0|=1 and b0 is large, the subtractive variant
           will be very expensive.
           We want b0 + k*a0 < hI, i.e., b0 - hI + 1 <= k*(-a0),
           i.e., k = ceil((b0-hI+1)/a0). */
	k /= a0; b0 -= k * a0; b1 -= k * a1;
    } else {
        if (UNLIKELY(!b0)) return 0;
        /* we switch to the classical algorithm here too */
	k /= b0; a0 += k * b0; a1 += k * b1;
    }
    pli->a0 = a0; pli->a1 = a1; pli->b0 = b0; pli->b1 = b1;
#if 0
    /* Too expensive checks */
    /*
    ASSERT (a1 > 0);
    ASSERT (b1 > 0);
    ASSERT ((a0 <= 0) && (a0 > -hI));
    ASSERT ((b0 >= 0) && (b0 <  hI));
    ASSERT (b0 - a0 >= hI);
    */
    int32_t J = si->J;
#if MOD2_CLASSES_BS
#if 1
#endif

#if 1 || !MOD2_CLASSES
    // left-shift everybody, since the following correspond to the
    // lattice 2p.
    a0 <<= 1; a1 <<= 1;
    b0 <<= 1; b1 <<= 1;
#endif
#endif

    pli->a = ((a1 << si->conf->logI) + a0);
    pli->c = ((b1 << si->conf->logI) + b0);
    if (a1 > J || (a1 == J && a0 > 0)) { pli->a = INT32_MAX/2; }
    if (b1 > J || (b1 == J && b0 > 0)) { pli->c = INT32_MAX/2; }

    /* It's difficult to encode everybody in 32 bits, and still keep
     * relevant information...
     */
    pli->bound0 = -a0;
    pli->bound1 = si->I - b0;
#endif
    return 1;
}

NOPROFILE_INLINE int
reduce_plattice1 (plattice_info_t *pli, const fbprime_t p, const fbprime_t r, sieve_info_srcptr si)
{
  int32_t a1 = 0, a0 = -((int32_t) p), b1 = 1, b0 = (int32_t) r, k;
#if MOD2_CLASSES_BS
  const int32_t hI = (int32_t) ((si->I) >> 1);
#else
  const int32_t hI = (int32_t) (si->I);
#endif
  const int32_t mhI = -hI;
  for(;;) {
    if (b0 < hI ) break; k = a0 / b0; a0 %= b0; a1 -= k * b1;
    if (a0 > mhI) break; k = b0 / a0; b0 %= a0; b1 -= k * a1;
    if (b0 < hI ) break; k = a0 / b0; a0 %= b0; a1 -= k * b1;
    if (a0 > mhI) break; k = b0 / a0; b0 %= a0; b1 -= k * a1;
    if (b0 < hI ) break; k = a0 / b0; a0 %= b0; a1 -= k * b1;
    if (a0 > mhI) break; k = b0 / a0; b0 %= a0; b1 -= k * a1;
    if (b0 < hI ) break; k = a0 / b0; a0 %= b0; a1 -= k * b1;
    if (a0 > mhI) break; k = b0 / a0; b0 %= a0; b1 -= k * a1;
    if (b0 < hI ) break; k = a0 / b0; a0 %= b0; a1 -= k * b1;
    if (a0 > mhI) break; k = b0 / a0; b0 %= a0; b1 -= k * a1;
    if (b0 < hI ) break; k = a0 / b0; a0 %= b0; a1 -= k * b1;
    if (a0 > mhI) break; k = b0 / a0; b0 %= a0; b1 -= k * a1;
    if (b0 < hI ) break; k = a0 / b0; a0 %= b0; a1 -= k * b1;
    if (a0 > mhI) break; k = b0 / a0; b0 %= a0; b1 -= k * a1;
    if (b0 < hI ) break; k = a0 / b0; a0 %= b0; a1 -= k * b1;
    if (a0 > mhI) break; k = b0 / a0; b0 %= a0; b1 -= k * a1;
  }
  k = b0 - hI - a0;
  if (b0 > -a0) {
    if (UNLIKELY(!a0)) return 0;
    k /= a0; b0 -= k * a0; b1 -= k * a1;
  } else {
    if (UNLIKELY(!b0)) return 0;
    k /= b0; a0 += k * b0; a1 += k * b1;
  }
  pli->a0 = (int32_t) a0; pli->a1 = (uint32_t) a1; pli->b0 = (int32_t) b0; pli->b1 = (uint32_t) b1;
  return 1;
}

NOPROFILE_INLINE int
reduce_plattice2 (plattice_info_t *pli, const fbprime_t p, const fbprime_t r, sieve_info_srcptr si)
{
  double a1 = 0., a0 = (double) -((int) p), b1 = 1., b0 = (double) r, k;
#if MOD2_CLASSES_BS
  const double hI = (double) ((si->I) >> 1);
#else
  const double hI = (double) (si->I);
#endif
  const double mhI = -hI;
  for(;;) {
    if (b0 < hI ) break; k = trunc(a0 / b0); a0 -= k * b0; a1 -= k * b1;
    if (a0 > mhI) break; k = trunc(b0 / a0); b0 -= k * a0; b1 -= k * a1;
    if (b0 < hI ) break; k = trunc(a0 / b0); a0 -= k * b0; a1 -= k * b1;
    if (a0 > mhI) break; k = trunc(b0 / a0); b0 -= k * a0; b1 -= k * a1;
    if (b0 < hI ) break; k = trunc(a0 / b0); a0 -= k * b0; a1 -= k * b1;
    if (a0 > mhI) break; k = trunc(b0 / a0); b0 -= k * a0; b1 -= k * a1;
    if (b0 < hI ) break; k = trunc(a0 / b0); a0 -= k * b0; a1 -= k * b1;
    if (a0 > mhI) break; k = trunc(b0 / a0); b0 -= k * a0; b1 -= k * a1;
    if (b0 < hI ) break; k = trunc(a0 / b0); a0 -= k * b0; a1 -= k * b1;
    if (a0 > mhI) break; k = trunc(b0 / a0); b0 -= k * a0; b1 -= k * a1;
    if (b0 < hI ) break; k = trunc(a0 / b0); a0 -= k * b0; a1 -= k * b1;
    if (a0 > mhI) break; k = trunc(b0 / a0); b0 -= k * a0; b1 -= k * a1;
    if (b0 < hI ) break; k = trunc(a0 / b0); a0 -= k * b0; a1 -= k * b1;
    if (a0 > mhI) break; k = trunc(b0 / a0); b0 -= k * a0; b1 -= k * a1;
    if (b0 < hI ) break; k = trunc(a0 / b0); a0 -= k * b0; a1 -= k * b1;
    if (a0 > mhI) break; k = trunc(b0 / a0); b0 -= k * a0; b1 -= k * a1;
  }
  k = b0 - hI - a0;
  if (b0 > -a0) {
    if (UNLIKELY(a0 == 0.)) return 0;
    k = trunc(k / a0); b0 -= k * a0; b1 -= k * a1;
  } else {
    if (UNLIKELY(b0 == 0.)) return 0;
    k = trunc(k / b0); a0 += k * b0; a1 += k * b1;
  }
  pli->a0 = (int32_t) a0; pli->a1 = (uint32_t) a1; pli->b0 = (int32_t) b0; pli->b1 = (uint32_t) b1;
  return 1;
}

/* I have tried to use the intrinsics only. But it's too boring, the code
   is too long & with many tests + loop. So, ASM x86 is your friend.
   But it's obvious the initialisation and the end must be done in intrinsics
   to avoid a waste of general registers (gcc could use the stack if possible).
   NB: To avoid a stupid gcc false warning, I have to clear hI for nothing. Grrr. */
NOPROFILE_INLINE int
reduce_plattice3 (plattice_info_t *pli, const fbprime_t p, const fbprime_t r, sieve_info_srcptr si)
{
  int ret;
  static const double one = 1.0;
  __m128d a, b, hI, mhI; /* args 1 to 4; ret is the 0th arg */
  a = b = hI = mhI = _mm_setzero_pd();
  b = _mm_cvtsi32_sd(b, p); a = _mm_sub_sd (a, b);     /* a   = 0     -p */
  b = _mm_loadh_pd(b, &one); b = _mm_cvtsi32_sd(b, r); /* b   = 1      r */ 
#if MOD2_CLASSES_BS
  hI = _mm_cvtsi32_sd(hI, (si->I) >> 1);               /* hI  = 0  (si->I)>>1 */
#else
  hI = _mm_cvtsi32_sd(hI, si->I);                      /* hI  = 0  si->I */
#endif
  mhI = _mm_sub_sd (mhI, hI);                          /* mhI = -hI */
  __asm__ ( "xor %0, %0\n"
    ".balign 8\n 0:\n"          /* Main loop */

#define INLOOP								\
    "ucomisd %3, %2\n"							\
      "movsd %1, %%xmm0\n"						\
      "jb 1f\n"                   /* if (b0 < hI ) break; */		\
      "divsd %2, %%xmm0\n"						\
      "unpcklpd %%xmm0, %%xmm0\n" /* a0/b0 a0/b0 */			\
      "cvttpd2dq %%xmm0, %%xmm0\n"					\
      "cvtdq2pd %%xmm0, %%xmm0\n" /* trunc(a0/b0) trunc(a0/b0) */	\
      "mulpd %2, %%xmm0\n"						\
      "subpd %%xmm0, %1\n"        /* a -= b * trunc(a0/b0) */		\
									\
      "ucomisd %1, %4\n"						\
      "movsd %2, %%xmm0\n"						\
      "jb 1f\n"                   /* if (mhI < a0) break; */		\
      "divsd %1, %%xmm0\n"						\
      "unpcklpd %%xmm0, %%xmm0\n" /* b0/a0 b0/a0 */			\
      "cvttpd2dq %%xmm0, %%xmm0\n"					\
      "cvtdq2pd %%xmm0, %%xmm0\n" /* trunc(b0/a0) trunc(b0/a0) */	\
      "mulpd %1, %%xmm0\n"						\
      "subpd %%xmm0, %2\n"        /* b -= a * trunc(b0/a0) */		\
    
    INLOOP INLOOP INLOOP INLOOP INLOOP INLOOP INLOOP INLOOP
#undef INLOOP 
    "jmp 0b\n"
    
    "1:\n"
    "xorpd %%xmm0, %%xmm0\n"
    "subsd %1, %%xmm0\n"        /* 0 -a0 */
    "ucomisd %%xmm0, %2\n"  
    "xorpd %%xmm1, %%xmm1\n"
    "addsd %2, %%xmm0\n"    
    "jb 2f\n"                   /* if (-a0 < b0) goto 2nd part*/
    
    "ucomisd %%xmm1, %1\n"      /* part 1 (then) : -a0 >= b0 */
    "subsd %3, %%xmm0\n"        /* 0 b0-a0-hI */
    "je 4f\n"                   /* if (a0 == 0) goto 4 (error) */
    
    "divsd %1, %%xmm0\n"
    "cvttpd2dq %%xmm0, %%xmm0\n"
    "cvtdq2pd %%xmm0, %%xmm0\n"
    "unpcklpd %%xmm0, %%xmm0\n" /* trunc((b0-hI-a0)/a0) trunc((b0-hI-a0)/a0) */
    "mulpd %1, %%xmm0\n"
    "subpd %%xmm0, %2\n"        /* b -= a * trunc((b0-hI-a0)/a0) */
    "jmp 3f\n"
    
    "2:\n"           
    "ucomisd %%xmm1, %2\n"      /* part 2 (else) : -a0 < b0 */
    "subsd %3, %%xmm0\n"        /* 0 b0-a0-hI */
    "je 4f\n"                   /* if (b0 == 0) goto return 1 (error) */
    
    "divsd %2, %%xmm0\n"
    "cvttpd2dq %%xmm0, %%xmm0\n"
    "cvtdq2pd %%xmm0, %%xmm0\n"
    "unpcklpd %%xmm0, %%xmm0\n" /* trunc((b0-hI-a0)/b0) trunc((b0-hI-a0)/b0) */
    "mulpd %2, %%xmm0\n"
    "addpd %%xmm0, %1\n"        /* a += b * trunc((b0-hI-a0)/b0) */
    
    "3:\n"
    "cvttpd2dq %1, %1\n"
    "cvttpd2dq %2, %2\n"
    "unpcklpd %2, %1\n"         /* final result : a = b1 b0 a1 a0 in (u)int32_t * 4 */
    "inc %0\n"
    
    "4:\n"
	    
    : "=&r"(ret), "+&x"(a), "+&x"(b) : "x"(hI), "x"(mhI) : "%xmm0", "%xmm1", "cc");
  _mm_storeu_pd((double *) pli, a);
  return(ret);
}

#if MOD2_CLASSES_BS
#define PLI_COEFF(pli, ab01) (pli->ab01 << 1)
#else
#define PLI_COEFF(pli, ab01) (pli->ab01)
#endif
static inline plattice_x_t plattice_a(const plattice_info_t * pli, sieve_info_srcptr si)
{
    int64_t a0 = PLI_COEFF(pli, a0);
    uint64_t a1 = PLI_COEFF(pli, a1);
    return (a1 << si->conf->logI) + a0;
}

static inline plattice_x_t plattice_c(const plattice_info_t * pli, sieve_info_srcptr si)
{
    int64_t b0 = PLI_COEFF(pli, b0);
    uint64_t b1 = PLI_COEFF(pli, b1);
    return (b1 << si->conf->logI) + b0;
}

static inline uint32_t plattice_bound0(const plattice_info_t * pli, sieve_info_srcptr si MAYBE_UNUSED)
{
    return - PLI_COEFF(pli, a0);
}

static inline uint32_t plattice_bound1(const plattice_info_t * pli, sieve_info_srcptr si)
{
    return si->I - PLI_COEFF(pli, b0);
}


/* This is for working with congruence classes only */
NOPROFILE_INLINE
plattice_x_t plattice_starting_vector(const plattice_info_t * pli, sieve_info_srcptr si, unsigned int par )
{
    /* With MOD2_CLASSES_BS set up, we have computed by the function
     * above an adapted basis for the band of total width I/2 (thus from
     * -I/4 to I/4). This adapted basis is in the coefficients a0 a1 b0
     *  b1 of the pli data structure.
     *
     * Now as per Proposition 1 of FrKl05 applied to I/2, any vector
     * whose i-coordinates are within ]-I/2,I/2[ (<ugly>We would like a
     * closed interval on the left. Read further on for that case</ugly>)
     * can actually be written as a combination with positive integer
     * coefficients of these basis vectors a and b.
     *
     * We also know that the basis (a,b) has determinant p, thus odd. The
     * congruence class mod 2 that we want to reach is thus accessible.
     * It is even possible to find coefficients (k,l) in {0,1} such that
     * ka+lb is within this congruence class. This means that we're going
     * to consider either a,b,or a+b as a starting vector. The
     * i-coordinates of these, as per the properties of Proposition 1, are
     * within ]-I/2,I/2[. Now all other vectors with i-coordinates in
     * ]-I/2,I/2[ which also belong to the same congruence class, can be
     * written as (2k'+k)a+(2l'+l)b, with k' and l' necessarily
     * nonnegative.
     *
     * The last ingredient is that (2a, 2b) forms an adapted basis for
     * the band of width I with respect to the lattice 2p. It's just an
     * homothety.
     *
     * To find (k,l), we proceed like this. First look at the (a,b)
     * matrix mod 2:
     *                 a0&1    a1&1
     *                 b0&1    b1&1
     * Its determinant is odd, thus the inverse mod 2 is:
     *                 b1&1    a1&1
     *                 b0&1    a0&1
     * Now the congruence class is given by the parity argument. The
     * vector is:
     *                par&1,   par>>1
     * Multiplying this vector by the inverse matrix above, we obtain the
     * coordinates k,l, which are:
     *            k = (b1&par&1)^(b0&(par>>1));
     *            l = (a1&par&1)^(a0&(par>>1));
     * Now our starting vector is ka+lb. Instead of multiplying by k and
     * l with values in {0,1}, we mask with -k and -l, which both are
     * either all zeroes or all ones in binary
     *
     */
    /* Now for the extra nightmare. Above, we do not have the guarantee
     * that a vector whose i-coordinate is precisely -I/2 has positive
     * coefficients in our favorite basis. It's annoying, because it may
     * well be that if such a vector also has positive j-coordinate, then
     * it is in fact the first vector we will meet. An example is given
     * by the following data:
     *
            f:=Polynomial(StringToIntegerSequence("
                -1286837891385482936880099433527136908899552
                55685111236629140623629786639929578
                13214494134209131997428776417
                -319664171270205889372
                -17633182261156
                40500"));

            q:=165017009; rho:=112690811;
            a0:=52326198; b0:=-1; a1:=60364613; b1:=2;
            lI:=13; I:=8192; J:=5088;
            p:=75583; r0:=54375;
            > M;
            [-2241    19]
            [ 1855    18]
            > M[1]-M[2];
            (-4096     1)

    * Clearly, above, for the congruence class (0,1), we must start with
    * this vector, not with the sum.
    */
    int32_t a0 = pli->a0;
    int32_t a1 = pli->a1;
    int32_t b0 = pli->b0;
    int32_t b1 = pli->b1;

    int k = -((b1&par&1)^(b0&(par>>1)));
    int l = -((a1&par&1)^(a0&(par>>1)));
    int32_t v[2]= { (a0&k)+(b0&l), (a1&k)+(b1&l)};

    /* handle exceptional case as described above */
    if (k && l && a0-b0 == -(1 << (si->conf->logI-1)) && a1 > b1) {
        v[0] = a0-b0;
        v[1] = a1-b1;
    }
    return (v[1] << si->conf->logI) | (v[0] + (1 << (si->conf->logI-1)));
}

/* }}} */


/***************************************************************************/
/********        Main bucket sieving functions                    **********/

/* {{{ thread-related defines */
/* All of this exists _for each thread_ */
struct thread_side_data_s {
  m_bucket_array_t mBA; /* Not used if not fill_in_m_buckets (3 passes sort) */
  k_bucket_array_t kBA; /* Ditto for fill_in_k_buckets (2 passes sort) */
  bucket_array_t BA;    /* Always used */
  factorbase_degn_t *fb_bucket; /* copied from sieve_info. Keep ? XXX */
  // double bucket_fill_ratio;     /* inverse sum of bucket-sieved primes */
  
  /* For small sieve */
  int * ssdpos;
  int * rsdpos;
};
typedef struct thread_side_data_s thread_side_data[1];
typedef struct thread_side_data_s * thread_side_data_ptr;
typedef const struct thread_side_data_s * thread_side_data_srcptr;

struct thread_data_s {
  int id;
  thread_side_data sides[2];
  las_info_ptr las;
  sieve_info_ptr si;
  las_report rep;
};
typedef struct thread_data_s thread_data[1];
typedef struct thread_data_s * thread_data_ptr;
typedef const struct thread_data_s * thread_data_srcptr;
/* }}} */

/**************************************************************************
 * Global DEFINEs for fill_in_buckets, fill_in_k_buckets, fill_in_m_buckets 
 **************************************************************************/
#ifdef SKIP_GCD3
#define FILL_BUCKET_SKIP_GCD3()				\
  && (!is_divisible_3_u32 (i + I) ||			\
      !is_divisible_3_u32 ((uint32_t) (x >> logI)))			
#else
#define FILL_BUCKET_SKIP_GCD3()
#endif

#ifdef TRACE_K								
#define FILL_BUCKET_TRACE_K(X) do {					\
    if (trace_on_spot_x(X)) {						\
      WHERE_AM_I_UPDATE(w, N, (X) >> 16);				\
      WHERE_AM_I_UPDATE(w, x, (uint16_t) (X));				\
      fprintf (stderr, "# Pushed (%u, %u) (%u, %s) to BA[%u]\n",	\
	       (unsigned int) (uint16_t) (X), logp, p, sidenames[side],	\
	       (unsigned int) ((X) >> 16));				\
      ASSERT(test_divisible(w));					\
    }									\
  } while(0)
#else
#define FILL_BUCKET_TRACE_K(X)
#endif

#ifdef HAVE_SSE2							
#define FILL_BUCKET_PREFETCH(PT) do {				\
    _mm_prefetch((char *)(PT), _MM_HINT_T0);			\
  } while (0)
#else
#define FILL_BUCKET_PREFETCH(PT)
#endif

#define FILL_BUCKET_INC_X() do {					\
    if (i >= bound1) x += inc_a;					\
    if (i < bound0)  x += inc_c;					\
  } while (0)
/************************************************************************/

/* {{{ */
void
fill_in_buckets(thread_data_ptr th, int side, where_am_I_ptr w MAYBE_UNUSED)
{
  WHERE_AM_I_UPDATE(w, side, side);
  sieve_info_srcptr si = th->si;
  bucket_array_t BA = th->sides[side]->BA;  /* local copy. Gain a register + use stack */
  // Loop over all primes in the factor base.
  //
  // Note that dispatch_fb already arranged so that all the primes
  // which appear here are >= bucket_thresh and <= pmax (the latter
  // being for the moment unconditionally set to FBPRIME_MAX by the
  // caller of dispatch_fb).
  
  fb_iterator t;
  fb_iterator_init_set_fb(t, th->sides[side]->fb_bucket);
  unsigned char last_logp = 0;
  for( ; !fb_iterator_over(t) ; fb_iterator_next(t)) {
    fbprime_t p = t->fb->p;
    unsigned char logp = t->fb->plog;
    ASSERT_ALWAYS (p & 1);
    WHERE_AM_I_UPDATE(w, p, p);
    /* Write new set of pointers if the logp value changed */
    if (UNLIKELY(last_logp != logp)) {
      aligned_medium_memcpy((uint8_t *)BA.logp_idx + BA.size_b_align * BA.nr_logp, BA.bucket_write, BA.size_b_align);
      BA.logp_val[BA.nr_logp++] = last_logp = logp;
    }
    /* If we sieve for special-q's smaller than the factor
       base bound, the prime p might equal the special-q prime q. */
    if (UNLIKELY(!mpz_cmp_ui(si->doing->p, p))) continue;
    fbprime_t R = fb_iterator_get_r(t), r = fb_root_in_qlattice(p, R, t->fb->invp, si);
    
#ifdef SKIP_GCD3
    const uint32_t I = si->I;
    const unsigned int logI = si->conf->logI;
#endif
    const uint32_t maskI = si->I-1;
    const uint64_t even_mask = (1ULL << si->conf->logI) | 1ULL;
    const uint64_t IJ = ((uint64_t) si->J) << si->conf->logI;

    /* Special cases */
    if (UNLIKELY((!r) || (r >= p))) {
      if (r > p) /* should only happen for lattice-sieved prime powers,
		    which is not possible currently since maxbits < I */
	continue;
      /* r == p or r == 0.
	 1. If r == 0 (mod p), this prime hits for i == 0 (mod p), but since p > I,
	 this implies i = 0 or i > I. We don't sieve i > I. Since gcd(i,j) |
	 gcd(a,b), for i = 0 we only need to sieve j = 1. 
	 So, x = j*I + (i + I/2) = I + I/2.
	 2. r == p means root at infinity, which hits for j == 0 (mod p). Since q > I > J,
	 this implies j = 0 or j > J. This means we sieve only (i,j) = (1,0) here.
	 FIXME: what about (-1,0)? It's the same (a,b) as (1,0) but which of these two
	 (if any) do we sieve? */
      uint64_t x = (r ? 1 : si->I) + (si->I >> 1);
      prime_hint_t prime = bucket_encode_prime (p);
      /*****************************************************************/
#define FILL_BUCKET_HEART() do {					\
	bucket_update_t **pbut = BA.bucket_write + (x >> 16);		\
	bucket_update_t *but = *pbut;					\
	FILL_BUCKET_TRACE_K(x);						\
	WHERE_AM_I_UPDATE(w, N, x >> 16);				\
	WHERE_AM_I_UPDATE(w, x, (uint16_t) x);				\
	but->p = prime;							\
	but->x = (uint16_t) x;						\
	*pbut = ++but;							\
	FILL_BUCKET_PREFETCH(but);					\
      } while (0)
      /*****************************************************************/
      FILL_BUCKET_HEART();
      continue;
    }
    /* If working with congruence classes, once the loop on the parity goes at the level
       above, this initialization should in fact either be done for each congruence class,
       or saved for later use within the factor base structure. */
    plattice_info_t pli;
    if (UNLIKELY(!reduce_plattice(&pli, p, r, si))) {
      pthread_mutex_lock(&io_mutex);
      fprintf (stderr, "# fill_in_buckets: reduce_plattice() returned 0 for p = "
	       FBPRIME_FORMAT ", r = " FBPRIME_FORMAT "\n", p, r);
	pthread_mutex_unlock(&io_mutex);
	continue; /* Simply don't consider that (p,r) for now.
		     FIXME: can we find the locations to sieve? */
      }
    /* OK, all special cases are done. */

    const uint32_t bound0 = plattice_bound0(&pli, si), bound1 = plattice_bound1(&pli, si);
#if !MOD2_CLASSES_BS
    const uint64_t inc_a = plattice_a(&pli, si), inc_c = plattice_c(&pli, si);
    uint64_t x = 1ULL << (si->conf->logI-1);
    uint32_t i = x;
    FILL_BUCKET_INC_X();
    if (x >= IJ) continue;
#else
    for(unsigned int parity = 1 ; parity < 4; parity++) {
      // The sieving point (0,0) is I/2 in x-coordinate
      uint64_t x = plattice_starting_vector(&pli, si, parity);
      if (x >= IJ) continue;
      const uint64_t inc_a = plattice_a(&pli, si), inc_c = plattice_c(&pli, si);
#endif
      const prime_hint_t prime = (prime_hint_t) bucket_encode_prime (p);
      
      /* Now, do the real work: the filling of the buckets */
      do {
	/***************************************************************/
#define FILL_BUCKET() do {						\
	  unsigned int i = x & maskI;					\
	  if (LIKELY(MOD2_CLASSES_BS || (x & even_mask)			\
		     FILL_BUCKET_SKIP_GCD3()				\
		     )) FILL_BUCKET_HEART();				\
	  FILL_BUCKET_INC_X();						\
	} while (0)
	/***************************************************************/
	FILL_BUCKET(); if (x >= IJ) break;
	FILL_BUCKET(); if (x >= IJ) break;
	FILL_BUCKET(); if (x >= IJ) break;
	FILL_BUCKET();
      } while (x < IJ);
#if MOD2_CLASSES_BS
    }
#endif
  }
  th->sides[side]->BA = BA;
}

/* Same than fill_in_buckets, but with 2 passes (k_buckets & buckets). */
void
fill_in_k_buckets(thread_data_ptr th, int side, where_am_I_ptr w MAYBE_UNUSED)
{
  WHERE_AM_I_UPDATE(w, side, side);
  sieve_info_srcptr si = th->si;
  bucket_array_t BA = th->sides[side]->BA;      /* local copy, gain a register + use stack */
  k_bucket_array_t kBA = th->sides[side]->kBA;
  // Loop over all primes in the factor base.
  //
  // Note that dispatch_fb already arranged so that all the primes
  // which appear here are >= bucket_thresh and <= pmax (the latter
  // being for the moment unconditionally set to FBPRIME_MAX by the
  // caller of dispatch_fb).
  
  fb_iterator t;
  fb_iterator_init_set_fb(t, th->sides[side]->fb_bucket);
  unsigned char last_logp = 0;
  for( ; !fb_iterator_over(t) ; fb_iterator_next(t)) {
    fbprime_t p = t->fb->p;
    unsigned char logp = t->fb->plog;
    ASSERT_ALWAYS (p & 1);
    WHERE_AM_I_UPDATE(w, p, p);
    
    /* Write new set of pointers if the logp value changed */
    if (UNLIKELY(last_logp != logp)) {
      aligned_medium_memcpy((uint8_t *)kBA.logp_idx + kBA.size_b_align * BA.nr_logp, kBA.bucket_write, kBA.size_b_align);
      BA.logp_val[BA.nr_logp++] = last_logp = logp;
    }
    
    /* If we sieve for special-q's smaller than the factor
       base bound, the prime p might equal the special-q prime q. */
    if (UNLIKELY(!mpz_cmp_ui(si->doing->p, p))) continue;
    fbprime_t R = fb_iterator_get_r(t), r = fb_root_in_qlattice(p, R, t->fb->invp, si);
    
#ifdef SKIP_GCD3
    const uint32_t I = si->I;
    const unsigned int logI = si->conf->logI;
#endif
    const uint32_t maskI = si->I-1;
    const uint64_t even_mask = (1ULL << si->conf->logI) | 1ULL;
    const uint64_t IJ = ((uint64_t) si->J) << si->conf->logI;

    /* Special cases */
    if (UNLIKELY((!r) || (r >= p))) {
      if (r > p) /* should only happen for lattice-sieved prime powers,
		    which is not possible currently since maxbits < I */
	continue;
      /* r == p or r == 0.
	 1. If r == 0 (mod p), this prime hits for i == 0 (mod p), but since p > I,
	 this implies i = 0 or i > I. We don't sieve i > I. Since gcd(i,j) |
	 gcd(a,b), for i = 0 we only need to sieve j = 1. 
	 So, x = j*I + (i + I/2) = I + I/2.
	 2. r == p means root at infinity, which hits for j == 0 (mod p). Since q > I > J,
	 this implies j = 0 or j > J. This means we sieve only (i,j) = (1,0) here.
	 FIXME: what about (-1,0)? It's the same (a,b) as (1,0) but which of these two
	 (if any) do we sieve? */
      uint64_t x = (r ? 1 : si->I) + (si->I >> 1);
      prime_hint_t prime = bucket_encode_prime(p);
      /* 1. pkbut must be volatile in BIG_ENDIAN: the write order of prime & x
	 is need because the last byte of x (always 0 because x <<= 8) must
	 be overwritten by the first byte of prime.
	 2. memcpy is good because it's its job to know if it's possible to
	 write more than 1 byte at an odd adress (ODD_ADDRESS_IO_INT).
	 gcc does a optimal job with memcpy & a little + constant length.
      */
      /**************************************************************************/
#ifdef CADO_LITTLE_ENDIAN
#define FILL_K_BUCKET_HEART() do {					\
	k_bucket_update_t **pkbut = kBA.bucket_write + (x >> 24);	\
	k_bucket_update_t *kbut = *pkbut;				\
	FILL_BUCKET_TRACE_K(x);						\
	WHERE_AM_I_UPDATE(w, N, x >> 16);				\
	WHERE_AM_I_UPDATE(w, x, (uint16_t) x);				\
	memcpy(kbut, &prime, sizeof(prime_hint_t));			\
	uint32_t i = (uint32_t) x;					\
	memcpy((uint8_t *) kbut + sizeof(prime_hint_t), &i, 4);		\
	*pkbut = ++kbut;						\
	FILL_BUCKET_PREFETCH(kbut);					\
      } while (0)
#else
#define FILL_K_BUCKET_HEART() do {					\
	k_bucket_update_t **pkbut = kBA.bucket_write + (x >> 24);	\
	volatile k_bucket_update_t *kbut = *pkbut;			\
	FILL_BUCKET_TRACE_K(x);						\
	WHERE_AM_I_UPDATE(w, N, x >> 16);				\
	WHERE_AM_I_UPDATE(w, x, (uint16_t) x);				\
	uint32_t i = (uint32_t) x << 8;					\
	memcpy(kbut, &i, 4);						\
	memcpy((uint8_t *) kbut + 3, &prime, sizeof(prime_hint_t));	\
	*pkbut = ++kbut;						\
	FILL_BUCKET_PREFETCH(kbut);					\
      } while(0)
#endif
      /**************************************************************************/
      FILL_K_BUCKET_HEART();
      continue;
    }
    /* If working with congruence classes, once the loop on the parity goes at the level
       above, this initialization should in fact either be done for each congruence class,
       or saved for later use within the factor base structure. */
    plattice_info_t pli;
    if (UNLIKELY(!reduce_plattice(&pli, p, r, si))) {
      pthread_mutex_lock(&io_mutex);
      fprintf (stderr, "# fill_in_buckets: reduce_plattice() returned 0 for p = "
	       FBPRIME_FORMAT ", r = " FBPRIME_FORMAT "\n", p, r);
	pthread_mutex_unlock(&io_mutex);
	continue; /* Simply don't consider that (p,r) for now.
		     FIXME: can we find the locations to sieve? */
      }
    /* OK, all special cases are done. */

    const uint32_t bound0 = plattice_bound0(&pli, si), bound1 = plattice_bound1(&pli, si);
#if !MOD2_CLASSES_BS
    const uint64_t inc_a = plattice_a(&pli, si), inc_c = plattice_c(&pli, si);
    uint64_t x = 1ULL << (si->conf->logI-1);
    uint32_t i = x;
    FILL_BUCKET_INC_X();
    if (x >= IJ) continue;
#else
    for(unsigned int parity = 1 ; parity < 4; parity++) {
      // The sieving point (0,0) is I/2 in x-coordinate
      uint64_t x = plattice_starting_vector(&pli, si, parity);
      if (x >= IJ) continue;
      const uint64_t inc_a = plattice_a(&pli, si), inc_c = plattice_c(&pli, si);
#endif
      const prime_hint_t prime = (prime_hint_t) bucket_encode_prime (p);
      
      /* Now, do the real work: the filling of the k-buckets */
      do { 
	/**********************************************************************/
#define FILL_K_BUCKET() do {						\
	  unsigned int i = x & maskI;					\
	  if (LIKELY(MOD2_CLASSES_BS || (x & even_mask)			\
		     FILL_BUCKET_SKIP_GCD3()				\
		     )) FILL_K_BUCKET_HEART();				\
	  FILL_BUCKET_INC_X();						\
	} while (0)
	/****************************************************************/
	FILL_K_BUCKET(); if (x >= IJ) break;
	FILL_K_BUCKET(); if (x >= IJ) break;
	FILL_K_BUCKET(); if (x >= IJ) break;
	FILL_K_BUCKET();
      } while (x < IJ);
#if MOD2_CLASSES_BS
    }
#endif
  }
  th->sides[side]->kBA = kBA;
  
  /* sort : 2nd pass; kBA -> BA */
  bucket_update_t **pbw = BA.bucket_write;
  for (uint32_t kb = 0; kb < kBA.n_bucket; ++kb) {
    uint8_t *kbs = (uint8_t *) (kBA.bucket_start[kb]);
    /* First part: we rewrite 1->256 buckets, and in the same time,
       we have to deal with the rewriting of logp_idx.
       I use a block in these part to see the range of variables */
    {
      size_t lg = (size_t) BA.bucket_write + BA.size_b_align - (size_t) pbw;
      if (LIKELY(lg > sizeof(bucket_update_t **) << 8)) lg = sizeof(bucket_update_t **) << 8;
      bucket_update_t **pbl = BA.logp_idx + (kb << 8);
      k_bucket_update_t **pkbl = kBA.logp_idx + kb;
      uint8_t *kbl = (uint8_t *) *pkbl;
      /* There are BA.nr_logp duplicates of all kBA.bucket_write in kBA.logp_idx. */
      for (uint8_t nr_logp = BA.nr_logp; nr_logp; --nr_logp) {
	/* Twelve kBA records in one time : it's the nearest of a cache line (60 bytes) */
	for (;
	     kbs +  sizeof(k_bucket_update_t)*12 <= kbl;
	     kbs += sizeof(k_bucket_update_t)*12) {
	  /*****************************************************************/
#ifdef CADO_LITTLE_ENDIAN
#define KBA_2_BA(A) do {						\
	    bucket_update_t **pbut, *but;				\
	    pbut = pbw + kbs[(A) + sizeof(bucket_update_t)];		\
	    but = *pbut;						\
	    memcpy(but, kbs + (A), optimal_move[sizeof(bucket_update_t)]); \
	    *pbut = ++but;						\
	  } while (0)
#else
#define KBA_2_BA(A) do {						\
	    bucket_update_t **pbut, *but;				\
	    pbut = pbw + kbs[A];					\
	    but = *pbut;						\
	    memcpy(but, kbs+(A)+1, optimal_move[sizeof(bucket_update_t)]); \
	    *pbut = ++but;						\
	  } while (0)
#endif
	  /*****************************************************************/
	  KBA_2_BA(0);                            KBA_2_BA(sizeof(k_bucket_update_t));
	  KBA_2_BA(sizeof(k_bucket_update_t)*2);  KBA_2_BA(sizeof(k_bucket_update_t)*3);
	  KBA_2_BA(sizeof(k_bucket_update_t)*4);  KBA_2_BA(sizeof(k_bucket_update_t)*5);
	  KBA_2_BA(sizeof(k_bucket_update_t)*6);  KBA_2_BA(sizeof(k_bucket_update_t)*7);
	  KBA_2_BA(sizeof(k_bucket_update_t)*8);  KBA_2_BA(sizeof(k_bucket_update_t)*9);
	  KBA_2_BA(sizeof(k_bucket_update_t)*10); KBA_2_BA(sizeof(k_bucket_update_t)*11);
	}
	for (; kbs < kbl; kbs += sizeof(k_bucket_update_t)) KBA_2_BA(0);
	/* OK, let's duplicate the current (at most) 256 pointers from
	   BA.bucket_write in BA.logp_idx */
	aligned_medium_memcpy(pbl, pbw, lg);
	pbl =    (bucket_update_t **) ((size_t)  pbl +  BA.size_b_align);
	pkbl = (k_bucket_update_t **) ((size_t) pkbl + kBA.size_b_align);
	kbl = (uint8_t *) *pkbl;
      }
    }
    /* 2nd part: BA.logp_idx is rewritten. We finish the rewrite of the current bucket */
    const uint8_t *kbw = (uint8_t *) (kBA.bucket_write[kb]);
    for (;
	 kbs +  sizeof(k_bucket_update_t)*12 <= kbw;
	 kbs += sizeof(k_bucket_update_t)*12) {
      KBA_2_BA(0);                            KBA_2_BA(sizeof(k_bucket_update_t));
      KBA_2_BA(sizeof(k_bucket_update_t)*2);  KBA_2_BA(sizeof(k_bucket_update_t)*3);
      KBA_2_BA(sizeof(k_bucket_update_t)*4);  KBA_2_BA(sizeof(k_bucket_update_t)*5);
      KBA_2_BA(sizeof(k_bucket_update_t)*6);  KBA_2_BA(sizeof(k_bucket_update_t)*7);
      KBA_2_BA(sizeof(k_bucket_update_t)*8);  KBA_2_BA(sizeof(k_bucket_update_t)*9);
      KBA_2_BA(sizeof(k_bucket_update_t)*10); KBA_2_BA(sizeof(k_bucket_update_t)*11);
    }
    for (; kbs < kbw; kbs += sizeof(k_bucket_update_t)) KBA_2_BA(0);
    pbw += 256;
  }
  th->sides[side]->BA = BA;
}

void
fill_in_m_buckets(thread_data_ptr th, int side, where_am_I_ptr w MAYBE_UNUSED)
{
  WHERE_AM_I_UPDATE(w, side, side);
  sieve_info_srcptr si = th->si;
  bucket_array_t BA = th->sides[side]->BA;  /* local copy; gain a register + use stack */
  k_bucket_array_t kBA = th->sides[side]->kBA;
  m_bucket_array_t mBA = th->sides[side]->mBA;
  // Loop over all primes in the factor base.
  //
  // Note that dispatch_fb already arranged so that all the primes
  // which appear here are >= bucket_thresh and <= pmax (the latter
  // being for the moment unconditionally set to FBPRIME_MAX by the
  // caller of dispatch_fb).
  
  fb_iterator t;
  fb_iterator_init_set_fb(t, th->sides[side]->fb_bucket);
  unsigned char last_logp = 0;
  for( ; !fb_iterator_over(t) ; fb_iterator_next(t)) {
    fbprime_t p = t->fb->p;
    unsigned char logp = t->fb->plog;
    ASSERT_ALWAYS (p & 1);
    WHERE_AM_I_UPDATE(w, p, p);
    
    /* Write new set of pointers if the logp value changed */
    if (UNLIKELY(last_logp != logp)) {
      aligned_medium_memcpy((uint8_t *)mBA.logp_idx + mBA.size_b_align * BA.nr_logp, mBA.bucket_write, mBA.size_b_align);
      BA.logp_val[BA.nr_logp++] = last_logp = logp;
    }
    
    /* If we sieve for special-q's smaller than the factor
       base bound, the prime p might equal the special-q prime q. */
    if (UNLIKELY(mpz_cmp_ui(si->doing->p, p) == 0)) continue;
    fbprime_t R = fb_iterator_get_r(t), r = fb_root_in_qlattice(p, R, t->fb->invp, si);
    
#ifdef SKIP_GCD3
    const uint32_t I = si->I;
    const unsigned int logI = si->conf->logI;
#endif
    const uint32_t maskI = si->I-1;
    const uint64_t even_mask = (1ULL << si->conf->logI) | 1ULL;
    const uint64_t IJ = ((uint64_t) si->J) << si->conf->logI;

    /* Special cases */
    if (UNLIKELY((!r) || (r >= p))) {
      if (r > p) /* should only happen for lattice-sieved prime powers,
		    which is not possible currently since maxbits < I */
	continue;
      /* r == p or r == 0.
	 1. If r == 0 (mod p), this prime hits for i == 0 (mod p), but since p > I,
	 this implies i = 0 or i > I. We don't sieve i > I. Since gcd(i,j) |
	 gcd(a,b), for i = 0 we only need to sieve j = 1. 
	 So, x = j*I + (i + I/2) = I + I/2.
	 2. r == p means root at infinity, which hits for j == 0 (mod p). Since q > I > J,
	 this implies j = 0 or j > J. This means we sieve only (i,j) = (1,0) here.
	 FIXME: what about (-1,0)? It's the same (a,b) as (1,0) but which of these two
	 (if any) do we sieve? */
      uint64_t x = (r ? 1 : si->I) + (si->I >> 1);
      prime_hint_t prime = bucket_encode_prime(p);
      /* memcpy is good because it's its job to know if it's possible to
	 write an int to all even addresses.
	 gcc does a optimal job with memcpy & a little + constant length.
      */
      /**************************************************************************/
#ifdef CADO_LITTLE_ENDIAN
#define FILL_M_BUCKET_HEART() do {					\
	m_bucket_update_t **pmbut = mBA.bucket_write + (x >> 32);	\
	m_bucket_update_t *mbut = *pmbut;				\
	FILL_BUCKET_TRACE_K(x);						\
	WHERE_AM_I_UPDATE(w, N, x >> 16);				\
	WHERE_AM_I_UPDATE(w, x, (uint16_t) x);				\
	memcpy(mbut, &prime, sizeof(prime_hint_t));			\
	uint32_t i = (uint32_t) x;					\
	memcpy((uint8_t *) mbut + sizeof(prime_hint_t), &i, 4);	\
	*pmbut = ++mbut;						\
	FILL_BUCKET_PREFETCH(mbut);					\
      } while (0)
#else
#define FILL_M_BUCKET_HEART() do {					\
	m_bucket_update_t **pmbut = mBA.bucket_write + (x >> 32);	\
	m_bucket_update_t *mbut = *pmbut;				\
	FILL_BUCKET_TRACE_K(x);						\
	WHERE_AM_I_UPDATE(w, N, x >> 16);				\
	WHERE_AM_I_UPDATE(w, x, (uint16_t) x);				\
	uint32_t i = (uint32_t) x;					\
	memcpy(mbut, &i, 4);						\
	memcpy((uint8_t *) mbut + 4, &prime, sizeof(prime_hint_t));	\
	*pmbut = ++mbut;						\
	FILL_BUCKET_PREFETCH(mbut);					\
      } while (0)
#endif
      /**************************************************************************/
      FILL_M_BUCKET_HEART();
      continue;
    }
    /* If working with congruence classes, once the loop on the parity goes at the level
       above, this initialization should in fact either be done for each congruence class,
       or saved for later use within the factor base structure. */
    plattice_info_t pli;
    if (UNLIKELY(!reduce_plattice(&pli, p, r, si))) {
      pthread_mutex_lock(&io_mutex);
      fprintf (stderr, "# fill_in_buckets: reduce_plattice() returned 0 for p = "
	       FBPRIME_FORMAT ", r = " FBPRIME_FORMAT "\n", p, r);
	pthread_mutex_unlock(&io_mutex);
	continue; /* Simply don't consider that (p,r) for now.
		     FIXME: can we find the locations to sieve? */
      }
    /* OK, all special cases are done. */

    const uint32_t bound0 = plattice_bound0(&pli, si), bound1 = plattice_bound1(&pli, si);
#if !MOD2_CLASSES_BS
    const uint64_t inc_a = plattice_a(&pli, si), inc_c = plattice_c(&pli, si);
    uint64_t x = 1ULL << (si->conf->logI-1);
    uint32_t i = x;
    FILL_BUCKET_INC_X();
    if (x >= IJ) continue;
#else
    for(unsigned int parity = 1 ; parity < 4; parity++) {
      // The sieving point (0,0) is I/2 in x-coordinate
      uint64_t x = plattice_starting_vector(&pli, si, parity);
      if (x >= IJ) continue;
      const uint64_t inc_a = plattice_a(&pli, si), inc_c = plattice_c(&pli, si);
#endif
      const prime_hint_t prime = bucket_encode_prime (p);

      /* Now, do the real work: the filling of the m-buckets */
      do {
	/***************************************************************/
#define FILL_M_BUCKET() do {						\
	  unsigned int i = x & maskI;					\
	  if (LIKELY(MOD2_CLASSES_BS || (x & even_mask)			\
		     FILL_BUCKET_SKIP_GCD3()				\
		     )) FILL_M_BUCKET_HEART();				\
	  FILL_BUCKET_INC_X();						\
	} while (0)
	/***************************************************************/
	FILL_M_BUCKET(); if (x >= IJ) break;
	FILL_M_BUCKET(); if (x >= IJ) break;
	FILL_M_BUCKET(); if (x >= IJ) break;
	FILL_M_BUCKET();
      } while (x < IJ);
#if MOD2_CLASSES_BS
    }
#endif
  }
  th->sides[side]->mBA = mBA;

  /* sort : 2nd pass; mBA -> kBA */
  k_bucket_update_t **pkbw = kBA.bucket_write;
  for (uint32_t mb = 0; mb < mBA.n_bucket; ++mb) {
    uint8_t *mbs = (uint8_t *) (mBA.bucket_start[mb]);
    /* First part: we rewrite 1->256 buckets, and in the same time,
       we have to deal with the rewriting of logp_idx.
       I use a block in these part to see the range of variables */
    {
      size_t lg = (size_t) kBA.bucket_write + kBA.size_b_align - (size_t) pkbw;
      if (LIKELY(lg > sizeof(k_bucket_update_t **) << 8)) lg = sizeof(k_bucket_update_t **) << 8;
      k_bucket_update_t **pkbl = kBA.logp_idx + (mb << 8);
      m_bucket_update_t **pmbl = mBA.logp_idx + mb;
      uint8_t *mbl = (uint8_t *) *pmbl;
      /* There are BA.nr_logp duplicates of all mBA.bucket_write in mBA.logp_idx. */
      for (uint8_t nr_logp = BA.nr_logp; nr_logp; --nr_logp) {
	/* Ten mBA records in one time : it's the nearest of a cache line (60 bytes) */
	for (;
	     mbs +  sizeof(m_bucket_update_t)*10 <= mbl;
	     mbs += sizeof(m_bucket_update_t)*10) {
	  /*****************************************************************/
#ifdef CADO_LITTLE_ENDIAN
#define MBA_2_KBA(A) do {						\
	    k_bucket_update_t **pkbut, *kbut;				\
	    pkbut = pkbw + mbs[(A)+sizeof(k_bucket_update_t)];		\
	    kbut = *pkbut;						\
	    memcpy(kbut, mbs+(A), optimal_move[sizeof(k_bucket_update_t)]); \
	    *pkbut = ++kbut;						\
	  } while (0)
#else
#define MBA_2_KBA(A) do {						\
	    k_bucket_update_t **pkbut, *kbut;				\
	    pkbut = pkbw + mbs[A];					\
	    kbut = *pkbut;						\
	    memcpy(kbut, mbs+(A)+1, optimal_move[sizeof(k_bucket_update_t)]); \
	    *pkbut = ++kbut;						\
	  } while (0)
#endif
	  /*****************************************************************/
	  MBA_2_KBA(0);	                          MBA_2_KBA(sizeof(m_bucket_update_t));
	  MBA_2_KBA(sizeof(m_bucket_update_t)*2); MBA_2_KBA(sizeof(m_bucket_update_t)*3);
	  MBA_2_KBA(sizeof(m_bucket_update_t)*4); MBA_2_KBA(sizeof(m_bucket_update_t)*5);
	  MBA_2_KBA(sizeof(m_bucket_update_t)*6); MBA_2_KBA(sizeof(m_bucket_update_t)*7);
	  MBA_2_KBA(sizeof(m_bucket_update_t)*8); MBA_2_KBA(sizeof(m_bucket_update_t)*9);
	}
	for (; mbs < mbl; mbs += sizeof(m_bucket_update_t)) MBA_2_KBA(0);
	/* OK, let's duplicate the current (at most) 256 pointers in
	   kBA.bucket_write in kBA.logp_idx */
	aligned_medium_memcpy(pkbl, pkbw, lg);
	pkbl = (k_bucket_update_t **) ((size_t) pkbl + kBA.size_b_align);
	pmbl = (m_bucket_update_t **) ((size_t) pmbl + mBA.size_b_align);
	mbl = (uint8_t *) *pmbl;
      }
    }
    /* 2nd part: kBA.logp_idx is rewritten. We finish the rewrite of the current bucket */
    const uint8_t *mbw = (uint8_t *) (mBA.bucket_write[mb]);
    for (;
	 mbs +  sizeof(m_bucket_update_t)*10 <= mbw;
	 mbs += sizeof(m_bucket_update_t)*10) {
      MBA_2_KBA(0);
      MBA_2_KBA(sizeof(m_bucket_update_t));
      MBA_2_KBA(sizeof(m_bucket_update_t)*2);
      MBA_2_KBA(sizeof(m_bucket_update_t)*3);
      MBA_2_KBA(sizeof(m_bucket_update_t)*4);
      MBA_2_KBA(sizeof(m_bucket_update_t)*5);
      MBA_2_KBA(sizeof(m_bucket_update_t)*6);
      MBA_2_KBA(sizeof(m_bucket_update_t)*7);
      MBA_2_KBA(sizeof(m_bucket_update_t)*8);
      MBA_2_KBA(sizeof(m_bucket_update_t)*9);
    }
    for (; mbs < mbw; mbs += sizeof(m_bucket_update_t)) MBA_2_KBA(0);
    pkbw += 256;
  }
  th->sides[side]->kBA = kBA;

  /* sort : 3th pass; kBA -> BA */
  bucket_update_t **pbw = BA.bucket_write;
  for (uint32_t kb = 0; kb < kBA.n_bucket; ++kb) {
    uint8_t *kbs = (uint8_t *) (kBA.bucket_start[kb]);
    /* First part: we rewrite 1->256 buckets, and in the same time,
       we have to deal with the rewriting of logp_idx.
       I use a block in these part to see the range of variables */
    {
      size_t lg = (size_t) BA.bucket_write + BA.size_b_align - (size_t) pbw;
      if (LIKELY(lg > sizeof(bucket_update_t **) << 8)) lg = sizeof(bucket_update_t **) << 8;
      bucket_update_t **pbl = BA.logp_idx + (kb << 8);
      k_bucket_update_t **pkbl = kBA.logp_idx + kb;
      uint8_t *kbl = (uint8_t *) *pkbl;
      /* There are BA.nr_logp duplicates of all kBA.bucket_write in kBA.logp_idx. */
      for (uint8_t nr_logp = BA.nr_logp; nr_logp; --nr_logp) {
	/* Twelve kBA records in one time : it's the nearest of a cache line (60 bytes) */
	for (;
	     kbs +  sizeof(k_bucket_update_t)*12 <= kbl;
	     kbs += sizeof(k_bucket_update_t)*12) {
	  /* See the define of KBA_2_BA in fill_in_k_buckets */
	  KBA_2_BA(0);                            KBA_2_BA(sizeof(k_bucket_update_t));
	  KBA_2_BA(sizeof(k_bucket_update_t)*2);  KBA_2_BA(sizeof(k_bucket_update_t)*3);
	  KBA_2_BA(sizeof(k_bucket_update_t)*4);  KBA_2_BA(sizeof(k_bucket_update_t)*5);
	  KBA_2_BA(sizeof(k_bucket_update_t)*6);  KBA_2_BA(sizeof(k_bucket_update_t)*7);
	  KBA_2_BA(sizeof(k_bucket_update_t)*8);  KBA_2_BA(sizeof(k_bucket_update_t)*9);
	  KBA_2_BA(sizeof(k_bucket_update_t)*10); KBA_2_BA(sizeof(k_bucket_update_t)*11);
	}
	for (; kbs < kbl; kbs += sizeof(k_bucket_update_t)) KBA_2_BA(0);
	/* OK, let's duplicate the current (at most) 256 pointers from
	   BA.bucket_write in BA.logp_idx */
	aligned_medium_memcpy(pbl, pbw, lg);
	pbl =    (bucket_update_t **) ((size_t)  pbl +  BA.size_b_align);
	pkbl = (k_bucket_update_t **) ((size_t) pkbl + kBA.size_b_align);
	kbl = (uint8_t *) *pkbl;
      }
    }
    /* 2nd part: BA.logp_idx is rewritten. We finish the rewrite of the current bucket */
    const uint8_t *kbw = (uint8_t *) (kBA.bucket_write[kb]);
    for (;
	 kbs +  sizeof(k_bucket_update_t)*12 <= kbw;
	 kbs += sizeof(k_bucket_update_t)*12) {
      KBA_2_BA(0);                            KBA_2_BA(sizeof(k_bucket_update_t));
      KBA_2_BA(sizeof(k_bucket_update_t)*2);  KBA_2_BA(sizeof(k_bucket_update_t)*3);
      KBA_2_BA(sizeof(k_bucket_update_t)*4);  KBA_2_BA(sizeof(k_bucket_update_t)*5);
      KBA_2_BA(sizeof(k_bucket_update_t)*6);  KBA_2_BA(sizeof(k_bucket_update_t)*7);
      KBA_2_BA(sizeof(k_bucket_update_t)*8);  KBA_2_BA(sizeof(k_bucket_update_t)*9);
      KBA_2_BA(sizeof(k_bucket_update_t)*10); KBA_2_BA(sizeof(k_bucket_update_t)*11);
    }
    for (; kbs < kbw; kbs += sizeof(k_bucket_update_t)) KBA_2_BA(0);
    pbw += 256;
  }
  th->sides[side]->BA = BA;
}

void * fill_in_buckets_both(thread_data_ptr th)
{
    where_am_I w;
    WHERE_AM_I_UPDATE(w, si, th->si);

    if (th->sides[ALGEBRAIC_SIDE]->BA.n_bucket < THRESHOLD_K_BUCKETS)
      fill_in_buckets(th, ALGEBRAIC_SIDE, w);
    else if (th->sides[ALGEBRAIC_SIDE]->BA.n_bucket < THRESHOLD_M_BUCKETS)
      fill_in_k_buckets(th, ALGEBRAIC_SIDE, w);
    else
      fill_in_m_buckets(th, ALGEBRAIC_SIDE, w);

    if (th->sides[RATIONAL_SIDE]->BA.n_bucket < THRESHOLD_K_BUCKETS)
      fill_in_buckets(th, RATIONAL_SIDE, w);
    else if (th->sides[RATIONAL_SIDE]->BA.n_bucket < THRESHOLD_M_BUCKETS)
      fill_in_k_buckets(th, RATIONAL_SIDE, w);
    else
      fill_in_m_buckets(th, RATIONAL_SIDE, w);
    return NULL;
}
/* }}} */

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
              pthread_mutex_lock(&io_mutex);
              fprintf (stderr,
                       "# Error, p = %lu does not divide at (N,x) = (%u,%d)\n",
                       (unsigned long) prime.p, N, x);
              pthread_mutex_unlock(&io_mutex);
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

    if (trial_div_very_verbose) {
        pthread_mutex_lock(&io_mutex);
        gmp_fprintf (stderr, "# trial_div() entry, x = %d, norm = %Zd\n",
                x, norm);
        pthread_mutex_unlock(&io_mutex);
    }

    // handle 2 separately, if it is in fb
    if (fb->p == 2) {
        int bit = mpz_scan1(norm, 0);
        int i;
        for (i = 0; i < bit; ++i) {
            fl->fac[fl->n] = 2;
            fl->n++;
        }
        if (trial_div_very_verbose) {
            pthread_mutex_lock(&io_mutex);
            gmp_fprintf (stderr, "# x = %d, dividing out 2^%d, norm = %Zd\n",
                    x, bit, norm);
            pthread_mutex_unlock(&io_mutex);
        }
        mpz_tdiv_q_2exp(norm, norm, bit);
        fb = fb_next (fb); // cannot do fb++, due to variable size !
    }

    // remove primes in "primes" that map to x
    divide_primes_from_bucket (fl, norm, N, x, primes, fbb);
#ifdef TRACE_K /* {{{ */
    if (trace_on_spot_ab(a,b) && fl->n) {
        fprintf(stderr, "# divided by 2 + primes from bucket that map to %u: ", x);
        if (!factor_list_fprint(stderr, *fl)) fprintf(stderr, "(none)");
        gmp_fprintf(stderr, ", remaining norm is %Zd\n", norm);
    }
#endif /* }}} */
    if (trial_div_very_verbose) {
        pthread_mutex_lock(&io_mutex);
        gmp_fprintf (stderr, "# x = %d, after dividing out bucket/resieved norm = %Zd\n",
                x, norm);
        pthread_mutex_unlock(&io_mutex);
    }

    do {
      /* Trial divide primes with precomputed tables */
#define TRIALDIV_MAX_FACTORS 32
      int i;
      unsigned long factors[TRIALDIV_MAX_FACTORS];
      if (trial_div_very_verbose)
      {
          pthread_mutex_lock(&io_mutex);
          fprintf (stderr, "# Trial division by ");
          for (i = 0; trialdiv_data[i].p != 1; i++)
              fprintf (stderr, " %lu", trialdiv_data[i].p);
          fprintf (stderr, "\n");
          pthread_mutex_unlock(&io_mutex);
      }

      nr_factors = trialdiv (factors, norm, trialdiv_data, TRIALDIV_MAX_FACTORS);

      for (i = 0; i < MIN(nr_factors, TRIALDIV_MAX_FACTORS); i++)
      {
          if (trial_div_very_verbose) {
              pthread_mutex_lock(&io_mutex);
              fprintf (stderr, " %lu", factors[i]);
              pthread_mutex_unlock(&io_mutex);
          }
          factor_list_add (fl, factors[i]);
      }
      if (trial_div_very_verbose) {
          pthread_mutex_lock(&io_mutex);
          gmp_fprintf (stderr, "\n# After trialdiv(): norm = %Zd\n", norm);
          pthread_mutex_unlock(&io_mutex);
      }
    } while (nr_factors == TRIALDIV_MAX_FACTORS + 1);
}
/* }}} */

/************************ cofactorization ********************************/

/* {{{ cofactoring area */

/* Return 0 if the leftover norm n cannot yield a relation.
   FIXME: need to check L^k < n < B^(k+1) too.
   XXX: In doing this, pay attention to the fact that for the descent,
   we might have B^2<L.

   Possible cases, where qj represents a prime in [B,L], and rj a prime > L:
   (0) n >= 2^mfb
   (a) n < L:           1 or q1
   (b) L < n < B^2:     r1 -> cannot yield a relation
   (c) B^2 < n < B*L:   r1 or q1*q2
   (d) B*L < n < L^2:   r1 or q1*q2 or q1*r2
   (e) L^2 < n < B^3:   r1 or q1*r2 or r1*r2 -> cannot yield a relation
   (f) B^3 < n < B^2*L: r1 or q1*r2 or r1*r2 or q1*q2*q3
   (g) B^2*L < n < L^3: r1 or q1*r2 or r1*r2
   (h) L^3 < n < B^4:   r1 or q1*r2, r1*r2 or q1*q2*r3 or q1*r2*r3 or r1*r2*r3
*/
static int
check_leftover_norm (mpz_t n, sieve_info_ptr si, int side)
{
  size_t s = mpz_sizeinbase (n, 2);
  unsigned int lpb = si->conf->sides[side]->lpb;
  unsigned int mfb = si->conf->sides[side]->mfb;

  if (s > mfb)
    return 0; /* n has more than mfb bits, which is the given limit */
  /* now n < 2^mfb */
  if (s <= lpb)
    return 1; /* case (a) */
    /* Note also that in the descent case where L > B^2, if we're below L
     * it's still fine of course, but we have no guarantee that our
     * cofactor is prime... */
  /* now n >= L=2^lpb */
  if (mpz_cmp (n, si->BB[side]) < 0)
    return 0; /* case (b) */
  /* now n >= B^2 */
  if (2 * lpb < s)
    {
      if (mpz_cmp (n, si->BBB[side]) < 0)
        return 0; /* case (e) */
      if (3 * lpb < s && mpz_cmp (n, si->BBBB[side]) < 0)
        return 0; /* case (h) */
    }
  if (mpz_probab_prime_p (n, 1))
    return 0; /* n is a pseudo-prime larger than L */
  return 1;
}

/* Adds the number of sieve reports to *survivors,
   number of survivors with coprime a, b to *coprimes */

NOPROFILE_STATIC int
factor_survivors (thread_data_ptr th, int N, unsigned char * S[2], where_am_I_ptr w MAYBE_UNUSED)
{
    las_info_ptr las = th->las;
    sieve_info_ptr si = th->si;
    cado_poly_ptr cpoly = si->cpoly;
    sieve_side_info_ptr rat = si->sides[RATIONAL_SIDE];
    sieve_side_info_ptr alg = si->sides[ALGEBRAIC_SIDE];
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

    for(int side = 0 ; side < 2 ; side++) {
        f[side] = alloc_mpz_array (8);
        m[side] = alloc_uint32_array (8);

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
        gmp_fprintf(stderr, "# Remaining norms which have not been accounted for in sieving: (%Zd, %Zd)\n", traced_norms[0], traced_norms[1]);
    }
#endif  /* }}} */

    /* XXX: Don't believe that resieve_start is easily changeable... */
    const int resieve_start = RATIONAL_SIDE;

    /* This is the one which gets the merged information in the end */
    unsigned char * SS = S[resieve_start];

#ifdef UNSIEVE_NOT_COPRIME
    unsieve_not_coprime (SS, N, si);
#endif

    for (int x = 0; x < BUCKET_REGION; ++x)
    {
#ifdef TRACE_K /* {{{ */
        if (trace_on_spot_Nx(N, x)) {
            fprintf(stderr, "# alg->Bound[%u]=%u, rat->Bound[%u]=%u\n",
                    alg_S[trace_Nx.x], alg_S[x] <= alg->bound ? 0 : alg->bound,
                    rat_S[trace_Nx.x], rat_S[x] <= rat->bound ? 0 : rat->bound);
        }
#endif /* }}} */
        unsigned int X;
        unsigned int i, j;

        if (!sieve_info_test_lognorm(alg->bound, rat->bound, alg_S[x], rat_S[x]))
        {
            SS[x] = 255;
            continue;
        }
        th->rep->survivor_sizes[rat_S[x]][alg_S[x]]++;
        surv++;

        X = x + (((uint64_t) N) << LOG_BUCKET_REGION);
        i = abs ((int) (X & (si->I - 1)) - si->I / 2);
        j = X >> si->conf->logI;
#ifndef UNSIEVE_NOT_COPRIME
        if (bin_gcd_int64_safe (i, j) != 1)
        {
#ifdef TRACE_K
            if (trace_on_spot_Nx(N, x)) {
                fprintf(stderr, "# Slot [%u] in bucket %u has non coprime (i,j)=(%d,%u)\n",
                        trace_Nx.x, trace_Nx.N, i, j);
            }
#endif
            SS[x] = 255;
            continue;
        }
#endif
    }

    /* Copy those bucket entries that belong to sieving survivors and
       store them with the complete prime */
    /* FIXME: choose a sensible size here */

    for(int z = 0 ; z < 2 ; z++) {
        int side = resieve_start ^ z;
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
    const int together = sizeof(__m128i);
    __m128i ones128 = (__m128i) {-1,-1};
#else
    const int together = sizeof(unsigned long);
#endif

    for (int xul = 0; xul < BUCKET_REGION; xul += together) {
#ifdef TRACE_K
        if ((unsigned int) N == trace_Nx.N && (unsigned int) xul <= trace_Nx.x && (unsigned int) xul + sizeof (unsigned long) > trace_Nx.x) {
            fprintf(stderr, "# Slot [%u] in bucket %u has value %u\n",
                    trace_Nx.x, trace_Nx.N, SS[trace_Nx.x]);
        }
#endif
#if defined(HAVE_SSE41) && defined(SSE_SURVIVOR_SEARCH)
        if (_mm_testc_si128(*(__m128i *)(SS + xul), ones128))
            continue;
#else
        if (*(unsigned long *)(SS + xul) == (unsigned long)(-1L)) 
            continue;
#endif
        for (int x = xul; x < xul + (int) together; ++x) {
            if (SS[x] == 255) continue;

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
                xul = BUCKET_REGION;
                fprintf(las->output, "# [descent] Aborting, deadline passed\n");
                break;
            }
#endif  /* DLP_DESCENT */

            // Compute algebraic and rational norms.
            NxToAB (&a, &b, N, x, si);

#ifdef TRACE_K
            if (trace_on_spot_ab(a, b)) {
                fprintf(stderr, "# about to start cofactorization for (%" PRId64 ",%" PRIu64 ")  %d %u\n",a,b, x, SS[x]);
            }
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
            if (b == 0 || (mpz_cmp_ui(si->doing->p, b) <= 0 && b % mpz_get_ui(si->doing->p) == 0))
                continue;

            copr++;

            int pass = 1;

            for(int z = 0 ; pass && z < 2 ; z++) {
                int side = RATIONAL_SIDE ^ z;   /* start with rational */
                mpz_t * f = si->sides[side]->fij;
                int deg = cpoly->pols[side]->degree;
                int lim = si->conf->sides[side]->lim;
                int i;
                unsigned int j;

                // Trial divide rational norm
                /* Compute the norms using the polynomials transformed to 
                   i,j-coordinates. The transformed polynomial on the 
                   special-q side is already divided by q */
                NxToIJ (&i, &j, N, x, si);
                mp_poly_homogeneous_eval_siui (norm[side], f, deg, i, j);
#ifdef TRACE_K
                if (trace_on_spot_ab(a, b)) {
                    gmp_fprintf(stderr, "# start trial division for norm=%Zd on %s side for (%" PRId64 ",%" PRIu64 ")\n",norm[side],sidenames[side],a,b);
                }
#endif
                trial_div (&factors[side], norm[side], N, x,
                        si->sides[side]->fb,
                        &primes[side], si->sides[side]->trialdiv_data,
                        lim, a, b);

                pass = check_leftover_norm (norm[side], si, side);
#ifdef TRACE_K
                if (trace_on_spot_ab(a, b)) {
                    gmp_fprintf(stderr, "# checked leftover norm=%Zd on %s side for (%" PRId64 ",%" PRIu64 "): %d\n",norm[side],sidenames[side],a,b,pass);
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

            unsigned int nbits[2];
            int first;

            for (int z = 0; z < 2; z++)
              nbits[z] = mpz_sizeinbase (norm[z], 2);

            if (cof_calls[0][nbits[0]] > 0 && cof_calls[1][nbits[1]] > 0)
              {
                if (cof_fails[0][nbits[0]] / cof_calls[0][nbits[0]] >
                    cof_fails[1][nbits[1]] / cof_calls[1][nbits[1]])
                  first = 0;
                else
                  first = 1;
              }
            else
              /* if norm[RATIONAL_SIDE] is above BLPrat, then it might not
               * be smooth. We factor it first. Otherwise we factor it last. */
              first = mpz_cmp (norm[RATIONAL_SIDE], BLPrat) > 0
                ? RATIONAL_SIDE : ALGEBRAIC_SIDE;

            for (int z = 0 ; pass > 0 && z < 2 ; z++)
              {
                int side = first ^ z;
                cof_calls[side][nbits[side]] ++;
                pass = factor_leftover_norm (norm[side], f[side], m[side],
                                             si, side);
                if (pass <= 0)
                  cof_fails[side][nbits[side]] ++;
#ifdef TRACE_K
                if (trace_on_spot_ab(a, b) && pass == 0)
                  gmp_fprintf (stderr, "# factor_leftover_norm failed on %s side for (%" PRId64 ",%" PRIu64 "), remains %Zd unfactored\n", sidenames[side], a, b, norm[side]);
#endif
              }

            if (pass <= 0) continue; /* a factor was > 2^lpb, or some
                                        factorization was incomplete */

            /* yippee: we found a relation! */

            if (stats == 1) /* learning phase */
                cof_succ[cof_rat_bitsize][cof_alg_bitsize] ++;

#ifdef UNSIEVE_NOT_COPRIME
            ASSERT (bin_gcd_int64_safe (a, b) == 1);
#endif

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

            relation_compress_rat_primes(rel);
            relation_compress_alg_primes(rel);

#ifdef TRACE_K
            if (trace_on_spot_ab(a, b)) {
                fprintf(stderr, "# Relation for (%" PRId64 ",%" PRIu64 ") printed\n", a, b);
            }
#endif

#if 0   /* incompatible with the todo list */
            if (!si->bench)
#else
            if (1)
#endif
            {
                pthread_mutex_lock(&io_mutex);
#if 0
                fprint_relation(las->output, rel);
#else
                /* This code will be dropped soon. The thing is
                 * that las is a moving target for the moment, and
                 * going through the fprint_relation code above changes
                 * the order of factors in printed relations. It's not so
                 * handy.
                 */
                if (create_descent_hints) {
                    fprintf (las->output, "(%1.4f) ", seconds() - tt_qstart);
                }
                fprintf (las->output, "%" PRId64 ",%" PRIu64, a, b);
                for(int z = 0 ; z < 2 ; z++) {
                    int side = RATIONAL_SIDE ^ z;
                    fprintf (las->output, ":");
                    int comma = factor_list_fprint (las->output, factors[side]);
                    for (unsigned int i = 0; i < f[side]->length; ++i) {
                        for (unsigned int j = 0; j < m[side]->data[i]; j++) {
                            if (comma++) fprintf (las->output, ",");
                            gmp_fprintf (las->output, "%Zx", f[side]->data[i]);
                        }
                    }
                    if (si->conf->side == side) {
                        if (comma++) fprintf (las->output, ",");
                        gmp_fprintf (las->output, "%Zx", si->doing->p);
                    }
                }
                fprintf (las->output, "\n");
                fflush (las->output);
#endif
                pthread_mutex_unlock(&io_mutex);
            }

            cpt++;
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
                        /* Currently we're assuming that the base of
                         * precomputed logs goes at least as far as 
                         * cpoly->pols[side]->lim, which is asserted to
                         * exceed si->conf->sides[side]->lim). Therefore
                         * the second if() is a no-op. */
                        if (mpz_cmp_ui(f[side]->data[i], cpoly->pols[side]->lim) <= 0)
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
                            // fprintf(las->output, "# [descent] Warning: cannot estimate refactoring time for relation involving %d%c\n", n, sidenames[side][0]);
                            continue;
                        }
                        descent_hint_ptr h = las->hint_table[k];
                        time_left += h->expected_time;
                    }
                }
                // fprintf(las->output, "# [descent] This relation entails an additional time of %.2f for the smoothing process\n", time_left);

                double new_deadline = seconds() + DESCENT_GRACE_TIME_RATIO * time_left;
                if (new_deadline < deadline) {
                    double delta = DBL_MAX;
                    if (deadline != DBL_MAX)
                        delta = (deadline-new_deadline)/DESCENT_GRACE_TIME_RATIO;
                    deadline = new_deadline;
                    relation_copy(winner, rel);
                    if (time_left == 0) {
                        fprintf(las->output, "# [descent] Yiippee, splitting done\n");
                        /* break both loops */
                        xul = BUCKET_REGION;
                        relation_clear(rel);
                        break;
                    } else if (delta != DBL_MAX) {
                        fprintf(las->output, "# [descent] Improved ETA by %.2f\n", delta);
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
            if (p <= cpoly->pols[side]->lim)
                continue;
            unsigned int n = ULONG_BITS - clzl(p);
            int k = las->hint_lookups[side][n];
            if (k < 0) continue;
            mpz_set_ui(q, p);
            /* Beware, cpoly->m mod p would be wrong ! */
            /* This can't work on 32-bits */
            mpz_set_ui(rho, findroot (winner->a, winner->b, p));
            gmp_fprintf(las->output, "# [descent] "HILIGHT_START"pushing %s (%Zd,%Zd) [%d%c]"HILIGHT_END" to todo list\n", sidenames[side], q, rho, mpz_sizeinbase(q, 2), sidenames[side][0]);
            las_todo_push_withdepth(&(las->todo), q, rho, side, si->doing->depth + 1);
        }
        computeroots(winner);
        for(int i = 0 ; i < winner->nb_ap ; i++) {
            int side = ALGEBRAIC_SIDE;
            unsigned long p = winner->ap[i].p;
            /* See comment above */
            if (p <= cpoly->pols[side]->lim)
                continue;
            unsigned int n = ULONG_BITS - clzl(p);
            int k = las->hint_lookups[side][n];
            if (k < 0) continue;
            mpz_set_ui(q, p);
            mpz_set_ui(rho, winner->ap[i].r);
            gmp_fprintf(las->output, "# [descent] "HILIGHT_START"pushing %s (%Zd,%Zd) [%d%c]"HILIGHT_END" to todo list\n", sidenames[side], q, rho, mpz_sizeinbase(q, 2), sidenames[side][0]);
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

/* {{{ factor_leftover_norm */

#define NMILLER_RABIN 1 /* in the worst case, what can happen is that a
                           composite number is declared as prime, thus
                           a relation might be missed, but this will not
                           affect correctness */
#define IS_PROBAB_PRIME(X) (0 != mpz_probab_prime_p((X), NMILLER_RABIN))
#define BITSIZE(X)      (mpz_sizeinbase((X), 2))
#define NFACTORS        8 /* maximal number of large primes */

/* This function was contributed by Jerome Milan (and bugs were introduced
   by Paul Zimmermann :-).
   Input: n - the number to be factored (leftover norm). Must be composite!
              Assumed to have no factor < B (factor base bound).
          L - large prime bound is L=2^l
   Assumes n > 0.
   Return value:
          -1 if n has a prime factor larger than L
          1 if all prime factors of n are < L
          0 if n could not be completely factored
   Output:
          the prime factors of n are factors->data[0..factors->length-1],
          with corresponding multiplicities multis[0..factors->length-1].
*/
int
factor_leftover_norm (mpz_t n, mpz_array_t* const factors,
                      uint32_array_t* const multis,
                      sieve_info_ptr si, int side)
{
  unsigned int lpb = si->conf->sides[side]->lpb;
  facul_strategy_t *strategy = si->sides[side]->strategy;
  uint32_t i, nr_factors;
  unsigned long ul_factors[16];
  int facul_code;

  /* For the moment this code can't cope with too large factors */
  ASSERT(lpb <= ULONG_BITS);

  factors->length = 0;
  multis->length = 0;

  /* If n < B^2, then n is prime, since all primes < B have been removed */
  if (mpz_cmp (n, si->BB[side]) < 0)
    {
      /* if n > L, return -1 */
      if (mpz_sizeinbase (n, 2) > lpb)
        return -1;

      if (mpz_cmp_ui (n, 1) > 0) /* 1 is special */
        {
          append_mpz_to_array (factors, n);
          append_uint32_to_array (multis, 1);
        }
      return 1;
    }

  /* use the facul library */
  //gmp_printf ("facul: %Zd\n", n);
  facul_code = facul (ul_factors, n, strategy);

  if (facul_code == FACUL_NOT_SMOOTH)
    return -1;

  ASSERT (facul_code == 0 || mpz_cmp_ui (n, ul_factors[0]) != 0);

  /* we use this mask to trap prime factors above bound */
  unsigned long oversize_mask = (-1UL) << lpb;
  if (lpb == ULONG_BITS) oversize_mask = 0;

  if (facul_code > 0)
    {
      nr_factors = facul_code;
      for (i = 0; i < nr_factors; i++)
	{
	  unsigned long r;
	  mpz_t t;
	  if (ul_factors[i] & oversize_mask) /* Larger than large prime bound? */
            return -1;
	  r = mpz_tdiv_q_ui (n, n, ul_factors[i]);
	  ASSERT_ALWAYS (r == 0UL);
	  mpz_init (t);
	  mpz_set_ui (t, ul_factors[i]);
	  append_mpz_to_array (factors, t);
	  mpz_clear (t);
	  append_uint32_to_array (multis, 1); /* FIXME, deal with repeated
						 factors correctly */
	}

      if (mpz_cmp (n, si->BB[side]) < 0)
        {
          if (mpz_sizeinbase (n, 2) > lpb)
            return -1;

          if (mpz_cmp_ui (n, 1) > 0) /* 1 is special */
            {
              append_mpz_to_array (factors, n);
              append_uint32_to_array (multis, 1);
            }
          return 1;
        }

      if (check_leftover_norm (n, si, side) == 0)
        return -1;
    }
  return 0; /* unable to completely factor n */
}
/*}}}*/

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
    __asm__
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
        S[side] = (unsigned char *) malloc_pagealigned(BUCKET_REGION + MEMSET_MIN);
    }
    unsigned char *SS = malloc_pagealigned(BUCKET_REGION);
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
            rep->tn[side] -= seconds ();
            init_rat_norms_bucket_region(S[side], i, si);
#if defined(TRACE_K) 
            if (trace_on_spot_N(w->N))
              fprintf (stderr, "# After init_rat_norms_bucket_region, N=%u S[%u]=%u\n", w->N, trace_Nx.x, S[side][trace_Nx.x]);
#endif
            rep->tn[side] += seconds ();

            /* Apply rational buckets */
            rep->ttbuckets_apply -= seconds();
            for (int j = 0; j < las->nb_threads; ++j)  {
                thread_data_ptr ot = th + j - th->id;
                apply_one_bucket(SS, ot->sides[side]->BA, i, w);
            }
	    SminusS(S[side], S[side] + BUCKET_REGION, SS);
            rep->ttbuckets_apply += seconds();

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
            rep->tn[side] -= seconds ();

            unsigned char * xS = S[side ^ 1];
            init_alg_norms_bucket_region(S[side], xS, i, si);
            rep->tn[side] += seconds ();
#if defined(TRACE_K) 
            if (trace_on_spot_N(w->N))
              fprintf (stderr, "# After init_alg_norms_bucket_region, N=%u S[%u]=%u\n", w->N, trace_Nx.x, S[side][trace_Nx.x]);
#endif

            /* Apply algebraic buckets */
            rep->ttbuckets_apply -= seconds();
            for (int j = 0; j < las->nb_threads; ++j) {
                thread_data_ptr ot = th + j - th->id;
                apply_one_bucket(SS, ot->sides[side]->BA, i, w);
            }
	    SminusS(S[side], S[side] + BUCKET_REGION, SS);
            rep->ttbuckets_apply += seconds();

            /* Sieve small algebraic primes */
            sieve_small_bucket_region(SS, i, s->ssd, ts->ssdpos, si, side, w);
	    SminusS(S[side], S[side] + BUCKET_REGION, SS);
#if defined(TRACE_K) 
            if (trace_on_spot_N(w->N))
              fprintf (stderr, "# Final value on algebraic side, N=%u alg_S[%u]=%u\n", w->N, trace_Nx.x, S[side][trace_Nx.x]);
#endif
        }

        /* Factor survivors */
        rep->ttf -= seconds ();
        rep->reports += factor_survivors (th, i, S, w);
        rep->ttf += seconds ();

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
    free_aligned(SS, BUCKET_REGION, pagesize());

    for(int side = 0 ; side < 2 ; side++) {
        thread_side_data_ptr ts = th->sides[side];
        free(ts->ssdpos);
        free(ts->rsdpos);
        free_aligned(S[side], BUCKET_REGION + MEMSET_MIN, 16);
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
    }
    return thrs;
}/*}}}*/

static void thread_data_free(thread_data * thrs, int n)/*{{{*/
{
    for (int i = 0; i < n ; ++i) {
        las_report_clear(thrs[i]->rep);
    }
    free(thrs);
}/*}}}*/

void thread_pickup_si(thread_data * thrs, sieve_info_ptr si, int n)/*{{{*/
{
    for (int i = 0; i < n ; ++i) {
        thrs[i]->si = si;
        for(int s = 0 ; s < 2 ; s++) {
            thrs[i]->sides[s]->fb_bucket = si->sides[s]->fb_bucket_threads[i];
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
    thread_data_ptr th = thrs[i];
    for(unsigned int side = 0 ; side < 2 ; side++) {
      thread_side_data_ptr ts = th->sides[side];
      uint32_t nb_buckets = thrs[i]->si->nb_buckets;
      uint64_t bucket_size = bucket_misalignment((uint64_t) (thrs[i]->si->sides[side]->max_bucket_fill_ratio * BUCKET_REGION), sizeof(bucket_update_t));
      /* The previous buckets are identical ? */
      if (ts->BA.n_bucket == nb_buckets && ts->BA.bucket_size == bucket_size) {
	/* Yes; so (bucket_write & bucket_read) = bucket_start; nr_logp = 0 */
	re_init_bucket_array(&(ts->BA), &(ts->kBA), &(ts->mBA));
	/* Buckets are ready to be filled */
	continue;
      }
      /* No. We free the buckets, if we have already malloc them. */
      if (ts->BA.n_bucket) clear_bucket_array(&(ts->BA), &(ts->kBA), &(ts->mBA));
      /* We (re)create the buckets */
      if (nb_buckets < THRESHOLD_K_BUCKETS)
	init_bucket_array   (nb_buckets, bucket_size, 255, &(ts->BA), &(ts->kBA), &(ts->mBA));
      else if (nb_buckets < THRESHOLD_M_BUCKETS)
	init_k_bucket_array (nb_buckets, bucket_size, 255, &(ts->BA), &(ts->kBA), &(ts->mBA));
      else
	init_m_bucket_array (nb_buckets, bucket_size, 255, &(ts->BA), &(ts->kBA), &(ts->mBA));
    }
  }
}/*}}}*/

static void thread_buckets_free(thread_data * thrs, unsigned int n)/*{{{*/
{
  for (unsigned int i = 0; i < n ; ++i) {
    thread_side_data_ptr ts;
    ts = thrs[i]->sides[RATIONAL_SIDE];
    clear_bucket_array(&(ts->BA), &(ts->kBA), &(ts->mBA));
    ts = thrs[i]->sides[ALGEBRAIC_SIDE];
    clear_bucket_array(&(ts->BA), &(ts->kBA), &(ts->mBA));
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
    for (int i = 0; i < las->nb_threads; ++i) {
        las_report_accumulate(rep, thrs[i]->rep);
    }
    if (las->verbose) {
        fprintf (las->output, "# ");
      /* fprintf (las->output, "%lu survivors after rational sieve,", rep->survivors0); */
        fprintf (las->output, "%lu survivors after algebraic sieve, ", rep->survivors1);
        fprintf (las->output, "coprime: %lu\n", rep->survivors2);
    }
    gmp_fprintf (las->output, "# %lu relation(s) for %s (%Zd,%Zd)\n", rep->reports, sidenames[si->conf->side], si->doing->p, si->doing->r);
    double qtts = qt0 - rep->tn[0] - rep->tn[1] - rep->ttf;
    if (rep->both_even) {
        fprintf (las->output, "# Warning: found %lu hits with i,j both even (not a bug, but should be very rare)\n", rep->both_even);
    }
    fprintf (las->output, "# Time for this special-q: %1.4fs [norm %1.4f+%1.4f, sieving %1.4f"
            " (%1.4f + %1.4f + %1.4f),"
            " factor %1.4f]\n", qt0,
            rep->tn[RATIONAL_SIDE],
            rep->tn[ALGEBRAIC_SIDE],
            qtts,
            rep->ttbuckets_fill,
            rep->ttbuckets_apply,
            qtts-rep->ttbuckets_fill-rep->ttbuckets_apply,
            rep->ttf);
#if 0   /* incompatible with the todo list */
    rep_bench += rep->reports;
#endif
    las_report_accumulate(report, rep);
    las_report_clear(rep);
}/*}}}*/


static FILE *stats_output = NULL;
/* Print statistics both to the passed-in file handle, and to stats_output 
   if non-NULL. This is to allow scripts to get stats from stderr rather
   than having to decompress the whole relation file. */
static int
print_stats (FILE *file, char* fmt, ...)
{
    int retval;
    char *msg;
    va_list ap;

    va_start(ap, fmt);
    retval = vasprintf(&msg, fmt, ap);
    if (retval != -1) {
        fprintf (file, "%s", msg);
        if (stats_output != NULL)
            fprintf (stats_output, "%s", msg);
        free (msg);
    }
    va_end(ap);

    return retval;
}

/*************************** main program ************************************/


static void declare_usage(param_list pl)
{
  param_list_decl_usage(pl, "poly", "polynomial file");
  param_list_decl_usage(pl, "fb",   "factor base file");
  param_list_decl_usage(pl, "q0",   "left bound of special-q range");
  param_list_decl_usage(pl, "q1",   "right bound of special-q range");
  param_list_decl_usage(pl, "rho",  "sieve only root r mod q0");
  param_list_decl_usage(pl, "v",    "(switch) verbose mode");
  param_list_decl_usage(pl, "out",  "filename where relations are written, instead of stdout");
  param_list_decl_usage(pl, "mt",   "number of threads to use");
  param_list_decl_usage(pl, "ratq", "(switch) use rational special-q");

  param_list_decl_usage(pl, "I",    "set sieving region to 2^I");
  param_list_decl_usage(pl, "skew", "(alias S) skewness");
  param_list_decl_usage(pl, "rlim", "rational factor base bound");
  param_list_decl_usage(pl, "alim", "algebraic factor base bound");
  param_list_decl_usage(pl, "lpbr", "set rational large prime bound to 2^lpbr");
  param_list_decl_usage(pl, "lpba", "set algebraic large prime bound to 2^lpba");
  param_list_decl_usage(pl, "mfbr", "set rational cofactor bound 2^mfbr");
  param_list_decl_usage(pl, "mfba", "set algebraic cofactor bound 2^mfba");
  param_list_decl_usage(pl, "rlambda", "rational lambda value");
  param_list_decl_usage(pl, "alambda", "algebraic lambda value");
  param_list_decl_usage(pl, "rpowlim", "limit on powers on rat side");
  param_list_decl_usage(pl, "apowlim", "limit on powers on alg side");
  param_list_decl_usage(pl, "tdthresh", "trial-divide primes p/r <= ththresh (r=number of roots)");
  param_list_decl_usage(pl, "bkthresh", "bucket-sieve primes p >= bkthresh");

  param_list_decl_usage(pl, "allow-largesq", "(switch) allows large special-q, e.g. for a DL descent");
  param_list_decl_usage(pl, "stats-stderr", "(switch) print stats to stderr in addition to stdout/out file");
  param_list_decl_usage(pl, "todo", "provide file with a list of special-q to sieve instead of qrange");
  param_list_decl_usage(pl, "descent-hint", "hint file ?????");
  param_list_decl_usage(pl, "no-prepare-hints", "(switch) ?????");
  param_list_decl_usage(pl, "mkhint", "(switch) ?????");
}

int main (int argc0, char *argv0[])/*{{{*/
{
    las_info las;
    double t0, tts;
    unsigned long sq = 0;
    int allow_largesq = 0;
    double totJ = 0.0;
    /* following command-line values override those in the polynomial file */
    int argc = argc0;
    char **argv = argv0;
    double max_full = 0.;
#if 0   /* incompatible with the todo list */
    int bench = 0;
    int bench2 = 0;
    double skip_factor = 1.01;  /* next_q = q*skip_factor in bench mode */
    double bench_percent = 1e-3; 
    long bench_tot_rep = 0;
    double bench_tot_time = 0.0;
    const char *statsfilename = NULL;
    const char *sievestatsfilename = NULL;
#endif

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
#if 0   /* incompatible with the todo list */
    param_list_configure_switch(pl, "-bench", &bench);
    param_list_configure_switch(pl, "-bench2", &bench2);
#endif
    param_list_configure_switch(pl, "-mkhint", &create_descent_hints);
    param_list_configure_alias(pl, "-skew", "-S");

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


    if (param_list_parse_switch(pl, "-stats-stderr"))
        stats_output = stderr;

    memset(las, 0, sizeof(las_info));
    las_info_init(las, pl);    /* side effects: prints cmdline and flags */

    /*
    int descent_bootstrap = param_list_lookup_string(pl, "target") != NULL;
    */
    int descent_lower = param_list_lookup_string(pl, "descent-hint") != NULL;

#if 0   /* incompatible with the todo list */
    statsfilename = param_list_lookup_string (pl, "stats");
    sievestatsfilename = param_list_lookup_string (pl, "sievestats");
    param_list_parse_double(pl, "skfact", &skip_factor);
    param_list_parse_double(pl, "percent", &bench_percent);
    param_list_parse_double (pl, "stats_prob", &stats_prob);
#endif

    /* We have the following dependency chain (not sure the account below
     * is exhaustive).
     *
     * q0 -> q0d (double) -> si->sides[*]->{scale,logmax}
     * q0 -> (I, lpb, lambda) for the descent
     * 
     * scales -> logp's in factor base.
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


#if 0   /* incompatible with the todo list */
    /* {{{ stats stuff */
    si->bench=bench + bench2;
    if (statsfilename != NULL) /* a file was given */
      {
        /* if the file exists, we open it in read-mode, otherwise we create
           it */
        stats_file = fopen (statsfilename, "r");
        if (stats_file != NULL)
          stats = 2;
        else
          {
            stats_file = fopen (statsfilename, "w");
            if (stats_file == NULL)
              {
                fprintf (stderr, "Error, cannot create file %s\n",
                         statsfilename);
                exit (EXIT_FAILURE);
              }
            stats = 1;
          }
      }

    if (sievestatsfilename != NULL) /* a file was given */
      {
        sievestats_file = fopen (sievestatsfilename, "w");
        if (sievestats_file == NULL)
          {
            fprintf (stderr, "Error, cannot create file %s\n",
                     sievestatsfilename);
            exit (EXIT_FAILURE);
          }
      }
    if (stats != 0)
      {
        cof_call = (uint32_t**) malloc ((si->cpoly->rat->mfb + 1) * sizeof(uint32_t*));
        cof_succ = (uint32_t**) malloc ((si->cpoly->rat->mfb + 1) * sizeof(uint32_t*));
        for (i = 0; i <= si->cpoly->rat->mfb; i++)
          {
            cof_call[i] = (uint32_t*) malloc ((si->cpoly->alg->mfb + 1)
                                              * sizeof(uint32_t));
            cof_succ[i] = (uint32_t*) malloc ((si->cpoly->alg->mfb + 1)
                                              * sizeof(uint32_t));
            for (j = 0; j <= si->cpoly->alg->mfb; j++)
              cof_call[i][j] = cof_succ[i][j] = 0;
          }
        if (stats == 2)
          {
            fprintf (las->output,
                    "# Use learning file %s with threshold %1.2e\n",
                     statsfilename, stats_prob);
            while (!feof (stats_file))
              {
                uint32_t c, s;
                if (fscanf (stats_file, "%u %u %u %u\n", &i, &j, &c, &s) != 4)
                  {
                    fprintf (stderr, "Error while reading file %s\n",
                             statsfilename);
                    exit (EXIT_FAILURE);
                  }
                if (i <= si->cpoly->rat->mfb && j <= si->cpoly->alg->mfb)
                  {
                    /* When s=0 and c>0, whatever STATS_PROB, we will always
                       have s/c < STATS_PROB, thus (i,j) will be discarded.
                       We allow a small error by considering (s+1)/(c+1)
                       instead. In case s=0, (i,j) is discarded only when
                       1/(c+1) < STATS_PROB (always discarded for c=0). */
                    cof_call[i][j] = c + 1;
                    cof_succ[i][j] = s + 1;
                  }
              }
          }
      }
    int rep_bench = 0;
    int nbq_bench = 0;
    double t_bench = seconds();
    /* }}} */
#endif
    
    thread_data * thrs = thread_data_alloc(las, las->nb_threads);

    las_report report;
    las_report_init(report);

    t0 = seconds ();
    fprintf (las->output, "#\n");

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
                fprintf(las->output, "# Using existing sieving parameters from hint list for q~2^%d on the %s side [%d%c]\n", sc->bitsize, sidenames[sc->side], sc->bitsize, sidenames[sc->side][0]);
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
        if (si->cpoly->pols[si->doing->pside]->degree == 1) {
            /* compute the good rho */
            int n;
            n = poly_roots(&si->doing->r, si->cpoly->pols[si->doing->pside]->f, 1, si->doing->p);
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

        if (SkewGauss (si, si->cpoly->skew) != 0)
            continue;

        /* check |a0|, |a1| < 2^31 if we use fb_root_in_qlattice_31bits */
#ifndef SUPPORT_LARGE_Q
        if (si->a0 <= -2147483648L || 2147483648L <= si->a0 ||
            si->a1 <= -2147483648L || 2147483648L <= si->a1)
          {
            fprintf (stderr, "Error, too large special-q, define SUPPORT_LARGE_Q\n");
            exit (1);
          }
#endif

        /* FIXME: maybe we can discard some special q's if a1/a0 is too large,
         * see http://www.mersenneforum.org/showthread.php?p=130478 
         *
         * Just for correctness at the moment we're doing it in very
         * extreme cases, see bug 15617
         */
        if (sieve_info_adjust_IJ(si, si->cpoly->skew, las->nb_threads) == 0) {
            gmp_fprintf (las->output, "# "HILIGHT_START"Discarding %s q=%Zd; rho=%Zd;"HILIGHT_END" a0=%" PRId64 "; b0=%" PRId64 "; a1=%" PRId64 "; b1=%" PRId64 "; raw_J=%u;\n",
                    sidenames[si->conf->side],
                    si->doing->p, si->doing->r, si->a0, si->b0, si->a1, si->b1,
                    si->J);
            continue;
        }


        gmp_fprintf (las->output, "# "HILIGHT_START"Sieving %s q=%Zd; rho=%Zd;"HILIGHT_END" a0=%" PRId64 "; b0=%" PRId64 "; a1=%" PRId64 "; b1=%" PRId64 ";",
                sidenames[si->conf->side],
                si->doing->p, si->doing->r, si->a0, si->b0, si->a1, si->b1);
        if (si->doing->depth) {
            fprintf(las->output, " # within descent, currently at depth %d", si->doing->depth);
        }
        fprintf(las->output, "\n");
        sq ++;

        /* checks the value of J,
         * precompute the skewed polynomials of f(x) and g(x), and also
         * their floating-point versions */
        sieve_info_update (si, las->nb_threads);
        totJ += (double) si->J;
        if (las->verbose)
            fprintf (las->output, "# I=%u; J=%u\n", si->I, si->J);

        /* Only when tracing. This function gets called once per
         * special-q only. Here we compute the two norms corresponding to
         * the traced (a,b) pair, and start by dividing out the special-q
         * from the one where it should divide */
        trace_per_sq_init(si);

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
            small_sieve_info(las, "small sieve", side, s->ssd);

            small_sieve_extract_interval(s->rsd, s->ssd, s->fb_parts_x->rs);
            small_sieve_info(las, "resieve", side, s->rsd);
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

        /* thread_buckets_free(thrs, las->nb_threads); */

        trace_per_sq_clear(si);
#if 0   /* incompatible with the todo list */
        /* {{{ bench stats */
        if (bench) {
            mpz_t newq0, savq0;
            mpz_init_set(savq0, q0);
            mpz_init_set(newq0, q0);
            {
                mpq_t nq;
                mpq_init(nq);
                mpq_set_d(nq, skip_factor);
                mpz_mul(mpq_numref(nq), mpq_numref(nq), q0);
                mpz_fdiv_q(newq0, mpq_numref(nq), mpq_denref(nq));
                mpq_clear(nq);
            }
            // print some estimates for special-q's between q0 and the next
            int nb_q = 1;
            do {
                mpz_nextprime(q0, q0);
                nb_q ++;
            } while (mpz_cmp(q0, newq0) < 0);
            mpz_set(q0, newq0);
            nroots=0;
            t_bench = seconds() - t_bench;
            gmp_fprintf(las->output,
                    "# Stats for q=%Zd: %d reports in %1.1f s\n",
                    savq0, rep_bench, t0);
            fprintf(las->output,
                    "# Estimates for next %d q's: %d reports in %1.0f s, %1.2f s/r\n",
                    nb_q, nb_q*rep_bench, t_bench*nb_q, t_bench/((double)rep_bench));
            bench_tot_time += t_bench*nb_q;
            bench_tot_rep += nb_q*rep_bench;
            rep_bench = 0;
            fprintf(las->output, "# Cumulative (estimated): %lu reports in %1.0f s, %1.2f s/r\n",
                    bench_tot_rep, bench_tot_time,
                    (double) bench_tot_time / (double) bench_tot_rep);
            t_bench = seconds();
            mpz_clear(newq0);
            mpz_clear(savq0);
        }
        /* }}} */
        /* {{{ bench stats */
        if (bench2) {
            nbq_bench++;
            const int BENCH2 = 50;
            if (rep_bench >= BENCH2) {
                t_bench = seconds() - t_bench;
                fprintf(las->output,
                        "# Got %d reports in %1.1f s using %d specialQ\n",
                        rep_bench, t_bench, nbq_bench);
                double relperq = (double)rep_bench / (double)nbq_bench;
                double est_rep = (double)rep_bench;
                do {
                    mpz_nextprime(q0, q0);
                    est_rep += relperq;
                } while (est_rep <= BENCH2 / bench_percent);
                gmp_fprintf(las->output,
                        "# Extrapolate to %ld reports up to q = %Zd\n",
                        (long) est_rep, q0);
                bench_tot_time += t_bench / bench_percent;
                bench_tot_rep += BENCH2 / bench_percent;
                fprintf(las->output,
                        "# Cumulative (estimated): %lu reports in %1.0f s, %1.2f s/r\n",
                        bench_tot_rep, bench_tot_time,
                        (double) bench_tot_time / (double) bench_tot_rep);
                // reinit for next slice of bench:
                t_bench = seconds();
                nbq_bench = 0;
                rep_bench = 0;
                nroots=0;
            }
        }
        /* }}} */
#endif
      } // end of loop over special q ideals.

    thread_buckets_free(thrs, las->nb_threads);

    if (descent_lower) {
        fprintf(las->output, "# Now displaying again the results of all descents\n");
        las_descent_helper_display_all_trees(las->descent_helper, las->output);
        las_descent_helper_free(las->descent_helper);
    }

    t0 = seconds () - t0;
    print_stats (las->output, "# Average J=%1.0f for %lu special-q's, max bucket fill %f\n",
            totJ / (double) sq, sq, max_full);
    tts = t0;
    tts -= report->tn[0];
    tts -= report->tn[1];
    tts -= report->ttf;
    if (las->verbose)
        facul_print_stats (las->output);
#if 0
    /* {{{ stats */
    if (sievestats_file != NULL)
    {
        fprintf (sievestats_file, "# Number of sieve survivors and relations by sieve residue pair\n");
        fprintf (sievestats_file, "# Format: S1 S2 #relations #survivors ratio\n");
        fprintf (sievestats_file, "# where S1 is the sieve residue on the rational side, S2 algebraic side\n");
        fprintf (sievestats_file, "# Make a pretty graph with gnuplot:\n");
        fprintf (sievestats_file, "# splot \"sievestatsfile\" using 1:2:3 with pm3d\n");
        fprintf (sievestats_file, "# plots histogram for relations, 1:2:4 for survivors, 1:2:($3/$4) for ratio\n");
        for(int i1 = 0 ; i1 < 256 ; i1++) {
            for (int i2 = 0; i2 < 256; i2++) {
                unsigned long r1 = report->report_sizes[i1][i2];
                unsigned long r2 = report->survivor_sizes[i1][i2];
                if (r1 > r2) {
                    fprintf(stderr, "Error, statistics report more relations (%lu) than "
                            "sieve survivors (%lu) for (%d,%d)\n", r1, r2, i1, i2);
                }
                if (r2 > 0)
                    fprintf (sievestats_file, "%d %d %lu %lu\n", 
                            i1, i2, r1, r2);
            }
            fprintf (sievestats_file, "\n");
        }
        fprintf (sievestats_file, "# ");
        fclose(sievestats_file);
        sievestats_file = NULL;
    }
    /* }}} */
#endif
    /*{{{ Display tally */
    if (las->nb_threads > 1) 
        print_stats (las->output, "# Total wct time %1.1fs [precise timings available only for mono-thread]\n", t0);
    else
        print_stats (las->output, "# Total time %1.1fs [norm %1.2f+%1.1f, sieving %1.1f"
                " (%1.1f + %1.1f + %1.1f),"
                " factor %1.1f]\n", t0,
                report->tn[RATIONAL_SIDE],
                report->tn[ALGEBRAIC_SIDE],
                tts,
                report->ttbuckets_fill,
                report->ttbuckets_apply,
                tts-report->ttbuckets_fill-report->ttbuckets_apply,
                report->ttf);

    print_stats (las->output, "# Total %lu reports [%1.3gs/r, %1.1fr/sq]\n",
            report->reports, t0 / (double) report->reports,
            (double) report->reports / (double) sq);
    
#if 0
    fprintf (stderr, "rat:");
    for (int i = 0; i < 256; i++)
      if (cof_calls[RATIONAL_SIDE][i] > 0)
        fprintf (stderr, " %d:%1.0f/%1.0f", i, cof_fails[RATIONAL_SIDE][i],
                 cof_calls[RATIONAL_SIDE][i]);
    fprintf (stderr, "\n");
    fprintf (stderr, "alg:");
    for (int i = 0; i < 256; i++)
      if (cof_calls[ALGEBRAIC_SIDE][i] > 0)
        fprintf (stderr, " %d:%1.0f/%1.0f", i, cof_fails[ALGEBRAIC_SIDE][i],
                 cof_calls[ALGEBRAIC_SIDE][i]);
    fprintf (stderr, "\n");
#endif
    /*}}}*/
#if 0   /* incompatible with the todo list */
    /* {{{ stats */
    if (bench || bench2) {
        fprintf(las->output, "# Total (estimated): %lu reports in %1.1f s\n",
                bench_tot_rep, bench_tot_time);
    }
    if (bucket_prime_stats) 
    {
        printf ("# Number of bucket primes: %ld\n", nr_bucket_primes);
        printf ("# Number of divisibility tests of bucket primes: %ld\n", 
                nr_div_tests);
        printf ("# Number of compositeness tests of bucket primes: %ld\n", 
                nr_composite_tests);
        printf ("# Number of wrapped composite values while dividing out "
                "bucket primes: %ld\n", nr_wrap_was_composite);
      }
    if (stats == 2)
      fprintf (las->output, "# Rejected %u cofactorizations out of %u due to stats file\n", cof_call[0][0] - cof_succ[0][0], cof_call[0][0]);
    /* }}} */
#endif

    thread_data_free(thrs, las->nb_threads);

    las_report_clear(report);

    las_info_clear(las);

    param_list_clear(pl);

#if 0   /* incompatible with the todo list */
    if (stats != 0) /* {{{ */
      {
        for (i = 0; i <= si->cpoly->rat->mfb; i++)
          {
            if (stats == 1)
              for (j = 0; j <= si->cpoly->alg->mfb; j++)
                fprintf (stats_file, "%u %u %u %u\n", i, j, cof_call[i][j],
                         cof_succ[i][j]);
            free (cof_call[i]);
            free (cof_succ[i]);
          }
        free (cof_call);
        free (cof_succ);
        fclose (stats_file);
      }
    /* }}} */
#endif
    return 0;
}/*}}}*/
