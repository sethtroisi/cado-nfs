#define _XOPEN_SOURCE 600       // should be defined before others
                                // in order to have posix_memalign
                                // on old systems.
                                // This is used in bucket.h
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <inttypes.h> /* for PRIx64 macro and strtoumax */
#include <math.h>   // for ceiling, floor in cfrac
#include "cado.h"
#include "../utils/mod_ul.h"
#include "fb.h"
#include "../utils/utils.h"
#include "../utils/manu.h"   /* for ctzl */
#include "basicnt.h" /* for gcd_ul */
#include <tifa.h>
#include "bucket.h"

/* comment to use classical sieving */
#define USE_BUCKETS

/* Trick to discard lognorms that will probably lead to non L-smooth
   cofactors. Disabled for now since it requires accurate sieving
   (in particular by prime powers and bad primes). */
// #define COFACTOR_TRICK

/* default sieve region side is 2^DEFAULT_I */
#define DEFAULT_I 12

/* optimal bit-size of the product of norms in Bernstein's algorithm */
#define BERNSTEIN_THRESHOLD 1000000

#define LOG_SCALE 1.4426950408889634 /* 1/log(2) to 17 digits, rounded to
                                        nearest. This is enough to uniquely
                                        identify the corresponding IEEE 754
                                        double precision number */

/* Guard for the logarithms of norms, so that the value does not wrap around
   zero due to roundoff errors. */
#define GUARD 4

/* GUARD+LOG_MAX should be as near as possible from 256, to get more accuracy
   in the norm computations, but not too much, otherwise a norm might be
   rounded to zero. */
#define LOG_MAX (255.9 - (double) GUARD)

/* define TRACE_K to -I/2+i+I*j to trace the sieve array entry corresponding
   to (i,j), i.e., a = a0*i+a1*j, b = b0*i+b1*j */
// #define TRACE_K 5418470

#define MAX_BUF_SIZE 512 /* maximal length of the factorization on each side */

/* concatenate an integer q at the end of s (which might be empty) */
#define sprintf_cat(s,q)                                \
  do                                                    \
    {                                                   \
      if (strlen (s) == 0)                              \
        sprintf (s, "%" PRIx64, q);                     \
      else                                              \
        sprintf (s + strlen (s), ",%" PRIx64, q);       \
    }                                                   \
  while (0)

/* same with gmp_sprintf */
#define gmp_sprintf_cat(s,q)                            \
  do                                                    \
    {                                                   \
      if (strlen (s) == 0)                              \
        gmp_sprintf (s, "%Zx", q);                      \
      else                                              \
        gmp_sprintf (s + strlen (s), ",%Zx", q);        \
    }                                                   \
  while (0)

/* uintmax_t is guaranteed to be larger or equal to uint64_t */
#define strtouint64(nptr,endptr,base) (uint64_t) strtoumax(nptr,endptr,base)

/* check is a "large" prime is really large, i.e., not in the factor base */
#define CHECK_LARGE_PRIME(p,lim,s)                                      \
  if (mpz_cmp_ui (p, lim) <= 0)                                         \
    gmp_fprintf (stderr, "# WARNING: factor base prime %Zd was missed " \
                 "on %s. side\n", p, s)

// General information about the siever
typedef struct {
    // sieving area
    uint32_t I;
    uint32_t J;
    int logI; // such that I = 1<<logI
    // description of the q-lattice
    uint64_t q;
    uint64_t rho;
    int32_t a0, b0, a1, b1;
    // parameters for bucket sieving
    int bucket_region;    // should be around L1 cache size, a power of 2
                        // and a multiple of I.
    int log_bucket_region;    
    int nb_buckets;
    int bucket_limit;   // maximal number of bucket_reports allowed in one bucket.
    unsigned int degree;   /* polynomial degree */
    double scale_alg; /* LOG_MAX / log (max |F(a0*i+a1*j, b0*i+b1*j)|) */
    double scale_rat; /* LOG_MAX / log (max |G(a0*i+a1*j, b0*i+b1*j)|) */
    mpz_t *fij;       /* coefficients of F(a0*i+a1*j, b0*i+b1*j) */
    double B;         /* bound for the norm computation */
    unsigned char alg_Bound[256]; /* non-zero for good lognorms */
    unsigned char rat_Bound[256]; /* non-zero for good lognorms */
    int checknorms;          /* if non-zero, completely factor the potential
                                relations */
    uint32_t *abadprimes;     /* primes which may appear on the algebraic side
                                 but are not in the factor base (end of list
                                 is 0) */
    uint32_t *rbadprimes;     /* primes which may appear on the rational side
                                 but are not in the factor base (end of list
                                 is 0) */
} sieve_info_t;

/************************** sieve info stuff *********************************/

static double get_maxnorm (cado_poly, sieve_info_t *, uint64_t);

/* initialize array C[0..255]: C[i] is non-zero whenever the log-norm i
   is considered as a potential report.
   The large prime bound is L = 2^l.
*/
static void
sieve_info_init_lognorm (unsigned char *C, unsigned char threshold
#ifdef COFACTOR_TRICK
                         , unsigned long B, unsigned long l, double scale
#endif
)
{
  unsigned long k;

  for (k = 0; k < 256; k++)
    C[k] = (k <= threshold);

#ifdef COFACTOR_TRICK
  {
    unsigned char k0, k1;
    /* for L < R <= B^2, R cannot be L-smooth (see lattice.tex) */
    k0 = (unsigned char) ((double) l * scale + 0.5) + GUARD;
    k1 = (unsigned char) (2.0 * log((double) B) * LOG_SCALE * scale + 0.5)
      + GUARD;
    for (k = k0 + GUARD; k + GUARD <= k1 && k < 256; k++)
      C[k] = 0;
  }
#endif
}

static void
sieve_info_init (sieve_info_t *si, cado_poly cpoly, int I, uint64_t q0)
{
  unsigned int d = cpoly->degree;
  unsigned int k;
  double r, scale;
  unsigned char alg_bound, rat_bound;

  si->degree = d;
  si->fij = malloc ((d + 1) * sizeof (mpz_t));
  for (k = 0; k <= d; k++)
    mpz_init (si->fij[k]);
  si->logI = I;
  si->I = 1 << si->logI;
  si->J = 1 << (si->logI - 1);

  /* initialize bounds for the norm computation, see lattice.tex */
  si->B = sqrt (2.0 * (double) q0 / (cpoly->skew * sqrt (3.0)));
  si->scale_alg = get_maxnorm (cpoly, si, q0); /* log(max norm) */
  fprintf (stderr, "# Alg. side: log(maxnorm)=%1.6f logbase=%1.6f",
           si->scale_alg, exp (si->scale_alg / LOG_MAX));
  si->scale_alg = LOG_MAX / si->scale_alg;
  scale = si->scale_alg / LOG_SCALE;
  alg_bound = (unsigned char) (cpoly->alambda * (double) cpoly->lpba * scale)
            + GUARD;
  fprintf (stderr, " bound=%u\n", alg_bound);
  sieve_info_init_lognorm (si->alg_Bound, alg_bound
#ifdef COFACTOR_TRICK
                           , cpoly->alim, cpoly->lpba, scale
#endif
                           );

  /* similar bound on the rational size: |a| <= s*I*B and |b| <= I*B */
  si->scale_rat = fabs (mpz_get_d (cpoly->g[1])) * cpoly->skew
                + fabs (mpz_get_d (cpoly->g[0]));
  si->scale_rat *= si->B * (double) si->I;
  si->scale_rat = log (si->scale_rat);
  /* on the rational side, we want that the non-reports on the algebraic
     side, which are set to 255, remain over the report bound R, even if
     the rational norm is totally smooth. For this, we simply add R to the
     maximal lognorm to compute the log base */
  r = cpoly->rlambda * (double) cpoly->lpbr; /* base-2 logarithm of the
                                                report bound */
  fprintf (stderr, "# Rat. side: log(maxnorm)=%1.6f ", si->scale_rat);
  si->scale_rat += r / LOG_SCALE;   /* add natural logarithm */
  fprintf (stderr, "logbase=%1.6f", exp (si->scale_rat / (LOG_MAX - GUARD)));
  /* we subtract again GUARD to avoid that non-reports overlap the report
     region due to roundoff errors */
  si->scale_rat = (LOG_MAX - GUARD) / si->scale_rat;
  scale = si->scale_rat / LOG_SCALE;
  rat_bound = (unsigned char) (r * scale) + GUARD;
  fprintf (stderr, " bound=%u\n", rat_bound);
  sieve_info_init_lognorm (si->rat_Bound, rat_bound
#ifdef COFACTOR_TRICK
                           , cpoly->rlim, cpoly->lpbr, scale
#endif
                           );

#ifdef USE_BUCKETS
  // bucket info
  // TODO: be more clever, here.
  si->log_bucket_region = 15;
  si->bucket_region = 1<<si->log_bucket_region;
  si->nb_buckets = 1 + ((si->I*si->J) / si->bucket_region);
  si->bucket_limit = 5*si->bucket_region; // This factor of 5 is huge!!!
  fprintf(stderr, "# log_bucket_region = %u\n", si->log_bucket_region);
  fprintf(stderr, "# bucket_region = %u\n", si->bucket_region);
  fprintf(stderr, "# nb_buckets = %u\n", si->nb_buckets);
  fprintf(stderr, "# bucket_limit = %u\n", si->bucket_limit);
#endif
}

static void
compute_badprimes (sieve_info_t *si, cado_poly cpoly,
                   factorbase_degn_t *fb_alg, factorbase_degn_t *fb_rat)
{
  unsigned long p;
  unsigned long l; /* number of bad primes */

  l = 0;
  si->abadprimes = (uint32_t*) malloc (sizeof (uint32_t));
  fprintf (stderr, "# Alg. bad primes:");
  for (p = 2; p <= cpoly->alim; p = getprime (p))
    {
      if (fb_alg->p != FB_END && fb_alg->p < p)
        fb_alg = fb_next (fb_alg);
      /* invariant: p <= fb_alg->p */
      if (mpz_divisible_ui_p (cpoly->f[cpoly->degree], p))
        {
          l ++;
          si->abadprimes = (uint32_t*) realloc (si->abadprimes,
                                                (l + 1) * sizeof (uint32_t));
          si->abadprimes[l - 1] = p;
          fprintf (stderr, " %lu", p);
        }
    }
  si->abadprimes[l] = 0; /* end of list marker */
  fprintf (stderr, "\n");
  getprime (0);

  l = 0;
  si->rbadprimes = (uint32_t*) malloc (sizeof (uint32_t));
  fprintf (stderr, "# Rat. bad primes:");
  for (p = 2; p <= cpoly->rlim; p = getprime (p))
    {
      if (fb_rat->p != FB_END && fb_rat->p < p)
        fb_rat = fb_next (fb_rat);
      /* invariant: p <= fb_rat->p */
      if (p != fb_rat->p && mpz_divisible_ui_p (cpoly->g[1], p))
        {
          l ++;
          si->rbadprimes = (uint32_t*) realloc (si->rbadprimes,
                                                (l + 1) * sizeof (uint32_t));
          si->rbadprimes[l - 1] = p;
          fprintf (stderr, " %lu", p);
        }
    }
  si->rbadprimes[l] = 0; /* end of list marker */
  fprintf (stderr, "\n");
  getprime (0);
}

static void
sieve_info_update (sieve_info_t *si, double skew)
{
  double s_over_a1, one_over_b1;

  /* check J */
  s_over_a1 = fabs (skew / (double) si->a1);
  one_over_b1 = fabs (1.0 / (double) si->b1);
  if (one_over_b1 < s_over_a1)
    s_over_a1 = one_over_b1; /* min(s/|a1|, 1/|b1|) */
  s_over_a1 *= si->B;
  if (s_over_a1 > 1.0)
    s_over_a1 = 1.0;
  si->J = (uint32_t) (s_over_a1 * (double) (si->I >> 1));
  fprintf (stderr, "# I=%u; J=%u\n", si->I, si->J);
  
}

static void
sieve_info_clear (sieve_info_t *si)
{
  unsigned int d = si->degree;
  unsigned int k;

  for (k = 0; k <= d; k++)
    mpz_clear (si->fij[k]);
  free (si->fij);
  free (si->abadprimes);
}

/************************** sieve reports stuff ******************************/

static void eval_fij (mpz_t, mpz_t *, unsigned int, long, unsigned long);

/* sieve reports */
typedef struct {
  int64_t a;
  uint64_t b;
  uint64_t q;          /* special-q (algebraic side for now) */
  mpz_t normr;         /* norm on rational side */
  mpz_t norma;         /* norm on algebraic side (after dividing by q) */
} sieve_report_t;

static void
sieve_report_init (sieve_report_t *r)
{
  mpz_init (r->normr);
  mpz_init (r->norma);
}

static void
sieve_report_clear (sieve_report_t *r)
{
  mpz_clear (r->normr);
  mpz_clear (r->norma);
}

static void
sieve_report_set (sieve_report_t *r, int64_t a, uint64_t b, uint64_t q,
                  cado_poly cpoly)
{
  r->a = a;
  r->b = b;
  r->q = q;
  eval_fij (r->norma, cpoly->f, cpoly->degree, a, b);
  if (q <= (uint64_t) ULONG_MAX)
    mpz_divexact_ui (r->norma, r->norma, q);
  else
    {
      mpz_set_uint64 (r->normr, q);
      mpz_divexact (r->norma, r->norma, r->normr);
    }
  eval_fij (r->normr, cpoly->g, 1, a, b);
}

typedef struct {
  unsigned long alloc;        /* number of allocated entries */
  unsigned long size;         /* number of sieve reports (size <= alloc) */
  unsigned long r_bitsize;    /* cumulated bit-size on rational side */
  unsigned long a_bitsize;    /* cumulated bit-size on algebraic side */
  sieve_report_t *report;     /* pointer to sieve reports */
} sieve_report_list_t;

static void
sieve_report_list_init (sieve_report_list_t *L)
{
  L->report = NULL;
  L->alloc = 0;
  L->size = 0;
  L->r_bitsize = 0;
  L->a_bitsize = 0;
}

static void
sieve_report_list_clear (sieve_report_list_t *L)
{
  unsigned long i;

  for (i = 0; i < L->size; i++)
    sieve_report_clear (L->report + i);
  free (L->report);
}

static void
sieve_report_list_reset (sieve_report_list_t *L)
{
  L->size = 0;
  L->r_bitsize = 0;
  L->a_bitsize = 0;
}

static void
sieve_report_list_add (sieve_report_list_t *L, int64_t a, uint64_t b,
                       uint64_t q, cado_poly cpoly)
{
  sieve_report_t *r;

  if (L->size == L->alloc)
    {
      L->alloc ++;
      L->report = (sieve_report_t*) realloc (L->report,
                                           L->alloc * sizeof (sieve_report_t));
      sieve_report_init (L->report + L->size);
    }
  r = L->report + L->size;
  sieve_report_set (r, a, b, q, cpoly);
  L->r_bitsize += mpz_sizeinbase (r->normr, 2);
  L->a_bitsize += mpz_sizeinbase (r->norma, 2);
  L->size ++;
}

/*****************************************************************************/


unsigned long invmod(unsigned long a, unsigned long b) {
  unsigned long u, v, fix;
  int t, lsh;

  fix = (b+1)>>1;

  assert (a < b);

  u = 1; v = 0; t = 0;
  do {
    do {
      b -= a; v += u;
      lsh = ctzl(b);
      b >>= lsh;
      t += lsh;
      u <<= lsh;
    } while (a<b);
    if (a == b)
      break;
    do {
      a -= b; u += v;
      lsh = ctzl(a);
      a >>= lsh;
      t += lsh;
      v <<= lsh;
    } while (b < a);
  } while (a != b);

  while (t>0) {
    unsigned long sig = u & 1UL;
    u >>= 1;
    if (sig)
      u += fix;
    --t;
  }
  return u;
}

// Compute the root r describing the lattice inside the q-lattice
// corresponding to the factor base prime (p,R).
// Formula: r = - (a1-R*b1)/(a0-R*b0) mod p
// In the case where denominator is zero, returns p.
// Otherwise r in [0,p-1]
fbprime_t
fb_root_in_qlattice(const fbprime_t p, const fbprime_t R, const sieve_info_t * si)
{
    int64_t aux;
    uint64_t num, den;

    // numerator
    aux = ((int64_t)R)*((int64_t)si->b1) - ((int64_t)si->a1); 
    if (aux >= 0)
        num = ((uint64_t)aux) % ((uint64_t)p);
    else {
        num = ((uint64_t)(-aux)) % ((uint64_t)p);
        num = (uint64_t)p - num;
    }
    if (num == 0)
        return 0;

    // denominator
    aux = ((int64_t)si->a0) - ((int64_t)R)*((int64_t)si->b0); 
    if (aux >= 0) {
        den = ((uint64_t)aux) % ((uint64_t)p);
        if (den == 0)
            return p;
    } else {
        den = ((uint64_t)(-aux)) % ((uint64_t)p);
        if (den == 0)
            return p;
        den = (uint64_t)p - den;
    }

    // divide
    den = invmod(den, p);
    num = num*den;
    return (fbprime_t)(num % ((uint64_t) p));
}


// Here, we have i/j == 0 mod p, so that j can be anything, and i is 0
// mod p.
static void
special_case_0(unsigned char *S, fbprime_t p, unsigned char logp, const sieve_info_t *si)
{
    long i, i0;
    unsigned long j;
    const uint32_t I = si->I;
    const long Is2 = (long)(I>>1);
    unsigned char *S_ptr;

#ifdef VERBSOE
    fprintf(stderr, "# Entering special_case_0(%u)\n", p);
#endif
    i0 = -(long)(I>>1) + ((unsigned long)(Is2) % p);
    S_ptr = S + (I>>1);
    for (j = 0; j < si->J; ++j) {
        for (i = i0; i < Is2; i += p)
          {
            S_ptr[i] -= logp;
#ifdef TRACE_K
            if ((S_ptr - S) + i == TRACE_K)
              fprintf (stderr, "Subtracted %u (%u) to S[%ld], remains %u\n",
                       logp, p, (S_ptr - S) + i, S_ptr[i]);
#endif
          }
        S_ptr += I;
    }
}

// Here, we have i/j == oo mod p, so that j is 0 mod p and i can be anything
static void
special_case_p(unsigned char *S, fbprime_t p, unsigned char logp, const sieve_info_t *si)
{
    long i;
    unsigned long j;
    const long Is2 = (long)(si->I>>1);
    const uint32_t J = si->J;
    unsigned char *S_ptr;

#ifdef VERBOSE
    fprintf(stderr, "# Entering special_case_p(%u)\n", p);
#endif
    S_ptr = S + Is2;
    for (j = 0; j < J; j += p) {
        for (i = -Is2; i < Is2; ++i)
          {
            S_ptr[i] -= logp;
#ifdef TRACE_K
            if ((S_ptr - S) + i  == TRACE_K)
              fprintf (stderr, "Subtracted %u (%u) to S[%ld], remains %u\n",
                       logp, p, (S_ptr - S) + i, S_ptr[i]);
#endif
          }
        S_ptr += p*si->I;
    }
}

#if 0
static void
sieve_slow (unsigned char *S, const factorbase_degn_t *fb,
        const sieve_info_t *si)
{
    fbprime_t p;
    fbprime_t r, R;
    unsigned char logp;
    unsigned char *S_ptr;

    // Loop over all primes in the factor base.
    while (fb->p != FB_END) {
        unsigned char nr;
        p = fb->p;
        logp = fb->plog;
        unsigned long Is2modp = (unsigned long)(si->I>>1) % p;

        for (nr = 0; nr < fb->nr_roots; ++nr) {
            R = fb->roots[nr];
            r = fb_root_in_qlattice(p, R, si);
            if (r == 0) {
                special_case_0 (p);
                continue;
            } 
            if (r == p) {
                special_case_p (p);
                continue;
            } 

            unsigned long j;
            for (j = 0; j < si->J; ++j) {
                const uint32_t I = si->I;
                long i0;
                // init i0 = -I/2 + ( (rj+I/2) mod p)
                {
                    modulus_t m;
                    residueul_t x, rr, jj;

                    // TODO:
                    // Assuming that p < 2^32 and we are on 64bit
                    // machine, this can be improved a lot! No need to
                    // use mod_ul lib.
                    mod_initmod_ul(m, p);
                    modul_initmod_ul(rr, r);
                    if (j >= p) 
                        modul_initmod_ul(jj, j % p);
                    else
                        modul_initmod_ul(jj, j);
                    modul_mul(x, rr, jj, m);            // here there is a modulo p
                    modul_add_ul(x, x, Is2modp, m);
                    i0 = (long)modul_get_ul(x, m) - (long)(I>>1) ;
                }

                S_ptr = S + (j*I + (I>>1));
                // sieve one j-line
                for (; i0 < (I>>1); i0 += p)
                    S_ptr[i0] -= logp;
            }
        }
        fb = fb_next (fb); // cannot do fb++, due to variable size !
    }
}
#endif


/*
 * Algorithm by Franke and Kleinjung for lattice sieving of largish
 * primes.
 */

typedef struct {
    int32_t alpha, beta, gamma, delta;  // coordinates of the basis
    uint32_t a, c;       // sieving offsets in x-coordinate
    uint32_t b0, b1;     // thresholds for the branch inside siever
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
//
// Return value:
// * non-zero if everything worked ok
// * zero when the algorithm failed. This can happen when p is a prime power,
//   and g, gcd(p,r) >= I, since then the subtractive Euclidean algorithm will
//   yield (a0=g, b0=0) at some point --- or the converse --- and the loop
//   while (|a0| >= I) a0 += b0 will loop forever.
//
static int
reduce_plattice (plattice_info_t *pli, const fbprime_t p, const fbprime_t r,
                 const sieve_info_t * si)
{
    int64_t a0, a1, b0, b1, I, J;
    I = si->I;
    J = si->J;
    a0 = -((int64_t)p); a1 = 0;
    b0 = r;  b1 = 1;

    /* subtractive variant of Euclid's algorithm */
    while ( b0 >= I )
    {
      /* a0 < 0 < b0 < -a0 */
        do {
            a0 += b0;
            a1 += b1;
        } while (a0 + b0 <= 0);
        /* -b0 < a0 <= 0 < b0 */
        if (-a0 < I)
          {
            if (a0 == 0)
              return 0;
            while (b0 >= I)
              {
                b0 += a0;
                b1 += a1;
              }
            goto case_even;
          }
        /* -b0 < a0 < 0 < b0 */
        do {
            b0 += a0;
            b1 += a1;
        } while (b0 + a0 >= 0);
        /* a0 < 0 <= b0 < -a0 */
    }
    if (b0 == 0)
      return 0;
    while (a0 <= -I)
      {
        a0 += b0;
        a1 += b1;
      }
 case_even:
    pli->alpha = (int32_t) a0;
    pli->beta = (int32_t) a1;
    pli->gamma = (int32_t) b0;
    pli->delta = (int32_t) b1;

    assert (pli->beta > 0);
    assert (pli->delta > 0);
    assert ((pli->alpha <= 0) && (pli->alpha > -I));
    assert ((pli->gamma >= 0) && (pli->gamma < I));
    assert (pli->gamma-pli->alpha >= I);

    // WARNING: Here, we assume a lot on a bound on I,J
    // TODO: clean these bound problems
    int64_t aa = ((int64_t)pli->beta)*I + (int64_t)(pli->alpha);
    if (aa > I*J)
        pli->a = (uint32_t)(INT32_MAX/2);
    else
        pli->a = (uint32_t)aa;
    int64_t cc = ((int64_t)pli->delta)*I + (int64_t)(pli->gamma);
    if (cc > I*J)
        pli->c = (uint32_t)(INT32_MAX/2);
    else
        pli->c = (uint32_t)cc;
    pli->b0 = -pli->alpha;
    pli->b1 = I - pli->gamma;

    return 1; /* algorithm was ok */
}

static void line_sieve(unsigned char *S, factorbase_degn_t **fb_ptr, 
        const sieve_info_t *si)
{
    const uint32_t I = si->I;
    double tm = seconds();
    unsigned char *S_ptr;
    fbprime_t p, r, R;
    unsigned char logp;
    while ((*fb_ptr)->p != FB_END && (*fb_ptr)->p <= I) {
        unsigned char nr;
        p = (*fb_ptr)->p;
        logp = (*fb_ptr)->plog;
        unsigned long Is2modp = (unsigned long)(I>>1) % p;

        for (nr = 0; nr < (*fb_ptr)->nr_roots; ++nr) {
            R = (*fb_ptr)->roots[nr];
            r = fb_root_in_qlattice(p, R, si);
            if (r == 0) {
                special_case_0(S, p, logp, si);
                continue;
            } 
            if (r == p) {
                special_case_p(S, p, logp, si);
                continue;
            } 
            
            { 
                // line sieving
                unsigned long j;
                long ii0 = Is2modp; // this will be (rj+I/2) mod p
                for (j = 0; j < si->J; ++j) {
                    long i0;
                    // init i0 = -I/2 + ( (rj+I/2) mod p)
                    i0 = ii0 - (long)(I>>1);
                    S_ptr = S + (j*I + (I>>1));
                    // sieve one j-line
                    for (; i0 < (long) (I>>1); i0 += p)
                      {
                        S_ptr[i0] -= logp;
#ifdef TRACE_K
                        if ((S_ptr - S) + i0 == TRACE_K)
                          fprintf (stderr, "Subtracted %u (%u) to S[%ld], remains %u\n",
                                   logp, p, (S_ptr - S) + i0, S_ptr[i0]);
#endif
                      }
                    // update starting point
                    ii0 += r;
                    if (ii0 >= (long) p)
                        ii0 -= p;

                }
            }
        }
        (*fb_ptr) = fb_next ((*fb_ptr)); // cannot do fb++, due to variable size !
    }
    fprintf (stderr, "# small primes sieved in %f sec\n", seconds()-tm);
}

#ifdef TRY_SSE
typedef uint32_t v4si __attribute__ ((vector_size (16)));

// r = (x < y) ? a : b
static inline void
xmm_cmovlt(v4si *r, v4si x, v4si y, v4si a, v4si b) {
    v4si tmp;
    tmp = __builtin_ia32_pcmpgtd128(y, x); // tmp = x<y ? 0xffff..fff. : 0
    a = __builtin_ia32_pand128(a, tmp);  // a = x<y ? a : 0
    tmp = __builtin_ia32_pandn128(tmp, b); // tmp = x<y ? 0 : b
    *r = __builtin_ia32_por128(a, tmp); 
}

// r = (x >= y) ? a : b
static inline void
xmm_cmovge(v4si *r, v4si x, v4si y, v4si a, v4si b) {
    v4si tmp0, tmp1;
    tmp0 = __builtin_ia32_pcmpgtd128(x, y); // tmp0 = x > y ? 0xffff..fff. : 0
    tmp1 = __builtin_ia32_pcmpeqd128(x, y); // tmp1 = x == y ? 0xffff..fff. : 0
    tmp0 |= tmp1;  // tmp0 = x >= y ? 0xffff..fff. : 0
    a = __builtin_ia32_pand128(a, tmp0);  
    tmp0 = __builtin_ia32_pandn128(tmp0, b); 
    *r = __builtin_ia32_por128(a, tmp0); 
}
#endif

#ifndef USE_BUCKETS
// This version of the sieve implements the following:
//   - naive line-sieving for small p
//   - Franke-Kleinjung lattice-sieving for large p
//   - write directly the reports in the sieve-array (no buckets)
// WARNING: STILL AT WORK!
static void
sieve_random_access (unsigned char *S, factorbase_degn_t *fb,
        const sieve_info_t *si)
{
    // Sieve primes <= I.
    line_sieve(S, &fb, si);

    double tm = seconds();
#ifdef TRY_SSE
    if (fb->p == FB_END)
        return;
    int finished = 0;
    int cur_root = 0;
    fbprime_t p[4];
    fbprime_t r[4];
    fbprime_t R[4];
    unsigned char logp[4];
    int how_many = 0;
    // when we enter this loop, we assume that 
    //   fb->p > 0
    //   cur_root < fb->nb_roots
    // describe the next prime ideal to process.
    // (except if finished is set, of course)
    while (!finished) {
        // Get next 4 primes in fb.
        how_many = 0;
        while (how_many < 4) {
            //printf("Dealing with p = %d, r = %d\n", fb->p, fb->roots[cur_root]);
            p[how_many] = fb->p;
            R[how_many] = fb->roots[cur_root];
            logp[how_many] = fb->plog;
            r[how_many] = fb_root_in_qlattice(p[how_many], R[how_many], si);
            how_many++;
            if (r[how_many-1] == 0) {
                special_case_0 (fb->p);
                how_many--;
            } else if (r[how_many-1] == p[how_many-1]) {
                special_case_p (fb->p);
                how_many--;
            }
            if (cur_root < fb->nr_roots-1)
                cur_root++;
            else {
                cur_root = 0;
                fb = fb_next (fb);
                if (fb->p == FB_END) {
                    finished = 1;
                    break;
                }
            }
        }
        if (how_many != 4)
            break;

        // Here we have 4 consecutives prime ideals to sieve.

        const uint32_t I = si->I;
        const uint32_t maskI = I-1;
        plattice_info_t pli[4];
        // precomupte appropriate basis for p-lattices.
        {
            int i;
            for (i = 0; i < 4; ++i) {
                reduce_plattice(pli + i, p[i], r[i], si);
            }
        }

        // Start sieving from (0,0) which is I/2 in x-coordinate
        uint32_t x0, x1, x2, x3;
        x0 = (I>>1);
        x1 = (I>>1);
        x2 = (I>>1);
        x3 = (I>>1);
        const uint32_t IJ = I*si->J;

        v4si vec_b0 = { pli[0].b0, pli[1].b0, pli[2].b0, pli[3].b0 };
        v4si vec_b1 = { pli[0].b1, pli[1].b1, pli[2].b1, pli[3].b1 };
        v4si vec_a = { pli[0].a, pli[1].a, pli[2].a, pli[3].a };
        v4si vec_c = { pli[0].c, pli[1].c, pli[2].c, pli[3].c };
        v4si vec_x = { x0, x1, x2, x3 };
        v4si vec_maskI = { maskI, maskI, maskI, maskI };
        v4si vec_IJ = { IJ, IJ, IJ, IJ };
        v4si vec_zero = { 0, 0, 0, 0};
        __asm__("## Inner sieving routine starts here!!!\n");
        while (x0 < IJ) {
            uint32_t i;
            v4si vec_i, vec_tmp;
            vec_i = vec_x & vec_maskI;
            S[x0] -= logp[0];
            S[x1] -= logp[1];
            S[x2] -= logp[2];
            S[x3] -= logp[3];
            xmm_cmovge(&vec_tmp, vec_i, vec_b1, vec_a, vec_zero);
            vec_x += vec_tmp;
            xmm_cmovlt(&vec_tmp, vec_i, vec_b0, vec_c, vec_zero);
            vec_x += vec_tmp;
            // x = ( IJ < x ) ? IJ : x
            xmm_cmovlt(&vec_x, vec_IJ, vec_x, vec_IJ, vec_x);
            x0 = ((uint32_t *)(&vec_x))[0];
            x1 = ((uint32_t *)(&vec_x))[1];
            x2 = ((uint32_t *)(&vec_x))[2];
            x3 = ((uint32_t *)(&vec_x))[3];
        }
        __asm__("## Inner sieving routine stops here!!!\n");
    }
#endif

    // Loop over all primes in the factor base > I
    while (fb->p != FB_END) {
        unsigned char nr;
        fbprime_t p = fb->p;
        unsigned char logp = fb->plog;

        for (nr = 0; nr < fb->nr_roots; ++nr) {
            fbprime_t r, R;
            R = fb->roots[nr];
            r = fb_root_in_qlattice(p, R, si);
            if (r == 0) {
                special_case_0 (S, p, logp, si);
                continue;
            } 
            if (r == p) {
                special_case_p (S, p, logp, si);
                continue;
            } 
            
            const uint32_t I = si->I;
            const uint32_t maskI = I-1;
 
            plattice_info_t pli;
            reduce_plattice(&pli, p, r, si);

            // Start sieving from (0,0) which is I/2 in x-coordinate
            uint32_t x;
            x = (I>>1);
            // TODO: to gain speed, by aligning the start of the
            // loop, in assembly.
            // Besides this small trick, gcc does a good job with
            // this loop: no branching for the if(), and a dozen of
            // instructions inside the while().
            __asm__("## Inner sieving routine starts here!!!\n");
            while (x < I*si->J) {
                uint32_t i;
                i = x & maskI;   // x mod I
                S[x] -= logp;
#ifdef TRACE_K
                if (x == TRACE_K)
                  fprintf (stderr, "Subtracted %u (%u) to S[%u], remains %u\n",
                           logp, p, x, S[x]);
#endif
                if (i >= pli.b1)
                    x += pli.a;
                if (i < pli.b0)
                    x += pli.c;
            }
            __asm__("## Inner sieving routine stops here!!!\n");
        }
        fb = fb_next (fb); // cannot do fb++, due to variable size !
    }
    fprintf (stderr, "# large primes sieved in %f sec\n", seconds()-tm);
}
#endif

#ifdef USE_BUCKETS
// This version of the sieve implements the following:
//   - naive line-sieving for small p
//   - Franke-Kleinjung lattice-sieving for large p
//   - write reports to buckets and then apply them to the sieve-array
// WARNING: STILL AT WORK!
static void
sieve_buckets (unsigned char *S, factorbase_degn_t *fb,
        const sieve_info_t *si)
{
    // Sieve primes <= I.
    line_sieve(S, &fb, si);

    double tm = seconds();
    // Init buckets
    bucket_array_t BA;
    BA = init_bucket_array(si->nb_buckets, si->bucket_limit);

    // Loop over all primes in the factor base > I
    while (fb->p != FB_END) {
        unsigned char nr;
        fbprime_t p = fb->p;
        unsigned char logp = fb->plog;

        for (nr = 0; nr < fb->nr_roots; ++nr) {
            fbprime_t r, R;
            R = fb->roots[nr];
            r = fb_root_in_qlattice(p, R, si);
            // Special cases should be quite rare for primes > I,
            // so we handle them naively.
            if (r == 0) {
                special_case_0 (S, p, logp, si);
                continue;
            } 
            if (r == p) {
                special_case_p (S, p, logp, si);
                continue;
            } 
            
            const uint32_t I = si->I;
            const uint32_t maskI = I-1;
            const uint32_t maskbucket = si->bucket_region - 1;
            const int shiftbucket = si->log_bucket_region;
 
            plattice_info_t pli;
            if (reduce_plattice(&pli, p, r, si) == 0)
              continue; /* Simply don't consider that (p,r) for now.
                           FIXME: can we find the locations to sieve? */
              
            // Start sieving from (0,0) which is I/2 in x-coordinate
            uint32_t x;
            x = (I>>1);
            // Skip (0,0), since this can not be a valid report.
            {
                uint32_t i;
                i = x & maskI;
                if (i >= pli.b1)
                    x += pli.a;
                if (i < pli.b0)
                    x += pli.c;
            }
            // TODO: to gain speed, by aligning the start of the
            // loop, in assembly.
            // Besides this small trick, gcc does a good job with
            // this loop: no branching for the if(), and a dozen of
            // instructions inside the while().
            __asm__("## Inner sieving routine starts here!!!\n");
            while (x < I*si->J) {
                uint32_t i;
                bucket_update_t update;
                i = x & maskI;   // x mod I
                // S[x] -= logp;
                update.x = (uint16_t) (x & maskbucket);
                update.logp = logp;
                update.p_low = 0;
                push_bucket_update(BA, x >> shiftbucket, update);
#ifdef TRACE_K
                if (x == TRACE_K)
                  fprintf (stderr, "Pushed (%u, %u) (%u) to BA[%u]\n",
                           x & maskbucket, logp, p, x >> shiftbucket);
#endif
                if (i >= pli.b1)
                    x += pli.a;
                if (i < pli.b0)
                    x += pli.c;
            }
            __asm__("## Inner sieving routine stops here!!!\n");
        }
        fb = fb_next (fb); // cannot do fb++, due to variable size !
    }
    fprintf (stderr, "# large primes pushed to buckets in %f sec\n", seconds()-tm);

    // Apply buckets to sieve array
    tm = seconds();
    int i;
    unsigned char * S_ptr = S;
    for (i = 0; i < si->nb_buckets; ++i) {
        int j = nb_of_updates(BA, i);
        for ( ; j > 0 ; --j) {
            bucket_update_t update = get_next_bucket_update(BA, i);
            S_ptr[update.x] -= update.logp;
#ifdef TRACE_K
            if ((update.x + i*si->bucket_region) == TRACE_K)
                fprintf (stderr, "Subtract %u to S[%u], from BA[%u]\n",
                        update.logp, TRACE_K, i);
#endif
        }
        S_ptr += si->bucket_region;
    }
    clear_bucket_array(BA);
    fprintf (stderr, "# buckets applied to sieve array in %f sec\n", seconds()-tm);
}
#endif

// Conversions between different representations for sieve locations:
//   x          is the index in the sieving array. x in [0,I*J[
//   (i,j)      is the coordinates in the q-lattice. i in [-I/2,I/2[
//                                                   j in [0,J[
//   (a,b)      is the original coordinates. a is signed, b is unsigned.

void 
xToIJ(int *i, unsigned int *j, int x, sieve_info_t * si)
{
    *i = (x % (si->I)) - (si->I >> 1);
    *j = x / si->I;
}

void
IJTox(int *x, int i, int j, sieve_info_t * si)
{
    *x = i + (si->I)*j + (si->I>>1);
}

void
IJToAB(int64_t *a, uint64_t *b, int i, int j, sieve_info_t * si)
{
    *a = i*si->a0 + j*si->a1;
    *b = i*si->b0 + j*si->b1;
}

/* Warning: b might be negative, in which case we return (-a,-b) */
void
xToAB(int64_t *a, uint64_t *b, int x, sieve_info_t * si)
{
    int i, j;
    int64_t c;

    i = (x % (si->I)) - (si->I >> 1);
    j = x / si->I;
    *a = (int64_t) i * (int64_t) si->a0 + (int64_t) j * (int64_t) si->a1;
    c =  (int64_t) i * (int64_t) si->b0 + (int64_t) j * (int64_t) si->b1;
    if (c >= 0)
      *b = c;
    else
      {
        *a = -*a;
        *b = -c;
      }
}

/*********************** norm computation ************************************/

/* Puts in fij[] the coefficients of f'(i) = F(a0*i+a1, b0*i+b1).
   Assumes the coefficients of fij[] are initialized. */
static void
fij_from_f (mpz_t *fij, mpz_t *f, int d, int32_t a0, int32_t a1,
            int32_t b0, int32_t b1)
{
  int k, l;
  mpz_t *g; /* will contain the coefficients of (b0*i+b1)^l */
  mpz_t f0;

  for (k = 0; k <= d; k++)
    mpz_set (fij[k], f[k]);

  g = malloc ((d + 1) * sizeof (mpz_t));
  for (k = 0; k <= d; k++)
    mpz_init (g[k]);
  mpz_init (f0);

  /* Let h(x) = quo(f(x), x), then F(x,y) = H(x,y)*x + f0*y^d, thus
     F(a0*i+a1, b0*i+b1) = H(a0*i+a1, b0*i+b1)*(a0*i+a1) + f0*(b0*i+b1)^d.
     We use that formula recursively. */

  mpz_set_ui (g[0], 1); /* g = 1 */

  for (k = d - 1; k >= 0; k--)
    {
      /* invariant: we have already translated coefficients of degree > k,
         in f[k+1..d], and g = (b0*i+b1)^(d - (k+1)), with coefficients in
         g[0..d - (k+1)]:
         f[k]   <- f[k] + a1*f[k+1]
         ...
         f[l] <- a0*f[l]+a1*f[l+1] for k < l < d
         ...
         f[d] <- a0*f[d] */
      mpz_swap (f0, fij[k]); /* save the new constant coefficient */
      mpz_mul_si (fij[k], fij[k + 1], a1);
      for (l = k + 1; l < d; l++)
        {
          mpz_mul_si (fij[l], fij[l], a0);
          mpz_addmul_si (fij[l], fij[l + 1], a1);
        }
      mpz_mul_si (fij[d], fij[d], a0);

      /* now compute (b0*i+b1)^(d-k) from the previous (b0*i+b1)^(d-k-1):
         g[d-k] = b0*g[d-k-1]
         ...
         g[l] = b1*g[l]+b0*g[l-1] for 0 < l < d-k
         ...
         g[0] = b1*g[0]
      */
      mpz_mul_si (g[d - k], g[d - k - 1], b0);
      for (l = d - k - 1; l > 0; l--)
        {
          mpz_mul_si (g[l], g[l], b1);
          mpz_addmul_si (g[l], g[l-1], b0);
        }
      mpz_mul_si (g[0], g[0], b1);

      /* now g has degree d-k, and we add f0*g */
      for (l = k; l <= d; l++)
        mpz_addmul (fij[l], g[l - k], f0);
    }

  mpz_clear (f0);
  for (k = 0; k <= d; k++)
    mpz_clear (g[k]);
  free (g);
}

/* return max |g(x)| for 0 <= x <= s,
   where g(x) = g[d]*x^d + ... + g[1]*x + g[0] */
static double
get_maxnorm_aux (double *g, const unsigned int d, double s)
{
  unsigned int k, l, sign_change, new_sign_change;
  double **dg;    /* derivatives of g */
  double a, va, b, vb;
  double *roots, gmax;
  
  dg = (double**) malloc (d * sizeof (double*));
  dg[0] = g;
  for (k = 1; k < d; k++) /* dg[k] is the k-th derivative, thus has
                             degree d-k, i.e., d-k+1 coefficients */
    dg[k] = (double*) malloc ((d - k + 1) * sizeof (double));
  roots = (double*) malloc (d * sizeof (double));
  for (k = 1; k < d; k++)
    for (l = 0; l <= d - k; l++)
      dg[k][l] = (l + 1) * dg[k - 1][l + 1];
  /* now dg[d-1][0]+x*dg[d-1][1] is the (d-1)-th derivative: it can have at
     most one sign change, iff dg[d-1][0] and dg[d-1][0]+dg[d-1][1] have
     different signs */
  if (dg[d-1][0] * (dg[d-1][0] + dg[d-1][1]) < 0)
    {
      sign_change = 1;
      roots[0] = - dg[d-1][0] / dg[d-1][1]; /* root of (d-1)-th derivative */
    }
  else
    sign_change = 0;
  roots[sign_change] = s; /* end of interval */
  for (k = d - 1; k-- > 1;)
    {
      /* invariant: sign_change is the number of sign changes of the
         (k+1)-th derivative, with corresponding roots in roots[0]...
         roots[sign_change-1], and roots[sign_change] = s. */
      a = 0.0;
      va = dg[k][0]; /* value of dg[k] at x=0 */
      new_sign_change = 0;
      for (l = 0; l <= sign_change; l++)
        {
          b = roots[l]; /* root of dg[k+1], or end of interval */
          vb = fpoly_eval (dg[k], d - k, b);
          if (va * vb < 0) /* root in interval */
            roots[new_sign_change++] = fpoly_dichotomy (dg[k], d - k,
                                                        a, b, va, 20);
          a = b;
          va = vb;
        }
      roots[new_sign_change] = s; /* end of interval */
      sign_change = new_sign_change;
    }
  /* now all extrema of g are 0, roots[0], ..., roots[sign_change] = s */
  gmax = fabs (g[0]);
  for (k = 0; k <= sign_change; k++)
    {
      va = fabs (fpoly_eval (g, d, roots[k]));
      if (va > gmax)
        gmax = va;
    }
  free (roots);
  for (k = 1; k < d; k++)
    free (dg[k]);
  free (dg);
  return gmax;
}

/* returns the maximal value of log |F(a,b)/q| for
   a = a0 * i + a1 * j, b = b0 * i + b1 * j and q >= q0,
   -I/2 <= i <= I/2, 0 <= j <= I/2*min(s*B/|a1|,B/|b1|)
   where B >= sqrt(2*q/s/sqrt(3)) for all special-q in the current range
   (s is the skewness, and B = si.B, see lattice.tex).

   Since |a0| <= s*B and |b0| <= B, then
   |a0 * i + a1 * j| <= s*B*I and |b0 * i + b1 * j| <= B*I,
   thus it suffices to compute M = max |F(x,y)| in the rectangle
   -s <= x <= s, 0 <= y <= 1, and to multiply M by (B*I)^deg(F).

   Since F is homogeneous, we know M = max |F(x,y)| is attained on the border
   of the rectangle, i.e.:
   (a) either on F(s,y) for 0 <= y <= 1
   (b) either on F(x,1) for -s <= x <= s
   (c) either on F(-s,y) for 0 <= y <= 1
   (d) or on F(x,0) for -s <= x <= s, but this maximum is f[d]*s^d,
       and is attained in (a) or (c).
*/
static double
get_maxnorm (cado_poly cpoly, sieve_info_t *si, uint64_t q0)
{
  unsigned int d = cpoly->degree, k;
  double *fd; /* double-precision coefficients of f */
  double norm, max_norm, pows, tmp;
  
  fd = (double*) malloc ((d + 1) * sizeof (double));
  for (k = 0; k <= d; k++)
    fd[k] = mpz_get_d (cpoly->f[k]);

  /* (b1) determine the maximum of |f(x)| for 0 <= x <= s */
  max_norm = get_maxnorm_aux (fd, d, cpoly->skew);

  /* (b2) determine the maximum of |f(-x)| for 0 <= x <= s */
  norm = get_maxnorm_aux (fd, d, -cpoly->skew);
  if (norm > max_norm)
    max_norm = norm;

  for (pows = 1.0, k = 0; k <= d; k++)
    {
      fd[k] *= pows;
      pows *= cpoly->skew;
    }
  /* swap coefficients; if d is odd, we need to go up to k = floor(d/2) */
  for (k = 0; k <= d / 2; k++)
    {
      tmp = fd[k];
      fd[k] = fd[d - k];
      fd[d - k] = tmp;
    }

  /* (a) determine the maximum of |g(y)| for 0 <= y <= 1 */
  norm = get_maxnorm_aux (fd, d, 1.0);
  if (norm > max_norm)
    max_norm = norm;

  /* (c) determine the maximum of |g(-y)| for 0 <= y <= 1 */
  norm = get_maxnorm_aux (fd, d, -1.0);
  if (norm > max_norm)
    max_norm = norm;
  
  free (fd);

  /* multiply by (B*I)^d and divide by q0 */
  return log (max_norm * pow (si->B * (double) si->I, (double) d)
              / (double) q0);
}

/* v <- f(i,j), where f is of degree d */
static void
eval_fij (mpz_t v, mpz_t *f, unsigned int d, long i, unsigned long j)
{
  unsigned int k;
  mpz_t jpow;

  mpz_init_set_ui (jpow, 1);
  mpz_set (v, f[d]);
  for (k = d; k-- > 0;)
    {
      mpz_mul_si (v, v, i);
      mpz_mul_ui (jpow, jpow, j);
      mpz_addmul (v, f[k], jpow);
    }
  mpz_clear (jpow);
}

/* initialize S[i+I/2+j*I] to round(log(|F_q(i,j)|/q)/log(rho))
   where log(maxnorm)/log(rho) = 255, i.e., rho = maxnorm^(1/255) */
static void
init_norms_alg (unsigned char *S, cado_poly cpoly, sieve_info_t *si)
{
  int i, halfI;
  unsigned int j, k, d = cpoly->degree;
  double *t, *u, invq, powj, norm;
  unsigned char *S_ptr;

  /* si->scale_alg = LOG_MAX / log(max |F(a0*i+a1*j, b0*i+b1*j)|) */
  
  /* si->fij is the f(i,j) polynomial at (0, 0) */

  halfI = si->I / 2;
  invq = 1.0 / (double) si->q;

  t = (double*) malloc ((d + 1) * sizeof (double));
  u = (double*) malloc ((d + 1) * sizeof (double));
  for (k = 0; k <= d; k++)
    t[k] = mpz_get_d (si->fij[k]) * invq;

  S_ptr = S + halfI;
  for (j = 0; j < si->J; j++)
    {
      /* scale by j^(d-k) the coefficients of fij */
      for (k = 0, powj = 1.0; k <= d; k++, powj *= (double) j)
        u[d - k] = t[d - k] * powj;
      
      /* now compute norms */
      for (i = -halfI; i < halfI; i++)
        {
          norm = fpoly_eval (u, d, (double) i);
          norm = fabs (norm);
          norm = log (norm);
          norm = norm * si->scale_alg;
          S_ptr[i] = GUARD + (unsigned char) (norm);
#ifdef TRACE_K
          if ((S_ptr - S) + i == TRACE_K)
            fprintf (stderr, "S[%ld] initialized to %u\n", (S_ptr - S) + i,
                     S_ptr[i]);
#endif          
        }
      S_ptr += si->I;
    }

  free (t);
  free (u);
}

/* check that the double x fits into an int32_t */
#define fits_int32_t(x) \
  ((double) INT32_MIN <= (x)) && ((x) <= (double) INT32_MAX)

void
SkewGauss (sieve_info_t *si, double skewness)
{
  double a[2], b[2], q;

  a[0] = (double) si->q;
  ASSERT_ALWAYS(a[0] < 9007199254740992.0); /* si.q should be less than 2^53
                                               so that a[0] is exact */
  b[0] = 0.0;
  a[1] = (double) si->rho;
  skewness = rint (skewness);
  b[1] = skewness;
  while (1)
    {
      /* reduce vector (a[0], b[0]) with respect to (a[1], b[1]) */
      q = (a[0] * a[1] + b[0] * b[1]) / (a[1] * a[1] + b[1] * b[1]);
      q = rint (q);
      if (q == 0.0)
        break;
      a[0] -= q * a[1];
      b[0] -= q * b[1];
      /* reduce vector (a[1], b[1]) with respect to (a[0], b[0]) */
      q = (a[0] * a[1] + b[0] * b[1]) / (a[0] * a[0] + b[0] * b[0]);
      q = rint (q);
      if (q == 0.0)
        break;
      a[1] -= q * a[0];
      b[1] -= q * b[0];
    }
  /* put the smallest vector in (a0,b0) */
  ASSERT(fits_int32_t(a[0]));
  ASSERT(fits_int32_t(b[0] / skewness));
  ASSERT(fits_int32_t(a[1]));
  ASSERT(fits_int32_t(b[1] / skewness));
  if (a[0] * a[0] + b[0] * b[0] < a[1] * a[1] + b[1] * b[1])
    {
      si->a0 = (int32_t) a[0];
      si->b0 = (int32_t) (b[0] / skewness);
      si->a1 = (int32_t) a[1];
      si->b1 = (int32_t) (b[1] / skewness);
    }
  else
    {
      si->a0 = (int32_t) a[1];
      si->b0 = (int32_t) (b[1] / skewness);
      si->a1 = (int32_t) a[0];
      si->b1 = (int32_t) (b[0] / skewness);
    }
}

/************************ factoring with TIFA ********************************/

/* FIXME: the value of 20 seems large. Normally, a few Miller-Rabin passes
   should be enough. See also http://www.trnicely.net/misc/mpzspsp.html */
#define NMILLER_RABIN 20
#define IS_PRIME(X)     (0 != mpz_probab_prime_p((X), NMILLER_RABIN))
#define BITSIZE(X)      (mpz_sizeinbase((X), 2))
#define NFACTORS        8 /* maximal number of large primes */

/* This function was contributed by Jerome Milan (and bugs were introduced
   by Paul Zimmermann :-).
   Input: n - the number to be factored (leftover norm)
          b - (large) prime bit size bound
   Return value:
          0 if n has a prime factor larger than 2^b
          1 if all prime factors of n are < 2^b
   Output:
          the prime factors of n are factors->data[0..factors->length-1],
          with corresponding multiplicities multis[0..factors->length-1].
*/          
int
factor_leftover_norm (mpz_t n, unsigned int b,
                      mpz_array_t* const factors, uint32_array_t* const multis)
{
  uint32_t i;
  ecode_t ecode;

  factors->length = 0;
  multis->length = 0;

  if (mpz_sgn (n) < 0)
    mpz_neg (n, n);

  /* it seems tifa_factor does not like 1 */
  if (mpz_cmp_ui (n, 1) == 0)
    return 1;

  if (IS_PRIME(n))
    {
      if (BITSIZE(n) > b)
        return 0;
      else
        {
          append_mpz_to_array (factors, n);
          append_uint32_to_array (multis, 1);
          return 1;
        }
    } 

  //  gmp_fprintf (stderr, "enter tifa_factor, n=%Zd\n", n);
  ecode = tifa_factor (factors, multis, n, FIND_COMPLETE_FACTORIZATION);
  // gmp_fprintf (stderr, "   exit tifa_factor\n");

  switch (ecode)
    {
    case COMPLETE_FACTORIZATION_FOUND:
      for (i = 0; i < factors->length; i++)
        {
          if (BITSIZE(factors->data[i]) > b)
            return 0;
        }
      return 1;

    case PARTIAL_FACTORIZATION_FOUND:
    case NO_FACTOR_FOUND:
    case FATAL_INTERNAL_ERROR:
    default:
        //
        // Should be rare but I cannot give any warranties here... We could
        // try to use another factoring library...
        //
        return 0;
    }
  return 0;
}

/************************ Bernstein's algorithm ******************************/

typedef struct {
  mpz_t *part;       /* each part[i] is a product of primes */
  unsigned long len; /* number of elements in part[] */
} acc_prime_t;

/* put in P all primes in fb and bp (badprimes), in products of size about 
   BERNSTEIN_THRESHOLD bits */
static void
acc_primes_init (acc_prime_t *P, factorbase_degn_t *fb, uint32_t *bp)
{
  unsigned long p, i;

  P->part = (mpz_t*) malloc (sizeof (mpz_t));
  P->len = 1;
  mpz_init_set_ui (P->part[0], 1);
  for (; bp[0] != 0; bp++)
    mpz_mul_ui (P->part[0], P->part[0], bp[0]);

  i = 0;
  p = fb->p;
  while (p != FB_END)
    {
      mpz_mul_ui (P->part[i], P->part[i], p);
      /* the threshold BT / 12 was determined by trial and error */
      if (mpz_sizeinbase (P->part[i], 2) >= BERNSTEIN_THRESHOLD / 12)
        {
          P->len = P->len + 1;
          P->part = (mpz_t*) realloc (P->part, P->len * sizeof (mpz_t));
          i = i + 1;
          mpz_init_set_ui (P->part[i], 1);
        }
      fb = fb_next (fb);
      p = fb->p;
    }
  if (mpz_cmp_ui (P->part[i], 1) == 0)
    {
      mpz_clear (P->part[i]);
      P->len = P->len - 1;
      P->part = (mpz_t*) realloc (P->part, P->len * sizeof (mpz_t));
    }
}

static void
acc_primes_clear (acc_prime_t *P)
{
  unsigned long i;

  for (i = 0; i < P->len; i++)
    mpz_clear (P->part[i]);
  free (P->part);
}

/* return ceil(log(k)/log(2)), assumes k >= 1 */
static unsigned int
ceil_log2 (unsigned long k)
{
  unsigned long l = 1;

  while (k > 1)
    {
      l ++;
      k = (k + 1) / 2;
    }
  return l;
}

/* we want a to be multiple of all primes in b anc c */
static void
combine_norms (mpz_t a, mpz_t b, mpz_t c)
{
  mpz_lcm (a, b, c);
}

static mpz_t**
BuildProductTree (sieve_report_t *report, unsigned int **lT0,
                  unsigned int *h0, unsigned long nreports)
{
  mpz_t **T;
  unsigned int *lT, h, l, k;

  h = ceil_log2 (nreports);
  T = (mpz_t**) malloc ((h + 1) * sizeof (mpz_t*));
  lT = (unsigned int*) malloc ((h + 1) * sizeof (unsigned int));
  for (l = 0; l <= h; l++)
    {
      /* T[0] is the level of the leaves, and T[h] is the root of the tree */
      lT[l] = 1 + ((nreports - 1) >> l);
      T[l] = (mpz_t*) malloc (lT[l] * sizeof (mpz_t));
      for (k = 0; k < lT[l]; k++)
        {
          mpz_init (T[l][k]);
          if (l == 0) /* compute rational norms */
            mpz_set (T[l][k], (report + k)->normr);
          else /* combine two subtrees */
            {
              if (2 * k + 1 < lT[l-1])
                combine_norms (T[l][k], T[l-1][2*k], T[l-1][2*k+1]);
              else
                mpz_set (T[l][k], T[l-1][2*k]);
            }
        }
    }
  *lT0 = lT;
  *h0 = h;
  return T;
}

/********************* factoring on rational side ****************************/

/* Those routines are a quick-and-dirty implementation which factors remaining
   norms on the rational side, after sieving on the algebraic side. It will
   probably become obsolete once sieving is done on the rational side.
*/

/* Compute gcd(P,t) for all elements of the tree T.
   T[h][0] is the root of the tree. */
void
GcdTree (mpz_t **T, unsigned int *lT, unsigned long h, mpz_t P)
{
  mpz_t *G;
  unsigned long l, m;
#ifdef VERBOSE
  double tm = cputime ();
#endif
  
  G = (mpz_t*) malloc (lT[0] * sizeof (mpz_t));
  for (l = 0; l < lT[0]; l++)
    mpz_init (G[l]);

  mpz_gcd (G[0], T[h][0], P); /* gcd of the root */
  /* dividing T[h][0] by the gcd does not seem to speed up things,
     similarly below for gcd(T[0][l], G[l]) for l > 0 */
#ifdef VERBOSE
  fprintf (stderr, "Tree root has %lu bits, prime product has %lu bits\n",
           mpz_sizeinbase (T[h][0], 2), mpz_sizeinbase (P, 2));
  fprintf (stderr, "Root gcd has size %lu bits\n", mpz_sizeinbase (G[0], 2));
  fprintf (stderr, "Root gcd took %1.0fms\n", cputime () - tm);
#endif
  for (l = h; l-- > 0;)
    { /* level l, we have lT[l] elements, the gcds of the upper level are
         in G[0...lT[l+1]-1] = G[0..ceil(lT[l]/2)-1] */
      for (m = lT[l]; m-- > 0;)
        mpz_gcd (G[m], T[l][m], G[m/2]);
    }
#ifdef VERBOSE
  fprintf (stderr, "Gcd tree took %1.0fms\n", cputime () - tm);
#endif

  /* now G[l] should divide T[0][l] for 0 <= l < lT[0] */
  for (l = 0; l < lT[0]; l++)
    {
      ASSERT (mpz_divisible_p (T[0][l], G[l]));
      do
        {
          mpz_divexact (T[0][l], T[0][l], G[l]);
          /* primes in G[l] might divide several times T[0][l] */
          mpz_gcd (G[l], T[0][l], G[l]);
        }
      while (mpz_cmp_ui (G[l], 1) != 0);
    }

  /* FIXME: now that the T[0][l] are smaller, should we recompute the
     product tree? Maybe from time to time. */

  for (l = 0; l < lT[0]; l++)
    mpz_clear (G[l]);
  free (G);
}

/* Divide norm by all primes in l, and print them in buf.
   Assumes end of list is 0 in l.
*/
static void
trial_divide_list (char *buf, mpz_t norm, uint32_t *l)
{
  unsigned long p = *l;

  while (p != 0)
    {
      while (mpz_divisible_ui_p (norm, p))
        {
          mpz_divexact_ui (norm, norm, p);
          buf += sprintf (buf, "%lx,", p);
        }
      p = *l++;
    }
  buf[0] = '\0';
}

/* Divide norm by all primes <= B, and print them in buf.
   If fb <> NULL, use only primes in fb.
 */
static void
trial_divide (char *buf, mpz_t norm, const unsigned long B,
              factorbase_degn_t *fb)
{
  unsigned long p, n = 0, P = B;
  double d;

  mpz_abs (norm, norm);
  /* we can stop as soon as the current prime is larger than sqrt(norm) */
  d = sqrt (mpz_get_d (norm));
  if (d < (double) P)
    P = (unsigned long) d;
  p = (fb == NULL) ? 2 : fb->p;
  for (;p <= P;)
    {
      while (mpz_divisible_ui_p (norm, p))
        {
          mpz_divexact_ui (norm, norm, p);
          buf += sprintf (buf, "%lx,", p);
          n ++;
          d = sqrt (mpz_get_d (norm));
          if (d < (double) P)
            P = (unsigned long) d;
        }
      /* if norm <= p^2 then norm can only have one remaining prime,
         thus it suffices to check p <= sqrt(norm) */
      if (fb == NULL)
        p = getprime (p);
      else
        {
          fb = fb_next (fb);
          p = fb->p;
          if (p == FB_END)
            break;
        }
    }
  if (mpz_cmp_ui (norm, B) <= 0)
    {
      buf += sprintf (buf, "%lx,", mpz_get_ui (norm));
      mpz_set_ui (norm, 1);
      n ++;
    }
  if (n == 0)
    buf[0] = '\0';
  else
    buf[-1] = '\0'; /* replace last ',' by end of string */
  if (fb == NULL)
    getprime (0);
}

/* Build a product tree with the rational norms,
   of height h = ceil(log2(reports)).
   Return the number of reports.
*/
static unsigned long
build_norms_rational (cado_poly cpoly, sieve_info_t *si,
                      sieve_report_list_t *Reports, acc_prime_t *P_rat,
                      acc_prime_t *P_alg, factorbase_degn_t *fb_alg)
{
  mpz_t **T, normr, norma;
  unsigned int k, l, h, delta;
  unsigned int *lT; /* lT[j] = length of T[j] */
  unsigned int nreports = Reports->size;
  mpz_t *P = NULL;         /* prime accumulator */
  unsigned long aP = 0;    /* allocated cells in P[] */
  double tm, ttrialr = 0, ttriala = 0, tt, tg = 0;
  unsigned long final_reports = 0, rat_reports = 0, ngcds = 0;
  int64_t a;
  uint64_t b;
  char bufr[MAX_BUF_SIZE], bufa[MAX_BUF_SIZE];
  size_t st; /* size in bits of root of product tree */
  mpz_array_t *f_r = NULL, *f_a = NULL;    /* large prime factors */
  uint32_array_t *m_r = NULL, *m_a = NULL; /* corresponding multiplicities */

  if (si->checknorms)
    {
      f_r = alloc_mpz_array (NFACTORS);
      f_a = alloc_mpz_array (NFACTORS);
      m_r = alloc_uint32_array (NFACTORS);
      m_a = alloc_uint32_array (NFACTORS);
    }

  mpz_init (normr);
  mpz_init (norma);

  /* we consider only the last reports whose product of rational norms
     has about BERNSTEIN_THRESHOLD bits */
  for (l = 0, k = nreports; l < BERNSTEIN_THRESHOLD && k > 0;)
    l += mpz_sizeinbase ((Reports->report + --k)->normr, 2);
  nreports = nreports - k;
  delta = k;               /* the first delta reports are ignored */

  tt = seconds ();
  T = BuildProductTree (Reports->report + delta, &lT, &h, nreports);
  st = mpz_sizeinbase (T[h][0], 2);
  tt = seconds () - tt;

  tg -= seconds ();
  for (k = 0; k < P_rat->len; k++, ngcds ++)
    GcdTree (T, lT, h, P_rat->part[k]);
  tg += seconds ();

  fprintf (stderr, "# Bernstein's algo on rat. side took %1.1fs\n", tt + tg);
  fprintf (stderr, "#   (norm tree %1.1f, %lu gcds %1.1f)\n", tt, ngcds, tg);

  tm = seconds ();
  /* check leftover norms that are smaller than the given bound */
  for (l = 0; l < nreports; l++)
    {
      int large_size = mpz_sizeinbase (T[0][l], 2);

      /* if T[0][l] is a prime > 2^lpbr, we can ignore the relation */
      if (large_size <= cpoly->mfbr &&
          !(large_size > cpoly->lpbr && mpz_probab_prime_p (T[0][l], 1)))
        {
          rat_reports ++; /* reports such that the algebraic side is
                             approximately smooth (up to rounding errors),
                             and the rational leftover norm is < 2^mfbr
                             and is not a prime > 2^lpbr */
          /* report_list[l] is I/2+i+j*I */

          /* since we know the norm is anyway smooth on the rational side, 
             we delay the trial division after testing the algebraic side */

          /* we know that q divides the norm on the algebraic side */
          mpz_swap (norma, (Reports->report + delta + l)->norma);
          /* on the algebraic side, we restrict to the factor base primes */
          ttriala -= seconds ();
          /* first divide out the bad primes */
          trial_divide_list (bufa, norma, si->abadprimes);
          /* then divide the primes in the algebraic factor base */
          trial_divide (bufa + strlen (bufa), norma, cpoly->alim, fb_alg);
          ttriala += seconds ();
          /* we do not take q into account for the mfb bound: if we want
             two large primes (lambda = 2, mfb = 2*lpb), then those large
             primes will come in addition to q */
          large_size = mpz_sizeinbase (norma, 2);
          if (large_size <= cpoly->mfba &&
              !(large_size > cpoly->lpba && mpz_probab_prime_p (norma, 1)))
            {
              /* now trial divide the rational norm */
              /* we know the product of large primes is T[0][l], thus we can
                 divide them out */
              mpz_divexact (normr, (Reports->report + delta + l)->normr,
                            T[0][l]);
              ttrialr -= seconds ();
              trial_divide (bufr, normr, cpoly->rlim, NULL);
              ttrialr += seconds ();
              /* since we divided out large primes, norm should be 1 here */
              ASSERT(mpz_cmp_ui (normr, 1) == 0);

              if (si->checknorms == 0 ||
                  (factor_leftover_norm (T[0][l], cpoly->lpbr, f_r, m_r) &&
                   factor_leftover_norm (norma, cpoly->lpba, f_a, m_a)))
                                                               
                {
                  final_reports ++;
                  /* add large prime on alg. side */
                  sprintf_cat (bufa, (Reports->report + delta + l)->q);
                  if (si->checknorms)
                    {
                      uint32_t i, j;
                      for (i = 0; i < f_r->length; i++)
                        {
                          CHECK_LARGE_PRIME (f_r->data[i], cpoly->rlim, "rat");
                          for (j = 0; j < m_r->data[i]; j++)
                            gmp_sprintf_cat (bufr, f_r->data[i]);
                        }
                      for (i = 0; i < f_a->length; i++)
                        {
                          CHECK_LARGE_PRIME (f_a->data[i], cpoly->alim, "alg");
                          for (j = 0; j < m_a->data[i]; j++)
                            gmp_sprintf_cat (bufa, f_a->data[i]);
                        }
                    }
                  a = (Reports->report + delta + l)->a;
                  b = (Reports->report + delta + l)->b;
                  printf ("%" PRId64 ",%" PRIu64 ":%s:%s\n", a, b, bufr, bufa);
                  fflush (stdout);
                }
            }
        }
    }
  fprintf (stderr, "# Trial division took %1.1fs", seconds () - tm);
  fprintf (stderr, " (rat. %1.1fs, alg. %1.1fs)\n", ttrialr, ttriala);

  for (l = 0; l < aP; l++)
    mpz_clear (P[l]);
  free (P);
  for (l = 0; l <= h; l++)
    {
      for (k = 0; k < lT[l]; k++)
        mpz_clear (T[l][k]);
      free (T[l]);
    }
  free (T);
  free (lT);
  mpz_clear (normr);
  mpz_clear (norma);

  if (si->checknorms)
    {
      clear_mpz_array (f_r);
      clear_mpz_array (f_a);
      clear_uint32_array (m_r);
      clear_uint32_array (m_a);
    }

  fprintf (stderr, "# Remaining reports on rational side: %lu\n", rat_reports);

  /* update number of reports */
  Reports->size -= nreports;

  return final_reports;
}

/* Add to report_list L values of k = i + I/2 + I*j where gcd(i,j)=1
   and norm is small on rational side.
*/
static void
sieve_report_list_accumulate (sieve_report_list_t *L, cado_poly cpoly,
                              sieve_info_t *si, unsigned char *S)
{
    /* sieve reports are those for which the remaining bit-size is less than
       lambda * lpb */
    unsigned char min = UCHAR_MAX;
    int kmin = 0;
    int i;
    unsigned int j, k;
    int64_t a;
    uint64_t b;

    /* for 0 <= k < I, we have j=0, thus a=a0*i and b=b0*i,
       and it suffices to consider i=1, i.e., k = I/2+1 */
    k = si->I / 2 + 1;
    if (S[k] < min)
      {
        kmin = k;
        min = S[k];
      }
    if (si->rat_Bound[S[k]] != 0)
      {
        a = si->a0;
        b = si->b0;
        ASSERT_ALWAYS(b != 0);
        /* gcd(a0,b0)=1, since it divides q */
        sieve_report_list_add (L, a, b, si->q, cpoly);
      }
    for (k = si->I; k < si->I * si->J; ++k)
      {
        if (S[k] < min)
          {
            kmin = k;
            min = S[k];
          }
        if (si->rat_Bound[S[k]] != 0)
          {
            /* we already checked gcd(i,j)=1 before, thus necessarily a and
               b are coprime here */
            xToAB (&a, &b, k, si);
            sieve_report_list_add (L, a, b, si->q, cpoly);
          }
      }
    xToIJ(&i, &j, kmin, si);
    fprintf (stderr, "# Min of %d is obtained for (i,j) = (%d,%d)\n", min, i, j);
}

/* Initialize lognorms on the rational side. We only initialize cells which
   correspond to sieve reports on the algebraic side, and for which gcd(a,b)=1,
   and put the maximal value (255) elsewhere. */
static void
init_norms_rat (unsigned char *S, cado_poly cpoly, sieve_info_t *si)
{
  double g1, g0, gi, gj, norm;
  int i, halfI = si->I / 2, kmin = 0;
  unsigned int j, k, reports = 0;
  unsigned char c, min = UCHAR_MAX;

  /* G_q(i,j) = g1*(a0*i+a1*j)+g0*(b0*i+b1*j)
              = (g1*a0+g0*b0)*i + (g1*a1+g0*b1)*j */

  g1 = mpz_get_d (cpoly->g[1]);
  g0 = mpz_get_d (cpoly->g[0]);
  gi = g1 * (double) si->a0 + g0 * (double) si->b0;
  gj = g1 * (double) si->a1 + g0 * (double) si->b1;

  /* for j=0, consider only i=1, i.e., k = I/2+1 */
  k = 0;
  c = S[si->I/2 + 1]; /* save value */
  for (i = -halfI; i < halfI; i++, k++)
    {
      if (S[k] < min)
        {
          kmin = k;
          min = S[k];
        }
      S[k] = UCHAR_MAX;
    }
  if (si->alg_Bound[c] != 0)
    {
      norm = gi;
      norm = fabs (norm);
      norm = log (norm);
      norm = norm * si->scale_rat;
      S[si->I/2 + 1] = GUARD + (unsigned char) (norm);
      reports ++;
    }
  for (j = 1; j < si->J; j++)
    {
      for (i = -halfI; i < halfI; i++, k++)
        {
          /* FIXME: we could process i and -i simultaneously, to perform
             only one gcd_ul */
          if (S[k] < min)
            {
              kmin = k;
              min = S[k];
            }
          if (si->alg_Bound[S[k]] != 0 && j != 0 &&
              gcd_ul ((i > 0) ? i : -i, j) == 1)
            {
              norm = gi * (double) i + gj * (double) j;
              norm = fabs (norm);
              norm = log (norm);
              norm = norm * si->scale_rat;
              S[k] = GUARD + (unsigned char) (norm);
              reports ++;
            }
          else
            S[k] = UCHAR_MAX;
        }
    }
  xToIJ(&i, &j, kmin, si);
  fprintf (stderr, "# Min of %d is obtained for (i,j) = (%d,%d)\n", min, i, j);
  fprintf (stderr, "# Number of half-reports: %u\n", reports);
}

/*************************** main program ************************************/

static void
usage (char *argv0)
{
  fprintf (stderr, "Usage: %s [-checknorms] [-I I] -poly xxx.poly -fb xxx.roots -q0 q0 [-q1 q1] [-rho rho]\n",
           argv0);
  fprintf (stderr, "          -checknorms     factor leftover norms\n");
  fprintf (stderr, "          -I i            sieving region has side 2^i [default %u]\n", DEFAULT_I);
  fprintf (stderr, "          -poly xxx.poly  use polynomial xxx.poly\n");
  fprintf (stderr, "          -fb xxx.roots   use factor base xxx.roots\n");
  fprintf (stderr, "          -q0 nnn         left bound of special-q range\n");
  fprintf (stderr, "          -q1 nnn         right bound of special-q range\n");
  fprintf (stderr, "          -rho r          sieve only algebraic root r mod q0\n");
  exit (1);
}

int
main (int argc, char *argv[])
{
    char *argv0 = argv[0];
    sieve_info_t si;
    char *fbfilename = NULL, *polyfilename = NULL;
    cado_poly cpoly;
    double t0, tnorma, tnormr, tfb, ts, tf, tq, ttn, tts, ttf, tp;
    uint64_t q0 = 0, q1 = 0, rho = 0;
    unsigned long *roots, nroots, reports, tot_reports = 0, i;
    unsigned char * S;
    factorbase_degn_t * fb_alg, * fb_rat;
    int checknorms = 0; /* factor or not the remaining norms */
    int I = DEFAULT_I;
    unsigned long sq = 0;
    double totJ = 0.0;
    sieve_report_list_t Reports[1];
    acc_prime_t P_rat, P_alg;

    fprintf (stderr, "# %s.r%s", argv[0], REV);
    for (i = 1; i < (unsigned int) argc; i++)
      fprintf (stderr, " %s", argv[i]);
    fprintf (stderr, "\n");

    while (argc > 1 && argv[1][0] == '-')
      {
        if (strcmp (argv[1], "-checknorms") == 0)
          {
            checknorms = 1;
            argc -= 1;
            argv += 1;
          }
        else if (argc > 2 && strcmp (argv[1], "-I") == 0)
          {
            I = atoi (argv[2]);
            argc -= 2;
            argv += 2;
          }
        else if (argc > 2 && strcmp (argv[1], "-fb") == 0)
          {
            fbfilename = argv[2];
            argc -= 2;
            argv += 2;
          }
        else if (argc > 2 && strcmp (argv[1], "-poly") == 0)
          {
            polyfilename = argv[2];
            argc -= 2;
            argv += 2;
          }
        else if (argc > 2 && strcmp (argv[1], "-q0") == 0)
          {
            /* uintmax_t is guaranteed to be larger or equal to uint64_t */
            q0 = strtouint64 (argv[2], NULL, 10);
            argc -= 2;
            argv += 2;
          }
        else if (argc > 2 && strcmp (argv[1], "-q1") == 0)
          {
            q1 = strtouint64 (argv[2], NULL, 10);
            argc -= 2;
            argv += 2;
          }
        else if (argc > 2 && strcmp (argv[1], "-rho") == 0)
          {
            rho = strtouint64 (argv[2], NULL, 10);
            argc -= 2;
            argv += 2;
          }
        else
          usage (argv0);
      }

    if (polyfilename == NULL || fbfilename == NULL || q0 == 0)
      usage (argv0);

    /* if -rho is given, we sieve only for q0, thus -q1 is not allowed */
    if (rho != 0 && q1 != 0)
      {
        fprintf (stderr, "Error, -q1 and -rho are mutually exclusive\n");
        exit (1);
      }

    /* if -q1 is not given, sieve only for q0 */
    if (q1 == 0)
      q1 = q0 + 1;

    /* check that q1 fits into an unsigned long */
    if (q1 > (uint64_t) ULONG_MAX)
      {
        fprintf (stderr, "Error, q1=%" PRIu64 " exceeds ULONG_MAX\n", q1);
        exit (1);
      }

    if (!read_polynomial (cpoly, polyfilename))
      {
        fprintf (stderr, "Error reading polynomial file\n");
        exit (EXIT_FAILURE);
      }

    /* this does not depend on the special-q */
    sieve_info_init (&si, cpoly, I, q0);
    si.checknorms = checknorms;

    // One more is for the garbage during SSE
    S = (unsigned char *) malloc ((1 + si.I * si.J)*sizeof(unsigned char));
    ASSERT_ALWAYS(S != NULL);

    /* the logarithm scale is LOG_MAX / log(max norm) */
    tfb = seconds ();
    fb_alg = fb_read (fbfilename, si.scale_alg, 0);
    ASSERT_ALWAYS(fb_alg != NULL);
    tfb = seconds () - tfb;
    fprintf (stderr, "# Reading algebraic factor base took %1.1fs\n", tfb);

    /* Prepare rational factor base */
    tfb = seconds ();
    fb_rat = fb_make_linear (cpoly->g, (fbprime_t) cpoly->rlim, si.scale_rat,
                             0);
    tfb = seconds () - tfb;
    fprintf (stderr, "# Creating rational factor base took %1.1fs\n", tfb);

    compute_badprimes (&si, cpoly, fb_alg, fb_rat);

    tp = seconds ();
    acc_primes_init (&P_rat, fb_rat, si.rbadprimes);
    acc_primes_init (&P_alg, fb_alg, si.abadprimes);
    tp = seconds () - tp;
    fprintf (stderr, "# Initializing %lu+%lu prime products took %1.1fs\n",
             P_rat.len, P_alg.len, tp);

    /* initialize sieve report list */
    sieve_report_list_init (Reports);

    /* special q (and root rho) */
    t0 = seconds ();
    roots = (unsigned long*) malloc (cpoly->degree * sizeof (unsigned long));
    q0 --; /* so that nextprime gives q0 if q0 is prime */
    nroots = 0;
    ttn = tts = ttf = 0.0;
    while (q0 < q1)
      {
        while (nroots == 0) /* go to next prime and generate roots */
          {
            unsigned long q;

            q0 = uint64_nextprime (q0);
            if (q0 >= q1)
              goto end;
            si.q = q0;
            q = q0; /* modul_roots_mod_long works on an unsigned long */
            nroots = modul_roots_mod_long (roots, cpoly->f, cpoly->degree, &q);
            fprintf (stderr, "# q=%" PRIu64 ": root(s)", q0);
            for (i = 1; i <= nroots; i++)
              fprintf (stderr, " %lu", roots[nroots-i]);
            fprintf (stderr, "\n");
          }
        tq = seconds ();

        /* computes a0, b0, a1, b1 from q, rho, and the skewness */
        si.rho = roots[--nroots];
        if (rho != 0 && si.rho != rho) /* if -rho, wait for wanted root */
          continue;
        SkewGauss (&si, cpoly->skew);
        fprintf (stderr, "# Sieving q=%" PRIu64 "; rho=%" PRIu64
                 "; a0=%d; b0=%d; a1=%d; b1=%d\n",
                 si.q, si.rho, si.a0, si.b0, si.a1, si.b1);
        sq ++;

        /* checks the value of J */
        sieve_info_update (&si, cpoly->skew);
        totJ += (double) si.J;

        /* norm computation */
        tnorma = seconds ();
        fij_from_f (si.fij, cpoly->f, cpoly->degree, si.a0, si.a1, si.b0,
                    si.b1);
        init_norms_alg (S, cpoly, &si);
        tnorma = seconds () - tnorma;
        fprintf (stderr, "# Initializing alg. norms took %1.1fs\n", tnorma);
#ifdef TRACE_K
        {
            int __i; unsigned int __j;
            xToIJ(&__i, &__j, TRACE_K, &si);
            fprintf(stderr, "Trace: %d, i.e. (i,j) = (%d,%d)\n", TRACE_K, __i, __j);
        }
#endif

        /* sieving on the algebraic side */
        ts = seconds ();
        //sieve_slow(S, fb_alg, &si);
#ifndef USE_BUCKETS
        sieve_random_access(S, fb_alg, &si);
#else
        sieve_buckets(S, fb_alg, &si);
#endif
        fprintf (stderr, "# Alg. sieving in %1.1fs\n", seconds () - ts);

        /* Sieving on the rational side */
        tnormr = seconds ();
        init_norms_rat (S, cpoly, &si); /* update the sieve array */
        tnormr = seconds () - tnormr;
        fprintf (stderr, "# Initializing rat. norms took %1.1fs\n", tnormr);

        /* call the siever */
#ifndef USE_BUCKETS
        sieve_random_access(S, fb_rat, &si);
#else
        sieve_buckets(S, fb_rat, &si);
#endif
        ts = seconds () - ts;
        fprintf (stderr, "# Total sieving in %1.1fs\n", ts);

        /* add reports to list */
        tf = seconds ();
        sieve_report_list_accumulate (Reports, cpoly, &si, S);
        tf = seconds () - tf;

        fprintf (stderr, "# Number of reports: %lu", Reports->size);
        fprintf (stderr, " (rat. bit size %lu, alg. bit size %lu)\n",
                 Reports->r_bitsize, Reports->a_bitsize);

        /* factor survivors */
        tf = seconds ();
        reports = 0;
        while (Reports->size > 0)
          reports += build_norms_rational (cpoly, &si, Reports, &P_rat, &P_alg,
                                           fb_alg);
        sieve_report_list_reset (Reports);
        tf = seconds () - tf;
        tot_reports += reports;
        tq = seconds() - tq;
        fprintf (stderr, "# Time for this (q,rho): %1.1fs [norm %1.1f,"
                 " sieving %1.1f, factor %1.1f]\n", tq, tnorma + tnormr, ts,
                 tf);
        ttn = ttn + tnorma + tnormr;
        tts = tts + ts;
        ttf = ttf + tf;
        if (tot_reports != 0)
          fprintf (stderr, "# Reports for this (q,rho): %lu, total %lu,"
                   " rate %1.2fs/r\n#\n", reports, tot_reports,
                   (seconds () - t0) / (double) tot_reports);
}

 end:
    t0 = seconds () - t0;
    fprintf (stderr, "# Average J = %1.0f\n", totJ / (double) sq);
    fprintf (stderr, "# Total time %1.1fs [norm %1.1f, sieving %1.1f,"
             " factor %1.1f]\n", t0, ttn, tts, ttf);
    fprintf (stderr, "# Total %lu reports [%1.2fs/r]\n",
             tot_reports, t0 / (double) tot_reports);

    free (fb_alg);
    free (fb_rat);
    sieve_info_clear (&si);
    clear_polynomial (cpoly);
    free (roots);
    sieve_report_list_clear (Reports);
    acc_primes_clear (&P_rat);
    acc_primes_clear (&P_alg);

    return 0;
}
