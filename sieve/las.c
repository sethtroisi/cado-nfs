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

/* As its name says, this is a ugly hack that initializes all lognorms to the
   maximal value (255) on the rational side. But it seems to work well, and to
   miss only about 5% relations wrt a more accurate estimation, thus we enable
   it until we have a fast lognorm computation. */
#define UGLY_HACK

/* Trick to discard lognorms that will probably lead to non L-smooth
   cofactors. Disabled for now since it requires accurate sieving
   (in particular by prime powers and bad primes). */
// #define COFACTOR_TRICK

/* default sieve region side is 2^DEFAULT_I */
#define DEFAULT_I 12

/* default bucket region: 2^15 = 32K == close to L1 size */
#ifndef BUCKET_REGION
#define BUCKET_REGION 15
#endif

/* This parameter controls the size of the buckets, i.e. the number of
 * updates that a bucket can accumulate.
 * This should be something like loglog(lpb) - loglog(I)
 */
#ifndef BUCKET_LIMIT_FACTOR
#define BUCKET_LIMIT_FACTOR 0.8
#endif

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

/* uintmax_t is guaranteed to be larger or equal to uint64_t */
#define strtouint64(nptr,endptr,base) (uint64_t) strtoumax(nptr,endptr,base)

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
sieve_info_init_lognorm (unsigned char *C, unsigned char threshold,
                         unsigned long B MAYBE_UNUSED,
                         unsigned long l MAYBE_UNUSED,
                         double scale MAYBE_UNUSED)
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
  
  /* We want some margin, (see below), so that we can set 255 to discard
   * non-survivors.*/
  si->scale_alg += (cpoly->alambda * (double) cpoly->lpba) / LOG_SCALE;

  fprintf (stderr, "# Alg. side: log(maxnorm)=%1.6f logbase=%1.6f",
           si->scale_alg, exp (si->scale_alg / LOG_MAX));
  // second guard, due to the 255 trick!
  si->scale_alg = (LOG_MAX-GUARD) / si->scale_alg;
  scale = si->scale_alg / LOG_SCALE;
  alg_bound = (unsigned char) (cpoly->alambda * (double) cpoly->lpba * scale)
            + GUARD;
  fprintf (stderr, " bound=%u\n", alg_bound);
  sieve_info_init_lognorm (si->alg_Bound, alg_bound, cpoly->alim, cpoly->lpba,
                           scale);

  //TODO: put back a clean, simple bound for rational side

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
#if 0  // desactivate 255 trick
  si->scale_rat += r / LOG_SCALE;   /* add natural logarithm */
#endif
  fprintf (stderr, "logbase=%1.6f", exp (si->scale_rat / LOG_MAX ));
  /* we subtract again GUARD to avoid that non-reports overlap the report
     region due to roundoff errors */
#if 0
  si->scale_rat = (LOG_MAX - GUARD) / si->scale_rat;
#else
  si->scale_rat = LOG_MAX / si->scale_rat;
#endif
  scale = si->scale_rat / LOG_SCALE;
  rat_bound = (unsigned char) (r*scale) + GUARD;
  fprintf (stderr, " bound=%u\n", rat_bound);
  sieve_info_init_lognorm (si->rat_Bound, rat_bound, cpoly->rlim, cpoly->lpbr,
                           scale);

  // bucket info
  // TODO: be more clever, here.
  si->log_bucket_region = BUCKET_REGION;
  ASSERT_ALWAYS(si->log_bucket_region >= si->logI);
  si->bucket_region = 1<<si->log_bucket_region;
  si->nb_buckets = 1 + ((si->I*si->J) / si->bucket_region);
  si->bucket_limit = (BUCKET_LIMIT_FACTOR)*si->bucket_region;
  fprintf(stderr, "# log_bucket_region = %u\n", si->log_bucket_region);
  fprintf(stderr, "# bucket_region = %u\n", si->bucket_region);
  fprintf(stderr, "# nb_buckets = %u\n", si->nb_buckets);
  fprintf(stderr, "# bucket_limit = %u\n", si->bucket_limit);
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

/*****************************************************************************/


// Compute the inverse of a modulo b, by binary xgcd.
// a must be less than b.
// a is modified in place
// return 1 on succes, 0 on failure
static inline int
invmod(unsigned long *pa, unsigned long b) {
  unsigned long a, u, v, fix;
  int t, lsh;
  a = *pa;

  // FIXME: here, we rely on internal representation of the modul module.
  if ((b & 1UL)==0)
      return modul_inv(pa, pa, &b);

  fix = (b+1)>>1;

  ASSERT (a < b);
  ASSERT (b & 1UL);

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
  if (a != 1)
    return 0;

  // Here, the inverse of a is u/2^t mod b.
  while (t>0) {
    unsigned long sig = u & 1UL;
    u >>= 1;
    if (sig)
      u += fix;
    --t;
  }
  *pa = u;
  return 1;
}

// Compute the root r describing the lattice inside the q-lattice
// corresponding to the factor base prime (p,R).
// Formula: r = - (a1-R*b1)/(a0-R*b0) mod p
// In the case where denominator is zero, returns p.
// In the case where denominator is non-invertible mod p, returns p+1.
// Otherwise r in [0,p-1]
static inline fbprime_t
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
    if (!invmod(&den, p))
        return p+1;
    num = num*den;
    return (fbprime_t)(num % ((uint64_t) p));
}


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

    ASSERT (pli->beta > 0);
    ASSERT (pli->delta > 0);
    ASSERT ((pli->alpha <= 0) && (pli->alpha > -I));
    ASSERT ((pli->gamma >= 0) && (pli->gamma < I));
    ASSERT (pli->gamma-pli->alpha >= I);

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

/***************************************************************************/

/********        Main bucket sieving functions                    **********/

void
fill_in_buckets(bucket_array_t BA, factorbase_degn_t *fb, const sieve_info_t * si) {
    double tm = seconds();
    // Loop over all primes in the factor base > I
    while (fb->p < si->I)     
        fb = fb_next (fb); // cannot do fb++, due to variable size !
    while (fb->p != FB_END) {
        unsigned char nr;
        fbprime_t p = fb->p;
        unsigned char logp = fb->plog;

        for (nr = 0; nr < fb->nr_roots; ++nr) {
            fbprime_t r, R;
            R = fb->roots[nr];
            r = fb_root_in_qlattice(p, R, si);
            // Special cases should be quite rare for primes > I,
            // so we ignore them.
            // TODO: should be line sieved in the non-bucket phase?
            // Or should we have a bucket line siever?
            if (r == 0 || r == p)
                continue;
            
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
            // TODO: check the generated assembly, in particular, the
            // push function should be reduced to a very simple step.
            asm("## Inner bucket sieving loop starts here!!!\n");
            while (x < I*si->J) {
                uint32_t i;
                bucket_update_t update;
                update.p = p;
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
            asm("## Inner bucket sieving loop stops here!!!\n");
        }
        fb = fb_next (fb); // cannot do fb++, due to variable size !
    }
    fprintf (stderr, "# medium primes pushed to buckets in %f sec\n", seconds()-tm);
}

static void
apply_one_bucket (unsigned char *S, bucket_array_t BA, int i,
                     sieve_info_t * si MAYBE_UNUSED)
{
    int j = nb_of_updates(BA, i);
    for (; j > 0; --j) {
       bucket_update_t update = get_next_bucket_update(BA, i);
       S[update.x] -= update.logp;
#ifdef TRACE_K
       if ((update.x + i * si->bucket_region) == TRACE_K)
           fprintf(stderr, "Subtract %u to S[%u], from BA[%u]\n",
                   update.logp, TRACE_K, i);
#endif
    }
}
/* Initialize lognorms on the rational side for the bucket_region
 * number N.
 * For the moment, nothing clever, wrt discarding (a,b) pairs that are
 * not coprime, except for the line j=0.
 */
static void
init_rat_norms_bucket_region (unsigned char *S, int N, cado_poly cpoly, sieve_info_t *si)
{
    double g1, g0, gi, gj, norm MAYBE_UNUSED;
    int i MAYBE_UNUSED, halfI MAYBE_UNUSED;
    unsigned int j, lastj;

    /* G_q(i,j) = g1*(a0*i+a1*j)+g0*(b0*i+b1*j)
       = (g1*a0+g0*b0)*i + (g1*a1+g0*b1)*j */

    g1 = mpz_get_d (cpoly->g[1]);
    g0 = mpz_get_d (cpoly->g[0]);
    gi = g1 * (double) si->a0 + g0 * (double) si->b0;
    gj = g1 * (double) si->a1 + g0 * (double) si->b1;

    /* bucket_region is a multiple of I. Let's find the starting j
     * corresponding to N and the last j.
     */
    halfI = si->I / 2;
    j = (N*si->bucket_region) >> si->logI;
    lastj = j + (si->bucket_region >> si->logI);
#ifdef UGLY_HACK
    memset (S, 255, (lastj - j) * si->I);
#else
    for (; j < lastj; ++j) {
        if (j == 0) {
            // compute only the norm for i = 1. Everybody else is 255.
            norm = gi;
            norm = fabs (norm);
            norm = log (norm);
            norm = norm * si->scale_rat;
            for (i = -halfI; i < halfI; i++) 
                *S++ = UCHAR_MAX;
            S[-halfI + 1] = GUARD + (unsigned char) (norm);
        } else {
            for (i = -halfI; i < halfI; i++) {
                norm = gi * (double) i + gj * (double) j;
                norm = fabs (norm);
                norm = log (norm);
                norm = norm * si->scale_rat;
                *S++ = GUARD + (unsigned char) (norm);
            }
        }
    }
#endif
}

/* Initialize lognorms on the algebraic side for the bucket_region
 * number N.
 * Only the survivors of the rational sieve will be initialized, the
 * others are set to 255. Case GCD(i,j)!=1 also gets 255.
 * return the number of reports (= number of norm initialisations)
 */
static int
init_alg_norms_bucket_region (unsigned char *S, int N, cado_poly cpoly, sieve_info_t *si)
{
  int i, halfI;
  unsigned int j, lastj, k, d = cpoly->degree;
  double *t, *u, invq, powj, norm;
  int report = 0;

  /* si->scale = LOG_MAX / log(max |F(a0*i+a1*j, b0*i+b1*j)|) */
  
  /* si->fij is the f(i,j) polynomial at (0, 0) */

  halfI = si->I >> 1;
  invq = 1.0 / (double) si->q;

  t = (double*) malloc ((d + 1) * sizeof (double));
  u = (double*) malloc ((d + 1) * sizeof (double));
  for (k = 0; k <= d; k++)
    t[k] = mpz_get_d (si->fij[k]) * invq;
    
  /* bucket_region is a multiple of I. Let's find the starting j
   * corresponding to N and the last j.
   */
  j = (N*si->bucket_region) >> si->logI;
  lastj = j + (si->bucket_region >> si->logI);
  S += halfI;
  for (; j < lastj; j++) {
      /* scale by j^(d-k) the coefficients of fij */
      for (k = 0, powj = 1.0; k <= d; k++, powj *= (double) j)
          u[d - k] = t[d - k] * powj;
      
      /* now compute norms */
      for (i = -halfI; i < halfI; i++) {
        if (si->rat_Bound[S[i]]) {
              norm = fpoly_eval (u, d, (double) i);
              norm = fabs (norm);
              norm = log (norm);
              norm = norm * si->scale_alg;
              S[i] = GUARD + (unsigned char) (norm);
              report++;
          } else {
              S[i] = UCHAR_MAX;
          }
      }
      S += si->I;
    }

  free (t);
  free (u);
  return report;
}

typedef struct {
    fbprime_t p;
    fbprime_t r;        // in ] 0, p [
    int next_position;  // start of the sieve for next bucket_region 
    unsigned char logp;
} small_typical_prime_data_t;

typedef struct {
    // nice primes
    small_typical_prime_data_t *nice_p;
    int nb_nice_p;
    // primes to sieve horizontally (r = oo)
    // ...
    // primes to sieve vertically (r = 0)
    // ...
} small_sieve_data_t;

void clear_small_sieve(small_sieve_data_t ssd) {
    free(ssd.nice_p);
}

// Prepare sieving of small primes: initialize a small_sieve_data_t
// structure to be used thereafter during sieving each region.
void init_small_sieve(small_sieve_data_t *ssd, factorbase_degn_t *fb,
        sieve_info_t *si)
{
    factorbase_degn_t *fb_sav = fb;
    int n = 0;
    unsigned int I = si->I;

    // Count prime ideals.
    while (fb->p != FB_END && fb->p <= I) {
        n += fb->nr_roots;
        fb = fb_next (fb); // cannot do fb++, due to variable size !
    }
    fb = fb_sav;
    // allocate space for these. n is an upper bound, since some of the
    // ideals might become special ones.
    ssd->nice_p = (small_typical_prime_data_t *)malloc(n*sizeof(small_typical_prime_data_t));
    n = 0;
    // Do another pass on fb, to fill in the data
    while (fb->p != FB_END && fb->p <= I) {
        int nr;
        fbprime_t p = fb->p;
        for (nr = 0; nr < fb->nr_roots; ++nr) {
            fbprime_t r;
            r = fb_root_in_qlattice(p, fb->roots[nr], si);
            if (r == 0) {
                // TODO: doit!
                // special_case_0(S, p, logp, si);
                continue;
            } 
            if (r == p) {
                // TODO: doit!
                // special_case_p(S, p, logp, si);
                continue;
            } 
            if (r == p+1) {
                // TODO:
                // Do something with those???
                continue;
            }

            // Fill in data for this nice prime
            ssd->nice_p[n].p = p;
            ssd->nice_p[n].r = r;
            ssd->nice_p[n].logp = fb->plog;
            ssd->nice_p[n].next_position = (I>>1)%p;
            n++;
        }
        fb = fb_next (fb); // cannot do fb++, due to variable size !
    }
    ssd->nb_nice_p = n;
}

// Sieve small primes (up to p < I) of the factor base fb in the sieve
// region number N.
void sieve_small_bucket_region(unsigned char *S, factorbase_degn_t *fb,
        int N, small_sieve_data_t ssd, sieve_info_t *si)
{
    const uint32_t I = si->I;
    unsigned char *S_ptr;
    unsigned long j, nj;
    int n;

    nj = (si->bucket_region >> si->logI);

    for (n = 0; n < ssd.nb_nice_p; ++n) {
        fbprime_t p, r;
        unsigned char logp;
        int x;
        p = ssd.nice_p[n].p;
        r = ssd.nice_p[n].r;
        logp = ssd.nice_p[n].logp;
        x = ssd.nice_p[n].next_position;
        int i, i0;
        S_ptr = S;
        i0 = x & ((1<<si->log_bucket_region) - 1);
        ASSERT(i0 < p);
        for (j = 0; j < nj; ++j) {
            for (i = i0 ; i < I; i += p) {
                S_ptr[i] -= logp;
            }
            i0 += r;
            if (i0 >= p)
                i0 -= p;
            S_ptr += I;
        }
        ssd.nice_p[n].next_position = (N+1)*si->bucket_region + i0;
    }
}

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

// print a comma-separated list of factors.
void factor_list_fprint(FILE *f, factor_list_t fl) {
    int i;
    for (i = 0; i < fl.n-1; ++i)
        fprintf(f, "%" PRIx64 ",", fl.fac[i]);
    fprintf(f, "%" PRIx64, fl.fac[fl.n-1]);
}

/* L is a list of bad primes (ended with 0) */
static void
trial_div (factor_list_t *fl, mpz_t norm, bucket_array_t BA, int N, int x,
           factorbase_degn_t *fb, uint32_t I, uint32_t *L)
{
  /* remove bad primes */
  while (L[0] != 0)
    {
      while (mpz_divisible_ui_p (norm, L[0])) {
        fl->fac[fl->n] = L[0];
        fl->n++;
        ASSERT_ALWAYS(fl->n <= FL_MAX_SIZE);
        mpz_divexact_ui (norm, norm, L[0]);
      }
      L ++;
    }

    // remove primes in fb_alg that are less than I
    while (fb->p != FB_END && fb->p <= I) {
        while (mpz_divisible_ui_p (norm, fb->p)) {
            fl->fac[fl->n] = fb->p;
            fl->n++;
            ASSERT_ALWAYS(fl->n <= FL_MAX_SIZE);
            mpz_divexact_ui (norm, norm, fb->p);
        }
        fb = fb_next (fb); // cannot do fb++, due to variable size !
    }

    // remove primes in BA that map to x
    rewind_bucket(BA, N);
    int nb;
    bucket_update_t update;
    for (nb = nb_of_updates(BA, N); nb > 0; --nb) {
        update = get_next_bucket_update(BA, N);
        if (update.x == x) {
            unsigned long p = update.p;
            while (mpz_divisible_ui_p (norm, p)) {
                fl->fac[fl->n] = p;
                fl->n++;
                ASSERT_ALWAYS(fl->n <= FL_MAX_SIZE);
                mpz_divexact_ui (norm, norm, p);
            }
        }
    }
}

void xToAB(int64_t *a, uint64_t *b, int x, sieve_info_t * si);
int factor_leftover_norm (mpz_t n, unsigned int b,
        mpz_array_t* const factors, uint32_array_t* const multis);
static void
eval_fij (mpz_t v, mpz_t *f, unsigned int d, long i, unsigned long j);

static int
factor_survivors (unsigned char *S, int N, bucket_array_t rat_BA,
                  bucket_array_t alg_BA, factorbase_degn_t *fb_rat,
                  factorbase_degn_t *fb_alg, cado_poly cpoly, sieve_info_t *si,
                  unsigned long *survivors)
{
    int x;
    int64_t a;
    uint64_t b;
    int cpt = 0;
    int surv = 0;
    mpz_t alg_norm, rat_norm;
    factor_list_t alg_factors, rat_factors;
    mpz_array_t *f_r = NULL, *f_a = NULL;    /* large prime factors */
    uint32_array_t *m_r = NULL, *m_a = NULL; /* corresponding multiplicities */

    f_r = alloc_mpz_array (8);
    f_a = alloc_mpz_array (8);
    m_r = alloc_uint32_array (8);
    m_a = alloc_uint32_array (8);

    factor_list_init(&alg_factors);
    factor_list_init(&rat_factors);
    mpz_init(alg_norm);
    mpz_init(rat_norm);

    for (x = 0; x < si->bucket_region; ++x) {
        if (!(si->alg_Bound[S[x]]))
            continue;
        // Compute algebraic and rational norms.
        xToAB(&a, &b, x + N*si->bucket_region, si);
        if (b == 0 || gcd_ul ((a > 0) ? a : -a, b) != 1)
          continue;
        surv++;
        eval_fij(alg_norm, cpoly->f, cpoly->degree, a, b);
        mpz_divexact_ui(alg_norm, alg_norm, si->q);
        // Trial divide algebraic norm
        trial_div(&alg_factors, alg_norm, alg_BA, N, x, fb_alg, si->I,
                  si->abadprimes);
        if (factor_leftover_norm(alg_norm, cpoly->lpba, f_a, m_a)) {
            eval_fij(rat_norm, cpoly->g, 1, a, b);
            // Trial divide rational norm
            trial_div(&rat_factors, rat_norm, rat_BA, N, x, fb_rat, si->I,
                      si->rbadprimes);
            if (factor_leftover_norm(rat_norm, cpoly->lpbr, f_r, m_r)) {
                unsigned int i, j;
                printf ("%" PRId64 ",%" PRIu64 ":", a, b);
                factor_list_fprint(stdout, rat_factors);
                for (i = 0; i < f_r->length; ++i)
                    for (j = 0; j < m_r->data[i]; j++)
                        gmp_printf(",%Zx", f_r->data[i]);
                printf (":");
                factor_list_fprint(stdout, alg_factors);
                for (i = 0; i < f_a->length; ++i)
                    for (j = 0; j < m_a->data[i]; j++)
                        gmp_printf(",%Zx", f_a->data[i]);
                printf ("\n");
                fflush (stdout);
                cpt++;
            }
            rat_factors.n = 0;
        }
        alg_factors.n = 0;
    }
    survivors[0] += surv;
    mpz_clear(alg_norm);
    mpz_clear(rat_norm);
    factor_list_clear(&alg_factors);
    factor_list_clear(&rat_factors);
    clear_mpz_array (f_r);
    clear_mpz_array (f_a);
    clear_uint32_array (m_r);
    clear_uint32_array (m_a);
    return cpt;
}


/****************************************************************************/

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
    double t0, tfb, tq, ttn, tts, ttf;
    uint64_t q0 = 0, q1 = 0, rho = 0;
    unsigned long *roots, nroots, tot_reports = 0, survivors0, survivors1;
    factorbase_degn_t * fb_alg, * fb_rat;
    int checknorms = 0; /* factor or not the remaining norms */
    int I = DEFAULT_I, i;
    unsigned long sq = 0;
    double totJ = 0.0;

    fprintf (stderr, "# %s.r%s", argv[0], REV);
    for (i = 1; i < argc; i++)
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

    /* special q (and root rho) */
    roots = (unsigned long*) malloc (cpoly->degree * sizeof (unsigned long));
    q0 --; /* so that nextprime gives q0 if q0 is prime */
    nroots = 0;
    tot_reports = 0;
    ttn = tts = ttf = 0.0;
    t0 = seconds ();
    fprintf (stderr, "#\n");
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
            if (nroots > 0)
              {
                fprintf (stderr, "### q=%" PRIu64 ": root%s", q0,
                         (nroots == 1) ? "" : "s");
                for (i = 1; i <= (int) nroots; i++)
                  fprintf (stderr, " %lu", roots[nroots-i]);
                fprintf (stderr, "\n");
              }
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
        /* precompute the skewed polynomial */
        fij_from_f (si.fij, cpoly->f, cpoly->degree, si.a0, si.a1, si.b0, si.b1);

        /* checks the value of J */
        sieve_info_update (&si, cpoly->skew);
        totJ += (double) si.J;

        /* Allocate alg buckets */
        bucket_array_t alg_BA;
        alg_BA = init_bucket_array(si.nb_buckets, si.bucket_limit);

        /* Fill in alg buckets */
        fill_in_buckets(alg_BA, fb_alg, &si);

        /* Allocate rational buckets */
        bucket_array_t rat_BA;
        rat_BA = init_bucket_array(si.nb_buckets, si.bucket_limit);

        /* Fill in rational buckets */
        fill_in_buckets(rat_BA, fb_rat, &si);

        /* Initialize data for sieving small primes */
        small_sieve_data_t ssd_alg, ssd_rat;
        init_small_sieve(&ssd_rat, fb_rat, &si);
        init_small_sieve(&ssd_alg, fb_alg, &si);

        /* Process each bucket region, one after the other */
        unsigned char *S;
        S = (unsigned char *)malloc(si.bucket_region*sizeof(unsigned char));
        int reports = 0;
        survivors0 = survivors1 = 0;
        for (i = 0; i < si.nb_buckets; ++i) {
            /* Init rational norms */
            ttn -= seconds ();
            init_rat_norms_bucket_region(S, i, cpoly, &si);
            ttn += seconds ();
            /* Apply rational bucket */
            apply_one_bucket(S, rat_BA, i, &si);
            /* Sieve small rational primes */
            sieve_small_bucket_region(S, fb_rat, i, ssd_rat, &si);

            /* Init algebraic norms */
            ttn -= seconds ();
            survivors0 += init_alg_norms_bucket_region(S, i, cpoly, &si);
            ttn += seconds ();
            /* Apply algebraic bucket */
            apply_one_bucket(S, alg_BA, i, &si);
            /* Sieve small algebraic primes */
            sieve_small_bucket_region(S, fb_alg, i, ssd_alg, &si);

            /* Factor survivors */
            ttf -= seconds ();
            reports += factor_survivors (S, i, rat_BA, alg_BA, fb_rat,
                                            fb_alg, cpoly, &si, &survivors1);
            ttf += seconds ();
        }
        clear_small_sieve(ssd_rat);
        clear_small_sieve(ssd_alg);
        tot_reports += reports;
        fprintf (stderr, "# %lu survivors after rational side,", survivors0);
        fprintf (stderr, " %lu after algebraic side\n", survivors1);
        fprintf (stderr, "# %d relation(s) for (%" PRIu64 ",%" PRIu64
                 "), total %lu [%1.3fs/r]\n#\n", reports, si.q, si.rho,
                 tot_reports, (seconds () - t0) / (double) tot_reports);
        clear_bucket_array(alg_BA);
        clear_bucket_array(rat_BA);
        free(S);
      } // end of loop over special q ideals.

 end:
    t0 = seconds () - t0;
    fprintf (stderr, "# Average J = %1.0f\n", totJ / (double) sq);
    tts = t0 - (ttn + ttf);
    fprintf (stderr, "# Total time %1.1fs [norm %1.1f, sieving %1.1f,"
             " factor %1.1f]\n", t0, ttn, tts, ttf);
    fprintf (stderr, "# Total %lu reports [%1.3fs/r]\n",
             tot_reports, t0 / (double) tot_reports);

    free (fb_alg);
    free (fb_rat);
    sieve_info_clear (&si);
    clear_polynomial (cpoly);
    free (roots);

    return 0;
}
