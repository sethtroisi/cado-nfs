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
#include "basicnt.h"         /* for bin_gcd */
#include "ecm/facul.h"
#include "bucket.h"
#include "trialdiv.h"

#ifdef SSE_NORM_INIT
#include <emmintrin.h>
#endif


/* As its name says, this is a ugly hack that initializes all lognorms to the
   maximal value (255) on the rational side. But it seems to work well, and to
   miss only about 7% to 8% relations wrt a more accurate estimation. */
//#define UGLY_HACK

/* number of bits used to estimate the norms */
#define NORM_BITS 10

/* define PROFILE to keep certain function from being inlined, in order to
   make them show up on profiler output */
//#define PROFILE

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

/* Some functions should not be inlined when we profile or it's hard or
   impossible to tell them apart from the rest in the profiler output */
#ifdef PROFILE
#define NOPROFILE_INLINE
#define NOPROFILE_STATIC
#else
#define NOPROFILE_INLINE static inline
#define NOPROFILE_STATIC static
#endif

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
    int bucket_thresh;    // bucket sieve primes >= bucket_thresh
    int bucket_region;    // should be around L1 cache size, a power of 2
                          // and a multiple of I.
    int log_bucket_region;
    int nb_buckets;
    int bucket_limit;   // maximal number of bucket_reports allowed in one bucket.
    unsigned int degree;   /* polynomial degree */
    double scale_alg;      /* norm scale used on the algebraic side */
    double scale_rat;      /* norm scale used on the rational side */
    double logmax_alg;     /* norms on the alg. side are < 2^logmax_alg */
    double logmax_rat;     /* norms on the rat. side are < 2^logmax_rat */
    mpz_t *fij;       /* coefficients of F(a0*i+a1*j, b0*i+b1*j) */
    double *fijd;     /* coefficients of F_q/q */
    double *tmpd;     /* temporary array */
    double B;         /* bound for the norm computation */
    unsigned char alg_Bound[256]; /* non-zero for good lognorms */
    unsigned char rat_Bound[256]; /* non-zero for good lognorms */
    int checknorms;          /* if non-zero, completely factor the potential
                                relations */
    fbprime_t *trialdiv_primes_alg;
    fbprime_t *trialdiv_primes_rat;
    trialdiv_divisor_t *trialdiv_data_alg;
    trialdiv_divisor_t *trialdiv_data_rat;
    unsigned char S_rat[1 << NORM_BITS];
    unsigned char S_alg[1 << NORM_BITS];
    facul_strategy_t *strategy;
} sieve_info_t;

/* Forward declarations of functions */
static void xToAB(int64_t *a, uint64_t *b, const unsigned int x,
                  const sieve_info_t *si);
/* Test if entry x in bucket region n is divisible by p */
void test_divisible_x (const fbprime_t p, const unsigned long x, const int n,
		       const sieve_info_t *si);
int factor_leftover_norm (mpz_t n, unsigned int b, mpz_array_t* const factors,
			  uint32_array_t* const multis,
			  facul_strategy_t *strategy);
NOPROFILE_STATIC void
eval_fij (mpz_t v, const mpz_t *f, const unsigned int d, const long i,
	  const unsigned long j);


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
    k1 = (unsigned char) (2.0 * log2 ((double) B) * scale + 0.5)
      + GUARD;
    for (k = k0 + GUARD; k + GUARD <= k1 && k < 256; k++)
      C[k] = 0;
  }
#endif
}

static void
sieve_info_init (sieve_info_t *si, cado_poly cpoly, int logI, uint64_t q0)
{
  unsigned int d = cpoly->degree;
  unsigned int k;
  double r, scale;
  unsigned char alg_bound, rat_bound;

  si->degree = d;
  si->fij = (mpz_t*) malloc ((d + 1) * sizeof (mpz_t));
  si->fijd = (double*) malloc ((d + 1) * sizeof (double));
  si->tmpd = (double*) malloc ((d + 1) * sizeof (double));
  for (k = 0; k <= d; k++)
    mpz_init (si->fij[k]);
  si->logI = logI;
  si->I = 1 << si->logI;
  si->J = 1 << (si->logI - 1);

  /* initialize bounds for the norm computation, see lattice.tex */
  si->B = sqrt (2.0 * (double) q0 / (cpoly->skew * sqrt (3.0)));
  si->logmax_alg = get_maxnorm (cpoly, si, q0); /* log2(max norm) */

  /* We want some margin, (see below), so that we can set 255 to discard
   * non-survivors.*/
  scale = si->logmax_alg + cpoly->alambda * (double) cpoly->lpba;

  fprintf (stderr, "# Alg. side: log2(maxnorm)=%1.2f logbase=%1.6f",
           scale, exp2 (scale / LOG_MAX));
  // second guard, due to the 255 trick!
  scale = (LOG_MAX - GUARD) / scale;
  alg_bound = (unsigned char) (cpoly->alambda * (double) cpoly->lpba *  scale)
            + GUARD;
  fprintf (stderr, " bound=%u\n", alg_bound);
  sieve_info_init_lognorm (si->alg_Bound, alg_bound, cpoly->alim, cpoly->lpba,
                           scale);
  si->scale_alg = scale;

  /* similar bound on the rational size: |a| <= s*I*B and |b| <= I*B */
  scale = fabs (mpz_get_d (cpoly->g[1])) * cpoly->skew
        + fabs (mpz_get_d (cpoly->g[0]));
  scale *= si->B * (double) si->I;
  si->logmax_rat = scale = log2 (scale);
  /* on the rational side, we want that the non-reports on the algebraic
     side, which are set to 255, remain over the report bound R, even if
     the rational norm is totally smooth. For this, we simply add R to the
     maximal lognorm to compute the log base */
  r = cpoly->rlambda * (double) cpoly->lpbr; /* base-2 logarithm of the
                                                report bound */
  fprintf (stderr, "# Rat. side: log2(maxnorm)=%1.2f ", scale);
  fprintf (stderr, "logbase=%1.6f", exp2 (scale / LOG_MAX ));
  /* we subtract again GUARD to avoid that non-reports overlap the report
     region due to roundoff errors */
  si->scale_rat = LOG_MAX / scale;
  rat_bound = (unsigned char) (r * si->scale_rat) + GUARD;
  fprintf (stderr, " bound=%u\n", rat_bound);
  sieve_info_init_lognorm (si->rat_Bound, rat_bound, cpoly->rlim, cpoly->lpbr,
                           si->scale_rat);

  // bucket info
  // TODO: be more clever, here.
  si->log_bucket_region = BUCKET_REGION;
  si->bucket_thresh = si->I; /* Default value */
  ASSERT_ALWAYS(si->log_bucket_region > si->logI); /* so that we sieve an even
                                                      number of rows */
  si->bucket_region = 1<<si->log_bucket_region;
  si->nb_buckets = 1 + (si->I * si->J - 1) / si->bucket_region;
  si->bucket_limit = (BUCKET_LIMIT_FACTOR)*si->bucket_region;
  fprintf(stderr, "# log_bucket_region = %u\n", si->log_bucket_region);
  fprintf(stderr, "# bucket_region = %u\n", si->bucket_region);
  fprintf(stderr, "# nb_buckets = %u\n", si->nb_buckets);
  fprintf(stderr, "# bucket_limit = %u\n", si->bucket_limit);
}


/* Finds prime factors p < lim of n and returns a pointer to a zero-terminated
   list of those factors. Repeated factors are stored only once. */
static fbprime_t *
factor_small (mpz_t n, fbprime_t lim)
{
  unsigned long p;
  unsigned long l; /* number of bad primes */
  fbprime_t *f;

  l = 0;
  f = (fbprime_t*) malloc (sizeof (fbprime_t));
  for (p = 2; p <= lim; p = getprime (p))
    {
      if (mpz_divisible_ui_p (n, p))
        {
          l ++;
          f = (fbprime_t*) realloc (f, (l + 1) * sizeof (fbprime_t));
          f[l - 1] = p;
        }
    }
  f[l] = 0; /* end of list marker */
  getprime (0);
  return f;
}

static void
sieve_info_update (sieve_info_t *si, const int verbose)
{
  if (verbose)
    fprintf (stderr, "# I=%u; J=%u\n", si->I, si->J);

  /* update number of buckets */
  si->nb_buckets = 1 + (si->I * si->J - 1) / si->bucket_region;
}

static void
sieve_info_clear (sieve_info_t *si)
{
  unsigned int d = si->degree;
  unsigned int k;

  for (k = 0; k <= d; k++)
    mpz_clear (si->fij[k]);
  free (si->fij);
  free (si->fijd);
  free (si->tmpd);
}

/*****************************************************************************/

// Redc_32 based on 64-bit arithmetic
// Assume:
//   * p is an odd prime < 2^32. FIXME: p < 2^31? (see below)
//   * invp is -1/p mod 2^32.
//   * x is some integer in [0, 2^32*p[
// Compute:
//   * x/2^32 mod p as an integer in [0, p[
static inline uint32_t
redc_u32(const uint64_t x, const uint32_t p, const uint32_t invp)
{
    uint32_t t = x*invp;                            /* t = x * invp mod 2^32 */
    uint32_t y = (x + (uint64_t)t*(uint64_t)p) >> 32; /* x + t*p is bounded
            by 2^32*p-1+(2^32-1)*p < 2*2^32*p: we might want p < 2^31 so that
            there is no overflow */
    if (y >= p)
        y -= p;
    ASSERT(y<p);
    return y;
}

// Signed redc_32 based on 64-bit arithmetic
// Assume:
//   * p is an odd prime < 2^32.
//   * invp is -1/p mod 2^32.
//   * x is some signed integer in ]-2^32*p, 2^32*p[
// Compute:
//   * x/2^32 mod p as an integer in [0, p[
static inline uint32_t
redc_32(const int64_t x, const uint32_t p, const uint32_t invp)
{
    uint32_t t = ((uint64_t)x)*invp;
    uint32_t y = (x + (uint64_t)t*(uint64_t)p) >> 32;
    // might be too large by p, or too small by p.
    if ((int32_t)y < 0)
        y += p;
    else if (y >= p)
        y -= p;
#ifndef NDEBUG
    if (y >= p) {
        fprintf(stderr, "BUG in redc_32. x = %" PRId64
                " p = %u, invp = %u, y = %d\n", x, p, invp, y);
        ASSERT(0);
    }
#endif

    return y;
}

// binary gcd for unsigned long

#if 0
static inline unsigned long
adjust_right_drop1 (unsigned long a)
{
  __asm__ ("1: shr $1, %0 # adjust_right_drop1\n\t"
	   "jc 2f\n\t"
	   "shr $1, %0\n\t"
	   "jnc 1b\n\t"
	   "2:"
	   : "+r" (a)::"cc");
  return a;
}

static inline unsigned long
bingcd(unsigned long a, unsigned long b) {
    int t, lsh;
    unsigned long A, B;

    if (UNLIKELY(a == 0))
        return b;
    if (UNLIKELY(b == 0))
        return a;

    t = ctzl(a);
    A = a >> (t + 1);
    lsh = ctzl(b);
    B = b >> (lsh + 1);
    t = MIN(t, lsh);

    // now a>>t and b>>t are odd, and the 2-part of the gcd is 2^t.
    // A = (a-1)/2, B = (b-1)/2
    // make A < B.
    if (UNLIKELY(A == B))
      return (2 * A + 1)<<t;
    if (A > B) {
        unsigned long tmp = A;
        A = B;
        B = tmp;
    }
    if (UNLIKELY(A == 0))
      return 1<<t;

    do {
        do {
	    B = adjust_right_drop1 (B - A);
        } while (B > A);
        if (A == B)
	    break;
        do {
	    A = adjust_right_drop1 (A - B);
        } while (B < A);
    } while (A != B);
    return (2*A+1)<<t;
}
#endif


/* Define to 0 or 1. Table lookup seems to be slightly faster than
   ctzl() on Opteron and slightly slower on Core2, but doesn't really
   make much difference in either case. */
#define LOOKUP_TRAILING_ZEROS 1

// Compute 2^32/a mod b for b odd,
// and 1/a mod b for b even, by binary xgcd.
// a must be less than b.
// a is modified in place
// return 1 on succes, 0 on failure
NOPROFILE_INLINE int
invmod(uint64_t *pa, uint64_t b) {
#if LOOKUP_TRAILING_ZEROS
  static const unsigned char trailing_zeros[256] =
    {8,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
     5,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
     6,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
     5,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
     7,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
     5,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
     6,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
     5,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0};
#endif
  uint64_t a, u, v, fix, p = b;
  int t, lsh;
  a = *pa;

  if (UNLIKELY(*pa == 0))
    return 0; /* or we get infinite loop */

  if (UNLIKELY(((b & 1UL)==0)))
    {
      modulusul_t m;
      residueul_t r;
      int rc;
      modul_initmod_ul (m, b);
      modul_init (r, m); /* With mod reduction */
      modul_set_ul (r, *pa, m);
      if ((rc = modul_inv(r, r, m)))
        *pa = modul_get_ul (r, m);
      modul_clear (r, m);
      modul_clearmod (m);
      return rc;
    }

  fix = (b+1)>>1;

  ASSERT (a < b);
  ASSERT (b & 1UL);

  u = 1; v = 0; t = 0;

  // make a odd
  lsh = ctzl(a);
  a >>= lsh;
  t += lsh;
  /* v <<= lsh; ??? v is 0 here */

  // Here a and b are odd, and a < b
  do {
    do {
      b -= a; v += u;
#if LOOKUP_TRAILING_ZEROS
      do {
	lsh = trailing_zeros [(unsigned char) b];
	b >>= lsh;
	t += lsh;
	u <<= lsh;
      } while (lsh == 8);
#else
      lsh = ctzl(b);
      b >>= lsh;
      t += lsh;
      u <<= lsh;
#endif
    } while (a<b);
    if (UNLIKELY(a == b))
      break;
    do {
      a -= b; u += v;
#if LOOKUP_TRAILING_ZEROS
      do {
	lsh = trailing_zeros [(unsigned char) a];
	a >>= lsh;
	t += lsh;
	v <<= lsh;
      } while (lsh == 8);
#else
      lsh = ctzl(a);
      a >>= lsh;
      t += lsh;
      v <<= lsh;
#endif
    } while (b < a);
  } while (a != b);
  if (a != 1)
    return 0;

  // Here, the inverse of a is u/2^t mod b.
  while (t>32) {
    unsigned long sig = u & 1UL;
    u >>= 1;
    if (sig)
      u += fix;
    --t;
  }
  while (t < 32)
    {
      u <<= 1;
      if (u >= p)
        u -= p;
      t ++;
    }
  *pa = u;
  return 1;
}

MAYBE_UNUSED
static inline fbprime_t
invmod_po2 (fbprime_t n)
{
  fbprime_t r;

  ASSERT (n % 2 != 0);

  r = (3 * n) ^ 2;
  r = 2 * r - r * r * n;
  r = 2 * r - r * r * n;
  r = 2 * r - r * r * n;
  if (sizeof (fbprime_t) > 4)
    r = 2 * r - r * r * n;

  return r;
}

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

NOPROFILE_INLINE fbprime_t
fb_root_in_qlattice (const fbprime_t p, const fbprime_t R,
        const uint32_t invp, const sieve_info_t * si)
{
  int64_t aux1, aux2;
  uint64_t u, v;
  fbprime_t add;

    /* Handle powers of 2 separately, REDC doesn't like them */
    if (UNLIKELY(p % 2 == 0))
      {
	fbprime_t u, v;
	ASSERT(p == (p & -p)); /* Test that p is power of 2 */
	if (R < p) /* Root in a,b-plane is non-projective */
	  {
	    u = R * si->b1 - si->a1;
	    v = si->a0 - R * si->b0;
	  }
	else /* Root in a,b-plane is projective */
	  {
	    u = si->b1 - (R - p) * si->a1;
	    v = (R - p) * si->a0 - si->b0;
	  }

	if (v % 2 != 0)
	  {
	    /* root in i,j-plane is non-projective */
	    v = invmod_po2 (v);
	    return (u * v) & (p - 1);
	  }
	else
	  {
	    /* root in i,j-plane is projective */
	    u = invmod_po2 (u);
	    return p + ((u * v) & (p - 1));
	  }
      }

    // Use Signed Redc for the computation:
    // Numerator and denominator will get divided by 2^32, but this does
    // not matter, since we take their quotient.

    if (LIKELY(R < p)) /* Root in a,b-plane is affine */
      {
	aux1 = ((int64_t)R)*((int64_t)si->b1) - ((int64_t)si->a1);
	aux2 = ((int64_t)si->a0) - ((int64_t)R)*((int64_t)si->b0);
      }
    else /* Root in a,b-plane is projective */
      {
	aux1 = ((int64_t)si->b1) - ((int64_t)(R - p))*((int64_t)si->a1);
	aux2 = ((int64_t)(R - p))*((int64_t)si->a0) - ((int64_t)si->b0);
      }
    u = redc_32(aux1, p, invp); /* 0 <= u < p */
    v = redc_32(aux2, p, invp); /* 0 <= den < p */

    add = 0;
    if (UNLIKELY(!invmod(&v, p)))
      {
	/* root in i,j-plane is projective */
	if (UNLIKELY(!invmod(&u, p)))
          {
            fprintf (stderr, "Error, root in (i,j)-plane is projective\n");
            exit (1); /* Should never happen! */
          }
	add = p;
      }
    u *= v;
    return (fbprime_t) redc_u32 (u, p, invp) + add;
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
NOPROFILE_INLINE int
reduce_plattice (plattice_info_t *pli, const fbprime_t p, const fbprime_t r,
                 const sieve_info_t * si)
{
    int32_t a0, a1, b0, b1, I, IJ, k;
    int logI;

    I = si->I;
    a0 = -((int32_t)p); a1 = 0;
    b0 = r;  b1 = 1;

    /* subtractive variant of Euclid's algorithm */
    while (b0 >= I)
    {
      /* a0 < 0 < b0 < -a0 */
        do {
            a0 += b0;
            a1 += b1;
        } while (a0 + b0 <= 0);
        /* -b0 < a0 <= 0 < b0 */
        if (-a0 < I)
          {
            if (UNLIKELY(a0 == 0))
              return 0;
            /* Now that |a0| < I, we switch to classical division, since
               if say |a0|=1 and b0 is large, the subtractive variant
               will be very expensive.
               We want b0 + k*a0 < I, i.e., b0 - I + 1 <= k*(-a0),
               i.e., k = ceil((b0-I+1)/a0). */
            k = 1 + (b0 - I) / (-a0);
            b0 += k * a0;
            b1 += k * a1;
            goto case_even;
          }
        /* -b0 < a0 < 0 < b0 */
        do {
            b0 += a0;
            b1 += a1;
        } while (b0 + a0 >= 0);
        /* a0 < 0 <= b0 < -a0 */
    }
    if (UNLIKELY(b0 == 0))
      return 0;
    /* we switch to the classical algorithm here too */
    k = 1 + (-a0 - I) / b0;
    a0 += k * b0;
    a1 += k * b1;

 case_even:
    pli->alpha = a0;
    pli->beta =  a1;
    pli->gamma = b0;
    pli->delta = b1;

    ASSERT (pli->beta > 0);
    ASSERT (pli->delta > 0);
    ASSERT ((pli->alpha <= 0) && (pli->alpha > -I));
    ASSERT ((pli->gamma >= 0) && (pli->gamma < I));
    ASSERT (pli->gamma - pli->alpha >= I);

    // WARNING: Here, we assume a lot on a bound on I,J
    // TODO: clean these bound problems
    logI = si->logI;
    IJ = si->J << logI;
    int64_t aa = ((int64_t)pli->beta << logI) + (int64_t)(pli->alpha);
    if (aa > IJ)
        pli->a = (uint32_t)(INT32_MAX/2);
    else
        pli->a = (uint32_t)aa;
    int64_t cc = ((int64_t)pli->delta << logI) + (int64_t)(pli->gamma);
    if (cc > IJ)
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
fill_in_buckets(bucket_array_t BA, factorbase_degn_t *fb,
		const sieve_info_t * si) {
unsigned char prev_logp = 0;
    // Loop over all primes in the factor base >= bucket_thresh
    while (fb->p != FB_END && fb->p < (fbprime_t) si->bucket_thresh)
        fb = fb_next (fb); // cannot do fb++, due to variable size !
    if (fb->p != FB_END)
      {
        prev_logp = fb->plog;
        bucket_new_logp (BA, fb->plog);
      }
    for (;fb->p != FB_END; fb = fb_next (fb)) {
        unsigned char nr;
        fbprime_t p = fb->p;
        unsigned char logp = fb->plog;

        if (fb->plog > prev_logp)
          {
            prev_logp = fb->plog;
            bucket_new_logp (BA, fb->plog);
          }

        /* If we sieve for special-q's smaller than the algebraic factor
           base bound, the prime p might equal the special-q prime q.
           Note that usually, this doesn't happen on the rational side, since
           the prime q cannot divide both sides, unless q divides Res(f,g). */
        if (UNLIKELY(p == si->q))
          continue;

        for (nr = 0; nr < fb->nr_roots; ++nr) {
            const uint32_t I = si->I;
            const uint32_t logI = si->logI;
            const uint32_t maskI = I-1;
            const uint32_t maskbucket = si->bucket_region - 1;
            const int shiftbucket = si->log_bucket_region;
            const uint32_t IJ = si->I * si->J;
            fbprime_t r, R;

            R = fb->roots[nr];
            r = fb_root_in_qlattice(p, R, fb->invp, si);
            // TODO: should be line sieved in the non-bucket phase?
            // Or should we have a bucket line siever?
            if (UNLIKELY(r == 0))
              {
                /* If r == 0 (mod p), this prime hits for i == 0 (mod p),
                   but since p > I, this implies i = 0 or i > I. We don't
                   sieve i > I. Since gcd(i,j) | gcd(a,b), for i = 0 we
                   only need to sieve j = 1 */
                /* x = j*I + (i + I/2) = I + I/2 */
                bucket_update_t update;
                uint32_t x = si->I + si->I / 2;
                update.x = (uint16_t) (x & maskbucket);
                update.p = p;
                update.logp = logp;
                push_bucket_update(BA, x >> shiftbucket, update);
                continue;
              }
            if (UNLIKELY(r == p))
              {
                /* r == p means root at infinity, which hits for
                   j == 0 (mod p). Since q > I > J, this implies j = 0
                   or j > J. This means we sieve only (i,j) = (1,0) here.
                   Since I < bucket_region, this always goes in bucket 0.
                   FIXME: what about (-1,0)? It's the same (a,b) as (1,0)
                   but which of these two (if any) do we sieve? */
                bucket_update_t update;
                update.x = (uint16_t) si->I / 2 + 1;
                update.p = p;
                update.logp = logp;
                push_bucket_update(BA, 0, update);
                continue;
              }
            if (UNLIKELY(r > p))
              {
                continue;
              }

            plattice_info_t pli;
            if (reduce_plattice(&pli, p, r, si) == 0)
              {
                fprintf (stderr, "# fill_in_buckets: reduce_plattice() "
                         "returned 0 for p = " FBPRIME_FORMAT ", r = "
                         FBPRIME_FORMAT "\n", p, r);
                continue; /* Simply don't consider that (p,r) for now.
                             FIXME: can we find the locations to sieve? */
              }

            // Start sieving from (0,0) which is I/2 in x-coordinate
            uint32_t x;
            x = (I>>1);
            // Skip (0,0), since this can not be a valid report.
            {
                uint32_t i;
                /* i = x & maskI;  ??? x = I/2, maskI = I-1, so x < maskI */
		i = x;
                if (i >= pli.b1)
                    x += pli.a;
                if (i < pli.b0)
                    x += pli.c;
            }
            // TODO: check the generated assembly, in particular, the
            // push function should be reduced to a very simple step.
            bucket_update_t update;
            update.p = p;
            update.logp = logp;
            __asm__("## Inner bucket sieving loop starts here!!!\n");
             while (x < IJ) {
                uint32_t i;
                i = x & maskI;   // x mod I
                /* if both i = x % I and j = x / I are even, then
                   both a, b are even, thus we can't have a valid relation */
                if ((x | (x >> logI)) & 1)
                  {
                    update.x = (uint16_t) (x & maskbucket);
#ifdef PROFILE
		    /* To make it visible in profiler */
		    *(BA.bucket_write[x >> shiftbucket])++ = update;
#else
		    push_bucket_update(BA, x >> shiftbucket, update);
#endif
                  }
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
            __asm__("## Inner bucket sieving loop stops here!!!\n");
        }
    }
}

NOPROFILE_STATIC void
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

/* Knowing the norm on the rational side is bounded by 2^(2^k), compute
   lognorms approximations for k bits of exponent + NORM_BITS-k bits
   of mantissa */
NOPROFILE_STATIC void
init_rat_norms (sieve_info_t *si)
{
  double e, m, norm, h;
  int i, j, k, l, K, L;
  unsigned char *S = si->S_rat;

  k = (int) ceil (log2 (si->logmax_rat));
  K = 1 << k;
  l = NORM_BITS - k;
  L = 1 << l;

  /* extract k bits from the exponent, and l bits from the mantissa */
  h = 1.0 / (double) L;
  for (e = 1.0, i = 0; i < K; i++, e *= 2.0)
    {
      /* e = 2^i for 0 <= i < 2^k */
      for (m = 1.0, j = 0; j < L; j++, m += h)
        {
          /* m = 1 + j/2^l for 0 <= j < 2^l */
          norm = m * e;
          /* Warning: since si->logmax_rat does not usually correspond to
             a power a two, and we consider full binades here, we have to
             take care that values > si->logmax_rat do not wrap around to 0 */
          norm = log2 (norm);
          if (norm >= si->logmax_rat)
            S[(i << l) + j] = 255;
          else
            {
              norm = norm * si->scale_rat;
              S[(i << l) + j] = GUARD + (unsigned char) norm;
            }
        }
    }
}

/* idem as init_rat_norms, but for the algebraic side */
NOPROFILE_STATIC void
init_alg_norms (sieve_info_t *si)
{
  double e, m, norm, h;
  int i, j, k, l, K, L;
  unsigned char *S = si->S_alg;

  k = (int) ceil (log2 (si->logmax_alg));
  K = 1 << k;
  ASSERT_ALWAYS(NORM_BITS >= k);
  l = NORM_BITS - k;
  L = 1 << l;

  /* extract k bits from the exponent, and l bits from the mantissa */
  h = 1.0 / (double) L;
  for (e = 1.0, i = 0; i < K; i++, e *= 2.0)
    {
      for (m = 1.0, j = 0; j < L; j++, m += h)
        {
          norm = m * e;
          norm = log2 (norm);
          if (norm >= si->logmax_alg)
            S[(i << l) + j] = 255;
          else
            {
              norm = norm * si->scale_alg;
              S[(i << l) + j] = GUARD + (unsigned char) norm;
            }
        }
    }
}

/* Initialize lognorms on the rational side for the bucket_region
 * number N.
 * For the moment, nothing clever, wrt discarding (a,b) pairs that are
 * not coprime, except for the line j=0.
 */
NOPROFILE_STATIC void
init_rat_norms_bucket_region (unsigned char *S, int N, cado_poly cpoly,
                              sieve_info_t *si)
{
    double g1, g0, gi, gj, norm MAYBE_UNUSED, gjj;
    int i MAYBE_UNUSED, halfI MAYBE_UNUSED, l;
    unsigned int j, lastj;
    uint64_t mask = (1 << NORM_BITS) - 1;
    union { double z; uint64_t x; } zx[1];

    l = NORM_BITS - (int) ceil (log2 (si->logmax_rat));

    /* G_q(i,j) = g1 * (a0*i+a1*j) + g0 * (b0*i+b1*j)
       = (g1*a0+g0*b0) * i + (g1*a1+g0*b1) * j
       = gi * i + gj * j
    */

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
    /* if j = 0, it will be the first value */
    if (j == 0)
      {
        // compute only the norm for i = 1. Everybody else is 255.
        norm = log2 (fabs (gi));
        norm = norm * si->scale_rat;
        for (i = -halfI; i < halfI; i++)
          *S++ = UCHAR_MAX;
        S[-halfI + 1] = GUARD + (unsigned char) (norm);
        ++j;
      }
    for (; j < lastj; ++j)
      {
        unsigned char n;
        gjj = gj * (double) j;
        zx->z = gjj - gi * (double) halfI;
        __asm__("### Begin rational norm loop\n");
#ifndef SSE_NORM_INIT
        uint64_t y;
        for (i = 0; i < (int) si->I; i++) {
          /* the double precision number 1.0 has high bit 0 (sign),
             then 11-bit biased exponent 1023, and 52-bit mantissa 0 */
          /* the magic constant here is simply 1023*2^52, where
             1023 is the exponent bias in binary64 */
          y = (zx->x - (uint64_t) 0x3FF0000000000000) >> (52 - l);
          n = si->S_rat[y & mask];
          ASSERT (n > 0);
          *S++ = n;
          zx->z += gi;
        }
#else
        union { __v2df dble;
            __v2di intg;
            struct {uint64_t y0; uint64_t y1; } intpair;
        } y_vec;

        __v2df gi_vec = { 2*gi, 2*gi };
        __v2di mask_vec = { mask, mask };
        __v2di cst_vec = { (uint64_t) 0x3FF0000000000000,
            (uint64_t) 0x3FF0000000000000 };

        // in spite of the appearance, only the low word gives the shift
        // count. The high word is ignored.
        __v2di shift_value = { (uint64_t)(52 - l), (uint64_t)(52 - l) };
        __v2df z_vec = { zx->z, zx->z+gi };

        for (i = 0; i < halfI; ++i) {
            y_vec.dble = z_vec;
            y_vec.intg -= cst_vec;
            y_vec.intg = _mm_srl_epi64(y_vec.intg, shift_value);
            y_vec.intg &= mask_vec;
//            *S++ = si->S_rat[((uint32_t *)(&y_vec.intg))[0]];
//             *S++ = si->S_rat[((uint32_t *)(&y_vec.intg))[2]];
            *S++ = si->S_rat[y_vec.intpair.y0];
            *S++ = si->S_rat[y_vec.intpair.y1];
            z_vec += gi_vec;
        }
#endif
        __asm__("### End rational norm loop\n");
      }
#endif
}

#ifdef SSE_NORM_INIT
static inline void
init_fpoly_v2df(__v2df *F, const double *f, const int deg)
{
    int i;
    for (i = 0; i <= deg; ++i) {
        __v2df tmp = { f[i], f[i] };
        F[i] = tmp;
    }
}


static inline __v2df
fpoly_eval_v2df_deg5(const __v2df *f, const __v2df x)
{
    __v2df r;
    r = f[5];
    r = r * x + f[4];
    r = r * x + f[3];
    r = r * x + f[2];
    r = r * x + f[1];
    r = r * x + f[0];
    return r;
}
#endif

/* Initialize lognorms on the algebraic side for the bucket_region
 * number N.
 * Only the survivors of the rational sieve will be initialized, the
 * others are set to 255. Case GCD(i,j)!=1 also gets 255.
 * return the number of reports (= number of norm initialisations)
 */
NOPROFILE_STATIC int
init_alg_norms_bucket_region (unsigned char *S, int N, cado_poly cpoly,
                              sieve_info_t *si)
{
  unsigned int j, lastj, k, d = cpoly->degree;
  double *t, *u, powj, i, halfI;
  int report = 0, l;
  uint64_t mask = (1 << NORM_BITS) - 1;
  union { double z; uint64_t x; } zx[1];
  uint64_t y;
  unsigned char *T = si->S_alg;

  l = NORM_BITS - (int) ceil (log2 (si->logmax_alg));
  halfI = (double) (si->I >> 1);

  t = si->fijd;
  u = si->tmpd;

  /* bucket_region is a multiple of I. Let's find the starting j
   * corresponding to N and the last j.
   */
  j = (N*si->bucket_region) >> si->logI;
  lastj = j + (si->bucket_region >> si->logI);
  for (; j < lastj; j++)
    {
      /* scale by j^(d-k) the coefficients of fij */
      for (k = 0, powj = 1.0; k <= d; k++, powj *= (double) j)
        u[d - k] = t[d - k] * powj;

#ifdef SSE_NORM_INIT
      __v2df u_vec[d];
      init_fpoly_v2df(u_vec, u, d);
      double i0 = 0.0 ;
      unsigned char *S_ptr0 = NULL;
      int cpt = 0;
#endif
      /* now compute norms */
      for (i = -halfI; i < halfI; i += 1.0)
        {
          if (si->rat_Bound[*S])
            {
#ifndef SSE_NORM_INIT
              unsigned char n;
              zx->z = fpoly_eval (u, d, i);
              /* 4607182418800017408 = 1023*2^52 */
              y = (zx->x - (uint64_t) 0x3FF0000000000000) >> (52 - l);
              report++;
              n = T[y & mask];
              ASSERT (n > 0);
              *S++ = n;
#else
              ASSERT(d == 5);
              __v2di mask_vec = { mask, mask };
              __v2di cst_vec = { (uint64_t) 0x3FF0000000000000,
                  (uint64_t) 0x3FF0000000000000 };
              __v2di shift_value = { (uint64_t)(52 - l), (uint64_t)(52 - l) };
              report++;
              if (cpt == 0) {
                  i0 = i;
                  cpt++;
                  S_ptr0 = S++;
              } else if (cpt == 1) {
                  cpt--;
                  __v2df i_vec = {i0, i};
                  union { __v2df dble;
                      __v2di intg;
                      struct {uint64_t y0; uint64_t y1; } intpair;
                  } fi_vec;
                  fi_vec.dble = fpoly_eval_v2df_deg5(u_vec, i_vec);
                  fi_vec.intg -= cst_vec;
                  fi_vec.intg = _mm_srl_epi64(fi_vec.intg, shift_value);
                  fi_vec.intg &= mask_vec;
                  *S_ptr0 = T[fi_vec.intpair.y0];
                  *S++ = T[fi_vec.intpair.y1];
              }
#endif
            }
          else
            *S++ = UCHAR_MAX;
        }
#ifdef SSE_NORM_INIT
      if (cpt == 1) { // odd number of computations
          zx->z = fpoly_eval (u, d, i0);
          y = (zx->x - (uint64_t) 0x3FF0000000000000) >> (52 - l);
          *S_ptr0 = T[y & mask];
      }
#endif
    }

  return report;
}

typedef struct {
    fbprime_t p;
    fbprime_t r;        // in [ 0, p [
    int next_position;  // start of the sieve for next bucket_region
    unsigned char logp; // TODO: word-align?
} small_typical_prime_data_t;

/* Small primes or powers of small primes p^k with projective root.
   These hit at i*v == j*u (mod p^k), but gcd(v, p^k) > 1.
   We store g = gcd(v, p^k), q = p^k / g, and U = u * (v/g)^(-1) (mod q).
   Then, for lines j = k*g, these primes hit at positions i == k*r (mod q).
   We initialise next_position = -I/2 + i % q,
   then within a line sieve next_position + k*q < I/2,
   and update next_position += r + g*I to get the first position to sieve
   in the next suitable line. */
typedef struct {
    fbprime_t g, q, U;
    int next_position;
    unsigned char logp; // TODO: word-align?
} small_bad_prime_data_t;

typedef struct {
    // nice primes
    small_typical_prime_data_t *nice_p;
    int nb_nice_p;
    // primes to sieve horizontally (r = oo)
    small_bad_prime_data_t *bad_p;
    int nb_bad_p;
    // primes to sieve vertically (r = 0)
    // ...
} small_sieve_data_t;

void clear_small_sieve(small_sieve_data_t ssd) {
    free(ssd.nice_p);
    ssd.nice_p = NULL;
    if (ssd.nb_bad_p > 0)
      free(ssd.bad_p);
    ssd.bad_p = NULL;
}

static void
copy_small_sieve (small_sieve_data_t *r, const small_sieve_data_t *s)
{
    r->nb_nice_p = s->nb_nice_p;
    if (s->nb_nice_p > 0)
      {
        const size_t size = s->nb_nice_p * sizeof (small_typical_prime_data_t);
	r->nice_p = (small_typical_prime_data_t *) malloc (size);
	ASSERT (r->nice_p != NULL);
	memcpy (r->nice_p, s->nice_p, size);
      }
    else
      r->nice_p = NULL;

    r->nb_bad_p = s->nb_bad_p;
    if (s->nb_bad_p > 0)
      {
        const size_t size = s->nb_bad_p * sizeof (small_bad_prime_data_t);
        r->bad_p = (small_bad_prime_data_t *) malloc (size);
        ASSERT (r->bad_p != NULL);
        memcpy (r->bad_p, s->bad_p, size);
      }
    else
      r->bad_p = NULL;
}

// Prepare sieving of small primes: initialize a small_sieve_data_t
// structure to be used thereafter during sieving each region.
void init_small_sieve(small_sieve_data_t *ssd, const factorbase_degn_t *fb,
                      const sieve_info_t *si, const char side)
{
    const factorbase_degn_t *fb_sav = fb;
    int size = 0, n;
    const unsigned int thresh = si->bucket_thresh;
    const int verbose = 0;
    const int do_bad_primes = 1;

    ASSERT (side == 'a' || side == 'r');

    // Count prime ideals of factor base primes p < thresh
    while (fb->p != FB_END && fb->p < thresh) {
        size += fb->nr_roots;
        fb = fb_next (fb); // cannot do fb++, due to variable size !
    }
    fb = fb_sav;

    // allocate space for these. n is an upper bound, since some of the
    // ideals might become special ones.
    ssd->nice_p = (small_typical_prime_data_t *)
      malloc(size * sizeof(small_typical_prime_data_t));
    n = 0;
    ssd->nb_bad_p = 0;
    ssd->bad_p = NULL;
    // Do another pass on fb and badprimes, to fill in the data
    // while we have any regular primes or bad primes < thresh left
    while (fb->p != FB_END && fb->p < thresh) {
      const fbprime_t p = fb->p;
      int nr;
      fbprime_t r;

      for (nr = 0; nr < fb->nr_roots; nr++)
        {
	  r = fb_root_in_qlattice (p, fb->roots[nr], fb->invp, si);
	  /* If this root is somehow interesting (projective in (a,b) or
	     in  (i,j) plane), print a message */
	  if (verbose && (fb->roots[nr] >= p || r >= p))
	    fprintf (stderr, "# Side %c, prime " FBPRIME_FORMAT
		     " root " FBPRIME_FORMAT " -> %d\n",
		     side, p, fb->roots[nr], r);
          if (r >= p)
	    {
	      fbprime_t q, g;
	      g = bin_gcd (p, r - p);
	      if (do_bad_primes && g < si->J)
		{
		  // Fill in data for this bad prime
		  /* Shouldn't happen that often so a realloc is ok */
		  ssd->bad_p = (small_bad_prime_data_t *)
		    realloc (ssd->bad_p, (ssd->nb_bad_p + 1) *
			     sizeof (small_bad_prime_data_t));
		  q = p / g;
		  ssd->bad_p[ssd->nb_bad_p].g = g;
		  ssd->bad_p[ssd->nb_bad_p].q = q;
		  if (q == 1)
		    {
		      ssd->bad_p[ssd->nb_bad_p].U = 0;
		      ssd->bad_p[ssd->nb_bad_p].next_position = g * si->I;
		    }
		  else
		    {
		      fprintf (stderr, "# Projective roots in (i,j)-plane for "
			       "prime powers not implemented yet\n");
		      fprintf (stderr, "# p = " FBPRIME_FORMAT
			       ", r = " FBPRIME_FORMAT "\n", p, r);
		      abort (); /* Use a general enough modinv!
				   ssd->bad_p[ssd->nb_bad_p].U = modinv ((r-p) / g, q); */
		      ssd->bad_p[ssd->nb_bad_p].next_position = (si->I >> 1) % q;
		    }
		  ssd->bad_p[ssd->nb_bad_p].logp = fb->plog;
		  ssd->nb_bad_p++;
		}
	    }
	  else
	    {
            // Fill in data for this nice prime
	    ASSERT (n < size);
            ssd->nice_p[n].p = p;
            ssd->nice_p[n].r = r;
            ssd->nice_p[n].logp = fb->plog;
            ssd->nice_p[n].next_position = (si->I >> 1)%p;
            n++;
          }
        }
      fb = fb_next (fb); // cannot do fb++, due to variable size !
    }
    ssd->nb_nice_p = n;
}

// Sieve small primes (up to p < bucket_thresh) of the factor base fb in the
// next sieve region S.
// Information about where we are is in ssd.
void sieve_small_bucket_region(unsigned char *S, const int bucket_nr,
			       small_sieve_data_t ssd, sieve_info_t *si,
			       const unsigned char side)
{
    const uint32_t I = si->I;
    unsigned char *S_ptr;
    unsigned long j, nj;
    int n;
    const int test_divisibility = 0; /* very slow, but nice for debugging */

    nj = (si->bucket_region >> si->logI); /* Nr. of lines per bucket region */
    ASSERT ((nj & 1) == 0);

    /* Handle 2 separately. The general code below assumes odd p (although
       it works correctly, it updates some i,j with gcd(i,j) = 2 for p = 2).
       Also, 2 can be sieved more quickly with long word transfers.
       TODO: use SSE2 */
    for (n = 0 ; n < ssd.nb_nice_p && ssd.nice_p[n].p == 2; n++) {
      unsigned long *S_ptr;
      unsigned long logps;
      const int r = ssd.nice_p[n].r;
      /* If r == 0, sieve all i == 0, j == 1 (mod 2).
         If r == 1, sieve all i == 1, j == 1 (mod 2).
         So either way, skip even j.
         r == 1/0, i.e. i == 1, j == 0 (mod 2) will be sieved as a
         "bad" prime below. */
      ASSERT (r < 2);
      logps = 0UL;
      for (j = 0; j < sizeof (unsigned long); j += 2)
	((unsigned char *)&logps)[r + j] = ssd.nice_p[n].logp;
      ASSERT (I % (4 * sizeof (unsigned long)) == 0);
      S_ptr = (unsigned long *) (S + I); /* Sieve only odd lines */
      for (j = 1; j < nj; j += 2)
      {
        const unsigned long *end = (unsigned long *)((char *) S_ptr + I);

        while (S_ptr < end)
          {
            *(S_ptr) -= logps;
            *(S_ptr + 1) -= logps;
            *(S_ptr + 2) -= logps;
            *(S_ptr + 3) -= logps;
            S_ptr += 4;
          }
        S_ptr = (unsigned long *)((char *) S_ptr + I);
      }
    }

    for ( ; n < ssd.nb_nice_p; ++n) {
        fbprime_t p, r, twop;
        unsigned char logp;
        unsigned int i, i0;
        p = ssd.nice_p[n].p;
        twop = p + p;
        r = ssd.nice_p[n].r;
        logp = ssd.nice_p[n].logp;
        i0 = ssd.nice_p[n].next_position;
        S_ptr = S;
        ASSERT(i0 < p);
        for (j = 0; j < nj; j += 2)
          {
            /* for j even, we sieve only odd i */
            for (i = (i0 & 1) ? i0 : i0 + p; i < I; i += twop)
	      {
		if (test_divisibility && side == 'a')
		  test_divisible_x (p, S_ptr + i - S, bucket_nr, si);
		S_ptr[i] -= logp;
	      }
            i0 += r;
            if (i0 >= p)
              i0 -= p;
            S_ptr += I;
            /* j odd */
            for (i = i0 ; i < I; i += p)
	      {
		if (test_divisibility && side == 'a')
		  test_divisible_x (p, S_ptr + i - S, bucket_nr, si);
		S_ptr[i] -= logp;
	      }
            i0 += r;
            if (i0 >= p)
              i0 -= p;
            S_ptr += I;
          }
        ssd.nice_p[n].next_position = i0;
    }

    /* Sieve the bad primes. We have
       i * g == j * U (mod p^k) where g = p^l. This hits only for g|j,
       then j = j' * g, and i == j' * U (mod p^(k-l)).
       In ssd we have stored g, q = p^(k-l), U, and next_position so that
       S + next_position is the next sieve entry that needs to be sieved.
       So if S + next_position is in the current bucket region, we
       update all  S + next_position + k*q  where  next_position + k*q < I,
       then set next_position = ((next_position % I) + U) % q) + I * g.  */
    for (n = 0; n < ssd.nb_bad_p; ++n) {
      if (!test_divisibility && ssd.bad_p[n].q == 1)
        {
	  /* q = 1, therefore U = 0, and we sieve all entries in lines
	     with g|j, beginning with the line starting at S[next_position] */
          unsigned long logps;
          unsigned int i0 = ssd.bad_p[n].next_position;
	  ASSERT (ssd.bad_p[n].U == 0);
          ASSERT (i0 % I == 0);
          ASSERT (I % (4 * sizeof (unsigned long)) == 0);
	  for (j = 0; j < sizeof (unsigned long); j++)
	    ((unsigned char *)&logps)[j] = ssd.bad_p[n].logp;
          while (i0 < (unsigned int) si->bucket_region)
            {
              unsigned long *S_ptr = (unsigned long *) (S + i0);
              unsigned long *end = S_ptr + I / sizeof (unsigned long);
              unsigned long logps2 = logps;
              if ((i0 & I) == 0) /* Even j coordinate? */
                {
                  /* Yes, update only odd i-coordinates */
                  /* Use array indexing to avoid endianness issues. */
		  for (j = 0; j < sizeof (unsigned long); j += 2)
		    ((unsigned char *)&logps2)[j] = 0;
                }
              while (S_ptr < end)
                {
                  *(S_ptr) -= logps2;
                  *(S_ptr + 1) -= logps2;
                  *(S_ptr + 2) -= logps2;
                  *(S_ptr + 3) -= logps2;
                  S_ptr += 4;
                }
              i0 += ssd.bad_p[n].g * I;
            }
          ssd.bad_p[n].next_position = i0 - si->bucket_region;
        }
      else
	{
	  /* q > 1, more general sieving code. */
	  unsigned int i0 = ssd.bad_p[n].next_position;
	  unsigned int line_idx = i0 & (I - 1U); /* Index in line to sieve */
	  unsigned char *S_ptr = S + i0 - line_idx; /* Start of line to sieve */
	  const fbprime_t g = ssd.bad_p[n].g;
	  const fbprime_t q = ssd.bad_p[n].q;
	  const fbprime_t U = ssd.bad_p[n].U;
	  const unsigned char logp = ssd.bad_p[n].logp;
	  const fbprime_t evenq = (q % 2 == 0) ? q : 2 * q;
	  ASSERT (U < q);
	  while (i0 < (unsigned int) si->bucket_region)
	    {
	      unsigned int i = line_idx;
	      if ((i0 & I) == 0) /* Is j even? */
		{
		  /* Yes, sieve only odd i values */
		  if (i % 2 == 0) /* Make i odd */
		    {
		      ASSERT (q % 2 == 1);
		      i += q;
		    }
		  for ( ; i < I; i += evenq)
		    {
		      if (test_divisibility && side == 'a')
			test_divisible_x (g*q, S_ptr + i - S, bucket_nr, si);
		      S_ptr[i] -= logp;
		    }
		}
	      else
		{
		  for ( ; i < I; i += q)
		    {
		      if (test_divisibility && side == 'a')
			test_divisible_x (g*q, S_ptr + i - S, bucket_nr, si);
		      S_ptr[i] -= logp;
		    }
		}

	      line_idx += U;
	      S_ptr += ssd.bad_p[n].g * I;
	      i0 += g * I + U;
	      if (line_idx >= q)
		{
		  line_idx -= q;
		  i0 -= q;
		}
	    }
	  ssd.bad_p[n].next_position = i0 - si->bucket_region;
	}
    }
}

/* Sieve small primes (p < I, p not in trialdiv_primes list) of the factor
   base fb in the next sieve region S, and add primes and the x position
   where they divide and where there's a sieve report to a bucket (rather
   than subtracting the log norm from S, as during sieving).
   Information about where we are is in ssd.
   Primes in trialdiv_primes must be in increasing order. */
void
resieve_small_bucket_region (bucket_array_t BA, unsigned char *S,
			     small_sieve_data_t *ssd,
			     const sieve_info_t *si,
			     fbprime_t *trialdiv_primes)
{
  const uint32_t I = si->I;
  unsigned char *S_ptr;
  unsigned long j, nj;
  int n, td_idx;
  const int resieve_very_verbose = 0, resieve_very_verbose_bad = 0;

  nj = (si->bucket_region >> si->logI);
  ASSERT ((nj & 1) == 0);

  td_idx = 0;
  for (n = 0 ; n < ssd->nb_nice_p; ++n) {
    const fbprime_t p = ssd->nice_p[n].p;
    while (trialdiv_primes[td_idx] != FB_END && trialdiv_primes[td_idx] < p)
      td_idx++;
    if (trialdiv_primes[td_idx] == FB_END || trialdiv_primes[td_idx] != p)
    {
      fbprime_t r, twop;
      unsigned int i0;
      twop = p + p;
      r = ssd->nice_p[n].r;
      i0 = ssd->nice_p[n].next_position;
      S_ptr = S;
      ASSERT(i0 < p);
      for (j = 0; j < nj; j += 2)
        {
          unsigned int i;
          /* for j even, we sieve only odd i */
          for (i = (i0 & 1) ? i0 : i0 + p ; i < I; i += twop)
            {
              if (S_ptr[i] != 255)
                {
                  bucket_update_t update;
                  unsigned int x = (j << (si->logI)) + i;
                  if (resieve_very_verbose)
                    fprintf (stderr, "resieve_small_bucket_region: root "
                             FBPRIME_FORMAT ",%d divides at x = "
                             "%d = %lu * %u + %d\n",
                             p, r, x, j, 1 << si->logI, i);
                  update.p = p;
                  update.x = x;
                  update.logp = 0;
                  push_bucket_update (BA, 0, update);
                }
            }
          i0 += r;
          if (i0 >= p)
              i0 -= p;
          S_ptr += I;
          /* j odd */
          for (i = i0 ; i < I; i += p)
              if (S_ptr[i] != 255)
                {
                  bucket_update_t update;
                  unsigned int x = ((j + 1) << (si->logI)) + i;
                  if (resieve_very_verbose)
                    fprintf (stderr, "resieve_small_bucket_region: root "
                             FBPRIME_FORMAT ",%d divides at x = "
                             "%d = %lu * %d + %d\n",
                             p, r, x, j + 1, 1 << si->logI, i);
                  update.p = p;
                  update.x = x;
                  update.logp = 0;
                  push_bucket_update (BA, 0, update);
                }
          i0 += r;
          if (i0 >= p)
              i0 -= p;
          S_ptr += I;
        }
      ssd->nice_p[n].next_position = i0;
    }
  }


  /* Resieve bad primes */
  td_idx = 0;
  for (n = 0; n < ssd->nb_bad_p; ++n) {
    const fbprime_t g = ssd->bad_p[n].g;
    while (trialdiv_primes[td_idx] != FB_END && trialdiv_primes[td_idx] < g)
      td_idx++;
    if (trialdiv_primes[td_idx] == FB_END || trialdiv_primes[td_idx] != g)
    {
      /* Test every p-th line, starting at S[next_position] */
      unsigned int i, i0 = ssd->bad_p[n].next_position;
      ASSERT (n == 0 || ssd->bad_p[n - 1].g <= ssd->bad_p[n].g);
      ASSERT (i0 % I == 0); /* make sure next_position points at start
                               of line */
      if (resieve_very_verbose_bad)
        fprintf (stderr, "# resieving bad prime " FBPRIME_FORMAT
                 ", i0 = %u\n", g, i0);
      while (i0 < (unsigned int) si->bucket_region)
        {
          unsigned char *S_ptr = S + i0;
          if ((i0 >> si->logI) % 2 == 0) /* Even j coordinate? */
            {
              /* Yes, test only odd i-coordinates */
              for (i = 1; i < I; i += 2)
                {
                  if (S_ptr[i] != 255)
                    {
                      bucket_update_t update;
                      const unsigned int x = i0 + i;
                      if (resieve_very_verbose_bad)
                        fprintf (stderr, "resieve_small_bucket_region even j: root "
                                 FBPRIME_FORMAT ",inf divides at x = %u\n",
                                 g, x);
                      update.p = g;
                      update.x = x;
                      update.logp = 0;
                      push_bucket_update (BA, 0, update);
                    }
                }
            }
          else
            {
              /* No, test all i-coordinates */
              for (i = 0; i < I; i++)
                {
                  if (S_ptr[i] != 255)
                    {
                      bucket_update_t update;
                      const unsigned int x = i0 + i;
                      if (resieve_very_verbose_bad)
                        fprintf (stderr, "resieve_small_bucket_region odd j: root "
                                 FBPRIME_FORMAT ",inf divides at x = %u\n",
                                 g, x);
                      update.p = g;
                      update.x = x;
                      update.logp = 0;
                      push_bucket_update (BA, 0, update);
                    }
                }
            }
          i0 += g * I;
        }
      ssd->bad_p[n].next_position = i0 - si->bucket_region;
      if (resieve_very_verbose_bad)
        fprintf (stderr, "# resieving: new i0 = %u, bucket_region = %d, "
                         "new next_position = %d\n",
                 i0, si->bucket_region, ssd->bad_p[n].next_position);
    }
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
// assumes there is at least one factor
void factor_list_fprint(FILE *f, factor_list_t fl) {
    int i;
    for (i = 0; i < fl.n-1; ++i)
        fprintf(f, "%" PRIx64 ",", fl.fac[i]);
    fprintf(f, "%" PRIx64, fl.fac[fl.n-1]);
}


/* Function stolen to sieve.c (AK) */
/* Returns 1 if fbptr->p divides norm, 0 otherwise */
static inline int
trialdiv_with_norm (factorbase_degn_t *fbptr, const mpz_t norm)
{
  modulusul_t m;
  residueul_t r;
  size_t i;
  int j;

  modul_initmod_ul (m, (unsigned long) fbptr->p);
  modul_init (r, m);

  for (i = 0; i < mpz_size (norm); i++)
    modul_addredc_ul (r, r, norm->_mp_d[i], fbptr->invp, m);

  j = (modul_get_ul (r, m) == 0);

  modul_clear (r, m);
  modul_clearmod (m);
  return j;
}


/* The entries in BA must be sorted in order of increasing x */
static void
divide_primes_from_bucket (factor_list_t *fl, mpz_t norm, const int x,
                           bucket_array_t BA, const int N)
{
  bucket_update_t update;
  while (!is_end (BA, N)) {
      update = get_next_bucket_update(BA, N);
      if (update.x > x)
        {
          rewind_bucket_by_1 (BA, N);
          break;
        }
      if (update.x == x) {
          unsigned long p = update.p;
          if (!mpz_divisible_ui_p (norm, p)) {
            fprintf (stderr,
                     "# Error, p = %lu does not divide at x = %d\n",
                     p, x);
          } else {
            do {
              fl->fac[fl->n] = p;
              fl->n++;
              ASSERT_ALWAYS(fl->n <= FL_MAX_SIZE);
              mpz_divexact_ui (norm, norm, p);
            } while (mpz_divisible_ui_p (norm, p));
          }
      }
  }
}


NOPROFILE_STATIC void
trial_div (factor_list_t *fl, mpz_t norm, bucket_array_t BA, int N, int x,
           factorbase_degn_t *fb, bucket_array_t resieved,
	   trialdiv_divisor_t *trialdiv_data)
{
    const int trial_div_very_verbose = 0; // (x == 30878);
    fl->n = 0; /* reset factor list */

    if (trial_div_very_verbose)
      gmp_fprintf (stderr, "# trial_div() entry, x = %d, norm = %Zd\n",
                   x, norm);

    // handle 2 separately, if it is in fb
    if (fb->p == 2) {
        int bit = mpz_scan1(norm, 0);
        int i;
        for (i = 0; i < bit; ++i) {
            fl->fac[fl->n] = 2;
            fl->n++;
        }
        if (trial_div_very_verbose)
          gmp_fprintf (stderr, "# x = %d, dividing out 2^%d, norm = %Zd\n",
                       x, bit, norm);
        mpz_tdiv_q_2exp(norm, norm, bit);
        fb = fb_next (fb); // cannot do fb++, due to variable size !
    }

    // remove primes in BA that map to x
    divide_primes_from_bucket (fl, norm, x, BA, N);
    if (trial_div_very_verbose)
      gmp_fprintf (stderr, "# x = %d, after dividing out buckets norm = %Zd\n",
                   x, norm);

    // remove primes found by resieving
    divide_primes_from_bucket (fl, norm, x, resieved, 0);
    if (trial_div_very_verbose)
      gmp_fprintf (stderr, "# x = %d, after dividing out resieved norm = %Zd\n",
                   x, norm);

    {
      /* Trial divide primes with precomputed tables */
      int nr_factors, i;
      unsigned long factors[32];
      if (trial_div_very_verbose)
        {
          fprintf (stderr, "# Trial division by ");
          for (i = 0; trialdiv_data[i].p != 1; i++)
            fprintf (stderr, " %lu", trialdiv_data[i].p);
          fprintf (stderr, "\n");
        }

      nr_factors = trialdiv (factors, norm, trialdiv_data);
      ASSERT (nr_factors <= 32);

      for (i = 0; i < nr_factors; i++)
        {
          if (trial_div_very_verbose)
            fprintf (stderr, " %lu", factors[i]);
          fl->fac[fl->n++] = factors[i];
        }
      if (trial_div_very_verbose)
        gmp_fprintf (stderr, "\n# After trialdiv(): norm = %Zd\n", norm);
    }
}

/* Return 0 if the leftover norm n cannot yield a relation:
   (a) if n > 2^mfb
   (b) if L < n < B^2
   (c) if L^2 < n < B^3

   FIXME: need to check L^k < n < B^(k+1) too.
*/
static int
check_leftover_norm (mpz_t n, size_t lpb, mpz_t BB, mpz_t BBB, size_t mfb)
{
  size_t s = mpz_sizeinbase (n, 2);

  if (s > mfb)
    return 0;
  if (lpb < s && mpz_cmp (n, BB) < 0)
    return 0;
  if (2 * lpb < s && mpz_cmp (n, BBB) < 0)
    return 0;
  if (lpb < s && mpz_probab_prime_p (n, 1))
    return 0;
  return 1;
}

/* Adds the number of sieve reports to *survivors,
   number of survivors with coprime a, b to *coprimes */

NOPROFILE_STATIC int
factor_survivors (unsigned char *S, int N, bucket_array_t rat_BA,
                  bucket_array_t alg_BA, factorbase_degn_t *fb_rat,
                  factorbase_degn_t *fb_alg, cado_poly cpoly, sieve_info_t *si,
                  unsigned long *survivors, unsigned long *coprimes,
		  small_sieve_data_t *srsd_alg, small_sieve_data_t *srsd_rat,
		  unsigned long *report_sizes_a, unsigned long *report_sizes_r,
		  const unsigned char *rat_S)
{
    int x, xul;
    int64_t a;
    uint64_t b;
    int cpt = 0;
    int surv = 0, copr = 0;
    mpz_t alg_norm, rat_norm, BBalg, BBrat, BBBalg, BBBrat, BLPrat;
    factor_list_t alg_factors, rat_factors;
    mpz_array_t *f_r = NULL, *f_a = NULL;    /* large prime factors */
    uint32_array_t *m_r = NULL, *m_a = NULL; /* corresponding multiplicities */

    f_r = alloc_mpz_array (8);
    f_a = alloc_mpz_array (8);
    m_r = alloc_uint32_array (8);
    m_a = alloc_uint32_array (8);

    factor_list_init(&alg_factors);
    factor_list_init(&rat_factors);
    mpz_init (alg_norm);
    mpz_init (rat_norm);
    mpz_init (BBalg);
    mpz_init (BBrat);
    mpz_init (BBBalg);
    mpz_init (BBBrat);
    mpz_init (BLPrat);
    mpz_ui_pow_ui (BBalg, cpoly->alim, 2);
    mpz_ui_pow_ui (BBrat, cpoly->rlim, 2);
    mpz_mul_ui (BBBalg, BBalg, cpoly->alim);
    mpz_mul_ui (BBBrat, BBrat, cpoly->rlim);
    mpz_set_ui (BLPrat, cpoly->alim);
    mpz_mul_2exp (BLPrat, BLPrat, cpoly->lpba); /* fb bound * lp bound */

    for (x = 0; x < si->bucket_region; ++x)
      {
        unsigned int X;
        unsigned int i, j;

        if (!(si->alg_Bound[S[x]]))
          {
            S[x] = 255;
            continue;
          }
        surv++;

        X = x + (N << si->log_bucket_region);
        i = abs ((int) (X & (si->I - 1)) - si->I / 2);
        j = X >> si->logI;
        if (bin_gcd (i, j) != 1)
          {
            S[x] = 255;
            continue;
          }
      }

    purge_bucket (rat_BA, N, S);
    bucket_sortbucket (rat_BA, N);
    rewind_bucket (rat_BA, N);
    purge_bucket (alg_BA, N, S);
    bucket_sortbucket (alg_BA, N);
    rewind_bucket (alg_BA, N);

    /* Resieve small primes for this bucket region */
    bucket_array_t resieved_alg, resieved_rat;

    /* FIXME: choose a sensible size here */
    resieved_alg = init_bucket_array (1, si->bucket_region);
    resieved_rat = init_bucket_array (1, si->bucket_region);
    resieve_small_bucket_region (resieved_alg, S, srsd_alg, si,
				 si->trialdiv_primes_alg);
    /* fprintf (stderr, "# Resieving bucket %d on alg side found %d primes\n",
       N, nb_of_updates (resieved_alg, 0)); */
    bucket_sortbucket (resieved_alg, 0);
    rewind_bucket (resieved_alg, 0);
    resieve_small_bucket_region (resieved_rat, S, srsd_rat, si,
				 si->trialdiv_primes_rat);
    /* fprintf (stderr, "# Resieving bucket %d on rat side found %d primes\n",
       N, nb_of_updates (resieved_rat, 0)); */
    bucket_sortbucket (resieved_rat, 0);
    rewind_bucket (resieved_rat, 0);

    /* Scan array one long word at a time. If any byte is <255, i.e. if
       the long word is != 0xFFFF...FF, examine the bytes */
    for (xul = 0; xul < si->bucket_region; xul += sizeof (unsigned long))
      if (*(unsigned long *)(S + xul) != (unsigned long)(-1L))
        for (x = xul; x < xul + (int) sizeof (unsigned long); ++x)
          {
            unsigned int i, j;

            if (S[x] == 255)
              continue;

            // Compute algebraic and rational norms.
            xToAB (&a, &b, x + N*si->bucket_region, si);

            /* since a,b both even were not sieved, either a or b should be odd */
            // ASSERT((a | b) & 1);
            if (UNLIKELY(((a | b) & 1) == 0))
              {
                fprintf (stderr, "# Error: a and b both even for N = %d, x = %d,\n"
                                 "i = %d, j = %d, a = %ld, b = %lu\n",
                         N, x, ((x + N*si->bucket_region) & (si->I - 1))
                           - (si->I >> 1),
                         (x + N*si->bucket_region) >> si->logI,
                         (long) a, (unsigned long) b);
                continue;
              }

            /* Since the q-lattice is exactly those (a, b) with
               a == rho*b (mod q), q|b  ==>  q|a  ==>  q | gcd(a,b) */
            if (b == 0 || (b >= si->q && b % si->q == 0))
              continue;

            copr++;

            /* For hunting missed relations */
            // if (a == -1128230057 && b == 66615)
            //  fprintf (stderr, "# Have relation -1128230057,66615 at x = %d\n", x);

            // Trial divide rational norm
            eval_fij (rat_norm, (const mpz_t *) cpoly->g, 1, a, b);
            trial_div (&rat_factors, rat_norm, rat_BA, N, x, fb_rat,
                       resieved_rat, si->trialdiv_data_rat);

            if (!check_leftover_norm (rat_norm, cpoly->lpbr, BBrat, BBBrat,
                                      cpoly->mfbr))
              continue;

            // Trial divide algebraic norm
            eval_fij (alg_norm, (const mpz_t *) cpoly->f, cpoly->degree, a, b);
            mpz_divexact_ui (alg_norm, alg_norm, si->q);
            trial_div (&alg_factors, alg_norm, alg_BA, N, x, fb_alg,
                       resieved_alg, si->trialdiv_data_alg);

            if (!check_leftover_norm (alg_norm, cpoly->lpba, BBalg, BBBalg,
                                      cpoly->mfba))
              continue;

            if (mpz_cmp (rat_norm, BLPrat) > 0)
              {
                /* rat_norm might not be smooth, factor it first */
                if (factor_leftover_norm (rat_norm, cpoly->lpbr, f_r, m_r,
                                          si->strategy) == 0)
                  continue;

                if (factor_leftover_norm (alg_norm, cpoly->lpba, f_a, m_a,
                                          si->strategy) == 0)
                  continue;
              }
            else
              {
                /* rat_norm is definitely smooth, factor it last */
                if (factor_leftover_norm (alg_norm, cpoly->lpba, f_a, m_a,
                                          si->strategy) == 0)
                  continue;

                if (factor_leftover_norm (rat_norm, cpoly->lpbr, f_r, m_r,
                                          si->strategy) == 0)
                  continue;
              }

            printf ("%" PRId64 ",%" PRIu64 ":", a, b);
            factor_list_fprint (stdout, rat_factors);
            for (i = 0; i < f_r->length; ++i)
              for (j = 0; j < m_r->data[i]; j++)
                gmp_printf (",%Zx", f_r->data[i]);
            printf (":");
            if (alg_factors.n != 0)
              {
                factor_list_fprint (stdout, alg_factors);
                printf (",");
              }
            for (i = 0; i < f_a->length; ++i)
              for (j = 0; j < m_a->data[i]; j++)
                gmp_printf ("%Zx,", f_a->data[i]);
            /* print special q */
            printf ("%" PRIx64 "", si->q);
            printf ("\n");
            fflush (stdout);
            cpt++;
            report_sizes_a[S[x]]++; /* Build histogram of lucky S[x] values */
            if (rat_S != NULL)
              report_sizes_r[rat_S[x]]++;
          }

    survivors[0] += surv;
    coprimes[0] += copr;
    clear_bucket_array (resieved_alg);
    clear_bucket_array (resieved_rat);
    mpz_clear (BBalg);
    mpz_clear (BBrat);
    mpz_clear (BBBalg);
    mpz_clear (BBBrat);
    mpz_clear (BLPrat);
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

#if 0
static void
xToIJ(int *i, unsigned int *j, int x, sieve_info_t * si)
{
    *i = (x % (si->I)) - (si->I >> 1);
    *j = x / si->I;
}

static void
IJTox(int *x, int i, int j, sieve_info_t * si)
{
    *x = i + (si->I)*j + (si->I>>1);
}

static void
IJToAB(int64_t *a, uint64_t *b, int i, int j, sieve_info_t * si)
{
    *a = i*si->a0 + j*si->a1;
    *b = i*si->b0 + j*si->b1;
}
#endif

/* Warning: b might be negative, in which case we return (-a,-b) */
static void
xToAB(int64_t *a, uint64_t *b, const unsigned int x, const sieve_info_t *si)
{
    int i, j;
    int64_t c;
    uint32_t I = si->I;

    i = (x & (I - 1)) - (I >> 1);
    j = x >> si->logI;
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

void test_divisible_x (const fbprime_t p, const unsigned long x, const int n,
		       const sieve_info_t *si)
{
  const unsigned long X = x + n * si->bucket_region;
  const long i = (long) (X & ((unsigned long) si->I - 1UL)) - (long) si->I / 2;
  const unsigned long j = X >> si->logI;
  mpz_t v;

  mpz_init (v);
  eval_fij (v, (const mpz_t *) si->fij, si->degree, i, j);
  /* gmp_fprintf (stderr, "# test_divisible_x (" FBPRIME_FORMAT ", %lu, %d, ): "
     "i = %ld, j = %lu, v = %Zd\n", p, x, n, i, j, v); */
  ASSERT (mpz_divisible_ui_p (v, (unsigned long) p));
  mpz_clear (v);
}


/*********************** norm computation ************************************/

/* Put in fij[] the coefficients of f'(i) = F(a0*i+a1, b0*i+b1).
   Assumes the coefficients of fij[] are initialized.
   Put in fijd[] a double-precision approximation of fij[]/q.
*/
static void
fij_from_f (sieve_info_t *si, mpz_t *f, int d)
{
  int k, l;
  mpz_t *g; /* will contain the coefficients of (b0*i+b1)^l */
  mpz_t f0;
  mpz_t *fij = si->fij;
  double invq;

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
      mpz_mul_si (fij[k], fij[k + 1], si->a1);
      for (l = k + 1; l < d; l++)
        {
          mpz_mul_si (fij[l], fij[l], si->a0);
          mpz_addmul_si (fij[l], fij[l + 1], si->a1);
        }
      mpz_mul_si (fij[d], fij[d], si->a0);

      /* now compute (b0*i+b1)^(d-k) from the previous (b0*i+b1)^(d-k-1):
         g[d-k] = b0*g[d-k-1]
         ...
         g[l] = b1*g[l]+b0*g[l-1] for 0 < l < d-k
         ...
         g[0] = b1*g[0]
      */
      mpz_mul_si (g[d - k], g[d - k - 1], si->b0);
      for (l = d - k - 1; l > 0; l--)
        {
          mpz_mul_si (g[l], g[l], si->b1);
          mpz_addmul_si (g[l], g[l-1], si->b0);
        }
      mpz_mul_si (g[0], g[0], si->b1);

      /* now g has degree d-k, and we add f0*g */
      for (l = k; l <= d; l++)
        mpz_addmul (fij[l], g[l - k], f0);
    }

  invq = 1.0 / (double) si->q;
  for (k = 0; k <= d; k++)
    si->fijd[k] = mpz_get_d (fij[k]) * invq;

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

/* returns the maximal value of log2 |F(a,b)/q| for
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
  return log2 (max_norm * pow (si->B * (double) si->I, (double) d)
               / (double) q0);
}

/* v <- |f(i,j)|, where f is of degree d */
NOPROFILE_STATIC void
eval_fij (mpz_t v, const mpz_t *f, const unsigned int d, const long i,
	  const unsigned long j)
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
  mpz_abs (v, v); /* avoids problems with negative norms */
  mpz_clear (jpow);
}


/* check that the double x fits into an int32_t */
#define fits_int32_t(x) \
  ((double) INT32_MIN <= (x)) && ((x) <= (double) INT32_MAX)

void
SkewGauss (sieve_info_t *si, double skewness)
{
  double a[2], b[2], q, maxab0, maxab1;

  a[0] = (double) si->q;
  ASSERT_ALWAYS(a[0] < 9007199254740992.0); /* si.q should be less than 2^53
                                               so that a[0] is exact */
  b[0] = 0.0;
  a[1] = (double) si->rho;
  b[1] = skewness;
  ASSERT(b[1] != 0);
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
  ASSERT(fits_int32_t(a[0]));
  ASSERT(fits_int32_t(b[0] / skewness));
  ASSERT(fits_int32_t(a[1]));
  ASSERT(fits_int32_t(b[1] / skewness));
  /* now b[0], b[1] should be of the form i*skewness, but this might not be
     exact due to rounding errors, thus we round them to the nearest integer */
  maxab0 = fabs (a[0]) > fabs (b[0]) ? fabs (a[0]) : fabs (b[0]);
  maxab1 = fabs (a[1]) > fabs (b[1]) ? fabs (a[1]) : fabs (b[1]);
  if (maxab0 <= maxab1)
    {
      si->a0 = (int32_t) a[0];
      si->b0 = (int32_t) rint (b[0] / skewness);
      si->a1 = (int32_t) a[1];
      si->b1 = (int32_t) rint (b[1] / skewness);
    }
  else /* swap (a0,b0) and (a1,b1) */
    {
      si->a1 = (int32_t) a[0];
      si->b1 = (int32_t) rint (b[0] / skewness);
      si->a0 = (int32_t) a[1];
      si->b0 = (int32_t) rint (b[1] / skewness);
      maxab1 = maxab0;
    }

  /* make sure J does not exceed I/2 */
  if (maxab1 >= si->B)
    si->J = (uint32_t) (si->B * skewness / maxab1 * (double) (si->I >> 1));
  else
    si->J = si->I >> 1;
}

/* return max(|a0|,|a1|)/min(|a0|,|a1|) */
double
skewness (sieve_info_t *si)
{
  double a0 = fabs ((double) si->a0);
  double a1 = fabs ((double) si->a1);

  return (a0 > a1) ? a0 / a1 : a1 / a0;
}

/************************ cofactorization ********************************/

/* FIXME: the value of 20 seems large. Normally, a few Miller-Rabin passes
   should be enough. See also http://www.trnicely.net/misc/mpzspsp.html */
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
          l - (large) prime bit size bound is L=2^l
   Assumes n > 0.
   Return value:
          0 if n has a prime factor larger than 2^l
          1 if all prime factors of n are < 2^l
   Output:
          the prime factors of n are factors->data[0..factors->length-1],
          with corresponding multiplicities multis[0..factors->length-1].
*/
int
factor_leftover_norm (mpz_t n, unsigned int l,
                      mpz_array_t* const factors, uint32_array_t* const multis,
		      facul_strategy_t *strategy)
{
  uint32_t i, nr_factors;
  unsigned long ul_factors[16];
  int facul_code;

  factors->length = 0;
  multis->length = 0;

  /* factoring programs do not like 1 */
  if (mpz_cmp_ui (n, 1) == 0)
    return 1;

  /* If n < L, we know that n is prime, since all primes < B have been
     removed, and L < B^2 in general, where B is the factor base bound,
     thus we only need a primality test when n > L. */
  if (BITSIZE(n) <= l)
    {
      append_mpz_to_array (factors, n);
      append_uint32_to_array (multis, 1);
      return 1;
    }
/* Input is required to be composite!
  else if (IS_PROBAB_PRIME(n))
    {
      return 0;
    }
*/

  /* use the facul library */
  facul_code = facul (ul_factors, n, strategy);

  if (facul_code == FACUL_NOT_SMOOTH)
    return 0;

  ASSERT (facul_code == 0 || mpz_cmp_ui (n, ul_factors[0]) != 0);

  if (facul_code > 0)
    {
      nr_factors = facul_code;
      for (i = 0; i < nr_factors; i++)
	{
	  unsigned long r;
	  mpz_t t;
	  if (ul_factors[i] > (1UL << l)) /* Larger than large prime bound? */
	    return 0;
	  r = mpz_tdiv_q_ui (n, n, ul_factors[i]);
	  ASSERT_ALWAYS (r == 0UL);
	  mpz_init (t);
	  mpz_set_ui (t, ul_factors[i]);
	  append_mpz_to_array (factors, t);
	  mpz_clear (t);
	  append_uint32_to_array (multis, 1); /* FIXME, deal with repeated
						 factors correctly */
	}

      if (mpz_cmp_ui (n, 1UL) == 0)
	return 1;
      unsigned int s = BITSIZE(n);
      if (s <= l)
        {
          append_mpz_to_array (factors, n);
          append_uint32_to_array (multis, 1);
          return 1;
        }
      /* If we still have more than two primes (or something non-smooth),
         bail out */
      if (s > 2*l)
        return 0;
      /* We always abort below, so let's skip the prp test
      if (IS_PROBAB_PRIME(n))
        return 0; */
    }
  /* When sieving for 3 large primes, here are so many left over, non-smooth
     numbers here that factoring them all takes a long time, for few
     additional relations */
  return 0;
}


/*************************** main program ************************************/

static void
usage (char *argv0)
{
  fprintf (stderr, "Usage: %s [-no-checknorms] [-I I] -poly xxx.poly -fb xxx.roots -q0 q0 [-q1 q1] [-rho rho]\n",
           argv0);
  fprintf (stderr, "          -no-checknorms  don't factor leftover norms\n");
  fprintf (stderr, "          -I i            sieving region has side 2^i [default %u]\n", DEFAULT_I);
  fprintf (stderr, "          -poly xxx.poly  use polynomial xxx.poly\n");
  fprintf (stderr, "          -fb xxx.roots   use factor base xxx.roots\n");
  fprintf (stderr, "          -q0 nnn         left bound of special-q range\n");
  fprintf (stderr, "          -q1 nnn         right bound of special-q range\n");
  fprintf (stderr, "          -rho r          sieve only algebraic root r mod q0\n");
  fprintf (stderr, "          -tdthresh nnn   trial-divide primes p/r <= nnn (r=number of roots)\n");
  fprintf (stderr, "          -bkthresh nnn   bucket-sieve primes p >= nnn\n");
  fprintf (stderr, "          -rlim     nnn   rational factor base bound nnn\n");
  fprintf (stderr, "          -alim     nnn   algebraic factor base bound nnn\n");
  fprintf (stderr, "          -lpbr     nnn   rational large prime bound 2^nnn\n");
  fprintf (stderr, "          -lpba     nnn   algebraic large prime bound 2^nnn\n");
  fprintf (stderr, "          -mfbr     nnn   rational cofactor bound 2^nnn\n");
  fprintf (stderr, "          -mfba     nnn   algebraic cofactor bound 2^nnn\n");
  fprintf (stderr, "          -rlambda  nnn   rational lambda value is nnn\n");
  fprintf (stderr, "          -alambda  nnn   algebraic lambda value is nnn\n");
  fprintf (stderr, "          -v              be verbose (print some sieving statistics)\n");
  exit (EXIT_FAILURE);
}

int
main (int argc, char *argv[])
{
    char *argv0 = argv[0];
    sieve_info_t si;
    char *fbfilename = NULL, *polyfilename = NULL;
    cado_poly cpoly;
    double t0, tfb, tq, tn_rat, tn_alg, tts, ttsm, ttf;
    uint64_t q0 = 0, q1 = 0, rho = 0;
    uint64_t *roots;
    unsigned long nroots, tot_reports = 0, survivors0, survivors1, survivors2;
    factorbase_degn_t * fb_alg, * fb_rat;
    int checknorms = 1; /* factor or not the remaining norms */
    int td_thresh = 1024; /* cost threshold trialdiv/resieving */
    int bucket_thresh = 0;
    int I = DEFAULT_I, i;
    int verbose = 0;
    unsigned long sq = 0;
    double totJ = 0.0;
    unsigned long report_sizes_a[256], report_sizes_r[256];
    /* following command-line values override those in the polynomial file */
    int rlim = 0, alim = 0;
    int lpbr = 0, lpba = 0;
    int mfbr = 0, mfba = 0;
    double rlambda = 0.0, alambda = 0.0;

    fprintf (stderr, "# %s.r%s", argv[0], REV);
    for (i = 1; i < argc; i++)
      fprintf (stderr, " %s", argv[i]);
    fprintf (stderr, "\n");

    while (argc > 1 && argv[1][0] == '-')
      {
        if (strcmp (argv[1], "-no-checknorms") == 0)
          {
            checknorms = 0;
            argc -= 1;
            argv += 1;
          }
        if (strcmp (argv[1], "-v") == 0)
          {
            verbose = 1;
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
        else if (argc > 2 && strcmp (argv[1], "-tdthresh") == 0)
          {
            td_thresh = atoi (argv[2]);
            argc -= 2;
            argv += 2;
          }
        else if (argc > 2 && strcmp (argv[1], "-bkthresh") == 0)
          {
            bucket_thresh = atoi (argv[2]);
            argc -= 2;
            argv += 2;
          }
        else if (argc > 2 && strcmp (argv[1], "-rlim") == 0)
          {
            rlim = atoi (argv[2]);
            argc -= 2;
            argv += 2;
          }
        else if (argc > 2 && strcmp (argv[1], "-alim") == 0)
          {
            alim = atoi (argv[2]);
            argc -= 2;
            argv += 2;
          }
        else if (argc > 2 && strcmp (argv[1], "-lpbr") == 0)
          {
            lpbr = atoi (argv[2]);
            argc -= 2;
            argv += 2;
          }
        else if (argc > 2 && strcmp (argv[1], "-lpba") == 0)
          {
            lpba = atoi (argv[2]);
            argc -= 2;
            argv += 2;
          }
        else if (argc > 2 && strcmp (argv[1], "-mfbr") == 0)
          {
            mfbr = atoi (argv[2]);
            argc -= 2;
            argv += 2;
          }
        else if (argc > 2 && strcmp (argv[1], "-mfba") == 0)
          {
            mfba = atoi (argv[2]);
            argc -= 2;
            argv += 2;
          }
        else if (argc > 2 && strcmp (argv[1], "-rlambda") == 0)
          {
            rlambda = atof (argv[2]);
            argc -= 2;
            argv += 2;
          }
        else if (argc > 2 && strcmp (argv[1], "-alambda") == 0)
          {
            alambda = atof (argv[2]);
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
        exit (EXIT_FAILURE);
      }

    /* if -q1 is not given, sieve only for q0 */
    if (q1 == 0)
      q1 = q0 + 1;

    /* check that q1 fits into an unsigned long */
    if (q1 > (uint64_t) ULONG_MAX)
      {
        fprintf (stderr, "Error, q1=%" PRIu64 " exceeds ULONG_MAX\n", q1);
        exit (EXIT_FAILURE);
      }

    cado_poly_init(cpoly);
    if (!cado_poly_read (cpoly, polyfilename))
      {
        fprintf (stderr, "Error reading polynomial file\n");
        exit (EXIT_FAILURE);
      }

    /* override sieving parameters if needed */
    if (rlim != 0)
      cpoly->rlim = rlim;
    if (alim != 0)
      cpoly->alim = alim;
    if (lpbr != 0)
      cpoly->lpbr = lpbr;
    if (lpba != 0)
      cpoly->lpba = lpba;
    if (mfbr != 0)
      cpoly->mfbr = mfbr;
    if (mfba != 0)
      cpoly->mfba = mfba;
    if (rlambda != 0.0)
      cpoly->rlambda = rlambda;
    if (alambda != 0.0)
      cpoly->alambda = alambda;
    fprintf (stderr, "# Sieving parameters: rlim=%lu alim=%lu lpbr=%d lpba=%d\n",
             cpoly->rlim, cpoly->alim, cpoly->lpbr, cpoly->lpba);
    fprintf (stderr, "#                     mfbr=%d mfba=%d rlambda=%1.1f alambda=%1.1f\n",
             cpoly->mfbr, cpoly->mfba, cpoly->rlambda, cpoly->alambda);

    /* this does not depend on the special-q */
    sieve_info_init (&si, cpoly, I, q0);
    si.checknorms = checknorms;
    if (bucket_thresh > si.bucket_thresh)
      si.bucket_thresh = bucket_thresh;

    {
      fbprime_t *leading_div;
      tfb = seconds ();
      leading_div = factor_small (cpoly->f[cpoly->degree], cpoly->alim);
      fb_alg = fb_read_addproj (fbfilename, si.scale_alg * LOG_SCALE, 0,
				leading_div);
      ASSERT_ALWAYS(fb_alg != NULL);
      tfb = seconds () - tfb;
      fprintf (stderr, "# Reading algebraic factor base took %1.1fs\n", tfb);
      free (leading_div);
      // fb_fprint (stderr, fb_alg);
    }

    /* Prepare rational factor base */
    tfb = seconds ();
    fb_rat = fb_make_linear (cpoly->g, (fbprime_t) cpoly->rlim,
                             si.scale_rat * LOG_SCALE, verbose, 1);
    tfb = seconds () - tfb;
    fprintf (stderr, "# Creating rational factor base took %1.1fs\n", tfb);

    init_rat_norms (&si);
    init_alg_norms (&si);

    /* Init refactoring stuff */
    int skip2;
    si.trialdiv_primes_rat = fb_extract_bycost (fb_rat, si.bucket_thresh,
                                                td_thresh);
    si.trialdiv_primes_alg = fb_extract_bycost (fb_alg, si.bucket_thresh,
                                                td_thresh);
    for (i = 0; si.trialdiv_primes_alg[i] != FB_END; i++);
    skip2 = (i > 0 && si.trialdiv_primes_alg[0] == 2) ? 1 : 0;
    si.trialdiv_data_alg = trialdiv_init (si.trialdiv_primes_alg + skip2,
                                          i - skip2);
    for (i = 0; si.trialdiv_primes_rat[i] != FB_END; i++);
    skip2 = (i > 0 && si.trialdiv_primes_rat[0] == 2) ? 1 : 0;
    si.trialdiv_data_rat = trialdiv_init (si.trialdiv_primes_rat + skip2,
                                          i - skip2);
    si.strategy = facul_make_strategy (15, MIN(cpoly->rlim, cpoly->alim),
                                       1UL << MIN(cpoly->lpbr, cpoly->lpba));

    for (i = 0; i < 256; i++)
      {
        report_sizes_a[i] = 0;
        report_sizes_r[i] = 0;
      }

    /* special q (and root rho) */
    roots = (uint64_t *) malloc (cpoly->degree * sizeof (uint64_t));
    q0 --; /* so that nextprime gives q0 if q0 is prime */
    nroots = 0;
    tot_reports = 0;
    tn_rat = tn_alg = tts = ttsm = ttf = 0.0;
    t0 = seconds ();
    fprintf (stderr, "#\n");
    while (q0 < q1)
      {
        while (nroots == 0) /* go to next prime and generate roots */
          {
            q0 = uint64_nextprime (q0);
            if (q0 >= q1)
              goto end;
            si.q = q0;
            nroots = poly_roots_uint64 (roots, cpoly->f, cpoly->degree, q0);
            if (nroots > 0)
              {
                fprintf (stderr, "### q=%" PRIu64 ": root%s", q0,
                         (nroots == 1) ? "" : "s");
                for (i = 1; i <= (int) nroots; i++)
                  fprintf (stderr, " %" PRIu64, roots[nroots-i]);
                fprintf (stderr, "\n");
              }
          }
        tq = seconds ();

        /* computes a0, b0, a1, b1 from q, rho, and the skewness */
        si.rho = roots[--nroots];
        if (rho != 0 && si.rho != rho) /* if -rho, wait for wanted root */
          continue;
        SkewGauss (&si, cpoly->skew);
        /* FIXME: maybe we can discard some special q's if a1/a0 is too large,
           see http://www.mersenneforum.org/showthread.php?p=130478 */
        /* if (skewness (&si) > 32.0)
           continue; */

        fprintf (stderr, "# Sieving q=%" PRIu64 "; rho=%" PRIu64
                 "; a0=%d; b0=%d; a1=%d; b1=%d\n",
                 si.q, si.rho, si.a0, si.b0, si.a1, si.b1);
        sq ++;
        /* precompute the skewed polynomial */
        fij_from_f (&si, cpoly->f, cpoly->degree);

        /* checks the value of J */
        sieve_info_update (&si, verbose);
        totJ += (double) si.J;

        /* Allocate alg buckets */
        ttsm -= seconds();
        bucket_array_t alg_BA;
        alg_BA = init_bucket_array(si.nb_buckets, si.bucket_limit);

        /* Fill in alg buckets */
        fill_in_buckets(alg_BA, fb_alg, &si);

        /* Allocate rational buckets */
        bucket_array_t rat_BA;
        rat_BA = init_bucket_array(si.nb_buckets, si.bucket_limit);

        /* Fill in rational buckets */
        fill_in_buckets(rat_BA, fb_rat, &si);
        ttsm += seconds();

        /* Initialize data for sieving small primes */
        small_sieve_data_t ssd_alg, ssd_rat;
        init_small_sieve(&ssd_rat, fb_rat, &si, 'r');
        init_small_sieve(&ssd_alg, fb_alg, &si, 'a');

	/* Copy small sieve data for resieving small primes */
	small_sieve_data_t srsd_alg, srsd_rat;
	copy_small_sieve (&srsd_alg, &ssd_alg);
	copy_small_sieve (&srsd_rat, &ssd_rat);

        /* Process each bucket region, one after the other */
        unsigned char *S, *rat_S;
        S = (unsigned char *)malloc(si.bucket_region*sizeof(unsigned char));
        if (verbose)
          rat_S = (unsigned char *)
                  malloc(si.bucket_region*sizeof(unsigned char));
        else
          rat_S = NULL;
        int reports = 0;
        survivors0 = survivors1 = survivors2 = 0;
        for (i = 0; i < si.nb_buckets; ++i) {
            /* Init rational norms */
            tn_rat -= seconds ();
            init_rat_norms_bucket_region(S, i, cpoly, &si);
            tn_rat += seconds ();
            /* Apply rational bucket */
            ttsm -= seconds();
            apply_one_bucket(S, rat_BA, i, &si);
            ttsm += seconds();
            /* Sieve small rational primes */
            sieve_small_bucket_region(S, i, ssd_rat, &si, 'r');
            /* Make copy of sieve report values for histogram */
            if (rat_S != NULL)
              memcpy (rat_S, S, si.bucket_region*sizeof(unsigned char));

            /* Init algebraic norms */
            tn_alg -= seconds ();
            survivors0 += init_alg_norms_bucket_region(S, i, cpoly, &si);
            tn_alg += seconds ();
            /* Apply algebraic bucket */
            ttsm -= seconds();
            apply_one_bucket(S, alg_BA, i, &si);
            ttsm += seconds();
            /* Sieve small algebraic primes */
            sieve_small_bucket_region(S, i, ssd_alg, &si, 'a');

            /* Factor survivors */
            ttf -= seconds ();
            reports += factor_survivors (S, i, rat_BA, alg_BA, fb_rat,
                                         fb_alg, cpoly, &si, &survivors1,
					 &survivors2, &srsd_alg, &srsd_rat,
					 report_sizes_a, report_sizes_r,
					 rat_S);
            ttf += seconds ();
        }
        clear_small_sieve(ssd_rat);
        clear_small_sieve(ssd_alg);
        clear_small_sieve(srsd_rat);
        clear_small_sieve(srsd_alg);
        tot_reports += reports;
        if (verbose)
          {
            fprintf (stderr, "# %lu survivors after rational sieve,",
                     survivors0);
            fprintf (stderr, " %lu survivors after algebraic sieve, ",
                     survivors1);
            fprintf (stderr, "coprime: %lu\n", survivors2);
          }
        fprintf (stderr, "# %d relation(s) for (%" PRIu64 ",%" PRIu64
                 "), total %lu [%1.3fs/r]\n#\n", reports, si.q, si.rho,
                 tot_reports, (seconds () - t0) / (double) tot_reports);
        clear_bucket_array(alg_BA);
        clear_bucket_array(rat_BA);
        free(S);
        if (rat_S != NULL)
          free (rat_S);
      } // end of loop over special q ideals.

 end:
    t0 = seconds () - t0;
    fprintf (stderr, "# Average J=%1.0f for %lu special-q's\n",
             totJ / (double) sq, sq);
    tts = t0 - (tn_rat + tn_alg + ttf);
    if (verbose)
      facul_print_stats (stderr);
    if (verbose)
      {
        fprintf (stderr, "# Histogram of sieve report values that led to "
                 "relations:\n# Algebraic side: ");
        for (i = 0; i < 256; i++)
          if (report_sizes_a[i] > 0)
            fprintf (stderr, "%d:%lu ", i, report_sizes_a[i]);
        fprintf (stderr, "\n# Rational side: ");
        for (i = 0; i < 256; i++)
          if (report_sizes_r[i] > 0)
            fprintf (stderr, "%d:%lu ", i, report_sizes_r[i]);
        fprintf (stderr, "\n");
      }
    fprintf (stderr, "# Total time %1.1fs [norm %1.2f+%1.1f, sieving %1.1f"
            " (%1.1f + %1.1f),"
             " factor %1.1f]\n", t0, tn_rat, tn_alg, tts, ttsm, tts-ttsm, ttf);
    fprintf (stderr, "# Total %lu reports [%1.3fs/r, %1.1fr/sq]\n",
             tot_reports, t0 / (double) tot_reports,
             (double) tot_reports / (double) sq);

    facul_clear_strategy (si.strategy);
    si.strategy = NULL;
    trialdiv_clear (si.trialdiv_data_alg);
    trialdiv_clear (si.trialdiv_data_rat);
    free (si.trialdiv_primes_alg);
    free (si.trialdiv_primes_rat);
    free (fb_alg);
    free (fb_rat);
    sieve_info_clear (&si);
    cado_poly_clear (cpoly);
    free (roots);

    return 0;
}
