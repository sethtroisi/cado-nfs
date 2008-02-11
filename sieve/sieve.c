#include "config.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <limits.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <time.h>
#include "cado.h"
#include "sieve.h"
#include "../utils/mod_ul.h"
#include "fb.h"
#include "sieve_aux.h"
#include "basicnt.h"
#include "../utils/utils.h"

#define MIN(a,b) ((a)<(b)?(a):(b))
#define MAX(a,b) ((a)>(b)?(a):(b))
#define LOG2 0.69314718055994530941723212145817656808
#define INVLOG2 1.4426950408889634073599246810018921374

#define SIEVE_PERMISSIBLE_ERROR ((unsigned char) 7)
#define PRP_REPS 1 /* number of tests in mpz_probab_prime_p */

#ifndef TRIALDIV_SKIPFORWARD
#define TRIALDIV_SKIPFORWARD 1
#define TRIALDIV_SKIPFORWARD_DEFAULT
#endif

/* We'll use ul_rho() instead of trial division if the size in bits 
   of the missing factor base prime is > UL_RHO_LOGTHRES */
#ifndef UL_RHO_LOGTHRES
#define UL_RHO_LOGTHRES 16
#define UL_RHO_LOGTHRES_DEFAULT
#endif

/* We'll use mpz_rho() instead of trial division if the size in bits 
   of the missing factor base prime is > UL_RHO_LOGTHRES */
#ifndef MPZ_RHO_LOGTHRES
#define MPZ_RHO_LOGTHRES 20
#define MPZ_RHO_LOGTHRES_DEFAULT
#endif

#define uc_add(a,b) ((unsigned char) ((a)+(b)))
#define uc_sub(a,b) ((unsigned char) ((a)-(b)))
#define add_error(a) ((unsigned char) ((a) + SIEVE_PERMISSIBLE_ERROR))


unsigned long sumprimes, nrprimes; /* For largest non-report fb primes */
unsigned long sumprimes2, nrprimes2;  /* For 2nd largest non-report fb prim. */
unsigned long ul_rho_called = 0, ul_rho_called1 = 0, 
  mpz_rho_called = 0, mpz_rho_called1 = 0;
unsigned long ul_rho_toolarge_nr, ul_rho_toolarge_sum; 
long long ul_rho_time;

#ifdef TRACE_RELATION_A
static void 
TRACE_A (long a, const char *fn, const int line, const char *s, ...)
{
  va_list ap;

  va_start (ap, s);

  if (a == TRACE_RELATION_A)
    {
      printf ("# TRACE RELATION a = %ld in %s(%d): ", a, fn, line);
      gmp_vprintf (s, ap);
      fflush (stdout);
    }

  va_end (ap);
}
#else
static void
TRACE_A (__attribute__ ((unused)) long a, 
	 __attribute__ ((unused)) const char *fn, 
	 __attribute__ ((unused)) const int line,
	 __attribute__ ((unused)) const char *s, ...)
{
}
#endif

/* Given an index i to the sieve array, an amin value and whether 
   only odd $a$ are in the sieve (odd=1) or not (odd=0),
   return the value of $a$ at sieve location i */
static long 
sieveidx_to_a (unsigned long i, long amin, int odd) 
{
  ASSERT_EXPENSIVE (odd == 0 || odd == 1);
  ASSERT_EXPENSIVE (odd == 0 || (amin & 1) == 1);
  return amin + ((long)(i) << (odd));
}

static unsigned long 
a_to_sieveidx (long a, long amin, int odd)
{
  ASSERT_EXPENSIVE (a >= amin);
  ASSERT_EXPENSIVE (odd == 0 || odd == 1);
  ASSERT_EXPENSIVE (odd == 0 || (a & 1) == 1);
  ASSERT_EXPENSIVE (odd == 0 || (amin & 1) == 1);
  return (a - amin) >> odd;
}

/* Some multiple precision functions we'll need */

void
mp_poly_eval (mpz_t r, mpz_t *poly, int deg, long a)
{
  int i;

  mpz_set (r, poly[deg]);
  for (i = deg - 1; i >= 0; i--)
    {
      mpz_mul_si (r, r, a);
      mpz_add (r, r, poly[i]);
    }
}


/* Scale coefficient f_i by a^i (inv == 1) or by c^(deg-i) (inv == -1) */
void
mp_poly_scale (mpz_t *r, mpz_t *poly, const int deg, const long c, 
	       const int inv)
{
  int i;
  mpz_t t;

  ASSERT (inv == 1 || inv == -1);

  mpz_init (t);
  mpz_set_ui (t, 1UL);
  for (i = 0; i <= deg; i++)
    {
      const int j = (inv == 1) ? i : deg - i;
      mpz_mul (r[j], poly[j], t);
      mpz_mul_si (t, t, c);
    }
  
  mpz_clear (t);
}


void 
mp_poly_print (mpz_t *poly, int deg, const char *name, int homogeneous)
{
  int i;
  printf ("%s", name);

  if (!homogeneous)
    {
      for (i = deg; i >= 2; i--)
	{
	  if (mpz_sgn (poly[i]) != 0)
	    gmp_printf ("%s%Zd * x^%d ", 
			mpz_sgn(poly[i]) > 0 ? "+" : "", poly[i], i);
	}
      if (deg >= 1 && mpz_sgn (poly[1]) != 0)
	gmp_printf ("%s%Zd * x ", mpz_sgn(poly[1]) > 0 ? "+" : "", poly[1]);
      if (deg >= 0 && mpz_sgn (poly[0]) != 0)
	gmp_printf ("%s%Zd", mpz_sgn(poly[0]) > 0 ? "+" : "", poly[0]);
    } 
  else 
    {
      for (i = deg; i >= 0; i--)
	{
	  if (mpz_sgn (poly[i]) != 0)
	    {
	      gmp_printf (" %s%Zd", 
			  mpz_sgn(poly[i]) > 0 && i < deg ? "+" : "", 
			  poly[i]);
	      if (i > 1)
		printf (" *a^%d", i);
	      if (i == 1)
		printf ("*a");
	      if (i + 1 < deg)
		printf ("*b^%d", deg - i);
	      if (i + 1 == deg)
		printf ("*b");
	    }
	}
    }
}

/* Evaluate the poynomial f of degree deg at point x */
static double
fpoly_eval (const double *f, const int deg, const double x)
{
  double r;
  int i;

  r = f[deg];
  for (i = deg - 1; i >= 0; i--)
    r = r * x + f[i];
  
  return r;
}

#if 0
/* Print polynomial with floating point coefficients. Assumes f[deg] != 0
   if deg > 0. */
static void 
fpoly_print (const double *f, const int deg, char *name)
{
  int i;

  printf (name);

  if (deg == 0)
    printf ("%f", f[0]);

  if (deg == 1)
    printf ("%f*x", f[1]);

  if (deg > 1)
    printf ("%f*x^%d", f[deg], deg);

  for (i = deg - 1; i >= 0; i--)
    {
      if (f[i] == 0.)
	continue;
      if (i == 0)
	printf (" %s %f", (f[i] > 0) ? "+" : "-", fabs(f[i]));
      else if (i == 1)
	printf (" %s %f*x", (f[i] > 0) ? "+" : "-", fabs(f[i]));
      else 
	printf (" %s %f*x^%d", (f[i] > 0) ? "+" : "-", fabs(f[i]), i);
    }

  printf ("\n");
}
#endif

/* Used only in compute_norms() and trialdiv_one_side()*/
static unsigned char
log_norm (const double *f, const int deg, const double x, 
	  const double log_scale, const double log_proj_roots)
{
  double r;
  unsigned char n;

  r = fpoly_eval (f, deg, x);
  n = fb_log (fabs (r), log_scale, - log_proj_roots);
  /* printf ("Norm at x = %.0f is %.0f, rounded log is %hhu\n", x, r, n); */
  return n;
}


/* Very slow but thorough way of computing norms. Returns the maximum rounded
   log norm in sievearray. */

unsigned char
compute_norms (unsigned char *sievearray, const long amin, const long amax, 
	       const unsigned long b, const double *poly, const int deg, 
	       const double proj_roots, const double log_scale, const int odd,
	       const int verbose)
{
  double f[MAXDEGREE + 1]; /* Polynomial in $a$ for a given fixed $b$ */
  double df[MAXDEGREE]; /* The derivative of f(x) */
  double bpow;
  clock_t tsc1, tsc2;
  const double log_proj_roots = log(proj_roots) * log_scale;
  long a, a2;
  int i, fsign1, fsign2, dfsign1, dfsign2;
  unsigned char n1, n2, nmax;
  const int stride = 128;
  
  ASSERT (odd == 0 || odd == 1);
  ASSERT (!odd || (amin & 1) == 1);
  ASSERT (!odd || (amax & 1) == 1);
  ASSERT (amin <= amax);

  tsc1 = clock ();
  
  bpow = 1.;
  for (i = deg; i >= 0; i--)
    {
      f[i] = poly[i] * bpow;
      bpow *= (double) b;
    }
  /* fpoly_print (f, deg, "# compute_norms: f(x) = "); */

  for (i = 1; i <= deg; i++)
    df[i - 1] = f[i] * (double) i;
  /* fpoly_print (df, deg - 1, "# compute_norms: df(x) = "); */
  
  n1 = log_norm (f, deg, (double) amin, log_scale, log_proj_roots);
  fsign1 = (fpoly_eval (f, deg, (double) amin) < 0) ? -1 : 1;
  dfsign1 = (fpoly_eval (df, deg - 1, (double) amin) < 0) ? -1 : 1;

  nmax = n1;

  a = amin;
  a2 = amin;
  /* a == amin + (d << odd), d = (a - amin) >> odd */

  while (a <= amax)
    {
      double n;
      /* We'll cover [a ... a2[ now */
      ASSERT (a == a2);
      ASSERT_EXPENSIVE (n1 == log_norm (f, deg, (double) a, log_scale, log_proj_roots));
      a2 = MIN(a + stride, amax + (1 << odd));
      ASSERT (!odd || (a2 - a) % 2 == 0);
      
      n = fpoly_eval (f, deg, (double) a2);
      n2 = fb_log (fabs(n), log_scale, -log_proj_roots);
      fsign2 = (n < 0.) ? -1 : 1;
      dfsign2 = (fpoly_eval (df, deg - 1, (double) a2) < 0.) ? -1 : 1;

      if (n1 == n2 && fsign1 == fsign2 && dfsign1 == dfsign2)
	{
	  /* Let's assume the log norm is n1 everywhere in this interval */
#ifdef TRACE_RELATION_A
	  if (a <= TRACE_RELATION_A && TRACE_RELATION_A < a2)
	    {
	      double na = log(fpoly_eval (f, deg, (double) a2)) * log_scale 
		- log_proj_roots;
	      printf ("# TRACE RELATION a = %ld in compute_norms: rounded log "
		      "norm without projective divisors (range computed) for "
		      "degree %d is %d, exact log of norm is %f\n", 
		      (long) TRACE_RELATION_A, deg, (int) n1, na);
	    }
#endif
	  memset (sievearray + ((a - amin) >> odd), n1, 
	          (size_t) ((a2 - a) >> odd));
	}
      else
	{
	  /* n1 and n2 or the sign of the derivatives are different. Do each 
	     a up to (exclusive) a2 individually */
	  /* printf ("log_c(F(%ld, %lu)) == %d != log_c(F(%ld, %lu)) == %d\n",
	     a, b, n1, a2, b, n2); */
	  sievearray[(a - amin) >> odd] = n1;
	  a += 1 << odd;
	  for ( ; a < a2; a += 1 << odd)
	    {
	      unsigned char n = log_norm (f, deg, (double) a, log_scale, 
					  log_proj_roots);

	      TRACE_A (a, __func__, __LINE__, 
		       "log norm (individually computed) for degree %d is "
		       "%d, norm is %f\n", 
		       deg, (int) n, fpoly_eval (f, deg, (double) a));

	      sievearray[(a - amin) >> odd] = n;
	      if (n > nmax)
		nmax = n;
	    }
	  
	}
      a = a2;
      n1 = n2;
      fsign1 = fsign2;
      dfsign1 = dfsign2;
    }
  
  tsc2 = clock ();
  if (verbose)
    {
      printf ("# Computing norms took %ld clocks\n", (long int) (tsc2 - tsc1));
      printf ("# Maximum rounded log norm is %u\n", (unsigned int) nmax);
    }
  
#ifdef PARI
  for (a = amin; a <= amax; a += 1 << odd)
    printf ("if(log_c(%c(%ld, %lu), %f) != %d,"
	    "print (\"log_c(%c(\", %ld,\", \", %lu\"), %f) \", %d)) /* PARI */\n", 
	    deg > 1 ? 'F' : 'G', a, b, log_proj_roots, sievearray[(a - amin) >> odd],
	    deg > 1 ? 'F' : 'G', a, b, log_proj_roots, sievearray[(a - amin) >> odd]); 
#endif

  return nmax;
}


/* Add a sieve report (a, p, l) to the report buffer *reports. 
   If the buffer becomes full, this function reallocates for more space. */

static void
add_sieve_report (sieve_reportbuffer_t *reports, const long a, 
		  const fbprime_t p, const unsigned char l)
{
  ASSERT (reports != NULL);
  ASSERT (reports->alloc > (size_t) 0);
  ASSERT (reports->nr <= reports->alloc);

  if (reports->nr == reports->alloc)
    {
      /* List is full, need to reallocate. Double its size */
      reports->alloc *= 2;
      TRACE_A (a, __func__, __LINE__, "increasing report buffer size to %ld"
	       "entries\n", (long) reports->alloc);
#if 0
      printf ("# Increasing reports buffer size to %ld entries\n", 
	      (long) reports->alloc);
#endif
      reports->reports = (sieve_report_t *) 
	realloc (reports->reports, reports->alloc * sizeof (sieve_report_t));
      ASSERT_ALWAYS (reports->reports != NULL);
    }

  TRACE_A (a, __func__, __LINE__, "adding report, p = %" FBPRIME_FORMAT
	   ", remaining log norm = %hhu\n", p, l);
  reports->reports[reports->nr].a = a;
  reports->reports[reports->nr].p = p;
  reports->reports[reports->nr].l = l;
  reports->nr++;
}

static void 
sieve_small_slow (unsigned char *sievearray, factorbase_small_inited_t *fb,
		  const unsigned int arraylen)
{
  fbprime_t p;
  unsigned int d;
  unsigned char l;

  for (; fb->p > 0; fb++)
    {
      p = fb->p;
      d = fb->loc_and_log & 0xffffff; /* Mask low 24 bits */
      l = fb->loc_and_log >> 24;
      while (d < arraylen)
	{
	  sievearray[d] -= l;
	  d += p;
	}
      d -= arraylen;
      fb->loc_and_log = d + (l << 24);
    }
}


/* sievearray must be long enough to hold amax-amin+1 chars.
   proj_roots is the rounded log of the projective roots */

/* odd: if 1, the sieve array contains only locations for odd a, 
        amin <= a <= amax */

/* Eventually we'll do:
   1. Generate bucket-sorted sieve updates for large fb primes. One bucket
      should hold the updates to positions in the same L2 chunk
   2. For each L2 chunk in this line:
     2.1 For each L1 chunk in this L2 chunk:
       2.1.1 Compute norms, dividing out tiny primes, i.e. p < 10
       2.1.2 Sieve by small primes, p < L1SIZE
     2.2 Do the sieve updates from this L2 chunk's bucket, writing
         sieve reports with useful primes and approx. log
*/


void 
sieve (unsigned char *sievearray, factorbase_degn_t *fb, 
       const long amin, const long amax, const unsigned long b, 
       const unsigned char threshold, sieve_reportbuffer_t *reports,
       const int odd)
{
  /* The sievearray[0] entry corresponds to (amin, b), and
     sievearray[d] to (amin + d * (1 + odd), b) */
  const uint32_t l = (amax - amin) / (1 + odd) + 1;
  uint32_t i, amin_p, p, d;
  unsigned char plog;
  const unsigned char threshold_with_error = add_error (threshold);

  ASSERT (odd == 0 || odd == 1);
  ASSERT (!odd || (amin & 1) == 1);
  ASSERT (!odd || (amax & 1) == 1);
  ASSERT (!odd || (b & 1) == 0);
  ASSERT (amin <= amax);

  /* Do the sieving */

  while (fb->p > 0)
    {
      ASSERT (fb_entrysize (fb) <= fb->size);
      p = fb->p;
      plog = fb->plog;

      /* Compute amin % p for the a value of the sieve location in 
	 sievearray[0] */
      /* FIXME This modular reduction should be simplified somehow. Do it 
	 once and store it in fb? */

      amin_p = signed_mod_longto32 (amin, p);
      
      for (i = 0; i < fb->nr_roots; i++)
        {
	  d = first_sieve_loc (p, fb->roots[i], amin_p, b, odd);
	  
	  /* Now d is the first index into sievearray where p divides. */

          /* Now update the sieve array. There is some blocking, but no
             bucket sieving atm. */
          for (; d < l; d += p)
	    {
	      unsigned char k;
	      k = sievearray[d] - plog;

	      TRACE_A (sieveidx_to_a(d, amin, odd), __func__, __LINE__, 
		       "new approx norm after sieving out " FBPRIME_FORMAT 
		       " is %d\n", p, (int) k);

	      sievearray[d] = k;
	      if (add_error (k) <= threshold_with_error)
		add_sieve_report (reports, sieveidx_to_a(d, amin, odd), p, k);
	    }
        }
      
      /* Move on to the next factor base prime */
      fb = fb_next (fb);
    }
}


static unsigned long
find_sieve_reports (const unsigned char *sievearray, 
		    sieve_reportbuffer_t *reports, 
                    const unsigned char threshold, 
                    const long amin, const unsigned long l, 
		    const unsigned long b, const int odd,
		    sieve_reportbuffer_t *other_reports)
{
  unsigned long d;
  const unsigned char threshold_with_error = add_error (threshold);
  size_t other_idx = 0;
  
  ASSERT (odd == 0 || odd == 1);
  ASSERT (!odd || (amin & 1) == 1);
  ASSERT (b > 0);
  ASSERT (reports->nr <= reports->alloc);
  ASSERT (other_reports == NULL || other_reports->nr <= other_reports->alloc);
  
  for (d = 0; d < l; d++)
    {
      if (add_error(sievearray[d]) <= threshold_with_error)
	{
	  const long a = sieveidx_to_a (d, amin, odd);
	  
	  TRACE_A (a, __func__, __LINE__, "remaining log norm %hhu is small"
		   "enough\n", sievearray[d]);

	  /* Testing n%3==0 is quite fast, and eliminates useless reports */
	  if (a % 3 == 0 && b % 3 == 0)
	    {
	      TRACE_A (a, __func__, __LINE__, "both a and b are divisible "
		       "by 3, not adding to report buffer\n");
	      continue;
	    }

	  /* If other_reports was passed, store this report only if there
	     is a matching report in other_reports */
	  if (other_reports != NULL)
	    {
	      while (other_idx < other_reports->nr && 
		     other_reports->reports[other_idx].a < a)
		{
		  /* Test that reports in other_reports are sorted in order of
		     increasing a value */
		  ASSERT (other_idx == 0 || 
			  other_reports->reports[other_idx - 1].a <=
			  other_reports->reports[other_idx].a);
	      
		  other_idx++;
		}
	      if (other_idx == other_reports->nr || 
		  other_reports->reports[other_idx].a != a)
		{
		  TRACE_A (a, __func__, __LINE__,
			   "no matching report in other_reports\n");
		  continue;
		}
	      else
		{
		  TRACE_A (a, __func__, __LINE__, 
			   "found matching report in other_reports\n");
		}
	    }

	  TRACE_A (a, __func__, __LINE__,
		   "sieve report a = %ld, p = 1, l = %d added to buffer\n", 
		   a, (int) sievearray[d]);

	  add_sieve_report (reports, a, (fbprime_t) 1, sievearray[d]);
	}
    }
}


static inline void
swap_reports (sieve_report_t *r1, sieve_report_t *r2)
{
  sieve_report_t t;

  t.a = r1->a;
  t.p = r1->p;
  t.l = r1->l;
  r1->a = r2->a;
  r1->p = r2->p;
  r1->l = r2->l;
  r2->a = t.a;
  r2->p = t.p;
  r2->l = t.l;
}

/* Sort sieve reports by increasing a value */

static void 
sort_sieve_reports_recurse (sieve_report_t *r, const size_t l)
{
  long p; /* pivot element */
  size_t b, h;
  const size_t m = (l - 1) / 2; /* The midpoint */

  if (l < 2)
    return;

  if (l == 2)
    {
      if (r[0].a > r[1].a)
	swap_reports (r, r + 1);
      return;
    }

  if (r[0].a > r[m].a)
    swap_reports (r, r + m);
  if (r[m].a > r[l - 1].a)
    swap_reports (r + m, r + (l - 1));
  if (r[0].a > r[m].a)
    swap_reports (r, r + m);
  /* Now r[0], r[m] and  r[l-1] are sorted in ascending order */

  if (l == 3)
    return;

  p = r[m].a;
  b = 1; 
  h = l - 1;
  
  /* Here r[0].a <= p, r[m].a == p, r[l-1] >= p */

  while (b < h)
    {
      while (b <= h && r[b].a <= p)
	b++;  /* r[i].a <= p for all 0 <= i < b */
      while (r[h].a > p)
	h--;
      if (b < h)
	swap_reports (r + b, r + h);
    }

  h = b;
  while (b > 0 && r[b - 1].a == p)
    b--;

  /* Here, r[i].a <= p for 0 <= i < b, r[i].a == p for b <= i < h,
     r[i].a > p for h <= i < l */

#ifdef WANT_ASSERT
  {
    unsigned long i;
    for (i = 0; i < b; i++)
      {
	ASSERT (r[i].a <= p);
      }
    for (i = b; i < h; i++)
      {
	ASSERT (r[i].a == p);
      }
    for (i = h; i < l; i++)
      {
	ASSERT (r[i].a > p);
      }
  }
#endif

  sort_sieve_reports_recurse (r, b);
  sort_sieve_reports_recurse (r + h, l - h);

#ifdef WANT_ASSERT
  {
    unsigned long i;
    for (i = 0; i < l - 1; i++)
      {
	ASSERT (r[i].a <= r[i + 1].a);
      }
  }
#endif
}

static void
sort_sieve_reports (sieve_reportbuffer_t *reports)
{
  sort_sieve_reports_recurse (reports->reports, reports->nr);
}

static void
sieve_block (const int lvl, unsigned char *sievearray, 
	     const unsigned long arraylen, factorbase_t fb, 
	     long long *times)
{
  clock_t tsc1, tsc2;
  unsigned long blockstart;
  const unsigned long b = fb->fbsmallbound[lvl];

  for (blockstart = 0; blockstart < arraylen; blockstart += b)
    {
      const unsigned long blocklen = MIN(b, arraylen - blockstart);

      /* Sieve smaller blocks */
      if (lvl > 0)
	sieve_block (lvl - 1, sievearray + blockstart, blocklen, fb, times);

      tsc1 = clock ();
      sieve_small_slow (sievearray + blockstart, fb->fbinit[lvl], blocklen);
      tsc2 = clock ();
      times[lvl] += tsc2 - tsc1;
    }
}


static unsigned long
sieve_one_side (unsigned char *sievearray, factorbase_t fb,
		sieve_reportbuffer_t *reports, 
		const unsigned char reports_threshold,
		const long amin, const long amax, const unsigned long b,
		const unsigned long proj_roots, const double log_scale,
		const double *dpoly, const unsigned int deg, 
		sieve_reportbuffer_t *other_reports, const int verbose)
{
  clock_t tsc1, tsc2;
  long long times[SIEVE_BLOCKING];
  int lvl;
  const int odd = 1 - (b & 1); /* If odd=1, only odd $a$ are sieved */
  const long eff_amin = amin + ((odd && (amin & 1) == 0) ? 1 : 0);
  const long eff_amax = amax - ((odd && (amax & 1) == 0) ? 1 : 0);
  const unsigned long l = ((eff_amax - eff_amin) >> odd) + 1;
  const int find_smallprime_reports = 1;

  fb_disable_roots (fb->fblarge, b, verbose);
  
  compute_norms (sievearray, eff_amin, eff_amax, b, dpoly, deg, 
                 (double) proj_roots, log_scale, odd, verbose);

  if (other_reports != NULL)
    {
      /* Set sievearray[n] to 255 if that $n$ is not in the list of
	 reports on the other side - we know this $a$ can never become
	 a relation and we don't want it to produce unnecessary sieve
	 reports */
      unsigned long startidx = 0, endidx;
      size_t i;
      for (i = 0; i < other_reports->nr; i++)
	{
	  endidx = a_to_sieveidx (other_reports->reports[i].a, eff_amin, odd);
	  /* Set sievearray[startidx ... endidx - 1] to 0 */
	  if (endidx > startidx)
	    memset (sievearray + startidx, 255, endidx - startidx);
	  startidx = endidx + 1;
	}
    }

  /* Init small primes for sieving this line */
  tsc1 = clock ();
  for (lvl = 0; lvl < SIEVE_BLOCKING; lvl++)
      fb_initloc_small (fb->fbinit[lvl], fb->fbsmall[lvl], eff_amin, b, odd);
  tsc2 = clock ();
  if (SIEVE_BLOCKING > 0 && verbose)
    printf ("# Initing small primes fb took %ld clocks\n", 
	    (long int) (tsc2 - tsc1));

  for (lvl = 0; lvl < SIEVE_BLOCKING; lvl++)
    times[lvl] = 0;

  /* Sieve small primes with multi-level blocking */
  if (SIEVE_BLOCKING > 0)
    sieve_block (SIEVE_BLOCKING - 1, sievearray, l, fb, times);
  
  if (SIEVE_BLOCKING > 0 && verbose)
    {
      printf ("# Sieving small prime blocks took: ");
      for (lvl = 0; lvl < SIEVE_BLOCKING; lvl++)
	printf ("level %d, %lld clocks;", lvl + 1, times[lvl]);
      printf ("\n");
    }

#ifdef TRACE_RELATION_A
  if (eff_amin <= TRACE_RELATION_A && TRACE_RELATION_A <= eff_amax)
    printf ("# TRACE RELATION a = %ld in sieve_one_side: approx norm after "
	    "sieving small primes is %d\n", (long) TRACE_RELATION_A, 
	    sievearray[(TRACE_RELATION_A - eff_amin) >> odd]);
#endif


  /* Set the report buffer to empty */
  reports->nr = 0;

  /* If the factor base limit is very small compared to the 
     largest block length, we may get many relations that are 
     smooth over the small fb primes and will not produce sieve 
     reports during sieving the large fb primes.
     Those reports are found here instead. */
  if (find_smallprime_reports)
    {
      tsc1 = clock ();
      find_sieve_reports (sievearray, reports, reports_threshold, 
			  eff_amin, l, b, odd, other_reports);
      tsc2 = clock ();
      if (verbose)
	{
	  printf ("# Finding sieve reports took %ld clocks\n", 
		  (long int) (tsc2 - tsc1));
	  printf ("# There were %lu sieve reports after sieving small "
		  "primes\n", reports->nr);
	}
    }
  
  tsc1 = clock ();
  sieve (sievearray, fb->fblarge, eff_amin, eff_amax, b, reports_threshold, 
	 reports, odd);
  tsc2 = clock ();
  if (verbose)
    printf ("# Sieving large fb primes took %ld clocks\n", 
	    (long int) (tsc2 - tsc1));
  
  fb_restore_roots (fb->fblarge, b, verbose);
  
  /* Sort the sieve reports in order of increasing a value */
  sort_sieve_reports (reports);
  
  return reports->nr;
}

static inline void
add_fbprime_to_list (fbprime_t *list, unsigned int *cursize, 
		     const unsigned int maxsize, const fbprime_t toadd)
{
  if ((*cursize) < maxsize)
    list[(*cursize)++] = toadd;
}


/* Divides the prime q and its powers out of C, appends each q divided out 
   to primes_a, subject to the length restriction *nr_primes < max_nr_primes. 
   Returns the exponent of q that divided (0 if q didn't divide at all). */

static inline int
trialdiv_one_prime (const fbprime_t q, mpz_t C, unsigned int *nr_primes,
                    fbprime_t *primes, const unsigned int max_nr_primes,
		    const int do_powers)
{
  int nr_divide = 0;

  while (mpz_divisible_ui_p (C, (unsigned long) q))
    {
      unsigned long r;
      nr_divide++;
      if (primes != NULL)
	add_fbprime_to_list (primes, nr_primes, max_nr_primes, q);
      r = mpz_tdiv_q_ui (C, C, (unsigned long) q);
      ASSERT (r == 0);
      if (!do_powers)
	break;
    }
  return nr_divide;
}

/* Factor an unsigned long into fb primes, adds primes to list. Not fast. */
static inline void 
trialdiv_slow (unsigned long C, unsigned int *nr_primes, fbprime_t *primes, 
	       const unsigned int max_nr_primes)
{
  fbprime_t q;

  while (C > 1)
    {
      q = iscomposite (C);
      if (q == 0)
	{
	  q = (fbprime_t) C;
	  ASSERT ((unsigned long) q == C); /* Check for truncation */
	}
      C /= (unsigned long) q;
      if (primes != NULL)
	add_fbprime_to_list (primes, nr_primes, max_nr_primes, q);
    }
}


/* Functions for finding prime factor of the norm during refactoring */

static inline int
trialdiv_with_root (factorbase_degn_t *fbptr, const long a, 
		    const unsigned long b)
{
  int i;
  modulus_t m;
  residue_t r;
  unsigned long absa;

  mod_initmod_ul (m, (unsigned long) fbptr->p);
  mod_init_noset0 (r, m);

  absa = labs(a);

  for (i = 0; i < fbptr->nr_roots ; i++)
    {
      if (a <= 0) /* We compute rb - a */
	mod_set_ul_reduced (r, (unsigned long) fbptr->roots[i], m);
      else /* We compute (-r)b + a */
	mod_set_ul_reduced (r, (unsigned long) (fbptr->p - fbptr->roots[i]), 
	                    m);

      mod_muladdredc_ul (r, r, b, absa, (unsigned long) fbptr->invp, m);
      
      if (mod_get_ul (r, m) == 0UL)
	{
	  mod_clear (r, m);
	  mod_clearmod (m);
	  return 1;
	}
    }

  mod_clear (r, m);
  mod_clearmod (m);
  return 0;
}


/* Returns 1 if fbptr->p divides norm, 0 otherwise */

static inline int
trialdiv_with_norm (factorbase_degn_t *fbptr, const mpz_t norm)
{
  modulus_t m;
  residue_t r;
  size_t i;
  int j;

  ASSERT (mpz_sgn (norm) >= 0);

  mod_initmod_ul (m, (unsigned long) fbptr->p);
  mod_init (r, m);
  
  for (i = 0; i < mpz_size (norm); i++)
    modul_addredc_ul (r, r, mpz_getlimbn (norm, i), fbptr->invp, m);

  j = (mod_get_ul (r, m) == 0);

  mod_clear (r, m);
  mod_clearmod (m);
  return j;
}


/* Same, assuming norm has 1 limb (plus the add value) */
static inline int
trialdiv_with_norm1 (factorbase_degn_t *fbptr, const mpz_t norm, 
		     const unsigned long add)
{
  modulus_t m;
  residue_t r;
  int i;

  ASSERT (mpz_sgn (norm) > 0);

  mod_initmod_ul (m, (unsigned long) fbptr->p);
  mod_init_noset0 (r, m);

  modul_redcsemi_ul_not0 (r, mpz_getlimbn (norm, (mp_size_t) 0), 
                          fbptr->invp, m);

  i = (mod_get_ul (r, m) + add == 0UL || mod_get_ul (r, m) + add == m[0]);

  mod_clear (r, m);
  mod_clearmod (m);
  return i;
}


/* Same, assuming norm has 2 limbs (plus the add value) */
static inline int
trialdiv_with_norm2 (factorbase_degn_t *fbptr, const mpz_t norm,
		     const unsigned long add)
{
  modulus_t m;
  residue_t r;
  int i;

  ASSERT (mpz_sgn (norm) > 0);

  mod_initmod_ul (m, (unsigned long)fbptr->p);
  mod_init_noset0 (r, m);

  modul_redcsemi_ul_not0 (r, mpz_getlimbn (norm, (mp_size_t) 0), 
                          fbptr->invp, m);
  modul_addredcsemi_ul (r, r, mpz_getlimbn (norm, (mp_size_t) 1), 
                        fbptr->invp, m);

  i = (mod_get_ul (r, m) + add == 0UL || mod_get_ul (r, m) + add == m[0]);

  mod_clear (r, m);
  mod_clearmod (m);
  return i;
}


/* Same, with assuming norm has 3 limbs (plus the add value) */
static inline int
trialdiv_with_norm3 (factorbase_degn_t *fbptr, const mpz_t norm,
		     const unsigned long add)
{
  modulus_t m;
  residue_t r;
  int i;

  ASSERT (mpz_sgn (norm) > 0);

  mod_initmod_ul (m, (unsigned long) fbptr->p);
  mod_init_noset0 (r, m);

  modul_redcsemi_ul_not0 (r, mpz_getlimbn (norm, (mp_size_t) 0), 
                          fbptr->invp, m);
  modul_addredcsemi_ul (r, r, mpz_getlimbn (norm, (mp_size_t) 1), 
                        fbptr->invp, m);
  modul_addredcsemi_ul (r, r, mpz_getlimbn (norm, (mp_size_t) 2), 
                        fbptr->invp, m);

  i = (mod_get_ul (r, m) + add == 0UL || mod_get_ul (r, m) + add == m[0]);

  mod_clear (r, m);
  mod_clearmod (m);
  return i;
}


/* Same, with assuming norm has 4 limbs (plus the add value) */
static inline int
trialdiv_with_norm4 (factorbase_degn_t *fbptr, const mpz_t norm,
		     const unsigned long add)
{
  modulus_t m;
  residue_t r;
  int i;

  ASSERT (mpz_sgn (norm) > 0);

  mod_initmod_ul (m, (unsigned long) fbptr->p);
  mod_init_noset0 (r, m);

  modul_redcsemi_ul_not0 (r, mpz_getlimbn (norm, (mp_size_t) 0), 
                          fbptr->invp, m);
  modul_addredcsemi_ul (r, r, mpz_getlimbn (norm, (mp_size_t) 1), 
                        fbptr->invp, m);
  modul_addredcsemi_ul (r, r, mpz_getlimbn (norm, (mp_size_t) 2), 
                        fbptr->invp, m);
  modul_addredcsemi_ul (r, r, mpz_getlimbn (norm, (mp_size_t) 3), 
                        fbptr->invp, m);

  i = (mod_get_ul (r, m) + add == 0UL || mod_get_ul (r, m) + add == m[0]);

  mod_clear (r, m);
  mod_clearmod (m);
  return i;
}


static unsigned long
ul_rho (const unsigned long N, const unsigned long cparm)
{
  modulusul_t m;
  residueul_t r1, r2, c, diff, accu;
  unsigned long invm, g;
  int i;
  unsigned int iterations = 0;
  const int iterations_between_gcd = 32;
  clock_t tsc1, tsc2;

  ASSERT (cparm < N);

  tsc1 = clock ();

  /* printf ("ul_rho: factoring %lu\n", N); */

  modul_initmod_ul (m, N);
  modul_init_noset0 (r1, m);
  modul_init_noset0 (r2, m);
  modul_init_noset0 (c, m);
  modul_init_noset0 (accu, m);
  modul_init_noset0 (diff, m);
  invm = -modul_invmodlong (m);
  modul_set_ul_reduced (c, cparm, m); /* We set c to cparm and use that as 
					 the Montgomery representation of 
					 cparm*2^-w % m */
  modul_set_ul_reduced (r1, 2UL, m);
  modul_set (r2, r1, m);
  modul_set (accu, r1, m);

  do {
    for (i = 0; i < iterations_between_gcd; i++)
      {
	/* FIXME: use modul_muladdredc_ul()? How does that change the 
	   arithmetic? */
	modul_mulredc (r1, r1, r1, invm, m);
	modul_add (r1, r1, c, m);

      	modul_mulredc (r2, r2, r2, invm, m);
	modul_add (r2, r2, c, m);
	modul_mulredc (r2, r2, r2, invm, m);
	modul_add (r2, r2, c, m);

	modul_sub (diff, r1, r2, m);
	modul_mulredc (accu, accu, diff, invm, m);
      }
    iterations += iterations_between_gcd;
  } while ((g = modul_gcd (accu, m)) == 1);

  modul_clear (r1, m);
  modul_clear (r2, m);
  modul_clear (c, m);
  modul_clear (accu, m);
  modul_clear (diff, m);
  modul_clearmod (m);

  tsc2 = clock ();
  ul_rho_time += tsc2 - tsc1;

  /* printf ("ul_rho: took %u iterations to find %lu\n", iterations, g); */

  return g;
}


static unsigned long
mpz_rho (const mpz_t N, const unsigned long c)
{
  mpz_t m, r1, r2, t, accu;
  int i;
  unsigned int iterations = 0;
  const int iterations_between_gcd = 32;
  unsigned long f;

  ASSERT (mpz_cmp_ui(N, c) > 0);

#if 0 || defined(RHODEBUG)
  gmp_printf ("mpz_rho: factoring %Zd\n", N);
#endif

  mpz_init (m);
  mpz_init (r1);
  mpz_init (r2);
  mpz_init (t);
  mpz_init (accu);

  mpz_set (m, N);
  mpz_set_ui (r1, 2UL);
  mpz_set_ui (r2, 2UL);
  mpz_set_ui (accu, 2UL);

  do {
    for (i = 0; i < iterations_between_gcd; i++)
      {
	mpz_mul (t, r1, r1);
	mpz_add_ui (t, t, c);
	mpz_mod (r1, t, m);

	mpz_mul (t, r2, r2);
	mpz_add_ui (t, t, c);
	mpz_mod (r2, t, m);
	mpz_mul (t, r2, r2);
	mpz_add_ui (t, t, c);
	mpz_mod (r2, t, m);

	mpz_sub (t, r1, r2);
	mpz_mul (t, t, accu);
	mpz_mod (accu, t, m);
      }
    iterations += iterations_between_gcd;
    mpz_gcd (t, accu, m);
  } while (mpz_cmp_ui (t, 1UL) == 0);

  /* FIXME: deal with too large factors better than this! */
  if (!mpz_fits_ulong_p (t))
    f = 1;
  else
    f = mpz_get_ui (t);


  mpz_clear (accu);
  mpz_clear (t);
  mpz_clear (r2);
  mpz_clear (r1);
  mpz_clear (m);

#if 0 || defined(RHODEBUG)
  {
    static double ratio = 0., nr = 0.;
    double r = (double)iterations/sqrt((double)f);
    ratio += r;
    nr += 1.;
    printf ("mpz_rho: took %u iterations to find %lu, i/sqrt(p) = %.3f, "
	    "avg = %.3f\n", 
	    iterations, f, r, ratio/nr);
    /* Interestingly, the average comes out as almost exactly 1! */
  }
#endif

  return f;
}


/* If n <= 10^10 then this function proves primality/compositeness
   of n. For larger n, return value 0 means n is composite, return 
   value 1 means it is probably prime */

/* if n == 1 (mod 3), then the only pseudoprimes < 10^10 
   to bases 2, 5 and 7 are 3014101261, 3215031751, 7535192941,
   to bases 2, 7 and 61 are 4759123141, 8411807377 */

/* If n == 2 (mod 3), then the only pseudoprimes < 10^13 to bases
   2, 3 and 5 are 244970876021, 405439595861, 1566655993781, 3857382025841,
   4074652846961, 5783688565841 (i.e. there are none < 10^10) */

int
ul_proven_prime (unsigned long n)
{
  modulusul_t m;
  residueul_t b;
  unsigned long invn;
  int r = 0;
  
  if (n % 2UL == 0UL) /* We want odd modulus for REDC */
    {
      r = (n == 2UL);
      goto end;
    }

  modul_initmod_ul (m, n);
  modul_init_noset0 (b, m);
  invn = -mod_invmodlong (m);

  modul_set_ul (b, 2UL, m);
  if (!modul_sprp (b, invn, m))
    goto end;

  if (n < 2047UL)
    {
      r = 1;
      goto end;
    }

  if (n % 3UL == 1UL)
    {
      modul_set_ul (b, 7UL, m);
      if (modul_sprp (b, invn, m))
	{
	  modul_set_ul (b, 61UL, m);
	  if (modul_sprp (b, invn, m)
#if (ULONG_MAX > 4294967295UL)
	      && n != 4759123141UL && n != 8411807377UL
#endif
	      )
	    r = 1;
	}
    }
  else
    {
      /* Case n % 3 == 0, 2 */
      
      modul_set_ul (b, 3UL, m);
      if (modul_sprp (b, invn, m))
	{
	  modul_set_ul (b, 5UL, m);
	  if (modul_sprp (b, invn, m))
	    r = 1;
	}
    }
  
 end:
#if defined(PARI)
  printf ("isprime(%lu) == %d /* PARI */\n", n, r);
#endif
  modul_clearmod (m);
  return r;
}


/* Returns a pointer to factor base entry which is either the 
   end-of-factorbase marker, or the first entry whose p divides norm.
   Requires that norm >= fbptr->p */

static factorbase_degn_t *
trialdiv_find_next (factorbase_degn_t *fbptr, const mpz_t norm)
{
  size_t s = mpz_size (norm);
  unsigned long add = 0;

  ASSERT (mpz_cmp_ui (norm, (unsigned long) fbptr->p) >= 0);

  if (mpz_getlimbn (norm, (mp_size_t) (s - 1)) < fbptr -> p)
    {
      add = mpz_getlimbn (norm, (mp_size_t) (s - 1));
      s--;
    }

  /* Choose good trial division routine outside of loop over fb primes */
  if (s == 1)
    {
      for (; fbptr->p != 0; fbptr = fb_next (fbptr))
	if (trialdiv_with_norm1 (fbptr, norm, add))
	  break;
    }
  else if (s == 2)
    {
      for (; fbptr->p != 0; fbptr = fb_next (fbptr))
	if (trialdiv_with_norm2 (fbptr, norm, add))
	  break;
    }
  else if (s == 3)
    {
      for (; fbptr->p != 0; fbptr = fb_next (fbptr))
	if (trialdiv_with_norm3 (fbptr, norm, add))
	  break;
    }
  else if (s == 4)
    {
      for (; fbptr->p != 0; fbptr = fb_next (fbptr))
	if (trialdiv_with_norm4 (fbptr, norm, add))
	  break;
    }
  else
    {
      for (; fbptr->p != 0; fbptr = fb_next (fbptr))
	if (trialdiv_with_norm (fbptr, norm))
	  break;
    }

  return fbptr;
}

/* 1. Compute the norm
   2. Divide out primes with projective roots
   3. Divide out the report prime(s), remembering the smallest one and 
      remembering the smallest approximate log from the reports
   4. Trial divide up to the smallest report prime (or up to the factor base 
      limit if there was no report prime)
   5. Check if the cofactor is small enough
   6. Check if the cofactor is < lpb, if yes report relation, assuming that 
      the cofactor is a prime
   7. Check if the cofactor is a prp, if yes skip to next report
   
   Factoring composite cofactors into large primes is left to the 
   next step in the tool chain.
*/

static int
trialdiv_one_side (mpz_t norm, mpz_t *scaled_poly, int degree, 
		   const long a, const unsigned long b, 
		   fbprime_t *primes, unsigned int *nr_primes, 
		   const unsigned int max_nr_primes,
		   const unsigned long proj_divisor, 
		   const unsigned int nr_proj_primes, 
		   const fbprime_t *proj_primes,
		   factorbase_degn_t *fullfb,
		   const sieve_reportbuffer_t *reports,
		   size_t report_idx, 
		   unsigned int *cof_toolarge, unsigned int *lp_toolarge,
		   const unsigned long fbb, const int lpb, const int mfb, 
		   const double lambda, const double log_scale)
{
  unsigned int k;
  factorbase_degn_t *fbptr;
  size_t log2size;
  unsigned char sievelog, finallog;
  const double log_proj_divisor = log ((double) proj_divisor) * log_scale;
  double dpoly[MAXDEGREE + 1];
  int r;
  fbprime_t p;
  const unsigned char MAX_SIEVELOG_ERROR = 1;

  /* 1. Compute norm */
  mp_poly_eval (norm, scaled_poly, degree, a);
  mpz_abs (norm, norm);
  TRACE_A (a, __func__, __LINE__, "norm for degree %d is %Zd\n", degree, norm);
  
  /* Compute approx. log of norm as was used for initialising the 
     sieve array. We assume that the absolute difference between this 
     sievelog here the the value from initialisation the sieve array is no 
     greater than MAX_SIEVELOG_ERROR */
  for (k = 0; (int) k <= degree; k++)
    dpoly[k] = mpz_get_d (scaled_poly[k]);
  sievelog = log_norm (dpoly, degree, (double) a, log_scale, log_proj_divisor);
  TRACE_A (a, __func__, __LINE__, "sievelog = %hhu\n", sievelog);
  
  
  /* 2. Divide the primes with projective roots and their powers out of norm */
  for (k = 0; k < nr_proj_primes; k++)
    {
      /* proj_primes contains repeated prime factors */
      r = trialdiv_one_prime (proj_primes[k], norm, nr_primes, primes, 
			      max_nr_primes, 0);
      ASSERT_ALWAYS (r > 0);
      TRACE_A (a, __func__, __LINE__, "dividing out " FBPRIME_FORMAT "^%d " 
	       "with proj. root. New norm is %Zd.\n", proj_primes[k], r, norm);
    }
  for (k = 0; k < nr_proj_primes; k++)
    {
      /* These primes may divide norm in a higher power than proj_divisor */
      r = trialdiv_one_prime (proj_primes[k], norm, nr_primes, primes, 
			      max_nr_primes, 1);
      TRACE_A (a, __func__, __LINE__, "dividing out " FBPRIME_FORMAT "^%d " 
	       "with proj. root. New norm is %Zd.\n", proj_primes[k], r, norm);
    }
  
  
  /* 3. Divide the report primes out of this norm and find the largest 
     approximate log */
  /* There must be at least one valid report */
  ASSERT_ALWAYS (report_idx < reports->nr);
  ASSERT_ALWAYS (reports->reports[report_idx].a == a);
  finallog = reports->reports[report_idx].l;
  p = reports->reports[report_idx].p;
  while (report_idx < reports->nr && reports->reports[report_idx].a == a)
    {
      /* We allow reports without a report prime. In that case the 
	 sieve_report_t contains p == 1. */
      if (reports->reports[report_idx].p != (fbprime_t) 1)
	{
	  r = trialdiv_one_prime (reports->reports[report_idx].p, norm, 
				  nr_primes, primes, max_nr_primes, 1);
	  ASSERT_ALWAYS (r > 0); /* If a report prime is listed but does not
				    divide the norm, there is a serious bug
				    in the sieve code */
	  TRACE_A (a, __func__, __LINE__, "dividing out report prime " 
		   FBPRIME_FORMAT "^%d. New norm = %Zd\n",
		   reports->reports[report_idx].p, r, norm);
	}

      /* Find the largest approximate log. */
      if (add_error (reports->reports[report_idx].l) > add_error (finallog))
	{
	  finallog = reports->reports[report_idx].l;
	  p = reports->reports[report_idx].p;
	}
      
      report_idx++;
    }
  /* If there was a report prime in the report with the largest log,
     add this prime's log (we want to reach the log the sieve had before
     dividing out this prime). */
  if (p != 1)
    finallog += fb_log ((double) p, log_scale, 0.);
  TRACE_A (a, __func__, __LINE__, "finallog = %hhu\n", finallog);

  /* 
     When we divide out a prime, we subtract log(p) from sievelog. When 
     sievelog +- MAX_SIEVELOG_ERROR == finallog, we know we are done
     (assuming the remaining factor base primes have log > MAX_SIEVELOG_ERROR)
     When sievelog +- MAX_SIEVELOG_ERROR - finallog < 2*log(p), we know that 
     exactly one more factor base prime divides.
  */

  /* Treat factor base prime p == 2 separately. We want to be able to use
     REDC in what follows and it requires an odd modulus */
  fbptr = fullfb;
  TRACE_A (a, __func__, __LINE__, "first fb prime is " FBPRIME_FORMAT "\n", 
	   fbptr->p);
  if (fbptr->p == (fbprime_t) 2)
    {
      if (mpz_even_p (norm))
	{
	  sievelog -= fbptr->plog;
	  r = trialdiv_one_prime ((fbprime_t) 2, norm, nr_primes, primes, 
				  max_nr_primes, 1);
	  ASSERT (r > 0);
	  TRACE_A (a, __func__, __LINE__, "dividing out fb prime 2^%d." 
		   "New norm = %Zd, sievelog = %hhu\n", r, norm, sievelog);
	}
      fbptr = fb_next (fbptr);
    }
  else
    {
      ASSERT_ALWAYS (!mpz_even_p (norm));
    }
  
  /* 4. Go through the factor base until 
     sievelog +- MAX_SIEVELOG_ERROR == finallog or 
     sievelog +- MAX_SIEVELOG_ERROR - finallog < 2*log(next fb base prime).
  */

  /* If sievelog <= finallog + MAX_SIEVELOG_ERROR, dividing out another 
     fb prime would surely push sievelog < finallog - MAX_SIEVELOG_ERROR
     (assuming log(fbprime) > 2*MAX_SIEVELOG_ERROR) + 1, so we can stop. */

  /* We only abort early if we know for sure that exactly 1 fb prime is 
     missing. Hence we loop while
     sievelog + MAX_SIEVELOG_ERROR - finallog >= 2*log(next fb base prime)
  */

  while (add_error (sievelog) > add_error (finallog) + MAX_SIEVELOG_ERROR && 
	 uc_sub(sievelog, finallog) + MAX_SIEVELOG_ERROR >= 2*fbptr->plog)
    {
      fbptr = trialdiv_find_next (fbptr, norm);
      
      if (fbptr->p == 0)
	{
	  fprintf (stderr, "Warning, reached end of fb for (%ld, %lu). "
		   "sievelog = %hhu, finallog = %hhu\n", 
		   a, b, sievelog, finallog);
	  gmp_fprintf (stderr, "Remaining norm = %Zd\n", norm);
	  break;
	}
      
      r = trialdiv_one_prime (fbptr->p, norm, nr_primes, primes,
			      max_nr_primes, 1);
      ASSERT_ALWAYS (r > 0);
      sievelog -= fbptr->plog;
      TRACE_A (a, __func__, __LINE__, "dividing out fb prime " FBPRIME_FORMAT
	       "^%d. New norm = %Zd, sievelog = %hhu\n", 
	       fbptr->p, r, norm, sievelog);
      fbptr = fb_next (fbptr);
    }

  /* Now at least one of these conditions hold:
     1. add_error (sievelog) <= add_error (finallog)
     2. uc_sub(sievelog, finallog) < 2*fbptr->plog
     3. fbptr->p == 0 

     If 1. or 3. hold, trial division over factor base primes is finished.
     If only 2. holds, there is exactly one more fb prime left in norm.
  */

  if (add_error (sievelog) > add_error (finallog) &&
      fbptr->p != 0)
    {
      factorbase_degn_t *oldfbptr = fbptr;
      /* There is exactly one more fb prime to divide out. */
      
      ASSERT (uc_sub(sievelog, finallog) < 2*fbptr->plog);
      /* That prime should have log size equal to sievelog - finallog,
	 assuming the sieve was initialised correctly. In case it was 
	 initialised incorrectly, we allow an extra 1 */
      
      if (add_error (sievelog) > add_error (finallog) && 
	  uc_sub (sievelog, finallog) + 1 < fbptr->plog)
	{
	  gmp_fprintf (stderr, "Error, a = %ld, b = %lu, sievelog = %hhu, "
		       "finallog = %hhu, p = " FBPRIME_FORMAT ", plog = %d. "
		       "Remaining norm = %Zd\n", 
		       a, b, sievelog, finallog, fbptr->p, (int) fbptr->plog, 
		       norm);
	  abort ();
	}
      
      /* Jump ahead in the factor base to primes of approximately 
	 the correct size. */
      nrprimes2++;
      sumprimes2 += (unsigned long) fbptr->p;
      
#if TRIALDIV_SKIPFORWARD
      const unsigned char missinglog = uc_sub (sievelog, finallog);
      /* We know add_error (sievelog) > add_error (finallog), so
	 missinglog is positive */
      /* Skip forward to the factor base prime of the right size */
      TRACE_A (a, __func__, __LINE__, "skipping forward in fb to prime of "
	       "log norm %hhu\n", missinglog);
      while (fbptr->p != 0 && 
	     fbptr->plog < missinglog)
	fbptr = fb_next (fbptr);
      TRACE_A (a, __func__, __LINE__, "skipped forward in fb, next prime "
	       "is " FBPRIME_FORMAT "\n", fbptr->p);
#endif
      fbptr = trialdiv_find_next (fbptr, norm);
      if (fbptr->p == (fbprime_t) 0)
	{
	  /* Reached end of factor base? Maybe sieve was not initialised
	     correctly and finallog is one too small, causing us to skip
	     ahead too far. Try again without skipping */
	  fbptr = trialdiv_find_next (oldfbptr, norm);
	}
      /* This time is must have worked */
      ASSERT_ALWAYS (fbptr->p != 0);
      
      r = trialdiv_one_prime (fbptr->p, norm, nr_primes, primes,
			      max_nr_primes, 1);
      ASSERT_ALWAYS (r > 0);
      sievelog -= fbptr->plog;
      TRACE_A (a, __func__, __LINE__, "dividing out fb prime " 
	       FBPRIME_FORMAT "^%d. New norm = %Zd, sievelog = %hhu\n", 
	       fbptr->p, r, norm, sievelog);
      /* Now we should have sievelog == finallog +- 1 */
      ASSERT_ALWAYS (add_error(sievelog) <= add_error(finallog) + 1 &&
		     add_error(sievelog) + 1 >= add_error(finallog));
    }

  /* Now we should have all factor base primes < the smallest report prime,
     and if the report buffer did not overflow, all primes >= the smallest
     report prime as well. */


  /* Lets try to factor the rest, i.e. remaining fb primes (if report buffer
     overflowed) and large primes.
     We still apply the mfb restriction, so this may discard relations 
     where not all report primes were given */
  
  /* 5. Check if the cofactor is small enough */
  log2size = mpz_sizeinbase (norm, 2);
  TRACE_A (a, __func__, __LINE__, "log2size of cofactor %Zd is %d\n", 
	   norm, (int) log2size);

  if ((double) log2size > lpb * lambda + SIEVE_PERMISSIBLE_ERROR)
    {
      gmp_fprintf (stderr, "Sieve report (%ld, %lu) is not smooth for "
                   "degree %d poly, cofactor is %Zd with %d bits\n", 
		   a, b, degree, norm, mpz_sizeinbase (norm, 2));
      return 0;
    }
  if ((int) log2size > mfb)
    {
      (*cof_toolarge)++;
      TRACE_A (a, __func__, __LINE__, "log2size %d > mfb %d, discarding "
	       "relation\n", (int) log2size, (int) mfb);
      return 0;
    }

  /* If cofactor is 1, nothing left to do */
  if (mpz_cmp_ui (norm, 1UL) == 0)
    {
      TRACE_A (a, __func__, __LINE__, "Cofactor is 1, this side is smooth\n");
      return 1;
    }

  p = fbptr->p;
  if (p == 0)
    p = fbb;
  /* All remaining prime factors of norm should be > p now, hence if 
     norm < p^2, it should be prime. */
  ASSERT(mpz_cmp_ui (norm, (unsigned long) p) > 0);
  {
    mpz_t t;
    mpz_init (t);
    mpz_set_ui (t, (unsigned long) p);
    mpz_mul_ui (t, t, (unsigned long) p);
    r = mpz_cmp (norm, t);
    mpz_clear (t);
  }
  /* If r < 0, assert it's really a prime */
  ASSERT (r >= 0 || mpz_probab_prime_p (norm, PRP_REPS));
  if (r < 0 || mpz_probab_prime_p (norm, PRP_REPS))
    {
      /* It's a prime */
      /* sizeinbase2 (n) = k  <==>  2^(k-1) <= n < 2^k */
      if (mpz_sizeinbase (norm, 2) <= (size_t) lpb)
	{
	  /* Prime and <2^lpb: relation is smooth */
	  TRACE_A (a, __func__, __LINE__, "cofactor %Zd is prime and <2^lpb"
		   " = 2^%d. This side is smooth.\n", norm, lpb);
	  ASSERT (mpz_fits_ulong_p (norm));
	  add_fbprime_to_list (primes, nr_primes, max_nr_primes, 
			       (fbprime_t) mpz_get_ui (norm));
	  return 1;
	}
      else
	{
	  /* Prime and >2^lpb: relation is not smooth */
	  TRACE_A (a, __func__, __LINE__, "cofactor %Zd is prime and >2^lpb"
		   " = 2^%d. This side is not smooth.\n", norm, lpb);
	  (*lp_toolarge)++;
	  return 0;
	}
    }

  TRACE_A (a, __func__, __LINE__, "cofactor %Zd is not prime and < 2^mfb, "
           "keeping relation\n", norm);

  /* So the cofactor is not a prime. Later we'll try to factor it somehow,
     for now we just print the relation with the large primes missing. */
  return 1;
}


/* return the number of printed relations */
static unsigned int
trialdiv_and_print (cado_poly poly, const unsigned long b, 
                    const sieve_reportbuffer_t *reports_a, 
		    const sieve_reportbuffer_t *reports_r, 
		    factorbase_t fba, factorbase_t fbr, 
		    const double log_scale, const int verbose)
{
  clock_t tsc1, tsc2;
  unsigned long proj_divisor_a, proj_divisor_r;
  unsigned int i, j, k;
  const unsigned int max_nr_primes = 128, max_nr_proj_primes = 16;
  mpz_t Fab, Gab, scaled_poly_a[MAXDEGREE], scaled_poly_r[2];
  fbprime_t primes_a[max_nr_primes], primes_r[max_nr_primes];
  fbprime_t proj_primes_a[max_nr_primes], proj_primes_r[max_nr_primes];
  unsigned int nr_primes_a, nr_primes_r, nr_proj_primes_a, nr_proj_primes_r;
  unsigned int matching_reports = 0, relations_found = 0;
  unsigned int lp_a_toolarge = 0, lp_r_toolarge = 0, 
    cof_a_toolarge = 0, cof_r_toolarge = 0, not_coprime = 0;
  int ok;

  mpz_init (Fab);
  mpz_init (Gab);
  for (i = 0; i <= (unsigned) poly->degree; i++)
    mpz_init (scaled_poly_a[i]);
  mpz_init (scaled_poly_r[0]);
  mpz_init (scaled_poly_r[1]);

  ul_rho_time = 0LL;
  tsc1 = clock ();

  /* Multiply f_i (and g_i resp.) by b^(deg-i) and put in scaled_poly */
  mp_poly_scale (scaled_poly_a, poly->f, poly->degree, b, -1); 
  mp_poly_scale (scaled_poly_r, poly->g, 1, b, -1); 

  /* Factor primes with projective roots on algebraic side. Prime factors
     go in the proj_primes_a[] list, the number of prime factors in 
     nr_proj_primes_a, the product is kept in proj_divisor_a. */
  nr_proj_primes_a = 0;
  proj_divisor_a = mpz_gcd_ui (NULL, poly->f[poly->degree], b);
  trialdiv_slow (proj_divisor_a, &nr_proj_primes_a, proj_primes_a, 
		 max_nr_proj_primes);

  /* Same for rational side */
  nr_proj_primes_r = 0;
  proj_divisor_r = mpz_gcd_ui (NULL, poly->g[1], b);
  trialdiv_slow (proj_divisor_r, &nr_proj_primes_r, proj_primes_r, 
		 max_nr_proj_primes);

  sumprimes = nrprimes = 0;
  sumprimes2 = nrprimes2 = 0;
  ul_rho_toolarge_nr = ul_rho_toolarge_sum = 0;

  for (i = 0, j = 0; i < reports_a->nr && j < reports_r->nr;)
    {
      if (reports_a->reports[i].a == reports_r->reports[j].a)
	if (gcd_ul ((unsigned long) labs(reports_a->reports[i].a), b) > 1UL)
	  {
	    not_coprime++;
	  }
	else
	{
          const long a = reports_a->reports[i].a;

	  matching_reports++;
      
          /* Do the algebraic side */
	  nr_primes_a = 0;
	  ok = trialdiv_one_side (Fab, scaled_poly_a, poly->degree, a, b, 
				  primes_a, &nr_primes_a, max_nr_primes,
				  proj_divisor_a, nr_proj_primes_a, 
				  proj_primes_a, fba->fullfb, 
				  reports_a, i,
				  &cof_a_toolarge, &lp_a_toolarge, 
				  poly->alim, poly->lpba, poly->mfba, 
				  poly->alambda, log_scale);

	  if (!ok) 
	    goto nextreport;

          /* Now the rational side */
	  nr_primes_r = 0;
	  ok = trialdiv_one_side (Gab, scaled_poly_r, 1, a, b, 
				  primes_r, &nr_primes_r, max_nr_primes,
				  proj_divisor_r, nr_proj_primes_r, 
				  proj_primes_r, fbr->fullfb, 
				  reports_r, j,
				  &cof_r_toolarge, &lp_r_toolarge, 
				  poly->rlim, poly->lpbr, poly->mfbr, 
				  poly->rlambda, log_scale);

	  if (!ok) 
	    goto nextreport;

	  /* Now print the relations */
	  relations_found++;
	  fb_sortprimes (primes_r, nr_primes_r);
	  fb_sortprimes (primes_a, nr_primes_a);

	  printf ("%ld,%lu:", a, b);
	  for (k = 0; k < nr_primes_r; k++)
	    printf ("%x%s", primes_r[k], k+1==nr_primes_r?"":",");
	  printf (":");
	  for (k = 0; k < nr_primes_a; k++)
	    printf ("%x%s", primes_a[k], k+1==nr_primes_a?"":",");
	  printf ("\n");
	  fflush (stdout);

	nextreport:
	  /* Skip over duplicates that might cause relations to be
	     trial divided/output repeatedly */
	  while (i + 1 < reports_a->nr && 
		 reports_a->reports[i].a == reports_a->reports[i + 1].a)
	    i++;
	  while (j + 1 < reports_r->nr && 
		 reports_r->reports[j].a == reports_r->reports[j + 1].a)
	    j++;
	}

      /* Assumes values in reports_a are sorted, same for reports_r  */
      if (reports_a->reports[i].a < reports_r->reports[j].a)
	i++;
      else
	j++;
    }

  if (verbose)
    {
      printf ("# Number of matching reports with a,b not coprime: %u\n",
	      not_coprime);
      printf ("# Number of matching reports with a,b coprime: %u\n", 
	      matching_reports);
      printf ("# Number of relations found: %u\n", relations_found);
      printf ("# Sum and number of largest non-report fb primes: %lu, %lu, avg: %.0f\n",
	      sumprimes, nrprimes, (double)sumprimes / (double)nrprimes);
      printf ("# Sum and number of 2nd largest non-report fb primes: %lu, %lu, avg: %.0f\n",
	      sumprimes2, nrprimes2, (double)sumprimes2 / (double)nrprimes2);
      printf ("# Called ul_rho %lu times with %lu repeats, mpz_rho %lu times "
	      "with %lu repeats\n",
	      ul_rho_called, ul_rho_called - ul_rho_called1, 
	      mpz_rho_called, mpz_rho_called - mpz_rho_called1);
      printf ("# %lu cofactors too large for ul_rho and too small for "
	      "mpz_rho, average %.3f", ul_rho_toolarge_nr, 
	      (double) ul_rho_toolarge_sum / (double) ul_rho_toolarge_nr);
      printf ("# Calls to ul_rho() took %llu clocks\n", ul_rho_time);
    }
  
  for (i = 0; i <= (unsigned)poly->degree; i++)
    mpz_clear (scaled_poly_a[i]);
  mpz_clear (scaled_poly_r[0]);
  mpz_clear (scaled_poly_r[1]);
      
  mpz_clear (Fab);
  mpz_clear (Gab);
  
  tsc2 = clock ();
  if (verbose)
  {
    printf ("# Trial factoring/printing%s took %ld clocks\n", 
#ifdef REDC_ROOTS
	    " (with    REDC)",
#else
	    " (without REDC)",
#endif
	    (long int) (tsc2 - tsc1));
    printf ("# Too large cofactors (discarded in this order): "
	    "alg > mfba: %d, alg prp > lpba: %d, "
	    "rat > mfbr: %d, rat prp > lpbr: %d\n", 
	    cof_a_toolarge, lp_a_toolarge, cof_r_toolarge, lp_r_toolarge);
  }

  return relations_found;
}


void rho_timing()
{
  mpz_t m;
  unsigned long n, q;
  unsigned int i;
  clock_t tsc1, tsc2;
  const unsigned int iterations = 10000;

#if (ULONG_MAX > 4294967295UL)
  n = 404428732271UL;
#else
  n = 4292870399UL; /* Too small really for meaningful tests */
#endif
  mpz_init (m);
  mpz_set_ui (m, n);

  tsc1 = clock ();
  for (i = 0; i < 10000; i++)
    q = ul_rho (n, 2UL);

  tsc2 = clock ();
  printf ("%u iteratios of ul_rho took %ld clocks\n", 
	  iterations, (long int) (tsc2 - tsc1));

  tsc1 = clock ();
  for (i = 0; i < 10000; i++)
    q = mpz_rho (m, 2UL);
  tsc2 = clock ();
  printf ("%u iteratios of mpz_rho took %ld clocks\n", 
	  iterations, (long int) (tsc2 - tsc1));

  fflush (stdout);
  mpz_clear (m);
}

static void
usage (void)
{
  fprintf (stderr, "Usage: sieve [-v] [-fb file] [-poly file] [-reports_a_len int] [-reports_r_len int] [-rhotiming] amin amax bmin bmax\n");
  exit (1);
}

int
main (int argc, char **argv)
{
  long amin, amax;
  unsigned long bmin, bmax, b;
  sieve_reportbuffer_t reports_a, reports_r;
  factorbase_t fba, fbr;
  char *fbfilename = NULL, *polyfilename = NULL;
  unsigned char *sievearray;
  unsigned int reports_a_len = 0, reports_r_len = 0;
  int verbose = 0;
  unsigned int deg;
  unsigned int i;
  double dpoly_a[MAXDEGREE], dpoly_r[2];
  const double log_scale = 1. / log (2.); /* Lets use log_2() for a start */
  cado_poly cpoly;
  char report_a_threshold, report_r_threshold;
  unsigned long relations_found = 0;
  double init_time = seconds (), sieve_time;

  while (argc > 1 && argv[1][0] == '-')
    {
      if (argc > 1 && strcmp (argv[1], "-v") == 0)
	{
	  verbose++;
	  argc--;
	  argv++;
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
      else if (argc > 2 && strcmp (argv[1], "-reports_a_len") == 0)
	{
	  reports_a_len = atoi (argv[2]);
	  argc -= 2;
	  argv += 2;
	}
      else if (argc > 2 && strcmp (argv[1], "-reports_r_len") == 0)
	{
	  reports_r_len = atoi (argv[2]);
	  argc -= 2;
	  argv += 2;
	}
      else  if (argc > 1 && strcmp (argv[1], "-rhotiming") == 0)
	{
	  rho_timing();
	  argc--;
	  argv++;
	}
      else
	break;
    }

  if (argc < 5)
    {
      usage ();
      exit (EXIT_FAILURE);
    }


  amin = atol (argv[1]);
  amax = atol (argv[2]);
  bmin = strtoul (argv[3], NULL, 10);
  bmax = strtoul (argv[4], NULL, 10);

  /* Check command line parameters for sanity */

  if (amin >= amax)
    {
      fprintf (stderr, "amin must be less than amax\n");
      exit (EXIT_FAILURE);
    }

  if (bmin > bmax)
    {
      fprintf (stderr, "bmin must be less than or equal to bmax\n");
      exit (EXIT_FAILURE);
    }

  if (polyfilename == NULL)
    {
      fprintf (stderr, "Please specify a polynomial file with -poly\n");
      exit (EXIT_FAILURE);
    }

  if (fbfilename == NULL)
    {
      fprintf (stderr, 
	       "Please specify a factor base file with the -fb option\n");
      exit (EXIT_FAILURE);
    }

  /* Print info about program and parameters */
  if (verbose)
    {
      printf ("# CADO line siever " REV "\n");
      printf ("# Compiled in parameters:\n");
#ifdef WANT_ASSERT
      printf ("# Assertions enabled\n");
#else
      printf ("# Assertions disabled\n");
#endif

#ifdef WANT_ASSERT_EXPENSIVE
      printf ("# Expensive assertions enabled\n");
#else
      printf ("# Expensive assertions disabled\n");
#endif

#if SIEVE_BLOCKING > 0
      printf ("# Sieving small factor base primes in blocks of size%s ",
	      (SIEVE_BLOCKING > 1) ? "s" : "");
      printf ("%lu", CACHESIZES[0]);
      for (i = 1; i < SIEVE_BLOCKING; i++)
	printf (", %lu", CACHESIZES[i]);
      printf (".\n");
#else
      printf ("# No block sieving of small factor basee primes.\n");
#endif

#if TRIALDIV_SKIPFORWARD
      printf ("# Can skip forward in factor base during refactoring");
#else
      printf ("# Cannot skip forward in factor base during refactoring");
#endif
#ifdef TRIALDIV_SKIPFORWARD_DEFAULT
      printf (" (using default setting).\n");
#else
      printf (".\n");
#endif

#ifdef TRACE_RELATION_A
      printf ("# Tracing the fate of relations with a=%ld.\n", 
	      (long) TRACE_RELATION_A);
#else
      printf ("# Relation tracing disabled.\n");
#endif

    }

  /* Read polynomial from file */
  if (!read_polynomial (cpoly, polyfilename))
    {
      fprintf (stderr, "Error reading polynomial file\n");
      exit (EXIT_FAILURE);
    }
  if (verbose)
  {
      printf ("# Read polynomial file %s\n", polyfilename);
      printf ("# Polynomials are:\n");
      mp_poly_print (cpoly->f, cpoly->degree, "# f(x) =", 0);
      printf ("\n");
      mp_poly_print (cpoly->g, 1, "# g(x) =", 0);
      printf ("\n");
  }

#ifdef PARI
  mp_poly_print (cpoly->f, cpoly->degree, "f(x) =");
  printf (" /* PARI */\n");
  printf ("F(a,b) = f(a/b)*b^%d /* PARI */\n", cpoly->degree);
  mp_poly_print (cpoly->g, 1, "g(x) =");
  printf (" /* PARI */\n");
  printf ("G(a,b) = g(a/b)*b /* PARI */\n");
#endif


  /* Read/generate and split the factor bases */

  /* Read the factor base for the algebraic side from file */
  fba->fullfb = fb_read (fbfilename, log_scale, verbose);
  if (fba->fullfb == NULL)
    {
      fprintf (stderr, "Could not read factor base\n");
      exit (EXIT_FAILURE);
    }
  fba->fblarge = fba->fullfb;
  /* Extract the small primes into their small prime arrays */
  for (i = 0; i < SIEVE_BLOCKING; i++)
    fb_extract_small (fba, CACHESIZES[i], i, verbose);


  /* Generate rational fb */
  fbr->fullfb = fb_make_linear (cpoly->g, (fbprime_t) cpoly->rlim, 
				log_scale, verbose);
  if (fbr == NULL)
    {
      fprintf (stderr, 
	       "Could not generate factor base for linear polynomial\n");
      exit (EXIT_FAILURE);
    }
  fbr->fblarge = fbr->fullfb;
  /* Extract the small primes into their small prime arrays */
  for (i = 0; i < SIEVE_BLOCKING; i++)
    fb_extract_small (fbr, CACHESIZES[i], i, verbose);
  
#ifdef WANT_ASSERT
  /* Check the factor bases for correctness */
  
  if (fb_check (fba, cpoly, 0) != 0)
    {
      fprintf (stderr, "Error in algebraic factor base\n");
      exit (EXIT_FAILURE);
    }
  else if (verbose)
    printf ("# Algebraic factor base test passed\n");

  if (fb_check (fbr, cpoly, 1) != 0)
    {
      fprintf (stderr, "Error in ratinal factor base\n");
      exit (EXIT_FAILURE);
    }
  else if (verbose)
    printf ("# Rational factor base test passed\n");
#endif

  deg = cpoly->degree;
  for (i = 0; i <= deg; i++)
      dpoly_a[i] = mpz_get_d (cpoly->f[i]);
  dpoly_r[0] = mpz_get_d (cpoly->g[0]);
  dpoly_r[1] = mpz_get_d (cpoly->g[1]);
  report_a_threshold = (unsigned char) ((double)(cpoly->lpba) * log(2.0) * 
					cpoly->alambda * log_scale + 0.5);
  report_r_threshold = (unsigned char) ((double)(cpoly->lpbr) * log(2.0) * 
					cpoly->rlambda * log_scale + 0.5);
  if (verbose)
    {
      printf ("# Report threshold for algebraic side: %d, rational side: %d\n",
	      (int) report_a_threshold, (int) report_r_threshold);
    }

#if 0
  if (verbose)
    for (s = 0; fb_skip(fb, s)->p != 0; s += fb_entrysize (fb_skip (fb, s)))
      fb_print_entry (fb_skip(fb, s));
#endif

  if (verbose)
    {
      printf ("# Allocating %ld bytes for sievearray\n",
	      (amax - amin + 1) * (long int) sizeof (char));
      fflush(stdout);
    }
  sievearray = (unsigned char *) malloc ((amax - amin + 1) * sizeof (char));

  /* Allocate report buffers */

  if (reports_a_len == 0) /* if not given on command line */
    reports_a_len = ((amax - amin + 1)) / 256 + 1000;
  if (verbose)
    {
      printf ("# Allocating %lu entries (%lu bytes) for reports_a\n",
	      (unsigned long) reports_a_len, (unsigned long) 
	      (reports_a_len * (long int) sizeof (sieve_report_t)));
      fflush(stdout);
    }
  reports_a.alloc = reports_a_len;
  reports_a.nr = 0;
  reports_a.reports = (sieve_report_t *) malloc (reports_a_len * 
						 sizeof (sieve_report_t));
  ASSERT_ALWAYS (reports_a.reports != NULL);
  if (verbose)
    {
      printf ("# Allocating %lu entries (%lu bytes) for reports_r\n",
	      (unsigned long) reports_r_len, 
	      (unsigned long) (reports_r_len * sizeof (sieve_report_t)));
      fflush(stdout);
    }

  if (reports_r_len == 0) /* if not given on command line */
    reports_r_len = ((amax - amin + 1)) / 256 + 1000;
  reports_r.alloc = reports_r_len;
  reports_r.nr = 0;
  reports_r.reports = (sieve_report_t *) malloc (reports_r_len * 
						 sizeof (sieve_report_t));
  ASSERT_ALWAYS (reports_r.reports != NULL);

  /* Do the sieving */

  sieve_time = seconds ();
  init_time = sieve_time - init_time;

  /* we sieve over bmin <= b < bmax, to avoid duplicates when sieving b0..b1,
     then b1..b2, and to figure out the initialization time */
  for (b = bmin; b < bmax; b++)
    {
      unsigned long proj_roots; /* = gcd (b, leading coefficient) */
      
      if (verbose)
	printf ("# Sieving line b = %lu\n", b);

      proj_roots = mpz_gcd_ui (NULL, cpoly->f[deg], b);
      if (verbose)
	printf ("# Primes with projective roots for b = %lu on algebtaic side "
		"divide %lu\n", b, proj_roots);

      if (verbose)
	  printf ("# Sieving algebraic side\n");

#ifdef PARI
      mp_poly_print (cpoly->f, cpoly->degree, "P(a,b) = ", 1);
      printf (" /* PARI */\n");
#endif
      
      sieve_one_side (sievearray, fba, &reports_a, report_a_threshold, 
		      amin, amax, b, proj_roots, log_scale, 
		      dpoly_a, deg, NULL, verbose);

      if (verbose)
	printf ("# There were %lu sieve reports on the algebraic side\n",
		(unsigned long) reports_a.nr);

      proj_roots = mpz_gcd_ui (NULL, cpoly->g[1], b);
      if (verbose)
	printf ("# Primes with projective roots for b = %lu on rational side "
		"divide %lu\n", b, proj_roots);

      if (verbose)
	printf ("# Sieving rational side\n");

#ifdef PARI
      mp_poly_print (cpoly->g, 1, "P(a,b) = ", 1);
      printf (" /* PARI */\n");
#endif

      sieve_one_side (sievearray, fbr, &reports_r, report_r_threshold, 
		      amin, amax, b, proj_roots, log_scale, 
		      dpoly_r, 1, &reports_a, verbose);
      if (verbose)
	printf ("# There were %lu sieve reports on the rational side\n",
		(unsigned long) reports_r.nr);

      /* Not trial factor the candidate relations */
      relations_found += trialdiv_and_print (cpoly, b, &reports_a, &reports_r, 
					     fba, fbr, log_scale, verbose);
    }

  sieve_time = seconds () - sieve_time;
  fprintf (stderr, "Found %lu relations in %1.0f seconds (%1.2e s/r)",
           relations_found, sieve_time, sieve_time / (double) relations_found);
  fprintf (stderr, " [init time = %1.0f seconds]\n", init_time);

  clear_polynomial (cpoly);
  fb_clear (fba);
  fb_clear (fbr);

  free (reports_a.reports);
  free (reports_r.reports);
  free (sievearray);

  return 0;
}
