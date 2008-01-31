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
  mpz_set_ui (t, 1);
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

static inline unsigned char
uc_add (unsigned char a, unsigned char b)
{
  return (unsigned char) (a + b);
}

static inline unsigned char
uc_sub (unsigned char a, unsigned char b)
{
  return (unsigned char) (a - b);
}

static inline unsigned char 
add_error (unsigned char n)
{
  return uc_add (n, SIEVE_PERMISSIBLE_ERROR);
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
  /* printf ("Norm at x = %.0f is %.0f, rounded log is %hhd\n", x, r, n); */
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
      n2 = fb_log (n, log_scale, log_proj_roots);
      fsign2 = (n < 0.) ? -1 : 1;
      dfsign2 = (fpoly_eval (df, deg - 1, (double) a2) < 0.) ? -1 : 1;

      if (n1 == n2 && fsign1 == fsign2 && dfsign1 == dfsign2)
	{
	  /* Let's assume the log norm is n1 everywhere in this interval */
#ifdef TRACE_RELATION_A
	  if (a <= TRACE_RELATION_A && TRACE_RELATION_A < a2)
	    printf ("# TRACE RELATION a = %ld in compute_norms: log norm "
		    "(range computed) for degree %d is %d\n", 
		    (long) TRACE_RELATION_A, deg, (int) n1);
#endif
	  memset (sievearray + ((a - amin) >> odd), n1, (a2 - a) >> odd);
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
       const unsigned char threshold, sieve_report_t *reports,
       unsigned int reports_length, const int odd)
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

	      TRACE_A (amin + ((long)d << odd), __func__, __LINE__, 
		       "new approx norm after sieving out " FBPRIME_FORMAT 
		       " is %d\n", p, (int) k);

	      sievearray[d] = k;
	      if (add_error (k) <= threshold_with_error)
		{
		  if (reports_length > 0)
		    {
		      TRACE_A ((amin + ((long)d << odd)), __func__, __LINE__,
			       "writing sieve report a = %ld, p = %u, l = %d\n",
			       (amin + ((long)d << odd)), p, (int) k);
		      
		      reports->a = amin + ((long)d << odd);
		      reports->p = p;
		      reports->l = k;
		      reports++;
		      reports_length--;
		    }
		  else
		    {
		      TRACE_A ((amin + ((long)d << odd)), __func__, __LINE__,
			       "Cannot writing sieve report a = %ld, p = %u, "
			       "l = %d, report buffer full\n",
			       (amin + ((long)d << odd)), p, (int) k);
		    }
		}
	    }
        }
      
      /* Move on to the next factor base prime */
      fb = fb_next (fb);
    }

  if (reports_length > 0)
    reports->p = 0; /* Put end marker */
}


static unsigned long
find_sieve_reports (const unsigned char *sievearray, sieve_report_t *reports, 
                    const unsigned int reports_len, 
                    const unsigned char threshold, 
                    const long amin, const unsigned long l, 
		    const unsigned long b, const int odd,
		    sieve_report_t *other_reports)
{
  unsigned long reports_nr, d;
  const unsigned char threshold_with_error = add_error (threshold);
  
  ASSERT (odd == 0 || odd == 1);
  ASSERT (!odd || (amin & 1) == 1);
  ASSERT (b > 0);

  reports_nr = 0;
  if (reports_nr == reports_len)
    return reports_nr;

  for (d = 0; d < l; d++)
    {
      if (add_error(sievearray[d]) <= threshold_with_error)
	{
	  const long a = amin + ((long)d << odd);
	  
	  /* Testing n%3==0 is quite fast, and eliminates useless reports */
	  if (a % 3 == 0 && b % 3 == 0)
	    continue;

	  /* If other_reports was passed, store this report only if there
	     is a matching report in other_reports */
	  if (other_reports != NULL)
	    {
	      while (other_reports->p != 0 && other_reports->a < a)
		other_reports++;
	      if (other_reports->a != a)
		{
		  TRACE_A (a, __func__, __LINE__,
			   "no matching report in other_reports\n");
		  
		  continue;
		}
	      else
		{
		  TRACE_A (a, __func__, __LINE__, 
			   "matching report in other_reports\n");
		}
	    }

	  reports[reports_nr].a = a;
	  reports[reports_nr].p = 1;
	  reports[reports_nr].l = sievearray[d];
	  reports_nr++;
	  TRACE_A (a, __func__, __LINE__,
		   "sieve report a = %ld, p = 1, l = %d added to list\n", 
		   a, (int) sievearray[d]);

	  if (reports_nr == reports_len)
	    {
	      fprintf (stderr, "Warning: find_sieve_reports: reports list "
		       "full at a = %ld\n", a);
	      return reports_nr;
	    }
	}
    }
  
  return reports_nr;
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
sort_sieve_reports (sieve_report_t *r, const unsigned long l)
{
  long p; /* pivot element */
  unsigned long b, h;
  const unsigned long m = (l - 1) / 2; /* The midpoint */

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
      while (r[b].a <= p && b <= h)
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

  sort_sieve_reports (r, b);
  sort_sieve_reports (r + h, l - h);

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
		sieve_report_t *reports, const unsigned int reports_length,
		const unsigned char reports_threshold,
		const long amin, const long amax, const unsigned long b,
		const unsigned long proj_roots, const double log_scale,
		const double *dpoly, const unsigned int deg, 
		sieve_report_t *other_reports, const int verbose)
{
  clock_t tsc1, tsc2;
  long long times[SIEVE_BLOCKING];
  unsigned long reports_nr = 0;
  int lvl;
  const int odd = 1 - (b & 1); /* If odd=1, only odd $a$ are sieved */
  const long eff_amin = amin + ((odd && (amin & 1) == 0) ? 1 : 0);
  const long eff_amax = amax - ((odd && (amax & 1) == 0) ? 1 : 0);
  const unsigned long l = ((eff_amax - eff_amin) >> odd) + 1;
  const int find_smallprime_reports = 1;

  fb_disable_roots (fb->fblarge, b, verbose);
  
  compute_norms (sievearray, eff_amin, eff_amax, b, dpoly, deg, proj_roots, 
		 log_scale, odd, verbose);

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

  /* If the factor base limit is very small compared to the 
     largest block length, we may get many relations that are 
     smooth over the small fb primes and will not produce sieve 
     reports during sieving the large fb primes.
     Those reports are found here instead. */
  if (find_smallprime_reports)
    {
      tsc1 = clock ();
      reports_nr +=
	find_sieve_reports (sievearray, 
			    reports + reports_nr, 
			    reports_length - reports_nr, 
			    reports_threshold, 
			    eff_amin, l, b, odd, other_reports);
      tsc2 = clock ();
      if (verbose)
	{
	  printf ("# Finding sieve reports took %ld clocks\n", 
		  (long int) (tsc2 - tsc1));
	  printf ("# There were %lu sieve reports after sieving small "
		  "primes\n", reports_nr);
	}
    }
  
  tsc1 = clock ();
  sieve (sievearray, fb->fblarge, eff_amin, eff_amax, b, reports_threshold, 
	 reports + reports_nr, reports_length - reports_nr, odd);
  tsc2 = clock ();
  if (verbose)
    printf ("# Sieving large fb primes took %ld clocks\n", 
	    (long int) (tsc2 - tsc1));
  
  fb_restore_roots (fb->fblarge, b, verbose);
  
  /* Count the sieve reports and sort by increasing a */
  for (reports_nr = 0; reports_nr < reports_length; reports_nr++)
    {
      if (reports[reports_nr].p == 0)
	break;
    }

  sort_sieve_reports (reports, reports_nr);
  
  return reports_nr;
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

  mod_initmod_ul (m, fbptr->p);
  mod_init_noset0 (r, m);

  absa = labs(a);

  for (i = 0; i < fbptr->nr_roots ; i++)
    {
      if (a <= 0) /* We compute rb - a */
	mod_set_ul_reduced (r, fbptr->roots[i], m);
      else /* We compute (-r)b + a */
	mod_set_ul_reduced (r, fbptr->p - fbptr->roots[i], m);

      mod_muladdredc_ul (r, r, b, absa, fbptr->invp, m);
      
      if (mod_get_ul (r, m) == 0)
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
  unsigned int i;

  ASSERT (mpz_sgn (norm) >= 0);

  mod_initmod_ul (m, fbptr->p);
  mod_init (r, m);
  
  for (i = 0; i < mpz_size (norm); i++)
    modul_addredc_ul (r, r, mpz_getlimbn (norm, i), fbptr->invp, m);

  i = (mod_get_ul (r, m) == 0);

  mod_clear (r, m);
  mod_clearmod (m);
  return i;
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

  mod_initmod_ul (m, fbptr->p);
  mod_init_noset0 (r, m);

  modul_redcsemi_ul_not0 (r, mpz_getlimbn (norm, 0), fbptr->invp, m);

  i = (mod_get_ul (r, m) + add == 0 || mod_get_ul (r, m) + add == m[0]);

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

  mod_initmod_ul (m, fbptr->p);
  mod_init_noset0 (r, m);

  modul_redcsemi_ul_not0 (r, mpz_getlimbn (norm, 0), fbptr->invp, m);
  modul_addredcsemi_ul (r, r, mpz_getlimbn (norm, 1), fbptr->invp, m);

  i = (mod_get_ul (r, m) + add == 0 || mod_get_ul (r, m) + add == m[0]);

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

  mod_initmod_ul (m, fbptr->p);
  mod_init_noset0 (r, m);

  modul_redcsemi_ul_not0 (r, mpz_getlimbn (norm, 0), fbptr->invp, m);
  modul_addredcsemi_ul (r, r, mpz_getlimbn (norm, 1), fbptr->invp, m);
  modul_addredcsemi_ul (r, r, mpz_getlimbn (norm, 2), fbptr->invp, m);

  i = (mod_get_ul (r, m) + add == 0 || mod_get_ul (r, m) + add == m[0]);

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

  mod_initmod_ul (m, fbptr->p);
  mod_init_noset0 (r, m);

  modul_redcsemi_ul_not0 (r, mpz_getlimbn (norm, 0), fbptr->invp, m);
  modul_addredcsemi_ul (r, r, mpz_getlimbn (norm, 1), fbptr->invp, m);
  modul_addredcsemi_ul (r, r, mpz_getlimbn (norm, 2), fbptr->invp, m);
  modul_addredcsemi_ul (r, r, mpz_getlimbn (norm, 3), fbptr->invp, m);

  i = (mod_get_ul (r, m) + add == 0 || mod_get_ul (r, m) + add == m[0]);

  mod_clear (r, m);
  mod_clearmod (m);
  return i;
}


/* Square root. Returns the largest integer x so that x^2 <= n.
   If e != NULL, stores n - x^2 in *e. Should be compiled with 
   -funroll-loops for best performance. */

static inline unsigned long  
ul_sqrtint (const unsigned long n, unsigned long *e)
{
  int i;   
  unsigned long xs, c, d, s2;

  d = n; /* d = n - x^2 */
  xs = 0UL;
  s2 = 1UL << (sizeof (unsigned long) * 8 - 2);

  for (i = sizeof (unsigned long) * 4 - 1; i != 0; i--)
    {
      /* Here, s2 = 1 << (2*i) */
      /* xs = x << (i + 1), the value of x shifted left i+1 bits */

      c = xs + s2; /* c = (x + 2^i) ^ 2 - x^2 = 2^(i+1) * x + 2^(2*i) */
      xs >>= 1; /* Now xs is shifted only i positions */
      if (d >= c)
        {
          d -= c;
          xs |= s2; /* x |= 1UL << i <=> xs |= 1UL << (2*i) */
        }
      s2 >>= 2;
    }

  c = xs + s2; 
  xs >>= 1;
  if (d >= c)
    {
      d -= c;   
      xs |= s2;
    }
 
  if (e != NULL)
    *e = d;
  return xs;
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
  modul_set_ul_reduced (r1, 2, m);
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
  mpz_set_ui (r1, 2);
  mpz_set_ui (r2, 2);
  mpz_set_ui (accu, 2);

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
  } while (mpz_cmp_ui (t, 1) == 0);

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

  if (n < 2047)
    {
      r = 1;
      goto end;
    }

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
  
 end:
#if defined(PARI)
  printf ("isprime(%lu) == %d /* PARI */\n", n, r);
#endif
  modul_clearmod (m);
  return r;
}


/* Returns a pointer to factor base entry which is either the 
   end-of-factorbase marker, or the first entry whose p divides norm.
   Requires that norm > fbptr->p */

static factorbase_degn_t *
trialdiv_find_next (factorbase_degn_t *fbptr, const mpz_t norm)
{
  size_t s = mpz_size (norm);
  unsigned long add = 0;

  ASSERT (mpz_cmp_ui (norm, fbptr->p) >= 0);

  if (mpz_getlimbn (norm, s - 1) < fbptr -> p)
    {
      add = mpz_getlimbn (norm, s - 1);
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
   3. Divide out the report prime(s), remembering the smallest one and remember
      the smallest approximate log from the reports
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
		   const sieve_report_t *reports,
		   unsigned int *cof_toolarge,
		   unsigned int *lp_toolarge,
		   const unsigned long fbb, const int lpb, const int mfb, 
		   const double lambda, const double log_scale)
{
  unsigned int k;
  factorbase_degn_t *fbptr;
  size_t log2size;
  unsigned char reportlog, missinglog;
  const double log_proj_divisor = log ((double) proj_divisor) * log_scale;
  double dpoly[MAXDEGREE + 1];
  int nr_report_primes;

  /* 1. Compute norm */
  mp_poly_eval (norm, scaled_poly, degree, a);
  mpz_abs (norm, norm);
  TRACE_A (a, __func__, __LINE__, "norm for degree %d is %Zd\n", degree, norm);
  
  /* Compute approx. log of norm as was used for initialising the 
     sieve array */
  for (k = 0; (int) k <= degree; k++)
    dpoly[k] = mpz_get_d (scaled_poly[k]);
  missinglog = log_norm (dpoly, degree, (double) a, log_scale, 
			 log_proj_divisor);
  TRACE_A (a, __func__, __LINE__, "missinglog = %hhu\n", missinglog);
  
  
  /* 2. Divide out the primes with projective roots */
  if (proj_divisor > 1)
    {
      unsigned long r;
      ASSERT_ALWAYS (nr_proj_primes > 0);
      r = mpz_tdiv_q_ui (norm, norm, proj_divisor);
      ASSERT_ALWAYS (r == 0);
      
      TRACE_A (a, __func__, __LINE__, "dividing out proj. divisor %lu. "
	       "New norm is %Zd.\n", proj_divisor, norm);

      for (k = 0; k < nr_proj_primes; k++)
	{
	  ASSERT (proj_divisor % proj_primes[k] == 0);
	  add_fbprime_to_list (primes, nr_primes, max_nr_primes, 
			       proj_primes[k]);
	  while (mpz_divisible_ui_p (norm, proj_primes[k]))
	    {
	      mpz_tdiv_q_ui (norm, norm, proj_primes[k]);
	      add_fbprime_to_list (primes, nr_primes, max_nr_primes, 
				   proj_primes[k]);
	      TRACE_A (a, __func__, __LINE__, "dividing out extra power of "
		       "proj. ""prime %lu. New norm = %Zd\n", 
		       proj_primes[k], norm);
	    }
	}
    }
  
  
  /* 3. Divide the report primes out of this norm and find the smallest 
     approximate log */
  /* There must be at least one valid report */
  ASSERT_ALWAYS (reports->a == a && reports->p != 0);
  reportlog = reports->l;
  nr_report_primes = 0;
  while (reports->a == a && reports->p != 0)
    {
      /* We allow reports without a report prime. In that case the 
	 sieve_report_t may contain p == 1. */
      if (reports->p != 1)
	{
	  int r;

	  r = trialdiv_one_prime (reports->p, norm, nr_primes, 
				  primes, max_nr_primes, 0);
	  ASSERT_ALWAYS (r > 0); /* If a report prime is listed but does not
				    divide the norm, there is a serious bug
				    in the sieve code */
	  missinglog -= fb_log ((double) reports->p, log_scale, 0.);
	  nr_report_primes++;

	  TRACE_A (a, __func__, __LINE__, "dividing out report prime " 
		   FBPRIME_FORMAT ". New new norm = %Zd, missinglog = %hhd\n", 
		   reports->p, norm, missinglog);
	}

      /* Find the smallest approximate log. */
      if (add_error (reports->l) < add_error (reportlog))
	reportlog = reports->l;
	  
      reports++;
    }
  missinglog -= reportlog;
  TRACE_A (a, __func__, __LINE__, "reportlog is %hhd, new missinglog = %hhd\n", 
	   reportlog, missinglog);

  /* 
     When we divide out a prime, we subtract log(p) from missinglog. When 
     missinglog == 0, we know we are done. When missinglog < 2*log(p), 
     we know that exactly one more factor base prime divides.
  */

  /* Treat factor base prime p == 2 separately. We want to be able to use
     REDC in what follows and it requires an odd modulus */
  fbptr = fullfb;
  TRACE_A (a, __func__, __LINE__, "first fb prime is " FBPRIME_FORMAT "\n", 
	   fbptr->p);
  if (fbptr->p == 2)
    {
      if (mpz_even_p (norm))
	missinglog -= fbptr->plog;
      while (mpz_even_p (norm))
	{
	  mpz_tdiv_q_2exp (norm, norm, 1);
	  add_fbprime_to_list (primes, nr_primes, max_nr_primes, 2);
	  TRACE_A (a, __func__, __LINE__, "dividing out fb prime 2." 
		   "New norm = %Zd, missinglog = %hhd\n", norm, missinglog);
	}
      fbptr = fb_next (fbptr);
    }
  else
    {
      ASSERT_ALWAYS (!mpz_even_p (norm));
    }
  
  /* 4. Go through the factor base until missinglog == 0 */

  while (missinglog != 0)
    {
      fbptr = trialdiv_find_next (fbptr, norm);

      if (fbptr->p == 0)
	{
	  fprintf (stderr, "Warning, reached end of fb for (%ld, %lu). "
		   "missinglog = %d\n", a, b, (int) missinglog);
	  gmp_fprintf (stderr, "Remaining norm = %Zd\n", norm);
	  break;
	}
      
      ASSERT (mpz_divisible_ui_p (norm, fbptr->p));
      trialdiv_one_prime (fbptr->p, norm, nr_primes, primes,
			  max_nr_primes, 1);

      ASSERT_ALWAYS (missinglog >= fbptr->plog);
      missinglog -= fbptr->plog;
      
      TRACE_A (a, __func__, __LINE__, "dividing out fb prime " FBPRIME_FORMAT
	       ". New norm = %Zd, missinglog = %hhd\n", 
	       fbptr->p, norm, missinglog);
      
      /* If we aren't done yet, there must remain at least one fb prime
	 greater than the one we just processed. Check that the
	 difference of approximate logs allows for that prime. */
      if (missinglog != 0 && missinglog < fbptr->plog)
	{
	  gmp_fprintf (stderr, "Error, a = %ld, b = %lu, missinglog = %d, "
		       "p = " FBPRIME_FORMAT ", plog = %d. "
		       "Remaining norm = %Zd\n", 
		       a, b, missinglog, fbptr->p, (int) fbptr->plog, norm);
	  abort ();
	}

      if (missinglog < 2 * fbptr->plog)
	{
	  /* We know there's exactly one prime (possibly repeated) 
	     missing, and we know its approximate size. We can choose 
	     a more efficient algorithm here to find it. */
	  
	  nrprimes2++;
	  sumprimes2 += (unsigned long) fbptr->p;
	  
	  if (mpz_fits_ulong_p (norm))
	    {
	      /* We can use the fast ul_rho() function */
	      if (missinglog > UL_RHO_LOGTHRES)
		break;
	    }
	  else
	    {
	      /* We need to use the slower mpz_rho() function */
	      if (missinglog > MPZ_RHO_LOGTHRES)
		break;
	    }
	  
	  /* If we can't use Pollard rho, let's try skipping forward in 
	     the factor base to those primes that have the correct 
	     rounded log. Keep in mind that the for() loop advances 
	     fbptr once after this while loop exits! */
#if TRIALDIV_SKIPFORWARD
	  while (fb_next (fbptr)->p != 0 && 
		 fb_next (fbptr)->plog < missinglog)
	    fbptr = fb_next (fbptr);
	  TRACE_A (a, __func__, __LINE__, "skipping forward in fb, next prime "
		   "is " FBPRIME_FORMAT "\n", fb_next (fbptr)->p);
#endif
	}
      fbptr = fb_next (fbptr);
    }

  if (fbptr->p != 0)
    {
      sumprimes += (unsigned long) (fbptr->p);
      nrprimes++;
    }

  if (missinglog != 0)
    {
      /* There is one more factor base prime (possibly in a power) that we 
	 want to find with a more efficient algorithm. */

      unsigned long q, r, s = 5;

      ASSERT_ALWAYS (missinglog >= fbptr->plog && 
		     missinglog < 2 * fbptr->plog);


      /* While in this loop, norm contains a factor base prime (possibly
	 in a power > 1) whose log norm = missinglog. 
	 That is also the smallest prime factor in norm. */
      do
	{
	  ASSERT_ALWAYS (mpz_cmp_ui (norm, 1) > 0);
	  
	  /* Check if norm is itself a prime (power), i.e. if there are no 
	     large primes present. */
	  if (mpz_perfect_power_p (norm))
	    {
	      /* Will rarely happen, need not be that fast */
	      mpz_t t;
	      unsigned long i;

	      mpz_init (t);
	      for (i = 2UL; !mpz_root (t, norm, i); i++);
	      ASSERT_ALWAYS (mpz_fits_ulong_p (t));
	      q = mpz_get_ui (t);
	      mpz_clear (t);
	      
#ifdef RHODEBUG
	      gmp_printf ("# %Zd is an %lu-th power of %lu for a = %ld, "
			  "b = %lu, no rho necessary\n", norm, i, q, a, b);
	      fflush (stdout);
#endif

	      /* norm may be a power of a composite value! I.e. if the large
		 prime also appears in the same power. If q is composite, 
		 find the smallest prime dividing it. This will happen very 
		 rarely. */
	      if (!ul_proven_prime (q))
		q = iscomposite (q);

	      TRACE_A (a, __func__, __LINE__, "norm %Zd is a power, %lu "
		       "divides it\n", norm, q);
	      ASSERT_ALWAYS (q != 0);
	      trialdiv_one_prime (q, norm, nr_primes, primes, max_nr_primes, 
				  1);
	      /* q may not have been the only prime that divides norm.
	         If it isn't, we simply loop again */
	      if (mpz_cmp_ui (norm, 1UL) == 0)
		break;
	    }
	  else if (mpz_cmp_ui (norm, fbb) <= 0)
	    {
	      /* Smaller than factor base bound and not a prime power, 
		 so it should be prime. */
#ifdef RHODEBUG
	      gmp_printf ("# %Zd is prime for a = %ld, b = %lu, no rho "
			  "necessary\n", norm, a, b);
	      fflush (stdout);
#endif
	      TRACE_A (a, __func__, __LINE__, "norm %Zd is < fbb %lu, "
		       "assuming it is prime and adding to list of primes\n", 
		       norm, fbb);
	      q = mpz_get_ui (norm);
	      ASSERT (ul_proven_prime (q));
	      trialdiv_one_prime (q, norm, nr_primes, primes, max_nr_primes, 
				  1);
	      break;
	    }
	  else
	    {
	      /* Not <fbb or prime power. Since we know there is a fb prime
		 in norm, not being < fbb implies being composite. 
		 Try to factor it */
	      TRACE_A (a, __func__, __LINE__, "Trying Pollard rho on norm "
		       "%Zd\n", norm);
	      if (mpz_fits_ulong_p (norm))
		{
		  unsigned long n;
		  n = mpz_get_ui (norm);
		  ul_rho_called1++;
		  do {
		    q = ul_rho (n, s++);
		    ul_rho_called++;
		  } while (q == n);
		  ASSERT (n % q == 0);
		}
	      else
		{
		  mpz_rho_called1++;
		  do {
		    q = mpz_rho (norm, s++);
		    mpz_rho_called++;
		  } while (q == 1);
		  ASSERT (mpz_divisible_ui_p (norm, q));
		}
	      TRACE_A (a, __func__, __LINE__, "Pollard rho found %lu\n", q);

	      /* Pollard rho has a way of discovering powers of primes if
		 they divide the input number. Let's at least check for a
		 square here. */
	      {
		unsigned long e, t;
		t = ul_sqrtint (q, &e);
		if (e == 0)
		  q = t;
	      }

	      /* So we have a divisor of norm, q. See if it is the factor
	         base prime we were looking for. */

	      if (q < fbb) /* It must be the missing factor base prime,
			      maybe a power of it divides norm */
		{
		  ASSERT (ul_proven_prime (q));
		  add_fbprime_to_list (primes, nr_primes, max_nr_primes, q);
		  r = mpz_tdiv_q_ui (norm, norm, q);
		  ASSERT_ALWAYS (r == 0);
		  /* Maybe a power of it divides norm */
		  trialdiv_one_prime (q, norm, nr_primes, primes, 
				      max_nr_primes, 1);
		  break;
		}

	      /* So it's not the factor base prime, but possibly a large 
		 prime. If it is, divide it out and try again */
	      if (ul_proven_prime (q))
		{
		  if (q > (1UL << lpb))
		    return 0;

		  add_fbprime_to_list (primes, nr_primes, max_nr_primes, q);
		  r = mpz_tdiv_q_ui (norm, norm, q);
		  ASSERT_ALWAYS (r == 0);
		}
	      else
		{
		  /* Composite factor that is not a power of the factor 
		     base prime. Let's do it the hard way, should happen 
		     rarely enough. */
#ifdef RHODEBUG
		    fprintf (stderr, 
			     "# Composite factor %lu was found by rho, "
			     "factoring it slowly\n", q);
#endif
		  q = iscomposite (q);
		  ASSERT_ALWAYS (q != 0);
		  trialdiv_one_prime (q, norm, nr_primes, primes, 
				      max_nr_primes, 1);
		  if (q < fbb)
		    break;
		}
	    }
	} while (1);

      if (fb_log (q, log_scale, 0.) != missinglog)
	{
	  /* q was not the missing factor base prime ? */
	  fprintf (stderr, "Warning, expected to find fb prime of log size "
		   "%hhu, but found %lu. a = %ld, b = %lu\n", 
		   missinglog, q, a, b);
	}
    }

  /* 5. Check if the cofactor is small enough */
  log2size = mpz_sizeinbase (norm, 2);
  TRACE_A (a, __func__, __LINE__, "log2size of cofactor %Zd is %d\n", 
	   norm, log2size);

  if ((double) log2size > 
      lpb * lambda + SIEVE_PERMISSIBLE_ERROR)
    {
      gmp_fprintf (stderr, 
		   "Sieve report (%ld, %lu) is not smooth for degree %d poly, "
		   "cofactor is %Zd with %d bits\n", 
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
  
  /* 6. Check if the cofactor is < lpb */
  if (mpz_cmp_ui (norm, 1UL) == 0)
    {
      TRACE_A (a, __func__, __LINE__, "Cofactor is 1. This side is smooth.\n");
      return 1;    
    }


  if ((int) log2size <= lpb)
    {
      fbprime_t q;

      if (mpz_probab_prime_p (norm, PRP_REPS))
	{
	  /* If this cofactor is a prp, add it to the list of primes */
	  q = (fbprime_t) mpz_get_ui (norm);
	  add_fbprime_to_list (primes, nr_primes, max_nr_primes, q);
	  TRACE_A (a, __func__, __LINE__, "cofactor" FBPRIME_FORMAT 
		   "is <= 2^lpb and prp. This side is smooth.\n", q);
	  mpz_set_ui (norm, 1UL);
	  return 1;
	}
      else
	{
	  /* If not prp, just ignore it, print the relation and let the 
	     next program in the tool chain figure out the prime factors 
	     of this composite */
	  TRACE_A (a, __func__, __LINE__, "cofactor %Zd is <= 2^lpb and "
		   "composite. This side is smooth.\n", norm);

	  /* It is a bit odd though to have a composite factor <lpb.
	     Print a warning about it. */
	  gmp_fprintf (stderr, "Warning: cofactor %Zd <lpb, but not "
		       "prime for (%ld, %lu), degree %d polynomial\n", 
		       norm, a, b, degree);
	  return 1;
	}
    }
  else
    {
      /* If this cofactor is a prp, since it's > lpb, 
	 skip this report */
      if (mpz_probab_prime_p (norm, PRP_REPS))
	{
	  TRACE_A (a, __func__, __LINE__, "cofactor %Zd is > 2^lpb and prp. "
		   "Discarding relation.\n", norm);
	  (*lp_toolarge)++;
	  return 0;
	}
      else
	{
	  TRACE_A (a, __func__, __LINE__, "cofactor %Zd is > 2^lpb and "
		   "composite. This side may be smooth.\n", norm);
	  return 1;
	}
    }
}

/* return the number of printed relations */
static unsigned int
trialdiv_and_print (cado_poly poly, const unsigned long b, 
                    const sieve_report_t *reports_a, 
		    const unsigned int reports_a_nr, 
		    const sieve_report_t *reports_r, 
		    const unsigned int reports_r_nr, 
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
    cof_a_toolarge = 0, cof_r_toolarge = 0;
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

  for (i = 0, j = 0; i < reports_a_nr && j < reports_r_nr;)
    {
      if (reports_a[i].a == reports_r[j].a && 
	  gcd(labs(reports_a[i].a), b) == 1)
	{
          const long a = reports_a[i].a;

	  matching_reports++;
      
          /* Do the algebraic side */
	  nr_primes_a = 0;
	  ok = trialdiv_one_side (Fab, scaled_poly_a, poly->degree, a, b, 
				  primes_a, &nr_primes_a, max_nr_primes,
				  proj_divisor_a, nr_proj_primes_a, 
				  proj_primes_a, fba->fullfb, 
				  reports_a + i,
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
				  reports_r + j,
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
	  while (i < reports_a_nr && reports_a[i].a == reports_a[i + 1].a)
	    i++;
	  while (j < reports_r_nr && reports_r[j].a == reports_r[j + 1].a)
	    j++;
	}

      /* Assumes values in reports_a are sorted, same for reports_r  */
      if (reports_a[i].a < reports_r[j].a)
	i++;
      else
	j++;
    }

  if (verbose)
    {
      printf ("# Number of matching reports: %u\n", matching_reports);
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
  sieve_report_t *reports_a, *reports_r;
  factorbase_t fba, fbr;
  char *fbfilename = NULL, *polyfilename = NULL;
  unsigned char *sievearray;
  unsigned int reports_a_len = 0, reports_r_len = 0;
  unsigned int reports_a_nr, reports_r_nr;
  int verbose = 0;
  unsigned int deg;
  unsigned int i;
  double dpoly_a[MAXDEGREE], dpoly_r[2];
  const double log_scale = 1. / log (2.); /* Lets use log_2() for a start */
  cado_poly cpoly;
  char report_a_threshold, report_r_threshold;
  unsigned long relations_found = 0;
  double total_time = seconds ();

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
  if(!read_polynomial (cpoly, polyfilename))
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
  if (reports_a_len == 0) /* if not given on command line */
    reports_a_len = ((amax - amin + 1)) / 5 + 1000;
  if (reports_r_len == 0) /* if not given on command line */
    reports_r_len = ((amax - amin + 1)) / 5 + 1000;

  if (verbose)
    {
      printf ("# Allocating %ld bytes for reports_a\n",
	      reports_a_len * (long int) sizeof (sieve_report_t));
      fflush(stdout);
    }
  reports_a = (sieve_report_t *) malloc (reports_a_len * 
					 sizeof (sieve_report_t));
  ASSERT (reports_a != NULL);
  if (verbose)
    {
      printf ("# Allocating %ld bytes for reports_r\n",
	      reports_r_len * (long int) sizeof (sieve_report_t));
      fflush(stdout);
    }
  reports_r = (sieve_report_t *) malloc (reports_r_len * 
					 sizeof (sieve_report_t));
  ASSERT (reports_r != NULL);

  for (b = bmin; b <= bmax; b++)
    {
      unsigned long proj_roots; /* = gcd (b, leading coefficient) */
      
      if (verbose)
	printf ("# Sieving line b = %lu\n", b);

      proj_roots = mpz_gcd_ui (NULL, cpoly->f[deg], b);
      if (verbose)
	printf ("# Projective roots for b = %lu on algebtaic side are: %lu\n", 
	        b, proj_roots);

      if (verbose)
	  printf ("# Sieving algebraic side\n");

#ifdef PARI
      mp_poly_print (cpoly->f, cpoly->degree, "P(a,b) = ", 1);
      printf (" /* PARI */\n");
#endif

      reports_a_nr = 
	sieve_one_side (sievearray, fba,  
			reports_a, reports_a_len, report_a_threshold, 
			amin, amax, b, proj_roots, log_scale, 
			dpoly_a, deg, NULL, verbose);

      if (verbose)
	printf ("# There were %d sieve reports on the algebraic side\n",
		reports_a_nr);

      proj_roots = mpz_gcd_ui (NULL, cpoly->g[1], b);
      if (verbose)
	printf ("# Projective roots for b = %lu on rational side are: %lu\n", 
	        b, proj_roots);

      if (verbose)
	printf ("# Sieving rational side\n");

#ifdef PARI
      mp_poly_print (cpoly->g, 1, "P(a,b) = ", 1);
      printf (" /* PARI */\n");
#endif

      reports_r_nr = 
	sieve_one_side (sievearray, fbr, 
			reports_r, reports_r_len, report_r_threshold, 
			amin, amax, b, proj_roots, log_scale, 
			dpoly_r, 1, reports_a, verbose);
      if (verbose)
	printf ("# There were %d sieve reports on the rational side\n",
		reports_r_nr);

      if (reports_a_len == reports_a_nr)
	  fprintf (stderr, "Warning: sieve reports list on algebraic side "
		   "full with %u entries for b=%lu\n", reports_a_len, b);

      if (reports_r_len == reports_r_nr)
	  fprintf (stderr, "Warning: sieve reports list on rational side "
		   "full with %u entries for b=%lu\n", reports_r_len, b);

      /* Not trial factor the candidate relations */
      relations_found += trialdiv_and_print (cpoly, b, reports_a, reports_a_nr,
                                             reports_r, reports_r_nr, fba, fbr,
                                             log_scale, verbose);
    }

  total_time = seconds () - total_time;
  fprintf (stderr, "Found %lu relations in %1.0f seconds (%1.2e s/r)\n",
           relations_found, total_time, total_time / (double) relations_found);

#if 0
  /* FIXME implement this */
  fb_clear (fba);
  fb_clear (fbr);
#endif
  free (reports_a);
  free (reports_r);
  free (sievearray);

  return 0;
}
