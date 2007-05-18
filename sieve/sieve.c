/* Apply the sieve updates to sievearray */
/* Output: sievearray. Must have enough memory allocated */
/* Inputs: amin, amax, b the $a$ range and the line $b$ to sieve */
/* factorbase, polynomial */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <asm/msr.h>
#include "cado.h"
#include "mod_ul.c"
#define MAXDEGREE 10
#ifdef WANT_ASSERT
  #include <assert.h>
  #define ASSERT(x) assert(x)
#else
  #define ASSERT(x)
#endif


#define REFAC_PRP_THRES 1000
#define SIEVE_PERMISSIBLE_ERROR 5

/*****************************************************************
 *                      Functions for calculus                   *
 *****************************************************************/

/* Evaluate poly */

double 
calc_poly_eval (const double *poly, const int deg, const double x)
{
    double r;
    int i;

    ASSERT (deg <= MAXDEGREE);
    r = poly[deg];
    for (i = deg - 1; i >= 0; i--)
	r = r * x + poly[i];

    return r;
}

/* deriv == poly is ok */
void
calc_deriv (double *deriv, const double *poly, const int deg)
{
    int i;
    ASSERT (deg <= MAXDEGREE);
    for (i = 1; i <= deg; i++)
	deriv[i - 1] = (double) i * poly[i];
}

/* Newton root finding */
double
calc_newton_root (const double *poly, const int deg, const double start)
{
    double x, lastx, deriv[MAXDEGREE + 1];

    ASSERT (deg <= MAXDEGREE);
    calc_deriv (deriv, poly, deg);

    x = start;
    lastx = x + 1.;
    while (fabs (x - lastx) > 1.e-5) /* 1.e-5 should not take too long,
					even for repeated roots */
    {
	lastx = x;
	x -= calc_poly_eval(poly, deg, x) / calc_poly_eval(deriv, deg - 1, x);
    }

    /* Two more, with quadratic convergence that should give error < 1.e-20
       which is good enough even for long double. For repeated roots: FIXME.*/
    x -= calc_poly_eval(poly, deg, x) / calc_poly_eval(deriv, deg - 1, x);
    x -= calc_poly_eval(poly, deg, x) / calc_poly_eval(deriv, deg - 1, x);

    return x;
}

/* Divide a polynomial by (x - root). Returns the remainder, which should
   be small if root was indeed a root */
double
calc_divide_root (double *quotient, const double *poly, const int deg, 
		  const double root)
{
    double a, b;
    int i;

    a = poly[deg];
    for (i = deg; i > 0; i--)
    {
	b = a;
	a = poly[i - 1] + a * root;
	quotient[i - 1] = b;
    }
    
    return a;
}

int
calc_find_extrema (double *extrema, const double *poly, const int deg)
{
    double deriv[MAXDEGREE + 1];
    double root;
    double *derivptr = deriv;
    int i, j = 0;
    
    ASSERT(deg <= MAXDEGREE);
    ASSERT(poly[deg] != 0.);
    /* Constant polys have no extrema and would mess up code below */
    if (deg == 0)
	return 0;
    calc_deriv (deriv, poly, deg);
    /* Deal with multiple roots at 0 which may be plenty in SNFS polys */
    if (*derivptr == 0.)
    {
	extrema[j++] = 0.;
	while (*derivptr == 0.)
	    derivptr++; 
    }
    /* Find number of real roots. TBD. */
    /* FIXME for repeated roots */
    for (i = deg - 1; i > 0; i--)
    {
	root = calc_newton_root (deriv, i, 0.);
	extrema[j++] = root;
	printf ("calc_find_extrema: f´(%f) = %f\n", 
		root, calc_poly_eval(deriv, i, root));
	calc_divide_root (deriv, deriv, i, root);
    }
    return j;
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
mp_poly_print (mpz_t *poly, int deg, const char *name)
{
  int i;
  printf ("%s = ", name);
  for (i = deg; i >= 0; i--)
    {
      if (mpz_sgn (poly[i]) != 0)
	gmp_printf ("%s%Zd * x^%d ", 
		    mpz_sgn(poly[i]) > 0 ? "+" : "", poly[i], i);
    }
  printf ("\n");
}

/*****************************************************************
 *           Simple functions for (modular) arithmetic           *
 *****************************************************************/


/* Returns 0 if n is prime, otherwise the smallest prime factor of n */
static fbprime_t
iscomposite (const fbprime_t n)
{
  fbprime_t i, i2;

  if (n % 2 == 0)
    return (n == 2) ? 0 : 2;

  /* (i + 2)^2 = i^2 + 4*i + 4 */
  for (i = 3, i2 = 9; i2 <= n; i2 += (i+1) * 4, i += 2)
    if (n % i == 0)
	return i;

  return 0;
}

static unsigned long
gcd (unsigned long a, unsigned long b)
{
  unsigned long t;

  if (a == 0)
    return b;

  if (a >= b)
    a %= b;

  while (a > 0)
    {
      t = b % a;
      b = a;
      a = t;
    }

  return b;
}


/*****************************************************************
 *                Functions for the factor base                  *
 *****************************************************************/

/* Hack to get around C's automatic multiplying constants to be added to 
   pointers by the pointer's base data type size */

static inline factorbase_t *
fb_skip (const factorbase_t *fb, const size_t s)
{
  return (factorbase_t *)((char *)fb + s);
}

static inline factorbase_t *
fb_next (const factorbase_t *fb)
{
  return (factorbase_t *)((char *)fb + fb->size);
}


static size_t
fb_entrysize (const factorbase_t *fb)
{
  return (sizeof (factorbase_t) + fb->nr_roots * sizeof (fbroot_t));
}

static size_t
fb_entrysize_uc (const unsigned char n)
{
  return (sizeof (factorbase_t) + n * sizeof (fbroot_t));
}

void 
fb_print_entry (factorbase_t *fb)
{
  int i;
  printf ("Prime " FBPRIME_FORMAT " with rounded log %d and roots ", 
	  fb->p, (int) fb->plog);
  for (i = 0; i + 1 < fb->nr_roots; i++)
    printf (FBROOT_FORMAT ", ", fb->roots[i]);
  printf (FBROOT_FORMAT "\n", fb->roots[i]);
}


/* Add fb_add to (void *)fb + fbsize. If a realloc failed, returns NULL.
   fb_add->size need not be set by caller, this function does it */

factorbase_t *
fb_add_to (factorbase_t *fb, size_t *fbsize, size_t *fballoc,
	   const size_t allocblocksize, factorbase_t *fb_add)
{
  const size_t fb_addsize  = fb_entrysize (fb_add); 
  factorbase_t *newfb = fb;

  ASSERT(fb_addsize <= allocblocksize); /* Otherwise we still might not have
					   enough mem after the realloc */

  /* Do we need more memory for fb? */
  if (*fballoc < *fbsize + fb_addsize)
    {
      *fballoc += allocblocksize;
      newfb = (factorbase_t *) realloc (fb, *fballoc);
      if (newfb == NULL)
	{
	  fprintf (stderr, 
		   "Could not reallocate factor base to %lu bytes\n",
		   (unsigned long) *fballoc);
	  return NULL;
	}
    }
  memcpy (fb_skip(newfb, *fbsize), fb_add, fb_addsize);
  fb_skip(newfb, *fbsize)->size = fb_addsize;
  
  *fbsize += fb_addsize;

  return newfb;
}

static factorbase_t *
fb_find_p (factorbase_t *fb, const fbprime_t p)
{
  while (fb->p > 0 && fb->p < p) /* Assumes fb is sorted in asc. order */
    fb = fb_next (fb);

  if (fb->p == p)
    return fb;

  return NULL;
}

static unsigned char
fb_log (double n, double log_scale, double offset)
{
  return (unsigned char) floor (log (n) * log_scale + offset + 0.5);
}

/* Generate a factor base with primes <= bound for a linear polynomial */

factorbase_t *
fb_make_linear (mpz_t *poly, const fbprime_t bound, const double log_scale, 
		const int verbose)
{
  long long tsc1, tsc2;
  fbprime_t p;
  factorbase_t *fb = NULL, *fb_cur, *fb_new;
  size_t fbsize = 0, fballoc = 0;
  const size_t allocblocksize = 1<<20;

  rdtscll (tsc1);
  fb_cur = (factorbase_t *) malloc (fb_entrysize_uc (1));
  ASSERT (fb_cur != NULL);

  fb_cur->nr_roots = 1;
  fb_cur->size = fb_entrysize_uc (1);

  if (verbose)
    {
      printf ("# Making factor base for polynomial g(x)");
      mp_poly_print (poly, 1, "");
    }

  p = 2;
  while (p <= bound)
    {
      modulus m;
      residue r1, r2;

      fb_cur->p = p;
      fb_cur->plog = fb_log (fb_cur->p, log_scale, 0.);

      mod_initmod_ul (m, p);
      mod_init (r1, m);
      mod_init (r2, m);
      
      mod_set_ul_reduced (r2, mpz_fdiv_ui (poly[1], p), m);
      /* If p | g1 and !(p | g0), p|G(a,b) <=> p|b and that will be handeled
	 by the projective roots */
      if (mod_get_ul (r2, m) == 0)
	continue;
      mod_set_ul_reduced (r1, mpz_fdiv_ui (poly[0], p), m);

      /* We want g_1 * a + g_0 * b == 0 <=> a/b == - g0 / g1 */
      mod_inv (r2, r2, m); /* r2 = 1 / g1 */

      mod_mul (r2, r1, r2, m); /* r2 = g0 / g1 */
      mod_neg (r2, r2, m); /* r2 = - g0 / g1 */

      fb_cur->roots[0] = mod_get_ul (r2, m);

#ifdef WANT_ASSERT
      {
	residue r3;
	mod_init (r3, m);
	mod_set_ul_reduced (r3, mpz_fdiv_ui (poly[1], p), m);
	mod_mul (r3, r3, r2, m); /* g1 * (- g0 / g1) */
	mod_add (r3, r3, r1, m); /* g1 * (- g0 / g1) + g0 */
	ASSERT (mod_get_ul (r3, m) == 0);
	mod_clear (r3, m);
      }
#endif

      mod_clear (r1, m);
      mod_clear (r2, m);
      mod_clearmod (m);

      fb_new = fb_add_to (fb, &fbsize, &fballoc, allocblocksize, fb_cur);
      if (fb_new == NULL)
	{
	  free (fb);
	  fb = NULL;
	  break;
	}
      /* fb_print_entry (fb_cur); */
      fb = fb_new;
      
      /* FIXME: handle prime powers */

      /* Skip to next prime. FIXME: use a sieve */
      do {p++;} while (iscomposite (p));
    }

  if (fb != NULL) /* If nothing went wrong so far, put the end-of-fb mark */
    {
      fb_cur->p = 0;
      fb_cur->nr_roots = 0;
      fb_new = fb_add_to (fb, &fbsize, &fballoc, allocblocksize, fb_cur);
      if (fb_new == NULL)
	free (fb);
      fb = fb_new;
    }

  rdtscll (tsc2);
  if (verbose)
    printf ("# Generating rational factor base took %lld clocks\n", 
	    tsc2 - tsc1);  

  return fb;
}


factorbase_t *
fb_read (const char *filename, const double log_scale, int verbose)
{
  factorbase_t *fb = NULL, *fb_cur, *fb_new;
  FILE *fbfile;
  size_t fbsize = 0, fballoc = 0;
  const size_t linesize = 255;
  size_t linelen;
  char line[linesize];
  char *lineptr;
  unsigned int i;
  int ok;
  unsigned long linenr = 0;
  const size_t allocblocksize = 1<<20; /* Allocate in MB chunks */
  fbprime_t p;
  fbprime_t maxprime = 0;
  unsigned long nr_primes = 0;

  fbfile = fopen (filename, "r");
  if (fbfile == NULL)
    {
      fprintf (stderr, "Could not open file %s for reading\n", filename);
      return NULL;
    }

  fb_cur = (factorbase_t *) malloc (sizeof (factorbase_t) + 
				      MAXDEGREE * sizeof(fbroot_t));
  if (fb_cur == NULL)
    {
      fprintf (stderr, "Could not allocate memory for factor base\n");
      fclose (fbfile);
      return NULL;
    }

  while (!feof(fbfile))
    {
      if (fgets (line, linesize, fbfile) == NULL)
	break;
      linenr++;
      linelen = strlen (line);
      if (linelen > 0 && line[linelen - 1] == '\n')
	linelen--; /* Remove newline */
      for (i = 0; i < linelen; i++) /* Skip comments */
	if (line[i] == '#')
	  break;
      linelen = i;
      while (linelen > 0 && isspace (line[i - 1]))
	linelen--; /* Skip whitespace at end of line */
      if (linelen == 0) /* Skip empty/comment lines */
	continue;
      line[linelen] = '\0';

      /* Parse the line */
      lineptr = line;
      fb_cur->p = strtoul (lineptr, &lineptr, 10);
      ok = 0;
      if (fb_cur->p == 0)
	fprintf (stderr, "fb_read: prime is not an integer on line %lu\n", 
		 linenr);
      else if (*lineptr != ':')
	fprintf (stderr, "fb_read: prime is not followed by colon on line %lu",
		 linenr);
      else
	ok = 1;
      if (!ok)
	continue;

      if (fb_cur->p <= maxprime)
	{
	  fprintf (stderr, "Error, primes in factor base file are "
		   "not in increasing order\n");
	  free (fb);
	  fb = NULL;
	  break;
	}

      p = iscomposite (fb_cur->p);
      if (p == 0) /* It's a prime, do a normal log */
	fb_cur->plog = fb_log (fb_cur->p, log_scale, 0.);
      else
	{
	  /* Take into account the logs of the smaller powers of p that
	     have been added already. We should not just use 
	     fb_cur->plog = log(p) as that would allow rounding errors to 
	     accumulate. */
	  double oldlog;
	  oldlog = fb_log (fb_cur->p / p, log_scale, 0.);
	  fb_cur->plog = fb_log (fb_cur->p, log_scale, - oldlog);
	}

      lineptr++; /* Skip colon */
      ok = 1;
      fb_cur->nr_roots = 0;
      /* Read roots */
      while (ok && *lineptr != '\0' && fb_cur->nr_roots <= MAXDEGREE)
	{
	  fbroot_t r = strtoul (lineptr, &lineptr, 10);
	  fb_cur->roots[fb_cur->nr_roots++] = r;
	  if (r >= fb_cur->p) /* Check if root is properly reduced */
	    ok = 0;
	  if (*lineptr != '\0' && *lineptr != ',')
	    ok = 0;
	  if (*lineptr == ',')
	    lineptr++;
	}

      if (!ok)
	{
	  fprintf (stderr, "Incorrect format in factor base file line %lu\n",
		   linenr);
	  continue;
	}

      if (fb_cur->nr_roots > MAXDEGREE)
	{
	  printf ("Error, too many roots for prime " FBPRIME_FORMAT 
		  "in factor base line %lu\n", p, linenr);
	  continue;
	}
      
      fb_new = fb_add_to (fb, &fbsize, &fballoc, allocblocksize, fb_cur);
      if (fb_new == NULL)
	{
	  free (fb);
	  fb = NULL;
	  break;
	}
      /* fb_print_entry (fb_cur); */
      fb = fb_new;
      maxprime = fb_cur->p;
      nr_primes++;
    }      

  if (fb != NULL) /* If nothing went wrong so far, put the end-of-fb mark */
    {
      fb_cur->p = 0;
      fb_cur->nr_roots = 0;
      fb_new = fb_add_to (fb, &fbsize, &fballoc, allocblocksize, fb_cur);
      if (fb_new == NULL)
	free (fb);
      fb = fb_new;
    }

  if (fb != NULL && verbose)
    {
      printf ("Factor base sucessfully read, %lu primes, largest was "
	      FBPRIME_FORMAT "\n", nr_primes, maxprime);
    }
  
  fclose (fbfile);
  free (fb_cur);
  return fb;
}

void
fb_disable_roots (factorbase_t *fb, const unsigned long b, const int verbose)
{
  long long tsc1, tsc2;
  unsigned long t;
  
  /* Remove the roots of primes that divide b, and of the powers of those 
     primes, from factor base */
  rdtscll (tsc1);
  t = b;
  while (t > 1)
    {
      factorbase_t *fb_del;
      unsigned long p, ppow;
      p = iscomposite (t);
      if (p == 0)
	p = t;
      do t /= p; while (t % p == 0);
      ppow = p;
      while ((fb_del = fb_find_p (fb, ppow)) != NULL)
	{
	  if (verbose)
	    {
	      printf ("# Temporarily removing from factor base ");
	      fb_print_entry (fb_del);
	    }
	  fb_del->nr_roots = 0;
	  ppow *= p;
	}
    }
  rdtscll (tsc2);
  if (verbose)
    printf ("# Removing primes from fb took %lld clocks\n", tsc2 - tsc1);
}

void
fb_restore_roots (factorbase_t *fb, const unsigned long b, int verbose)
{
  long long tsc1, tsc2;
  unsigned long t;

  /* Put the roots back in */
  rdtscll (tsc1);
  t = b;
  while (t > 1)
    {
      factorbase_t *fb_restore;
      unsigned long p, ppow;
      p = iscomposite (t);
      if (p == 0)
	p = t;
      t /= p;
      ppow = p;
      while ((fb_restore = fb_find_p (fb, ppow)) != NULL)
	{
	  fb_restore->nr_roots = 
	    (fb_restore->size - sizeof (factorbase_t)) 
	    / sizeof (fbroot_t);
	  if (verbose)
	    {
	      printf ("# Restored to factor base ");
	      fb_print_entry (fb_restore);
	    }
	  ppow *= p;
	}
    }
  rdtscll (tsc2);
  if (verbose)
    printf ("# Restoring primes to fb took %lld clocks\n", tsc2 - tsc1);
}

fbprime_t
fb_maxprime (factorbase_t *fb)
{
    fbprime_t maxp;
    factorbase_t *fbptr = fb;
    while (fbptr->p != 0)
    {
	ASSERT (maxp <= fbptr->p);
	maxp = fbptr->p;
	fbptr = fb_next (fbptr);
    }
    
    return maxp;
}


static unsigned char
log_norm (const double *f, const int deg, const double x, 
	  const double log_scale, const double log_proj_roots)
{
  double r;
  int i;
  unsigned char n;

  r = f[deg];
  for (i = deg - 1; i >= 0; i--)
    r = r * x + f[i];

  n = fb_log (fabs (r), log_scale, - log_proj_roots);
  /* printf ("Norm at x = %.0f is %.0f, rounded log is %d\n", x, r, (int) n); */
  return n;
}


/* Very slow but thorough way of computing norms. Returns the maximum rounded
   log norm in sievearray. */

unsigned char
compute_norms (unsigned char *sievearray, const long amin, const long amax, 
	       const unsigned long b, const double *poly, const int deg, 
	       const double proj_roots, const double log_scale, 
	       const int verbose)
{
    double f[MAXDEGREE + 1]; /* Poly in a for a given fixed b */
    double bpow;
    const double log_proj_roots = log(proj_roots) * log_scale;
    unsigned long long tsc1, tsc2;
    long a, a2;
    int i;
    unsigned char n1, n2, nmax;
    const int stride = 256;

    rdtscll (tsc1);

    bpow = 1.;
    for (i = deg; i >= 0; i--)
    {
	f[i] = poly[i] * bpow;
	bpow *= (double) b;
    }

    n1 = log_norm (f, deg, (double) amin, log_scale, log_proj_roots);
    nmax = n1;

    a2 = amin;
    for (a = amin; a <= amax;)
    {
      ASSERT (a = a2);
      ASSERT (n1 = log_norm (f, deg, (double) a, log_scale, log_proj_roots));
      /* We'll cover [a ... a2 - 1] now */
      a2 = (a + stride <= amax + 1) ? a + stride : amax + 1;

      n2 = log_norm (f, deg, (double) a2, log_scale, log_proj_roots);

      if (n1 == n2) /* Let's assume the log norm is n1 everywhere in this
		       interval */
	{
	  /* printf ("log_c(F(%ld, %lu)) == log_c(F(%ld, %lu)) == %d\n",
	     a, b, a2, b, n1); */
	  memset (sievearray + a - amin, n1, a2 - a);
	}
      else
	{
	  /* n1 and n2 are different. Do each a up to (exclusive) a2 
	     individually */
	  /* printf ("log_c(F(%ld, %lu)) == %d != log_c(F(%ld, %lu)) == %d\n",
	     a, b, n1, a2, b, n2); */
	  sievearray[a - amin] = n1;
	  for ( ; a < a2; a++)
	    {
	      unsigned char n = log_norm (f, deg, (double) a, log_scale, 
					       log_proj_roots);
	      sievearray[a - amin] = n;
	      if (n > nmax)
		nmax = n;
	    }

	}
      a = a2;
      n1 = n2;
    }

    rdtscll (tsc2);
    if (verbose)
      {
	printf ("# Computing norms took %lld clocks\n", tsc2 - tsc1);
	printf ("# Maximum rounded log norm is %u\n", (unsigned int) nmax);
      }

    return nmax;
}

/* sievearray must be long enough to hold amax-amin+1 chars.
   proj_roots is the rounded log of the projective roots */

void 
sieve (unsigned char *sievearray, factorbase_t *fb, 
       const long amin, const long amax, const unsigned long b, 
       const unsigned char threshold, fbprime_t *useful_primes_par,
       const fbprime_t useful_threshold, const unsigned int useful_length_par)
{
  uint32_t i, amin_p, p, d;
  unsigned int useful_length = useful_length_par;
  fbprime_t *useful_primes = useful_primes_par;
  unsigned char plog;
  const unsigned long l = amax - amin + 1;

  ASSERT (amax > amin);

  if (useful_length > 0)
    useful_length--; /* Make sure we have space for a 0 mark at the end */

  /* Do the sieving */

  while (fb->p > 0)
    {
      ASSERT (fb_entrysize(fb) <= fb->size);
      p = fb->p;
      plog = fb->plog;

      /* This modular reduction should be simplified somehow. Do it once
         and store it in fb? */
      if (amin < 0)
	{
	  amin_p = p - ((unsigned long)(-amin) % p);
	  if (amin_p == p)
	    amin_p = 0; /* FIXME ugly */
	}
      else
        amin_p = (unsigned long) amin % p;

      for (i = 0; i < fb->nr_roots; i++)
        {
	  modulus m;
	  residue r1, r2;

          /* Find first index in sievearray where p on this root divides.
             So we want a/b == r (mod p) <=> a == br (mod p). Then we want
             d so that amin + d == a (mod p) <=> d == a - amin (mod p)
          */
	  ASSERT (fb->roots[i] < p);
	  ASSERT (b % p != 0);
	  mod_initmod_ul (m, p); /* Most of the mod_*() calls are no-ops */
	  mod_init (r1, m);
	  mod_init (r2, m);
	  mod_set_ul (r1, b, m); /* Modular reduction */
	  mod_set_ul_reduced (r2, fb->roots[i], m);
          mod_mul (r1, r1, r2, m); /* Multiply and mod reduction. If we keep 
	    r_i in Montgomery representation, a single mul/REDC will compute 
	    the product, reduce it mod p and return it as an integer. TBD. */
	  mod_sub_ul (r1, r1, amin_p, m);
          d = mod_get_ul (r1, m); 
	  mod_clear (r1, m);
	  mod_clear (r2, m);
	  mod_clearmod (m);
	  
	  /* Now d is the first index into sievearray where p divides. */

          /* Now update the sieve array. There is no partitioning, blocking, 
             bucket sieving or anything atm. */
          for (; d < l; d += p)
	    {
	      unsigned char k;
	      k = sievearray[d] - plog;
	      sievearray[d] = k;
	      if (k <= threshold)
		if (p >= useful_threshold && useful_length > 0)
		  {
		    *useful_primes++ = p;
		    useful_length--;
		  }
	    }
        }
      
      /* Move on to the next factor base prime */
      fb = fb_next (fb);
    }
  if (useful_length_par > 0)
    *useful_primes++ = 0;
}

#define TYPE_STRING 0
#define TYPE_MPZ 1
#define TYPE_INT 2
#define TYPE_ULONG 3
#define TYPE_DOUBLE 4
#define PARSE_MATCH -1
#define PARSE_ERROR 0
#define PARSE_NOMATCH 1

/* Return value: 1: tag does not match, 0: error occurred, -1 tag did match */

static int
parse_line (void *target, char *line, const char *tag, int *have, 
	    const int type)
{
  char *lineptr = line;

  if (strncmp (lineptr, tag, strlen (tag)) != 0)
    return PARSE_NOMATCH;

  if (have != NULL && *have != 0)
    {
      fprintf (stderr, "parse_line: %s appears twice\n", tag);
      return PARSE_ERROR;
    }
  if (have != NULL)
    *have = 1;
  
  lineptr += strlen (tag);
  if (type == TYPE_STRING) /* character string of length up to 256 */
    {
      strncpy ((char *)target, lineptr, 256);
      ((char *)target)[255] = '0';
    }
  else if (type == TYPE_MPZ)
    {
      mpz_init (*(mpz_t *)target);
      mpz_set_str (*(mpz_t *)target, lineptr, 0);
    }
  else if (type == TYPE_INT)
    *(int *)target = atoi (lineptr);
  else if (type == TYPE_ULONG)
    *(unsigned long *)target = strtoul (lineptr, NULL, 10);
  else if (type == TYPE_DOUBLE)
    *(double *)target = atof (lineptr);
  else return PARSE_ERROR;
  
  return PARSE_MATCH;
}

cado_poly *
read_polynomial (char *filename)
{
  FILE *file;
  const int linelen = 512;
  char line[linelen];
  cado_poly *poly;
  int have_name = 0, have_n = 0, have_Y0 = 0, have_Y1 = 0;
  int i, ok;

  file = fopen (filename, "r");
  if (file == NULL)
    {
      fprintf (stderr, "read_polynomial: could not open %s\n", filename);
      return NULL;
    }

  poly = (cado_poly *) malloc (sizeof(cado_poly));
  (*poly)->f = (mpz_t *) malloc (MAXDEGREE * sizeof (mpz_t));
  (*poly)->g = (mpz_t *) malloc (2 * sizeof (mpz_t));

  (*poly)->name[0] = '\0';
  (*poly)->degree = -1;
  (*poly)->type[0] = '\0';

  while (!feof (file))
    {
      ok = 1;
      if (fgets (line, linelen, file) == NULL)
	break;
      if (line[0] == '#')
	continue;

      ok *= parse_line (&((*poly)->name), line, "name: ", &have_name, TYPE_STRING);
      ok *= parse_line (&((*poly)->n), line, "n: ", &have_n, TYPE_MPZ);
      ok *= parse_line (&((*poly)->skew), line, "skew: ", NULL, TYPE_DOUBLE);
      ok *= parse_line (&((*poly)->g[0]), line, "Y0: ", &have_Y0, TYPE_MPZ);
      ok *= parse_line (&((*poly)->g[1]), line, "Y1: ", &have_Y1, TYPE_MPZ);
      ok *= parse_line (&((*poly)->type), line, "type: ", NULL, TYPE_STRING);
      ok *= parse_line (&((*poly)->rlim), line, "rlim: ", NULL, TYPE_ULONG);
      ok *= parse_line (&((*poly)->alim), line, "alim: ", NULL, TYPE_ULONG);
      ok *= parse_line (&((*poly)->lpbr), line, "lpbr: ", NULL, TYPE_INT);
      ok *= parse_line (&((*poly)->lpba), line, "lpba: ", NULL, TYPE_INT);
      ok *= parse_line (&((*poly)->mfbr), line, "mfbr: ", NULL, TYPE_INT);
      ok *= parse_line (&((*poly)->mfba), line, "mfba: ", NULL, TYPE_INT);
      ok *= parse_line (&((*poly)->rlambda), line, "rlambda: ", NULL, TYPE_DOUBLE);
      ok *= parse_line (&((*poly)->alambda), line, "alambda: ", NULL, TYPE_DOUBLE);
      ok *= parse_line (&((*poly)->qintsize), line, "qintsize: ", NULL, TYPE_INT);


      if (ok == 1 && line[0] == 'c' && isdigit (line[1]) && line[2] == ':' &&
	       line[3] == ' ')
	{
	  int index = line[1] - '0', i;
	  for (i = (*poly)->degree + 1; i <= index; i++)
	    mpz_init ((*poly)->f[i]);
	  if (index > (*poly)->degree)
	    (*poly)->degree = index;
	  mpz_set_str ((*poly)->f[index], line + 4, 0);
	  ok = -1;
	}

      if (ok == PARSE_NOMATCH)
	{
	  fprintf (stderr, 
		   "read_polynomial: Cannot parse line %s\nIgnoring.\n", 
		   line);
	  continue;
	}
      
      if (ok == PARSE_ERROR)
	break;
    }
  
  if (ok != PARSE_ERROR)
    {
      if (have_n == 0)
	{
	  fprintf (stderr, "n ");
	  ok = PARSE_ERROR;
	}
      if (have_Y0 == 0)
	{
	  fprintf (stderr, "Y0 ");
	  ok = PARSE_ERROR;
	}
      if (have_Y1 == 0)
	{
	  fprintf (stderr, "Y1 ");
	  ok = PARSE_ERROR;
	}
      if (ok == PARSE_ERROR)
	fprintf (stderr, "are missing in polynomial file\n");
    }

  if (ok == PARSE_ERROR)
    {
      for (i = 0; i <= (*poly)->degree; i++)
	mpz_clear ((*poly)->f[i]);
      if (have_n)
	mpz_clear ((*poly)->n);
      if (have_Y0)
	mpz_clear ((*poly)->g[0]);
      if (have_Y1)
	mpz_clear ((*poly)->g[1]);
      free ((*poly)->f);
      free ((*poly)->g);
      free (poly);
      poly = NULL;
    }
  
  return poly;
}

void
print_useful (const fbprime_t *useful_primes, 
	      const unsigned int useful_length)
{
  unsigned int useful_count = 0;

  if (useful_primes == NULL)
    return;

  if (*useful_primes == 0)
    {
      printf ("# There were no useful primes\n");
      return;
    }

  printf ("# Useful primes were: ");
  while (*useful_primes != 0)
    {
      useful_count ++;
      printf (FBPRIME_FORMAT " ", *useful_primes++);
    }
  printf ("\n");
  if (useful_count == useful_length - 1)
    printf ("#Storage for useful primes is full, "
	    "consider increasing useful_length\n");
}


unsigned long
sieve_one_side (unsigned char *sievearray, factorbase_t *fb,
		long *reports, const unsigned int reports_len, 
		const unsigned char reports_threshold,
		fbprime_t *useful_primes, const fbprime_t useful_threshold, 
		const unsigned int useful_length,
		const long amin, const long amax, const unsigned long b,
		const unsigned long proj_roots, const double log_scale,
		const double *dpoly, const unsigned int deg, 
		const int verbose)
{
  long long tsc1, tsc2;
  long a;
  unsigned long reports_nr;

  fb_disable_roots (fb, b, verbose);
  
  compute_norms (sievearray, amin, amax, b, dpoly, deg, proj_roots, 
		 log_scale, verbose);
  
  rdtscll (tsc1);
  sieve (sievearray, fb, amin, amax, b, reports_threshold, useful_primes, 
	 useful_threshold, useful_length);
  rdtscll (tsc2);
  if (verbose)
    printf ("# Sieving took %lld clocks\n", tsc2 - tsc1);
  
  fb_restore_roots (fb, b, verbose);
  
  if (verbose)
    print_useful (useful_primes, useful_length);
  
  /* Store the sieve reports */
  rdtscll (tsc1);
  reports_nr = 0;
  for (a = amin; reports_nr < reports_len && a <= amax; a++)
    {
      /* The + 10 is to deal with accumulated rounding that might have
	 cause the sieve value to drop below 0 and wrap around */
      if ((unsigned char) (sievearray[a - amin] + 10) <= reports_threshold + 10
	  && gcd(labs(a), b) == 1)
	reports[reports_nr++] = a;
    }
  rdtscll (tsc2);
  if (verbose)
    {
      printf ("# There were %lu sieve reports\n", reports_nr);
      printf ("# Finding sieve reports took %lld clocks\n", tsc2 - tsc1);
    }

  return reports_nr;
}


/* Divides the prime q and its powers out of C, appends each q divided out
   to primes_a. If anything was divided out, returns 1 if cofactor is 1 or
   a probable prime. Otherwise returns 0. */

static int
trialdiv_one_prime (const fbprime_t q, mpz_t C, unsigned int *nr_primes_a,
                    fbprime_t *primes_a, const unsigned int max_nr_primes)
{
  int is_prp = 0;

  if (mpz_divisible_ui_p (C, q))
    {
      do
	{
	  if (*nr_primes_a < max_nr_primes)
	    primes_a[(*nr_primes_a)++] = q;
	  mpz_tdiv_q_ui(C, C, q);
	} while (mpz_divisible_ui_p (C, q));

      if (mpz_cmp_ui (C, 1) == 0 || 
	  (q > REFAC_PRP_THRES && mpz_probab_prime_p (C, 1)))
        {
          is_prp = 1;
        }
    }
  return is_prp;
}

static void
trialdiv_and_print (cado_poly *poly, const unsigned long b, 
                    const long *reports_a, const unsigned int reports_a_nr, 
		    fbprime_t *useful_primes_a, const long *reports_r, 
		    const unsigned int reports_r_nr, 
		    fbprime_t *useful_primes_r, 
		    factorbase_t *fba, const int verbose)
{
  unsigned long long tsc1, tsc2;
  unsigned int i, j, k;
  const unsigned int max_nr_primes = 128;
  mpz_t Fab, Gab, scaled_poly_a[MAXDEGREE], scaled_poly_r[2];
  fbprime_t primes_a[max_nr_primes], primes_r[max_nr_primes];
  fbprime_t q, incrq;
  unsigned int nr_primes_a, nr_primes_r;

  mpz_init (Fab);
  mpz_init (Gab);
  for (i = 0; i <= (unsigned) (*poly)->degree; i++)
    mpz_init (scaled_poly_a[i]);
  mpz_init (scaled_poly_r[0]);
  mpz_init (scaled_poly_r[1]);

  /* Multiply f_i by b^(deg-i) and put in scaled_poly */
  mp_poly_scale (scaled_poly_a, (*poly)->f, (*poly)->degree, b, -1); 
  mp_poly_scale (scaled_poly_r, (*poly)->g, 1, b, -1); 

  rdtscll (tsc1);
  for (i = 0, j = 0; i < reports_a_nr && j < reports_r_nr;)
    {
      if (reports_a[i] == reports_r[j])
      {
          const long a = reports_a[i];
          factorbase_t *nextfb;
      
	  mp_poly_eval (Fab, scaled_poly_a, (*poly)->degree, a);
	  mpz_abs (Fab, Fab);
	  mp_poly_eval (Gab, scaled_poly_r, 1, a);
	  mpz_abs (Gab, Gab);
	  
          nr_primes_a = nr_primes_r = 0;

          /* Do the algebraic side */
	  /* See if any of the "useful primes" divide this norm */
	  for (k = 0; useful_primes_a[k] != 0; k++)
	  {
	      const fbprime_t q = useful_primes_a[k];
	      if (mpz_divisible_ui_p (Fab, q))
	        {
	          if (nr_primes_a < max_nr_primes)
	            primes_a[nr_primes_a++] = q;
		  mpz_tdiv_q_ui (Fab, Fab, q);
	        }
	  }
	  
	  /* Go through the entire factor base */ 
	  for (nextfb = fba; nextfb->p != 0; nextfb = fb_next (nextfb))
            if (trialdiv_one_prime (nextfb->p, Fab, &nr_primes_a, primes_a,
                                    max_nr_primes))
	      break;
	  
	  /* Do the rational side */

	  while (mpz_even_p (Gab))
	    {
	      if (nr_primes_r < max_nr_primes)
	        primes_r[nr_primes_r++] = 2;
              mpz_div_2exp (Gab, Gab, 1);
	    }

	  for (k = 0; useful_primes_r[k] != 0; k++)
	  {
	      const fbprime_t q = useful_primes_r[k];
	      if (mpz_divisible_ui_p (Gab, q))
	        {
	          if (nr_primes_r < max_nr_primes)
	            primes_r[nr_primes_r++] = q;
		  mpz_tdiv_q_ui (Gab, Gab, q);
	        }
	  }
	  
          trialdiv_one_prime (3, Gab, &nr_primes_r, primes_r, max_nr_primes);
          
          for (q = 5, incrq = 2; q <= (*poly)->rlim; 
               q += incrq, incrq = 4 - incrq)
            if (trialdiv_one_prime (q, Gab, &nr_primes_r, primes_r,
                                    max_nr_primes))
	      break;

	  /* Test if the cofactor after trial division is small */
	  if ((double) mpz_sizeinbase (Fab, 2) > 
	      (*poly)->lpba * (*poly)->alambda + SIEVE_PERMISSIBLE_ERROR)
	  {
	      gmp_fprintf (stderr, 
			   "Sieve report %ld, %lu is not smooth on "
			   "algebraic side, cofactor is %Zd with %d bits\n", 
			   a, b, Fab, mpz_sizeinbase (Fab, 2));
	  }

	  else if ((double) mpz_sizeinbase (Gab, 2) > 
	      (*poly)->lpbr * (*poly)->rlambda + SIEVE_PERMISSIBLE_ERROR)
	  {
	      gmp_fprintf (stderr, 
			   "Sieve report %ld, %lu is not smooth on "
			   "rational side, cofactor is %Zd with %d bits\n", 
			   a, b, Gab, mpz_sizeinbase (Gab, 2));
	  } else {

	      /* Now print the relations */
	      printf ("%ld,%lu:", a, b);
	      for (k = 0; k < nr_primes_r; k++)
		  printf ("%x%s", primes_r[k], k+1==nr_primes_r?":":",");
	      for (k = 0; k < nr_primes_a; k++)
		  printf ("%x%s", primes_a[k], k+1==nr_primes_a?"\n":",");
	  }
        }
        
        /* Assumes values in reports_a are sorted and unique, 
	   same for reports_r  */
        if (reports_a[i] < reports_r[j])
          i++;
        else
          j++;
    }
  
  for (i = 0; i <= (unsigned)(*poly)->degree; i++)
    mpz_clear (scaled_poly_a[i]);
  mpz_clear (scaled_poly_r[0]);
  mpz_clear (scaled_poly_r[1]);
      
  mpz_clear (Fab);
  mpz_clear (Gab);
  
  rdtscll (tsc2);
  if (verbose)
    printf ("# Trial factoring/printing took %lld clocks\n", tsc2 - tsc1);
}


int
main (int argc, char **argv)
{
  long amin, amax;
  unsigned long bmin, bmax, b;
  fbprime_t *useful_primes_a, *useful_primes_r, useful_threshold;
  factorbase_t *fba, *fbr;
  long *reports_a, *reports_r;
  char *fbfilename = NULL, *polyfilename = NULL;
  unsigned char *sievearray;
  unsigned int useful_length_a, useful_length_r, reports_a_len, reports_a_nr, 
      reports_r_len, reports_r_nr;
  int verbose = 0;
  unsigned int deg;
  unsigned int i;
  double dpoly_a[MAXDEGREE], dpoly_r[1];
  const double log_scale = 1. / log (2.); /* Lets use log_2() for a start */
  cado_poly *cpoly;
  char report_a_threshold, report_r_threshold;


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

      else 
	break;
    }

  if (argc < 5)
    {
      fprintf (stderr, "Please specify amin amax bmin bmax\n");
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

  cpoly = read_polynomial (polyfilename);
  if (cpoly == NULL)
    {
      fprintf (stderr, "Error reading polynomial file\n");
      exit (EXIT_FAILURE);
    }

  fba = fb_read (fbfilename, log_scale, verbose);
  if (fba == NULL)
    {
      fprintf (stderr, "Could not read factor base\n");
      exit (EXIT_FAILURE);
    }

  fbr = fb_make_linear ((*cpoly)->g, (fbprime_t) (*cpoly)->rlim, log_scale,
			verbose);
  if (fbr == NULL)
    {
      fprintf (stderr, 
	       "Could not generate factor base for linear polynomial\n");
      exit (EXIT_FAILURE);
    }
  
  if (0)
  {
      factorbase_t *fbptr;
      for (fbptr = fbr; fbptr->p != 0; fbptr = fb_next (fbptr))
	  fb_print_entry(fbptr);
  }

  deg = (*cpoly)->degree;
  for (i = 0; i <= deg; i++)
      dpoly_a[i] = mpz_get_d ((*cpoly)->f[i]);
  dpoly_r[0] = mpz_get_d ((*cpoly)->g[0]);
  dpoly_r[1] = mpz_get_d ((*cpoly)->g[1]);
  report_a_threshold = (double)((*cpoly)->lpba) * log(2.0) * 
      (*cpoly)->alambda * log_scale + 0.5;
  report_r_threshold = (double)((*cpoly)->lpbr) * log(2.0) * 
      (*cpoly)->rlambda * log_scale + 0.5;

#if 0
  if (verbose)
    for (s = 0; fb_skip(fb, s)->p != 0; s += fb_entrysize (fb_skip (fb, s)))
      fb_print_entry (fb_skip(fb, s));
#endif

  sievearray = (unsigned char *) malloc ((amax - amin + 1) * sizeof (char));
  useful_length_a = ((amax - amin + 1)) / 100 + 1000;
  useful_length_r = ((amax - amin + 1)) / 10 + 1000;
  useful_primes_a = (fbprime_t *) malloc (useful_length_a* sizeof (fbprime_t));
  useful_primes_r = (fbprime_t *) malloc (useful_length_r* sizeof (fbprime_t));
  ASSERT (useful_primes_a != NULL && useful_primes_r != NULL);

  useful_threshold = fb_maxprime (fba) / 10;
  if (verbose)
    printf ("# Threshold for useful primes is " FBPRIME_FORMAT "\n", 
	    useful_threshold);

  reports_a_len = useful_length_a;
  reports_a = (long *) malloc (reports_a_len * sizeof (long));
  ASSERT (reports_a != NULL);
  reports_r_len = useful_length_r;
  reports_r = (long *) malloc (reports_r_len * sizeof (long));
  ASSERT (reports_r != NULL);

  for (b = bmin; b <= bmax; b++)
    {
      unsigned long proj_roots; /* = gcd (b, leading coefficient) */
      
      proj_roots = mpz_gcd_ui (NULL, (*cpoly)->f[deg], b);
      if (verbose)
	printf ("# Projective roots for b = %lu on algebtaic side are: %lu\n", 
	        b, proj_roots);

      if (verbose)
	  printf ("#Sieving algebraic side\n");

      reports_a_nr = 
	sieve_one_side (sievearray, fba, reports_a, reports_a_len, 
			report_a_threshold, useful_primes_a, useful_threshold,
			useful_length_a, amin, amax, b, proj_roots, log_scale, 
			dpoly_a, deg, verbose);

      proj_roots = mpz_gcd_ui (NULL, (*cpoly)->f[deg], b);
      if (verbose)
	printf ("# Projective roots for b = %lu on rational side are: %lu\n", 
	        b, proj_roots);

      if (verbose)
	  printf ("#Sieving rational side\n");

      reports_r_nr = 
	sieve_one_side (sievearray, fbr, reports_r, reports_r_len, 
			report_r_threshold, useful_primes_r, useful_threshold,
			useful_length_r, amin, amax, b, proj_roots, log_scale, 
			dpoly_r, 1, verbose);

      if (reports_a_len == reports_a_nr)
	  fprintf (stderr, "Warning: sieve reports list on algebraic side "
		   "full with %u entries for b=%lu\n", reports_a_len, b);

      if (reports_r_len == reports_r_nr)
	  fprintf (stderr, "Warning: sieve reports list on rational side "
		   "full with %u entries for b=%lu\n", reports_r_len, b);

      /* Not trial factor the candidate relations */
      trialdiv_and_print (cpoly, b, 
			  reports_a, reports_a_nr, useful_primes_a, 
			  reports_r, reports_r_nr, useful_primes_r, 
			  fba, verbose);
    }


  free (fba);
  free (useful_primes_a);
  free (useful_primes_r);
  free (reports_a);
  free (reports_r);
  free (sievearray);

  return 0;
}
