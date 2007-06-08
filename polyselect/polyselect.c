/* polyselect - polynomial selection

Author: Paul Zimmermann

Copyright 2007 INRIA

This file is part of the CADO project.

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA

   Usage: polyselect < c80 > c80.poly

   References:

   [1] "On polynomial selection for the general number field sieve",
       Thorsten Kleinjung, Mathematics of Computation 75 (2006), p. 2037-2047.
   [2] "Integer factorization, part 4: polynomial selection", Dan Bernstein,
       invited lecture, Arizona Winter School, University of Arizona, Tucson,
       Arizona, 14 March 2006.
   [3] "Polynomial Selection for the Number Field Sieve Integer Factorisation
       Algorithm", Brian Antony Murphy, Australian National University, 1999.
*/

#define VERSION 155 /* try to match the svn version */

#define METHOD 3 /* algorithm used for polynomial selection:
		    1: basic method in B^{d(d+1)/(d-1)} to save a factor B
		    2: Murphy's algorithm in B^{d+1}
		    3: skew polynomials */
#define MARGIN 2.0 /* margin for mu: keep only logmu <= best_logmu + MARGIN */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <limits.h> /* for ULONG_MAX */
#include <values.h> /* for DBL_MAX */
#include <math.h>   /* for log, fabs */
#include "cado.h"
#include <sys/types.h>
#include <sys/resource.h>

#define REPS 1 /* number of Miller-Rabin tests in isprime */
#define NUM_EST 32 /* number of iterations to estimate average valuations
                      for ``bad-behaved'' primes */
#define SQR(x) (x) * (x)

#ifdef DEBUG
#define ASSERT(x) assert(x)
#else
#define ASSERT(x)
#endif

static void
usage ()
{
  fprintf (stderr, "Usage: polyselect [-v] [-raw] [-e nnn] [-m xxx]< in > out\n\n");
  fprintf (stderr, "       -v     - verbose\n");
  fprintf (stderr, "       -raw   - does not output factor base parameters\n");
  fprintf (stderr, "       -e nnn - use effort nnn\n");
  fprintf (stderr, "       -m xxx - use decomposition in base xxx\n");
  fprintf (stderr, "       in     - input file (number to factor)\n");
  fprintf (stderr, "       out    - output file (polynomials)\n");
  exit (1);
}

int
cputime ()
{
  struct rusage rus;

  getrusage (0, &rus);
  return rus.ru_utime.tv_sec * 1000 + rus.ru_utime.tv_usec / 1000;
}

/************************** polynomial arithmetic ****************************/

/* allocate a polynomial of degree d, and initialize it */
mpz_t*
alloc_poly (int d)
{
  mpz_t *f;
  int i;

  f = (mpz_t*) malloc ((d + 1) * sizeof (mpz_t));
  for (i = 0; i <= d; i++)
    mpz_init (f[i]);
  return f;
}

/* free a polynomial of degree d */
void
poly_clear (mpz_t *f, int d)
{
  int i;
  
  for (i = 0; i <= d; i++)
    mpz_clear (f[i]);
  free (f);
}

/* f <- g */
void
copy_poly (mpz_t *f, mpz_t *g, int d)
{
  int i;

  for (i = 0; i <= d; i++)
    mpz_set (f[i], g[i]);
}

/* f <- diff(g,x) */
void
diff_poly (mpz_t *f, mpz_t *g, int d)
{
  int i;

  for (i = 0; i < d; i++)
    mpz_mul_ui (f[i], g[i + 1], i + 1);
}

void
fprint_polynomial (FILE *fp, mpz_t *f, const int d)
{
  int i;
  for (i = 0; i <= d; i++)
    if (mpz_cmp_ui (f[i], 0) > 0)
      gmp_fprintf (fp, "+%Zd*x^%d", f[i], i);
    else if (mpz_cmp_ui (f[i], 0) < 0)
      gmp_fprintf (fp, "%Zd*x^%d", f[i], i);
  fprintf (fp, ";\n");
}

#define print_polynomial(f,d) fprint_polynomial(stdout,f,d)

/* D <- |discriminant (f)| = |resultant(f,diff(f,x))/lc(f)| */
void
discriminant (mpz_t D, mpz_t *f0, const int d)
{
  mpz_t *f, *g, num, den;
  int df, dg, i, s;

  f = alloc_poly (d);
  g = alloc_poly (d - 1);
  copy_poly (f, f0, d);
  df = d;
  diff_poly (g, f0, d);
  dg = d - 1;
  mpz_init_set_ui (num, 1);
  mpz_init_set_ui (den, 1);
  while (dg > 0)
    {
      s = df;
      while (df >= dg)
	{
	  /* f <- f * lc(g), except f[df] since we'll divide it afterwards */
	  for (i = 0; i < df; i++)
	    mpz_mul (f[i], f[i], g[dg]);
	  s -= dg;
	  for (i = 0; i < dg; i++)
	    mpz_submul (f[i + df - dg], g[i], f[df]);
	  df --;
	  /* normalize f */
	  while (df > 0 && mpz_cmp_ui (f[df], 0) == 0)
	    df --;
	}
      /* den <- den * lc(g)^deg(f) */
      s -= df;
      if (s > 0)
	for (i = 0; i < s; i++)
	  mpz_mul (num, num, g[dg]);
      else
	for (i = 0; i < -s; i++)
	  mpz_mul (den, den, g[dg]);
      /* swap f and g */
      for (i = 0; i <= dg; i++)
	mpz_swap (f[i], g[i]);
      i = df;
      df = dg;
      dg = i;
    }
  /* num/den*g^deg(f) */
  mpz_pow_ui (g[0], g[0], df);
  mpz_mul (g[0], g[0], num);
  mpz_divexact (g[0], g[0], den);
  /* divide by lc(f0) to get the discriminant */
  mpz_divexact (D, g[0], f0[d]);
  mpz_clear (num);
  mpz_clear (den);
  poly_clear (f, d);
  poly_clear (g, d - 1);
}

/*************************** root properties *********************************/

/* v <- f(x,y) */
void
eval_poly (mpz_t v, mpz_t *f, int d, unsigned long x, unsigned long y)
{
  int i;
  mpz_t Y;

  mpz_init_set_ui (Y, 1);
  mpz_set (v, f[d]);
  for (i = d - 1; i >= 0; i--)
    {
      mpz_mul_ui (v, v, x);
      mpz_mul_ui (Y, Y, y); /* invariant: Y = y^(d-i) */
      mpz_addmul (v, f[i], Y);
    }
  mpz_clear (Y);
}

unsigned long
igcd (unsigned long a, unsigned long b)
{
  unsigned long c;

  while (b != 0)
    {
      c = a % b;
      a = b;
      b = c;
    }
  return a;
}

/* estimate the average valuation of p in f(x) for 0 <= x, y < X */
double
estimate_valuation (mpz_t *f, int d, unsigned long X, unsigned long p)
{
  double s = 0.0, t = 0.0;
  unsigned long x, y;
  mpz_t v;

  mpz_init (v);
  for (x = 1; x <= X; x++)
    for (y = 1; y <= X; y++)
      if (igcd (x, y) == 1)
        {
          eval_poly (v, f, d, x, y);
          if (mpz_cmp_ui (v, 0) != 0)
            {
              t ++;
              while (mpz_divisible_ui_p (v, p))
                {
                  s ++;
                  mpz_divexact_ui (v, v, p);
                }
            }
        }
  mpz_clear (v);
  return s / t;
}

/* Returns the number of roots of f(x) mod p, where f has degree d.
   Assumes p does not divide disc(f). */
int
roots_mod (mpz_t *f, int d, const unsigned long p)
{
  mpz_t *fp, *g, *h;
  int i, j, k, df, dg, dh;

  /* the number of roots is the degree of gcd(x^p-x,f) */

  /* we first compute fp = f/lc(f) mod p */
  while (d >= 0 && mpz_divisible_ui_p (f[d], p))
    d --;
  ASSERT (d >= 0); /* f is 0 mod p: should not happen since otherwise p would
		      divide N, because f(m)=N */
  fp = alloc_poly (d);
  mpz_set_ui (fp[d], p);
  mpz_invert (fp[d], f[d], fp[d]); /* 1/f[d] mod p */
  for (i = 0; i < d; i++)
    {
      mpz_mul (fp[i], f[i], fp[d]);
      mpz_mod_ui (fp[i], fp[i], p);
    }
  mpz_set_ui (fp[d], 1); /* useless? */
  df = d;

  /* we first compute x^p mod fp; since fp has degree d, all operations can
     be done with polynomials of degree < 2d */
  g = alloc_poly (2 * d - 1);
  h = alloc_poly (2 * d - 1);

  /* initialize g to x */
  mpz_set_ui (g[1], 1);
  dg = 1;

  mpz_set_ui (g[2], p);
  k = mpz_sizeinbase (g[2], 2); /* number of bits of p: bit 0 to k-1 */

  for (k -= 2; k >= 0; k--)
    {
      /* square g into h: g has degree dg -> h has degree 2*dg */
      for (i = 0; i <= dg; i++)
	for (j = i + 1; j <= dg; j++)
	  if (i == 0 || j == dg)
	    mpz_mul (h[i + j], g[i], g[j]);
	  else
	    mpz_addmul (h[i + j], g[i], g[j]);
      for (i = 1; i < 2 * dg; i++)
	mpz_mul_2exp (h[i], h[i], 1);
      mpz_mul (h[0], g[0], g[0]);
      mpz_mul (h[2 * dg], g[dg], g[dg]);
      for (i = 1; i < dg; i++)
	mpz_addmul (h[2 * i], g[i], g[i]);
      dh = 2 * dg;

      /* reduce mod p */
      ASSERT (dh < 2 * d);
      for (i = 0; i <= dh; i++)
	mpz_mod_ui (h[i], h[i], p);
      
      /* multiply h by x if bit k of p is set */
      if (p & (1 << k))
	{
	  for (i = dh; i >= 0; i--)
	    mpz_swap (h[i+1], h[i]);
	  mpz_set_ui (h[0], 0);
	  dh ++;
	  ASSERT (dh < 2 * d);
	}

      /* reduce mod fp */
      while (dh >= d)
	{
	  /* subtract h[dh]*fp*x^(dh-d) from h */
	  mpz_mod_ui (h[dh], h[dh], p);
	  for (i = 0; i < d; i++)
	    mpz_submul (h[dh - d + i], h[dh], fp[i]);
	  /* it is not necessary to reduce h[j] for j < dh */
	  dh --;
          ASSERT (dh < 2 * d);
	}

      /* reduce h mod p and copy to g */
      for (i = 0; i <= dh; i++)
	mpz_mod_ui (g[i], h[i], p);
      dg = dh;
      ASSERT (dg < 2 * d);
    }

  /* subtract x */
  mpz_sub_ui (g[1], g[1], 1);
  mpz_mod_ui (g[1], g[1], p);

  while (dg >= 0 && mpz_cmp_ui (g[dg], 0) == 0)
    dg --;

  /* take the gcd with fp */
  while (dg > 0)
    {
      while (df >= dg)
	{
	  /* divide f by g */
	  mpz_set_ui (h[0], p);
	  mpz_invert (h[0], g[dg], h[0]); /* 1/g[dg] mod p */
	  mpz_mul (h[0], h[0], fp[df]);
	  mpz_mod_ui (h[0], h[0], p); /* fp[df]/g[dg] mod p */
	  for (i = 0; i < dg; i++)
	    {
	      mpz_submul (fp[df - dg + i], h[0], g[i]);
	      mpz_mod_ui (fp[df - dg + i], fp[df - dg + i], p);
	    }
	  df --;
	  while (df >= 0 && mpz_cmp_ui (fp[df], 0) == 0)
	    df --;
	}
      /* now d < dg: swap f and g */
      for (i = 0; i <= dg; i++)
	mpz_swap (fp[i], g[i]);
      i = df;
      df = dg;
      dg = i;
    }

  /* if g=0 now, gcd is f, otherwise if g<>0, gcd is 1 */
  if (mpz_cmp_ui (g[0], 0) != 0)
    df = 0;

  poly_clear (g, 2 * d - 1);
  poly_clear (h, 2 * d - 1);
  poly_clear (fp, d);

  return df;
}

int
isprime (unsigned long p)
{
  mpz_t P;
  int res;
  
  mpz_init_set_ui (P, p);
  res = mpz_probab_prime_p (P, REPS);
  mpz_clear (P);
  return res;
}

/* Compute the value alpha(F) from Murphy's thesis, page 49:
   alpha(F) = sum(prime p <= B, (1 - q_p*p/(p+1)) log(p)/(p-1))
   where q_p is the number of roots of F mod p.

   alpha(F) is an estimate of the average logarithm of the part removed
   from sieving, compared to a random integer.

   We want alpha as small as possible, i.e., alpha negative with a large
   absolute value. Typical good values are alpha=-4, -5, ...
*/
double
get_alpha (mpz_t *f, const int d, unsigned long B, int verbose)
{
  double alpha, e;
  unsigned long p, q;
  mpz_t disc;

  mpz_init (disc);
  discriminant (disc, f, d);
  if (mpz_divisible_ui_p (disc, 2) == 0)
    {
      q = roots_mod (f, d, 2);
      ASSERT (q <= 2);
      alpha = (1.0 - 2.0 * (double) q / 3.0) * log (2.0);
    }
  else
    {
      e = estimate_valuation (f, d, NUM_EST, 2);
      if (verbose >= 2)
	fprintf (stderr, "2 divides disc(f): %1.2f\n", e);
      alpha = (1.0 - e) * log (2.0);
    }
  for (p = 3; p <= B; p += 2)
    {
      if (isprime (p))
	{
	  if (mpz_divisible_ui_p (disc, p))
	    {
	      e = estimate_valuation (f, d, NUM_EST, p);
	      if (verbose >= 2)
		fprintf (stderr, "%lu divides disc(f): %1.2f\n", p, e);
	      alpha += (1.0 / (double) (p - 1) - e) * log ((double) p);
	    }
	  else
	    {
	      q = roots_mod (f, d, p);
	      ASSERT (q <= p);
	      alpha += (1.0 - (double) q * (double) p / (double) (p + 1)) *
		log ((double) p) / (double) (p - 1);
	    }
	}
    }
  mpz_clear (disc);
  return alpha;
}

/***************************** input/output **********************************/

static void
init (cado_input in)
{
  mpz_init (in->n);
  in->degree = 0;
}

static void
clear (cado_input in)
{
  mpz_clear (in->n);
}

struct sd {
  size_t s;
  unsigned long d;
};

static struct sd default_degrees[DEFAULT_DEGREES_LENGTH] = DEFAULT_DEGREES;

/* default degree, when not given by user */
static int
default_degree (mpz_t n)
{
  size_t s = mpz_sizeinbase (n, 10);
  int i;

  for (i = 0; s > default_degrees[i].s; i++);
  return default_degrees[i].d;
}

static struct sd default_rlim_table[DEFAULT_RLIM_LENGTH] = DEFAULT_RLIM;
static struct sd default_alim_table[DEFAULT_ALIM_LENGTH] = DEFAULT_ALIM;

/* default degree, when not given by user */
static unsigned long
default_rlim (mpz_t n)
{
  size_t s = mpz_sizeinbase (n, 10);
  int i;

  for (i = 0; s > default_rlim_table[i].s; i++);
  return default_rlim_table[i].d;
}

/* default degree, when not given by user */
static unsigned long
default_alim (mpz_t n)
{
  size_t s = mpz_sizeinbase (n, 10);
  int i;

  for (i = 0; s > default_alim_table[i].s; i++);
  return default_alim_table[i].d;
}

static struct sd default_lpbr_table[DEFAULT_LPBR_LENGTH] = DEFAULT_LPBR;
static struct sd default_lpba_table[DEFAULT_LPBA_LENGTH] = DEFAULT_LPBA;

/* default rational large prime bound, when not given by user */
static unsigned long
default_lpbr (mpz_t n)
{
  size_t s = mpz_sizeinbase (n, 10);
  int i;

  for (i = 0; s > default_lpbr_table[i].s; i++);
  return default_lpbr_table[i].d;
}

/* default algebraic large prime bound, when not given by user */
static unsigned long
default_lpba (mpz_t n)
{
  size_t s = mpz_sizeinbase (n, 10);
  int i;

  for (i = 0; s > default_lpba_table[i].s; i++);
  return default_lpba_table[i].d;
}

static struct sd default_mfbr_table[DEFAULT_MFBR_LENGTH] = DEFAULT_MFBR;
static struct sd default_mfba_table[DEFAULT_MFBA_LENGTH] = DEFAULT_MFBA;

/* default rational factor-residual bound, when not given by user */
static unsigned long
default_mfbr (mpz_t n)
{
  size_t s = mpz_sizeinbase (n, 10);
  int i;

  for (i = 0; s > default_mfbr_table[i].s; i++);
  return default_mfbr_table[i].d;
}

/* default algebraic factor-residual bound, when not given by user */
static unsigned long
default_mfba (mpz_t n)
{
  size_t s = mpz_sizeinbase (n, 10);
  int i;

  for (i = 0; s > default_mfba_table[i].s; i++);
  return default_mfba_table[i].d;
}

struct sf {
  size_t s;
  double f;
};

static struct sf default_rlambda_table[DEFAULT_RLAMBDA_LENGTH] = DEFAULT_RLAMBDA;
static struct sf default_alambda_table[DEFAULT_ALAMBDA_LENGTH] = DEFAULT_ALAMBDA;

/* default rational lambda, when not given by user */
static double
default_rlambda (mpz_t n)
{
  size_t s = mpz_sizeinbase (n, 10);
  int i;

  for (i = 0; s > default_rlambda_table[i].s; i++);
  return default_rlambda_table[i].f;
}

/* default algebraic lambda, when not given by user */
static double
default_alambda (mpz_t n)
{
  size_t s = mpz_sizeinbase (n, 10);
  int i;

  for (i = 0; s > default_alambda_table[i].s; i++);
  return default_alambda_table[i].f;
}

static struct sd default_qint_table[DEFAULT_QINT_LENGTH] = DEFAULT_QINT;

/* default rational factor-residual bound, when not given by user */
static unsigned long
default_qint (mpz_t n)
{
  size_t s = mpz_sizeinbase (n, 10);
  int i;

  for (i = 0; s > default_qint_table[i].s; i++);
  return default_qint_table[i].d;
}

static void
init_poly (cado_poly p, cado_input in)
{
  int i;

  strcpy (p->name, "");
  mpz_init_set (p->n, in->n);
  p->degree = in->degree;
  if (p->degree == 0)
    p->degree = default_degree (p->n);
  p->f = alloc_poly (p->degree);
  p->g = alloc_poly (2);
  mpz_init (p->m);
  strcpy (p->type, "gnfs");
}

static void
clear_poly (cado_poly p)
{
  int i;

  mpz_clear (p->n);
  for (i = 0; i <= p->degree; i++)
    mpz_clear (p->f[i]);
  free (p->f);
  mpz_clear (p->g[0]);
  mpz_clear (p->g[1]);
  free (p->g);
  mpz_clear (p->m);
}

static void
parse_error (int c)
{
  if (c == 0)
    {
      fprintf (stderr, "Error, no field read\n");
      exit (1);
    }
  if (c == EOF)
    {
      fprintf (stderr, "Error, end of file\n");
      exit (1);
    }
}

/* read a string of characters [a-z]+, returns the number of characters read */
static int
read_string (char *s)
{
  int c, l = 0;
  
  while (1)
    {
      c = getchar ();
      if (c == EOF)
	return c;
      if (islower (c) == 0)
	{
	  ungetc (c, stdin);
	  s[l] = '\0';
	  return l;
	}
      s[l++] = c;
    }
}

static void
parse_input (cado_input in)
{
  int c;
  char s[256];

  while ((c = getchar ()) != EOF)
    {
      /* skip blank characters and lines */
      while (isspace (c) && c != EOF)
	c = getchar ();

      if (c == EOF)
	break;
      
      if (c != '#') /* comment line */
	{
	  ungetc (c, stdin);
	  c = read_string (s);
	  if (c == 0 || c == EOF)
	    parse_error (c);

	  /* read possible spaces or tabs */
	  while (isblank (c = getchar ()))
	    printf ("read '%c'\n", c);

	  /* read ':' */
	  if (c != ':')
	    {
	      fprintf (stderr, "Error, ':' expected\n");
	      exit (1);
	    }

	  if (strcmp (s, "n") == 0)
	    {
	      if (mpz_inp_str (in->n, stdin, 10) == 0)
		{
		  fprintf (stderr, "Error after n:\n");
		  exit (1);
		}
	    }
	  else if (strcmp (s, "deg") == 0)
	    {
	      c = scanf ("%d", &(in->degree));
	      if (c == 0 || c == EOF)
		parse_error (c);
	    }
	  else
	    {
	      fprintf (stderr, "Error, unrecognized field: %s\n", s);
	      exit (1);
	    }
	}

      /* read to end of line */
      while (!feof (stdin) && c != '\n')
	c = getchar ();
    }

  if (mpz_cmp_ui (in->n, 0) == 0)
    {
      fprintf (stderr, "Error, no input number read\n");
      exit (1);
    }
}

/* round to nearest, assume d > 0 */
static void
mpz_ndiv_qr (mpz_t q, mpz_t r, mpz_t n, mpz_t d)
{
  int s;

  ASSERT (mpz_cmp_ui (d, 0) != 0);
  mpz_fdiv_qr (q, r, n, d); /* round towards -inf */
  mpz_mul_2exp (r, r, 1);
  s = mpz_cmp (r, d);
  mpz_div_2exp (r, r, 1);
  if (s > 0)
    {
      mpz_add_ui (q, q, 1);
      mpz_sub (r, r, d);
    }
}

/* return a value with the sign of sum(g[i]*s^(2i-d-1), i=0..d) */
double
eval_double_poly (double *g, int d, double s)
{
  double v;
  int i;

  /* sum(g[i]*s^(2*i-d-1), i=0..d) = sum(g[i]*s^(2i), i=0..d)/s^(d+1).
     Since we only need the sign, no need to divide by s^(d+1). */
  for (v = g[d], i = d - 1; i >= 0; i--)
    v = s * v + g[i];
  return v;
}

/* Find the argmin of 
   [|p[d]|*S^d + ... + |p[0]|*S^{-d}] * [m/S + S],
   i.e., a zero of g(S) = (d+1) |p[d]| S^d + ... - (d-1) |p[0]| S^{-d}
   + (d-1) m |p[d]| S^{d-2} + ... - (d+1) m |p[0]| S^{-d-2}.
   Since there is only one sign change in the coefficients of g, it has
   exactly one root on [0, +Inf[. We search it by dichotomy.
*/
double
skewness (mpz_t *p, int d, double m)
{
  double *g, sa, sb, s;
  int i;

  g = (double*) malloc ((d + 2) * sizeof (double));
  /* g[0] stores the coefficient of degree -d-2, ..., 
     g[d+1] that of degree d */
  g[0] = 0.0;
  for (i = 0; i <= d; i++)
    {
      s = fabs (mpz_get_d (p[i]));
      g[i + 1] = (double) (2 * i - d + 1) * s;
      g[i] += (double) (2 * i - d - 1) * (double) m * s;
    }
  if (eval_double_poly (g, d, 1.0) > 0) /* g(1.0) > 0 */
    {
      /* search sa < 1 such that g(sa) <= 0 */
      for (sa = 0.5; eval_double_poly (g, d, sa) > 0; sa *= 0.5);
      sb = 2.0 * sa;
    }
  else /* g(1.0) <= 0 */
    {
      /* search sb such that g(sb) >= 0 */
      for (sb = 2.0; eval_double_poly (g, d, sb) < 0; sb *= 2.0);
      sa = 0.5 * sb;
    }
  /* now g(sa) < 0, g(sb) > 0, and sb = 2*sa */
  for (i = 0; i < 10; i++)
    {
      s = (sa + sb) / 2.0;
      if (eval_double_poly (g, d, s) < 0)
	sa = s;
      else
	sb = s;
    }
  free (g);
  return sqrt ((sa + sb) / 2.0);
}

/* mu(m,S), as defined in reference [2]:
   
   mu(m,S) = (m/S+S)*(|a_d S^d| + ... + |a_0 S^{-d}|)

   It gives an estimate of *value* of the numbers to be sieved (to be
   multiplied by the sieving region).

   We want mu(m,S) as small as possible.
*/
double
get_mu (mpz_t *p, unsigned long degree, mpz_t m, double S)
{
  double mu;
  double S2 = SQR(S);
  int i;

  mu = fabs (mpz_get_d (p[degree]));
  for (i = degree - 1; i >= 0; i--)
    mu = mu * S2 + fabs (mpz_get_d (p[i]));
  for (i = 0; i < degree; i++)
    mu /= S;
  return (mpz_get_d (m) / S + S) * mu;
}

/* Decompose p->n in base m, rounding to nearest. */
void
generate_base_m (cado_poly p, mpz_t m)
{
  unsigned long i;
  int d = p->degree;

  mpz_ndiv_qr (p->f[1], p->f[0], p->n, m);
  for (i = 1; i < d; i++)
    mpz_ndiv_qr (p->f[i+1], p->f[i], p->f[i], m);
  mpz_neg (p->g[0], m);
  mpz_set_ui (p->g[1], 1);
  mpz_set (p->m, m);
  p->skew = skewness (p->f, d, mpz_get_d (m));
}

void
generate_sieving_parameters (cado_poly out)
{
  out->rlim = default_rlim (out->n);
  out->alim = default_alim (out->n);
  out->lpbr = default_lpbr (out->n);
  out->lpba = default_lpba (out->n);
  out->mfbr = default_mfbr (out->n);
  out->mfba = default_mfba (out->n);
  out->rlambda = default_rlambda (out->n);
  out->alambda = default_alambda (out->n);
  out->qintsize = default_qint (out->n);
}

/* Returns a bound to compute alpha(F). In principle we should use the
   algebraic factor base bound alim, but since it is too expensive,
   we use the heuristic of using sqrt(alim). */
unsigned long
get_alpha_bound (cado_poly out)
{
  return (unsigned long) sqrt ((double) default_alim (out->n));
}

#if (METHOD == 3)
/* Method described in reference [2], slide 8
   (Bernstein considered degree d=5, here we generalize to arbitrary d):
   Take m \approx B n^{1/(d+1)}
   Then f_d \approx B^{-d} n^{1/(d+1)}.
   Decompose n in base m, and control f_{d-1} as in Method 2.

   The number of tries T is related to the decrease B by:
   For d=5, T = B^4.5; for d=4, T = B^3.3.
*/
void
generate_poly (cado_poly out, unsigned long T, int verbose)
{
  unsigned long d = out->degree, alim;
  double B, mu, logmu, best_logmu = DBL_MAX, best_E = DBL_MAX, alpha, E;
  double mB;
  mpz_t best_m, k, t, r;
  /* value of T = B^eff[d] for d <= 7 */
  static double eff[] = {0.0, 0.0, 3.0, 3.0, 3.333, 4.5, 6.6, 9.0};
  unsigned long calls_alpha = 0;
  int given_m;

  mpz_init (best_m);
  mpz_init (k);
  mpz_init (t);
  mpz_init (r);
  if (d > 7)
    {
      fprintf (stderr, "Error, too large degree in generate_poly\n");
      exit (1);
    }
  B = pow ((double) T, 1.0 / eff[d]); /* B = T^{1/eff[d]} */
  given_m = mpz_cmp_ui (out->m, 0) != 0;
  if (given_m == 0) /* else use given value of m */
    {
      mpz_root (out->m, out->n, d + 1); /* n^{1/(d+1)} */
      mB = B * mpz_get_d (out->m);
      mpz_set_d (out->m, mB); /* B^{1/(d-1)} n^{1/(d+1)} */
    }
  alim = get_alpha_bound (out); /* needed for alpha */
  /* when going from m to m+k, f[d-1] goes to f[d-1] - k d f[d] */
  while (T-- > 0)
    {
      if (given_m == 0)
        {
          mpz_pow_ui (t, out->m, d - 1); /* m^(d-1) */
          mpz_ndiv_qr (t, r, out->n, t); /* t = round(N / m^(d-1) */
          mpz_fdiv_qr (out->f[d], out->f[d - 1], t, out->m); /* f[d]*m+f[d-1]*/
          if (mpz_cmp_ui (out->f[d], 0) == 0)
            {
              gmp_fprintf (stderr, "Stopping selection at m=%Zd since leading coefficient is zero\n", out->m);
              break;
            }
          mpz_mul_ui (t, out->f[d], d);
          mpz_ndiv_qr (k, r, out->f[d - 1], t);
          ASSERT (mpz_cmp_ui (k, 0) >= 0);
          mpz_add (out->m, out->m, k);
        }
      generate_base_m (out, out->m);
      mu = get_mu (out->f, d, out->m, out->skew);
      logmu = log (mu);
      /* we keep only the values of mu near optimal, since computing alpha
	 is much more expensive */
      if (logmu < best_logmu + MARGIN)
	{
          if (logmu < best_logmu)
            best_logmu = logmu;
	  alpha = get_alpha (out->f, d, alim, verbose);
	  calls_alpha ++;
	  E = logmu + alpha;
	  if (E < best_E)
	    {
	      best_E = E;
	      fprintf (stderr, "m=");
	      mpz_out_str (stderr, 10, out->m);
	      fprintf (stderr, " skew=%1.2f logmu=%1.2f alpha~%1.2f E=%1.2f\n",
		       SQR(out->skew), logmu, alpha, E);
	      mpz_set (best_m, out->m);
	    }
	}
      mpz_add_ui (out->m, out->m, 1);
    }
  if (verbose)
    fprintf (stderr, "calls to get_alpha: %lu\n", calls_alpha);
  generate_base_m (out, best_m);
  generate_sieving_parameters (out);
  mpz_clear (best_m);
  mpz_clear (k);
  mpz_clear (t);
  mpz_clear (r);
}
#endif

/* st is the initial value cputime() at the start of the program */
void
print_poly (FILE *fp, cado_poly p, int argc, char *argv[], int st, int raw)
{
  int i;
  double S, mu, alpha, logmu;

  if (strlen (p->name) != 0)
    fprintf (fp, "name: %s\n", p->name);
  fprintf (fp, "n: ");
  mpz_out_str (fp, 10, p->n);
  fprintf (fp, "\n");
  S = p->skew;
  fprintf (fp, "skew: %1.3f\n", SQR(S));
  mu = get_mu (p->f, p->degree, p->m, S);
  alpha = get_alpha (p->f, p->degree, (p->alim > 1000) ? 1000 : p->alim, 0);
  logmu = log (mu);
  fprintf (fp, "# logmu: %1.2f, alpha: %1.2f logmu+alpha=%1.2f\n", logmu,
	   alpha, logmu + alpha);
  for (i = p->degree; i >= 0; i--)
    {
      fprintf (fp, "c%d: ", i);
      mpz_out_str (fp, 10, p->f[i]);
      fprintf (fp, "\n");
    }
  fprintf (fp, "# ");
  fprint_polynomial (fp, p->f, p->degree);
  fprintf (fp, "Y1: ");
  mpz_out_str (fp, 10, p->g[1]);
  fprintf (fp, "\n");
  fprintf (fp, "Y0: ");
  mpz_out_str (fp, 10, p->g[0]);
  fprintf (fp, "\n");
  fprintf (fp, "# m: ");
  mpz_out_str (fp, 10, p->m);
  fprintf (fp, "\n");
  fprintf (fp, "type: %s\n", p->type);
  if (raw == 0)
    {
      fprintf (fp, "rlim: %lu\n", p->rlim);
      fprintf (fp, "alim: %lu\n", p->alim);
      fprintf (fp, "lpbr: %d\n", p->lpbr);
      fprintf (fp, "lpba: %d\n", p->lpba);
      fprintf (fp, "mfbr: %d\n", p->mfbr);
      fprintf (fp, "mfba: %d\n", p->mfba);
      fprintf (fp, "rlambda: %1.1f\n", p->rlambda);
      fprintf (fp, "alambda: %1.1f\n", p->alambda);
      fprintf (fp, "qintsize: %d\n", p->qintsize);
    }
  fprintf (fp, "# generated by polyselect.r%d", VERSION);
  for (i = 1; i < argc; i++)
    fprintf (fp, " %s", argv[i]);
  fprintf (fp, " in %ds\n", (cputime () - st) / 1000);
}

/* return L_{1/3}(n,c) = exp(c log(n)^{1/3} log(log(n))^{2/3})
   with c = (64/9)^{1/3} ~ 1.922999
*/
double
L (mpz_t n)
{
  double logn, e;
  double c = 1.9229994270765444764;
  
  logn = log (mpz_get_d (n));
  e = log (logn);
  e = logn * SQR(e);
  e = c * pow (e, 0.33333333333333333333);
  return exp (e);
}

/***************************** main program **********************************/

int
main (int argc, char *argv[])
{
  cado_input in;
  cado_poly out;
  int verbose = 0, raw = 0;
  unsigned long effort = 1;
  char **argv0 = argv;
  int argc0 = argc;
  int st = cputime ();
  mpz_t m;

  mpz_init (m);

  /* parse options */
  while (argc > 1 && argv[1][0] == '-')
    {
      if (strcmp (argv[1], "-v") == 0)
	{
	  verbose ++;
	  argv ++;
	  argc --;
	}
      else if (strcmp (argv[1], "-raw") == 0)
	{
	  raw = 1;
	  argv ++;
	  argc --;
	}
      else if (argc > 2 && strcmp (argv[1], "-e") == 0)
	{
	  effort = atoi (argv[2]);
	  argv += 2;
	  argc -= 2;
	}
      else if (argc > 2 && strcmp (argv[1], "-m") == 0)
	{
          mpz_set_str (m, argv[2], 0);
	  argv += 2;
	  argc -= 2;
	}
      else
	{
	  fprintf (stderr, "Unknown option: %s\n", argv[1]);
	  usage ();
	}
    }

  if (argc != 1)
    usage ();

  if (mpz_cmp_ui (m, 0) != 0 && effort != 1)
    {
      fprintf (stderr, "Warning: option -m xxx implies -e 1\n");
      effort = 1;
    }

  init (in);
  parse_input (in);
  if (verbose)
    {
      fprintf (stderr, "Using method %d\n", METHOD);
      fprintf (stderr, "L[1/3,(64/9)^(1/3)](n)=%1.2e\n", L (in->n));
    }
  init_poly (out, in);
  mpz_set (out->m, m); /* if non zero, use given value of m */
  clear (in);

  generate_poly (out, effort, verbose);
  print_poly (stdout, out, argc0, argv0, st, raw);
  clear_poly (out);

  mpz_clear (m);

  return 0;
}
