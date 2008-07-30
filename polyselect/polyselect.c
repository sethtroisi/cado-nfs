/* polyselect - polynomial selection

Author: Paul Zimmermann

Copyright 2007, 2008 INRIA

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
   [4] http://www.crypto-world.com/code/Murphy, Magma code from Scott Contini,
       May 2000.
*/

#include <stdio.h>
#include <stdlib.h>
#include <limits.h> /* for ULONG_MAX */
#include <math.h>   /* for log, fabs */
#include "cado.h"
#include "utils/utils.h" /* for cputime() */
#include "aux.h" /* for common routines with kleinjung.c */

static void
usage ()
{
  fprintf (stderr, "Usage: polyselect [-v] [-full] [-b b] [-degree d] [-e nnn] [-m xxx]< in > out\n\n");
  fprintf (stderr, "       -v        - verbose\n");
  fprintf (stderr, "       -full     - also output factor base parameters\n");
  fprintf (stderr, "       -b b      - linear polynomial bx-m (default 1)\n");
  fprintf (stderr, "       -degree d - use algebraic polynomial of degree d\n");
  fprintf (stderr, "       -e nnn    - use effort nnn\n");
  fprintf (stderr, "       -m xxx    - use decomposition in base xxx\n");
  fprintf (stderr, "       in        - input file (number to factor)\n");
  fprintf (stderr, "       out       - output file (polynomials)\n");
  exit (1);
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

/*****************************************************************************/

typedef struct {
  int alloc;    /* number of allocated coefficients */
  int degree;   /* degree < alloc */
  mpz_t *coeff; /* coefficient list */
} __mpz_poly_struct;
typedef __mpz_poly_struct mpz_poly_t[1];

/* initialize a polynomial with maximal degree d */
void
mpz_poly_init (mpz_poly_t f, int d)
{
  f->alloc = d + 1;
  f->degree = -1; /* initialize to 0 */
  f->coeff = alloc_mpz_array (d + 1);
}

/* clear a polynomial */
void
mpz_poly_clear (mpz_poly_t f)
{
  clear_mpz_array (f->coeff, f->degree + 1);
}

void
mpz_poly_out (FILE *fp, mpz_poly_t f)
{
  if (f->degree < 0)
    fprintf (fp, "0\n");
  else
    fprint_polynomial (fp, f->coeff, f->degree);
}

/* fp <- f/lc(f) mod p. Return degree of fp (-1 if fp=0). */
int
mpz_poly_set_mod (mpz_poly_t fp, mpz_t *f, int d, unsigned long p)
{
  int i;

  while (d >= 0 && mpz_divisible_ui_p (f[d], p))
    d --;
  ASSERT (d >= 0); /* f is 0 mod p: should not happen since otherwise p would
		      divide N, because f(m)=N */
  mpz_set_ui (fp->coeff[d], p);
  mpz_invert (fp->coeff[d], f[d], fp->coeff[d]); /* 1/f[d] mod p */
  for (i = 0; i < d; i++)
    {
      mpz_mul (fp->coeff[i], f[i], fp->coeff[d]);
      mpz_mod_ui (fp->coeff[i], fp->coeff[i], p);
    }
  mpz_set_ui (fp->coeff[d], 1); /* useless? */
  fp->degree = d;

  return d;
}

/* realloc f to n coefficients */
void
mpz_poly_realloc (mpz_poly_t f, int n)
{
  if (f->alloc < n)
    {
      f->coeff = realloc_mpz_array (f->coeff, f->alloc, n);
      f->alloc = n;
    }
}

/* f <- x */
void
mpz_poly_set_x (mpz_poly_t f)
{
  mpz_poly_realloc (f, 2);
  f->degree = 1;
  mpz_set_ui (f->coeff[0], 0);
  mpz_set_ui (f->coeff[1], 1);
}

/* swap f and g */
void
mpz_poly_swap (mpz_poly_t f, mpz_poly_t g)
{
  int i;
  mpz_t *t;

  i = f->alloc;
  f->alloc = g->alloc;
  g->alloc = i;
  i = f->degree;
  f->degree = g->degree;
  g->degree = i;
  t = f->coeff;
  f->coeff = g->coeff;
  g->coeff = t;
}

/* h <- g^2, g and h must differ */
void
mpz_poly_sqr (mpz_poly_t h, mpz_poly_t g)
{
  int i, j, dg = g->degree;
  mpz_t *gc = g->coeff, *hc = h->coeff;

  for (i = 0; i <= dg; i++)
    for (j = i + 1; j <= dg; j++)
      if (i == 0 || j == dg)
	mpz_mul (hc[i + j], gc[i], gc[j]);
      else
	mpz_addmul (hc[i + j], gc[i], gc[j]);
  for (i = 1; i < 2 * dg; i++)
    mpz_mul_2exp (hc[i], hc[i], 1);
  mpz_mul (hc[0], gc[0], gc[0]);
  mpz_mul (hc[2 * dg], gc[dg], gc[dg]);
  for (i = 1; i < dg; i++)
    mpz_addmul (hc[2 * i], gc[i], gc[i]);
  h->degree = 2 * dg;
}

/* normalize h so that h->coeff[deg(h)] <> 0 */
void
mpz_poly_normalize (mpz_poly_t h)
{
  int dh = h->degree;

  while (dh >= 0 && mpz_cmp_ui (h->coeff[dh], 0) == 0)
    dh --;
  h->degree = dh;
}

/* g <- h mod p */
void
mpz_poly_mod_ui (mpz_poly_t g, mpz_poly_t h, unsigned long p)
{
  int i;

  mpz_poly_realloc (g, h->degree + 1);
  for (i = 0; i <= h->degree; i++)
    mpz_mod_ui (g->coeff[i], h->coeff[i], p);
  g->degree = h->degree;
  mpz_poly_normalize (h);
}

/* h <- x*h */
void
mpz_poly_mul_x (mpz_poly_t h)
{
  int i;

  mpz_poly_realloc (h, h->degree + 2);
  for (i = h->degree; i >= 0; i--)
    mpz_swap (h->coeff[i+1], h->coeff[i]);
  mpz_set_ui (h->coeff[0], 0);
  h->degree ++;
}

/* h <- h - x mod p */
void
mpz_poly_sub_x (mpz_poly_t g, unsigned long p)
{
  mpz_poly_realloc (g, 2);
  while (g->degree < 1)
    mpz_set_ui (g->coeff[++(g->degree)], 0);
  mpz_sub_ui (g->coeff[1], g->coeff[1], 1);
  mpz_mod_ui (g->coeff[1], g->coeff[1], p);
  mpz_poly_normalize (g);
}

/* h <- rem(h, f) mod p, f not necessarily monic */
void
mpz_poly_div_r (mpz_poly_t h, mpz_poly_t f, unsigned long p)
{
  int i, d = f->degree, dh = h->degree, monic;
  mpz_t *hc = h->coeff, t;

  mpz_init (t);
  monic = mpz_cmp_ui (f->coeff[d], 1) == 0;
  if (!monic)
    {
      mpz_set_ui (t, p);
      mpz_invert (t, f->coeff[d], t); /* 1/f[d] mod p */
    }
  while (dh >= d)
    {
      /* subtract h[dh]/f[d]*x^(dh-d)*f from h */
      if (!monic)
	{
	  mpz_mul (hc[dh], hc[dh], t);
	  mpz_mod_ui (hc[dh], hc[dh], p);
	}
      for (i = 0; i < d; i++)
	mpz_submul (hc[dh - d + i], hc[dh], f->coeff[i]);
      /* it is not necessary to reduce h[j] mod p here for j < dh */
      do
	{
	  dh --;
	  if (dh >= 0)
	    mpz_mod_ui (hc[dh], hc[dh], p);
	}
      while (dh >= 0 && mpz_cmp_ui (hc[dh], 0) == 0);
    }
  h->degree = dh;
  /* reduce mod p */
  mpz_poly_mod_ui (h, h, p);
  mpz_clear (t);
}

/* fp <- gcd (fp, g) */
void
mpz_poly_gcd (mpz_poly_t fp, mpz_poly_t g, unsigned long p)
{
  while (g->degree >= 0)
    {
      mpz_poly_div_r (fp, g, p);
      mpz_poly_mod_ui (fp, fp, p);
      /* now deg(fp) < deg(g): swap f and g */
      mpz_poly_swap (fp, g);
    }
}

/* Return the number of roots of f(x) mod p, where f has degree d.
*/
int
roots_mod (mpz_t *f, int d, const long p)
{
  mpz_poly_t fp, g, h;
  int k, df;

  /* the number of roots is the degree of gcd(x^p-x,f) */

  mpz_poly_init (fp, d);
  mpz_poly_init (g, 2 * d - 1);
  mpz_poly_init (h, 2 * d - 1);

  d = mpz_poly_set_mod (fp, f, d, p);

  /* we first compute x^p mod fp; since fp has degree d, all operations can
     be done with polynomials of degree < 2d: a square give at most degree
     2d-2, and multiplication by x gives 2d-1. */

  /* initialize g to x */
  mpz_poly_set_x (g);

  k = nbits (p);

  for (k -= 2; k >= 0; k--)
    {
      mpz_poly_sqr (h, g);             /* h <- g^2 */
      if (p & (1 << k))
	mpz_poly_mul_x (h);            /* h <- x*h */

      mpz_poly_div_r (h, fp, p);       /* h -> rem(h, fp) */
      mpz_poly_mod_ui (g, h, p);       /* g <- h mod p */
    }

  /* subtract x */
  mpz_poly_sub_x (g, p);

  /* fp <- gcd (fp, g) */
  mpz_poly_gcd (fp, g, p);

  /* now fp contains gcd(x^p-x, f) */

  df = fp->degree;

  mpz_poly_clear (fp);
  mpz_poly_clear (g);
  mpz_poly_clear (h);

  printf ("p=%ld:%d\n", p, df);

  return df;
}

/*****************************************************************************/

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

/*****************************************************************************/

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
}

/* Returns k such that f(x+k) has the smallest 1-norm, with the corresponding
   skewness.
   The coefficients f[] are modified to those of f(x+k).

   The linear polynomial b*x-m is changed into b*(x+k)-m, thus m is
   changed into m-k*b.
*/
long
translate (mpz_t *f, int d, mpz_t *g, mpz_t m, mpz_t b)
{
  double logmu0, logmu;
  int i, j, dir;
  long k;

  logmu0 = LOGNORM (f, d, SKEWNESS (f, d, SKEWNESS_DEFAULT_PREC));

  /* first compute f(x+1) */
  /* define f_i(x) = f[i] + f[i+1]*x + ... + f[d]*x^(d-i)
     then f_i(x+1) = f[i] + (x+1)*f_{i+1}(x+1) */
  for (i = d - 1; i >= 0; i--)
    {
      /* invariant: f[i+1], ..., f[d] are the coefficients of f_{i+1}(x+1),
         thus we have to do: f[i] <- f[i] + f[i+1], f[i+1] <- f[i+1] + f[i+2],
         ..., f[d-1] <- f[d-1] + f[d] */
      for (j = i; j < d; j++)
        mpz_add (f[j], f[j], f[j+1]);
    }
  mpz_sub (m, m, b);
  k = 1;

  /* invariant: the coefficients are those of f(x+k) */
  logmu = LOGNORM (f, d, SKEWNESS (f, d, SKEWNESS_DEFAULT_PREC));

  if (logmu < logmu0)
    dir = 1;
  else
    dir = -1;
  logmu0 = logmu;

  while (1)
    {
      /* translate from f(x+k) to f(x+k+dir) */
      for (i = d - 1; i >= 0; i--)
        for (j = i; j < d; j++)
          if (dir == 1)
            mpz_add (f[j], f[j], f[j+1]);
          else
            mpz_sub (f[j], f[j], f[j+1]);
      if (dir == 1)
        mpz_sub (m, m, b);
      else
        mpz_add (m, m, b);
      k = k + dir;
      logmu = LOGNORM (f, d, SKEWNESS (f, d, SKEWNESS_DEFAULT_PREC));
      if (logmu < logmu0)
        logmu0 = logmu;
      else
        break;
    }

  /* go back one step */
  for (i = d - 1; i >= 0; i--)
    for (j = i; j < d; j++)
      if (dir == 1)
        mpz_sub (f[j], f[j], f[j+1]);
      else
        mpz_add (f[j], f[j], f[j+1]);
  if (dir == 1)
    mpz_add (m, m, b);
  else
    mpz_sub (m, m, b);
  k = k - dir;

  /* change linear polynomial */
  mpz_neg (g[0], m);

  return k;
}

/* Method described in reference [2], slide 8
   (Bernstein considered degree d=5, here we generalize to arbitrary d):
   Take m \approx B n^{1/(d+1)}
   Then f_d \approx B^{-d} n^{1/(d+1)}.
   Decompose n in base m, and control f_{d-1} as in Method 2.

   The number of tries T is related to the decrease B by:
   For d=5, T = B^4.5; for d=4, T = B^3.3.

   b is the leading coefficient of the linear polynomial (default is 1).
*/
void
generate_poly (cado_poly out, double T, mpz_t b)
{
  unsigned long d = out->degree;
  double B, logmu, E, mB;
  double best_E = DBL_MAX;
  mpz_t best_m, k, t, r;
  /* value of T = B^eff[d] for d <= 7 */
  static double eff[] = {0.0, 0.0, 3.0, 3.0, 3.333, 4.5, 6.6, 9.0};
  int given_m;
  m_logmu_t *M; /* stores values of m, b and log(mu) found in 1st phase */
  unsigned long Malloc; /* number of allocated entries in M */
  unsigned long Msize;  /* number of stored entries in M */
  unsigned long i;
  double st;
  unsigned long Malloc2, Msize2;
  
  ASSERT (T <= 9007199254740992.0); /* if T > 2^53, T-- will loop */

  mpz_init (best_m);

  given_m = mpz_cmp_ui (out->m, 0) != 0;
  if (given_m) /* use given value of m */
    {
      mpz_set (best_m, out->m);
      goto end;
    }
  
  Malloc = (unsigned long) pow (T, 0.4);
  M = m_logmu_init (Malloc);
  mpz_init (k);
  mpz_init (t);
  mpz_init (r);

  if (d > 7)
    {
      fprintf (stderr, "Error, too large degree in generate_poly\n");
      exit (1);
    }
  B = pow (T, 1.0 / eff[d]); /* B = T^{1/eff[d]} */
  mpz_root (out->m, out->n, d + 1); /* n^{1/(d+1)} */
  mB = B * mpz_get_d (out->m); /* B^{1/(d-1)} n^{1/(d+1)} */
  mpz_set_d (out->m, mB);

  /* First phase: search for m which give a base-m decomposition with a small
     norm. When going from m to m+k, f[d-1] goes to f[d-1] - k d f[d] */
  st = seconds ();
  Msize = 0;
  while (T > 0)
    {
      mpz_pow_ui (t, out->m, d - 1); /* m^(d-1) */
      mpz_ndiv_qr (t, r, out->n, t); /* t = round (N / m^(d-1)) */
      /* we use mpz_fdiv_qr here (round towards -infinity) to ensure
         f[d-1] >= 0, which in turn ensures f[d-1]/(d f[d]) >= 0 */
      mpz_fdiv_qr (out->f[d], out->f[d - 1], t, out->m); /* t=f[d]*m+f[d-1] */
      if (mpz_cmp_ui (out->f[d], 0) == 0)
        {
          gmp_fprintf (stderr, "# Stopping selection at m=%Zd since leading coefficient is zero\n", out->m);
          break;
        }
      mpz_mul_ui (t, out->f[d], d);
      mpz_fdiv_q (k, out->f[d - 1], t);
      mpz_add (out->m, out->m, k);
      T--; /* T counts the number of calls to generate_base_mb */
      generate_base_mb (out, out->m, b);
      logmu = LOGNORM (out->f, d, out->skew);
      m_logmu_insert (M, Malloc, &Msize, b, out->m, logmu, "lognorm=", 0);
      if (T > 0 && mpz_sgn (out->f[d - 1]) > 0)
        { /* try another small f[d-1], avoiding a mpz_pow_ui call */
          mpz_add_ui (out->m, out->m, 1);
          T--;
          generate_base_mb (out, out->m, b);
          logmu = LOGNORM (out->f, d, out->skew);
          m_logmu_insert (M, Malloc, &Msize, b, out->m, logmu, "lognorm=", 0);
        }
      mpz_add_ui (out->m, out->m, 1);
    }
  fprintf (stderr, "# First phase took %.2fs and kept %lu polynomial(s)\n",
           seconds () - st, Msize);
  fflush (stderr);

  /* Second/third phases: loop over entries in M database, and try to find the
     best rotation for each one. In principle we should compute the
     contribution to alpha(F) for primes up to the algebraic factor base bound,
     but since it is too expensive, we use a fixed bound.
     To speed up things we first compute an approximation of alpha(F) by
     considering small primes only, and keep only the best polynomials, for
     which we compute an approximation of alpha(F) up to a larger prime bound.
  */
  Malloc2 = (unsigned long) sqrt ((double) Malloc);
  Msize2 = 0;
  st = seconds ();
  for (i = 0; i < Msize; i++)
    {
      generate_base_mb (out, M[i].m, M[i].b);
      /* we do not use translation here, since it has little effect on the
         norm, and moreover it does not permute with the base-m generation,
         thus we would need to save M[i].m, and redo all steps in the same
         order below */
      E = rotate (out->f, d, ALPHA_BOUND_SMALL, M[i].m, M[i].b, 0);
      m_logmu_insert (M, Malloc2, &Msize2, M[i].b, M[i].m, E, "E~", 0);
    }
  fprintf (stderr, "# Second phase took %.2fs and kept %lu polynomial(s)\n",
           seconds () - st, Msize2);
  fflush (stderr);

  st = seconds ();
  for (i = 0; i < Msize2; i++)
    {
      generate_base_mb (out, M[i].m, M[i].b);
      E = rotate (out->f, d, ALPHA_BOUND, M[i].m, M[i].b, 0);
      if (E < best_E)
        {
          best_E = E;
          gmp_fprintf (stderr, "# p=1 m=%Zd E=%1.2f\n", M[i].m, E);
          mpz_set (best_m, M[i].m);
        }
    }
  fprintf (stderr, "# Third phase took %.2fs\n", seconds () - st);
  fflush (stderr);

  mpz_clear (k);
  mpz_clear (t);
  mpz_clear (r);
  m_logmu_clear (M, Malloc);

 end:
  /* Warning: we must restart from M[i].m here, and follow exactly the same
     steps as above, to find the same polynomial. In particular, the base-m
     generation does not permute with the translation. */
  generate_base_mb (out, best_m, b);
  /* rotate */
  rotate (out->f, d, ALPHA_BOUND, out->m, b, 1);
  /* translate */
  translate (out->f, d, out->g, out->m, b);
  /* recompute skewness, since it might differ from the original
     polynomial with no rotation or translation */
  out->skew = SKEWNESS (out->f, out->degree, 20);
  generate_sieving_parameters (out);
  mpz_clear (best_m);
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
  e = logn * (e * e);
  e = c * pow (e, 0.33333333333333333333);
  return exp (e);
}

/***************************** main program **********************************/

int
main (int argc, char *argv[])
{
  param_list pl;
  cado_poly poly;
  int verbose = 0, raw = 1, i;
  double effort = 1.0;
  char **argv0 = argv;
  int argc0 = argc;
  double st = seconds ();
  mpz_t b;
  FILE *f;

  fprintf (stderr, "# %s.r%s", argv[0], REV);
  for (i = 1; i < argc; i++)
    fprintf (stderr, " %s", argv[i]);
  fprintf (stderr, "\n");

  param_list_init(pl);
  cado_poly_init(poly);
  mpz_init_set_ui (b, 1); /* leading coefficient of linear polynomial */

  argv++, argc--;
  for( ; argc ; ) {
      /* knobs first */
      if (strcmp(argv[0], "-v") == 0) { verbose++; argv++,argc--; continue; }
      if (strcmp(argv[0], "-full") == 0) { raw=0; argv++,argc--; continue; }
      /* Then aliases */
      if (param_list_update_cmdline_alias(pl, "effort", "-e", &argc, &argv))
          continue;
      /* Pick just everything from the rest that looks like a parameter */
      if (param_list_update_cmdline(pl, NULL, &argc, &argv)) { continue; }

      /* Now last resort measures */
      if (strspn(argv[0], "0123456789") == strlen(argv[0])) {
          param_list_add_key(pl, "n", argv[0], PARAMETER_FROM_CMDLINE);
          argv++,argc--;
          continue;
      }
      /* If something remains, then it could be an input file */
      if ((f = fopen(argv[0], "r")) != NULL) {
          param_list_read_stream(pl, f);
          fclose(f);
          argv++,argc--;
          continue;
      }
      /* bail out */
      fprintf(stderr, "Unhandled parameter %s\n", argv[0]);
      usage();
  }

  int have_n = param_list_parse_mpz(pl, "n", poly->n);

  if (!have_n) {
      if (verbose) {
          fprintf(stderr, "Reading n from stdin\n");
      }
      param_list_read_stream(pl, stdin);
      have_n = param_list_parse_mpz(pl, "n", poly->n);
  }

  if (!have_n) {
      fprintf(stderr, "No n defined ; sorry.\n");
      exit(1);
  }

  param_list_parse_double(pl, "effort", &effort);
  param_list_parse_mpz(pl, "m", poly->m);
  param_list_parse_int(pl, "degree", &(poly->degree));
  param_list_parse_mpz(pl, "b", b);

  if (param_list_warn_unused(pl)) {
      usage();
  }
  param_list_clear(pl);

  if (mpz_cmp_ui (poly->m, 0) != 0 && effort != 1.0)
    {
      fprintf (stderr, "Warning: option -m xxx implies -e 1\n");
      effort = 1.0;
    }
  
  /* check b >= 1 */
  if (mpz_cmp_ui (b, 1) < 0)
    {
      fprintf (stderr, "Error, b should be greater or equal to 1\n");
      exit (1);
    }

  if (verbose)
    fprintf (stderr, "L[1/3,(64/9)^(1/3)](n)=%1.2e\n", L (poly->n));

  if (poly->degree <= 0)
    poly->degree = default_degree (poly->n);
  poly->name[0] = '\0';
  strncpy (poly->type, "gnfs", sizeof (poly->type));

  generate_poly (poly, effort, b);

  print_poly (stdout, poly, argc0, argv0, st, raw);

  cado_poly_clear (poly);
  mpz_clear (b);

  return 0;
}
