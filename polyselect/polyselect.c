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
#include <float.h>  /* for DBL_MAX */
#include <math.h>   /* for log, fabs */
#include "cado.h"
#include "utils/utils.h" /* for cputime() */

#include "aux.c"

/* define LOGNORM as l1_lognorm, SKEWNESS as l1_skewness to get the 1-norm
   on coefficients */
#define LOGNORM  L2_lognorm
#define SKEWNESS L2_skewness

#define SKEWNESS_DEFAULT_PREC 10

#define MAX_ROTATE 16 /* try rotation by k*(b*x-m) for |k| <= MAX_ROTATE.
                         FIXME: we should adjust that bound dynamically,
                         so that the norm of the modified polynomial does
                         not change much. */

#define mpz_add_si(a,b,c)                       \
  if (c >= 0) mpz_add_ui (a, b, c);             \
  else mpz_sub_ui (a, b, -(c))

#define mpz_submul_si(a,b,c)                    \
  if (c >= 0) mpz_submul_ui (a, b, c);          \
  else mpz_addmul_ui (a, b, -(c))
  

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

/************************** polynomial arithmetic ****************************/

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

/* g <- content(f) where deg(f)=d */
void
content_poly (mpz_t g, mpz_t *f, int d)
{
  int i;

  ASSERT(d >= 1);

  mpz_gcd (g, f[0], f[1]);
  for (i = 2; i <= d; i++)
    mpz_gcd (g, g, f[i]);
}

/* v <- f(r), where f is of degree d */
void
eval_poly_ui (mpz_t v, mpz_t *f, int d, unsigned long r)
{
  int i;

  mpz_set (v, f[d]);
  for (i = d - 1; i >= 0; i--)
    {
      mpz_mul_ui (v, v, r);
      mpz_add (v, v, f[i]);
    }
}

/* v <- f'(r), where f is of degree d */
void
eval_poly_diff_ui (mpz_t v, mpz_t *f, int d, unsigned long r)
{
  int i;

  mpz_mul_ui (v, f[d], d);
  for (i = d - 1; i >= 1; i--)
    {
      mpz_mul_ui (v, v, r);
      mpz_addmul_ui (v, f[i], i); /* v <- v + i*f[i] */
    }
}

/* h(x) <- h(x + r/p), where the coefficients of h(x + r/p) are known to
   be integers */
void
poly_shift_divp (mpz_t *h, int d, long r, unsigned long p)
{
  int i, k;
  mpz_t t;

  mpz_init (t);
  for (i = 1; i <= d; i++)
    for (k = d - i; k < d; k++)
      { /* h[k] <- h[k] + r/p h[k+1] */
        ASSERT (mpz_divisible_ui_p (h[k+1], p) != 0);
        mpz_divexact_ui (t, h[k+1], p);
	if (r >= 0)
	  mpz_addmul_ui (h[k], t, r);
	else
	  mpz_submul_ui (h[k], t, -r);
      }
  mpz_clear (t);
}

/* D <- |discriminant (f)| = |resultant(f,diff(f,x))/lc(f)| */
void
discriminant (mpz_t D, mpz_t *f0, const int d)
{
  mpz_t *f, *g, num, den;
  int df, dg, i, s;

  f = alloc_mpz_array (d + 1);
  g = alloc_mpz_array (d);
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
  clear_mpz_array (f, d + 1);
  clear_mpz_array (g, d);
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

/* auxiliary routine for special_valuation(), see below */
double
special_val0 (mpz_t *f, int d, unsigned long p)
{
  double v;
  mpz_t c, *g, *h;
  int alloc, i, r0, nroots;
  unsigned long *roots, r;

  mpz_init (c);
  content_poly (c, f, d);
  for (v = 0.0; mpz_divisible_ui_p (c, p); v++, mpz_divexact_ui (c, c, p));
  alloc = v != 0.0;

  /* g <- f/p^v */
  if (alloc != 0)
    {
      g = alloc_mpz_array (d + 1);
      mpz_ui_pow_ui (c, p, (unsigned long) v); /* p^v */
      for (i = 0; i <= d; i++)
        mpz_divexact (g[i], f[i], c);
    }
  else
    g = f;

  h = alloc_mpz_array (d + 1);
  /* first compute h(x) = g(px) */
  mpz_set_ui (c, 1);
  for (i = 0; i <= d; i++)
    {
      mpz_mul (h[i], g[i], c);
      mpz_mul_ui (c, c, p);
    }
  /* Search for roots of g mod p */
  ASSERT (d > 0);
  roots = (unsigned long*) malloc (d * sizeof (unsigned long));
  FATAL_ERROR_CHECK(roots == NULL, "not enough memory");


  nroots = poly_roots_ulong (roots, g, d, p);
  ASSERT (nroots <= d);
  for (r0 = 0, i = 0; i < nroots; i++)
    {
      r = roots[i];
      eval_poly_diff_ui (c, g, d, r);
      if (mpz_divisible_ui_p (c, p) == 0) /* g'(r) <> 0 mod p */
	v += 1.0 / (double) (p - 1);
      else /* hard case */
	{
	  /* g(px+r) = h(x + r/p), thus we can go from h0(x)=g(px+r0)
	     to h1(x)=g(px+r1) by computing h0(x + (r1-r0)/p) */
	  poly_shift_divp (h, d, r - r0, p);
	  r0 = r;
	  v += special_val0 (h, d, p) / (double) p;
	}
    }
  free (roots);
  clear_mpz_array (h, d + 1);

  if (alloc != 0)
    clear_mpz_array (g, d + 1);
  mpz_clear (c);

  return v;
}

/* Compute the average valuation of F(a,b) for gcd(a,b)=1, for a prime p
   dividing the discriminant of f, using the following algorithm from 
   Guillaume Hanrot (which is some kind of p-adic variant of Uspensky's
   algorithm):

   val(f, p)
     return val0(f, p) * p / (p+1) + val0(f(1/(p*x))*(p*x)^d, p) * 1/(p+1)

   val0(f, p). 
     v <- valuation (content(f), p); 
     f <- f/p^v  

     r <- roots mod p(f, p)

     for r_i in r do
         if f'(r_i) <> 0 mod p then v +=  1/(p-1). 
         else
              f2 <- f(p*x + r_i)
              v += val0(f2, p) / p. 
         endif
     endfor
     Return v. 

A special case when:
(a) p^2 does not divide disc(f),
(b) p does not divide lc(f),
then the average valuation is (p q_p - 1)/(p^2 - 1), where q_p is the number
of roots of f mod p. When q_p=1, we get 1/(p+1).

Note: when p does not divide lc(f), the val0(f(1/(p*x))*(p*x)^d, p) call
always returns 0 in val(f,p).

Assumes p divides disc = disc(f), d is the degree of f.
*/
double special_valuation(mpz_t * f, int d, unsigned long p, mpz_t disc)
{
    double v;
    int p_divides_lc;
    int pvaluation_disc = 0;
    double pd = (double) p;

    if (mpz_divisible_ui_p(disc, p)) {
	mpz_t t;
	pvaluation_disc++;
	mpz_init(t);
	mpz_divexact_ui(t, disc, p);
	if (mpz_divisible_ui_p(t, p))
	    pvaluation_disc++;
	mpz_clear(t);
    }

    p_divides_lc = mpz_divisible_ui_p(f[d], p);

    if (pvaluation_disc == 0) {
	/* easy ! */
	int e;
	e = poly_roots_ulong(NULL, f, d, p);
	if (p_divides_lc) {
	    /* Or the discriminant would have valuation 1 at least */
	    ASSERT(mpz_divisible_ui_p(f[d - 1], p) == 0);
	    e++;
	}
	return (pd * e) / (pd * pd - 1);
    } else if (pvaluation_disc == 1 && (p_divides_lc == 0)) {
	/* special case where p^2 does not divide disc and p does not
	   divide lc(f) */
	int e;
	e = poly_roots_ulong(NULL, f, d, p);
	/* something special here. */
	return (pd * e - 1) / (pd * pd - 1);
    } else {
	v = special_val0(f, d, p) * (double) p;
	if (p_divides_lc) {
	    /* compute g(x) = f(1/(px))*(px)^d, i.e., g[i] = f[d-i]*p^i */
	    /* IOW, the reciprocal polynomial evaluated at px */
	    mpz_t *g;
	    mpz_t t;
	    int i;

	    g = alloc_mpz_array(d + 1);
	    mpz_init_set_ui(t, 1);	/* will contains p^i */
	    for (i = 0; i <= d; i++) {
		mpz_mul(g[i], f[d - i], t);
		mpz_mul_ui(t, t, p);
	    }
	    v += special_val0(g, d, p);
	    clear_mpz_array(g, d + 1);
            mpz_clear(t);
	}
	v /= (double) (p + 1);
    }
    return v;
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

/* Compute the value alpha(F) from Murphy's thesis, page 49:
   alpha(F) = sum(prime p <= B, (1 - q_p*p/(p+1)) log(p)/(p-1))
   where q_p is the number of roots of F mod p, including the number of
   projective roots (i.e., the zeros of the reciprocal polynomial mod p).

   alpha(F) is an estimate of the average logarithm of the part removed
   from sieving, compared to a random integer.

   We want alpha as small as possible, i.e., alpha negative with a large
   absolute value. Typical good values are alpha=-4, -5, ...
*/
double
get_alpha (mpz_t *f, const int d, unsigned long B)
{
  double alpha, e;
  unsigned long p;
  mpz_t disc;

  mpz_init (disc);
  discriminant (disc, f, d);

  /* special_valuation returns the expected average exponent of p in F(a,b)
     for coprime a, b, i.e., e = q_p*p/(p^2-1), thus the contribution for p
     is (1/(p-1) - e) * log(p) */

  /* prime p=2 */
  e = special_valuation (f, d, 2, disc);
  alpha = (1.0 - e) * log (2.0);

  /* FIXME: generate all primes up to B and pass them to get_alpha */
  for (p = 3; p <= B; p += 2)
    if (isprime (p))
      {
	e = special_valuation (f, d, p, disc);
	alpha += (1.0 / (double) (p - 1) - e) * log ((double) p);
      }
  mpz_clear (disc);
  return alpha;
}

/********************* data structures for first phase ***********************/

typedef struct {
  mpz_t m;
  double logmu;
} m_logmu_t;

m_logmu_t*
m_logmu_init (unsigned long Malloc)
{
  unsigned long i;
  m_logmu_t* M;

  M = (m_logmu_t*) malloc (Malloc * sizeof (m_logmu_t));
  for (i = 0; i < Malloc; i++)
    mpz_init (M[i].m);
  return M;
}

void
m_logmu_clear (m_logmu_t* M, unsigned long Malloc)
{
  unsigned long i;

  for (i = 0; i < Malloc; i++)
    mpz_clear (M[i].m);
  free (M);
}

/* Insert (m, logmu) in the M database. Current implementation of M is
   a sorted list, with increasing values of logmu.
   alloc: number of allocated entries in M.
   size:  number of stored entries in M (size <= alloc).
   Returns the new size = min (size + 1, alloc).
*/
unsigned long
m_logmu_insert (m_logmu_t* M, unsigned long alloc, unsigned long size,
                mpz_t m, double logmu)
{
  ASSERT(size <= alloc);

  if (size < alloc)
    {
      mpz_set (M[size].m, m);
      M[size].logmu = logmu;
      return size + 1;
    }
  else /* size=alloc: database is full, remove entry with smallest logmu */
    {
      if (logmu < M[alloc - 1].logmu)
        {
          size --;
          while (size > 0 && logmu < M[size - 1].logmu)
            {
              mpz_swap (M[size - 1].m, M[size].m);
              M[size].logmu = M[size - 1].logmu;
              size --;
            }
          /* now either size = 0 or M[size - 1].logmu <= logmu */
          mpz_set (M[size].m, m);
          M[size].logmu = logmu;
          if (size == 0)
            gmp_fprintf (stderr, "# m=%Zd logmu=%1.2f\n", m, logmu);
        }
      return alloc;
    }
}

/*****************************************************************************/

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

static void
init_poly (cado_poly p, cado_input in)
{
  /* strcpy (p->name, ""); */
  p->name[0] = '\0';
  mpz_init_set (p->n, in->n);
  p->degree = in->degree;
  if (p->degree == 0)
    p->degree = default_degree (p->n);
  p->f = alloc_mpz_array (p->degree + 1);
  p->g = alloc_mpz_array (3);
  mpz_init (p->m);
  strncpy (p->type, "gnfs", sizeof(p->type));
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

/************************* norm and skewness *********************************/

/* Return an approximation of the value of S that minimizes the expression:

                |p[d]| * S^{d/2} + ... + |p[0]| * S^{-d/2}

   S is a zero of g(S) = d |p[d]| S^{d-1} + ... - d |p[0]| S^{-d-1},
   thus S^2 is a zero of d |p[d]| x^d + ... - d |p[0]|.

   This correspond to a sieving region |a| <= S^{1/2} R, and |b| <= R/S^{1/2}.
   Note: Bernstein [2] defines the 'skewness' as S, i.e., the square root of
   the returned value, with |a| <= SR and |b| <= R/S. We prefer here to follow
   the convention from [1], which gives sieving bound Amax = S * Bmax.

   Since there is only one sign change in the coefficients of g, it has
   exactly one root on [0, +Inf[. We search it by dichotomy.
   There is one difference with respect to Definition 3.1 in [1]:
   we consider the L1 norm instead of the Linf norm.

   prec is the relative precision of the result.
*/
double
l1_skewness (mpz_t *p, int d, int prec)
{
  double *g, sa, sb, s;
  int i;

  g = (double*) malloc ((d + 1) * sizeof (double));
  /* g[0] stores the coefficient of degree -d-1 of g(S), ..., 
     g[d] that of degree d-1 */
  for (i = 0; i <= d; i++)
    g[i] = (double) (2 * i - d) * fabs (mpz_get_d (p[i]));
  if (fpoly_eval (g, d, 1.0) > 0) /* g(1.0) > 0 */
    {
      /* search sa < 1 such that g(sa) <= 0 */
      for (sa = 0.5; fpoly_eval (g, d, sa) > 0; sa *= 0.5);
      sb = 2.0 * sa;
    }
  else /* g(1.0) <= 0 */
    {
      /* search sb such that g(sb) >= 0 */
      for (sb = 2.0; fpoly_eval (g, d, sb) < 0; sb *= 2.0);
      sa = 0.5 * sb;
    }
  /* now g(sa) < 0, g(sb) > 0, and sb = 2*sa.
     We use 10 iterations only, which gives a relative precision of
     2^(-10) ~ 0.001, which is enough for our needs. */
  s = fpoly_dichotomy (g, d, sa, sb, -1.0, prec);
  free (g);
  return s;
}

/* Returns the logarithm of the 1-norm related to the coefficients,
   log(mu(f,S)), with mu(f,S) defined as follows:
   
   mu(f,S) = |a_d S^{d/2}| + ... + |a_0 S^{-d/2}|

   It gives an estimate of *value* of the numbers to be sieved (to be
   multiplied by the sieving region).

   We want mu(f,S) as small as possible.
*/
double
l1_lognorm (mpz_t *p, unsigned long degree, double S)
{
  double mu;
  long i;

  mu = fabs (mpz_get_d (p[degree]));
  for (i = degree - 1; i >= 0; i--)
    mu = mu * S + fabs (mpz_get_d (p[i]));
  /* divide by S^(d/2) */
  for (i = 0; (unsigned long) i < degree; i += 2)
    mu /= S;
  /* If the degree is even, say d = 2k, we divided by S for i=0, 2, ..., 2k-2,
     thus k = d/2 times.
     If the degree is odd, say d = 2k+1, we divided by S for i=0, 2, ..., 2k,
     thus k+1 times = (d+1)/2, thus we have to multiply back by sqrt(S). */
  if ((unsigned long) i > degree)
    mu *= sqrt (S);
  return log (mu);
}

/* same as L2_lognorm, but takes 'double' instead of 'mpz_t' as coefficients */
double
L2_lognorm_d (double *a, unsigned long d, double s)
{
  double n;

  if (d == 4)
    {
      double a4, a3, a2, a1, a0;

      a4 = a[4] * s * s * s * s;
      a3 = a[3] * s * s * s;
      a2 = a[2] * s * s;
      a1 = a[1] * s;
      a0 = a[0];
      n = 4.0 / 9.0 * (a4 * a4 + a0 * a0) + 8.0 / 21.0 * (a4 * a2 + a2 * a0)
        + 4.0 / 21.0 * (a3 * a3 + a1 * a1) + 8.0 / 25.0 * (a4 * a0 + a3 * a1)
        + 4.0 / 25.0 * a2 * a2;
      return 0.5 * log(n / (s * s * s * s));
    }
  else if (d == 5)
    {
      double a5, a4, a3, a2, a1, a0;

      a5 = a[5] * s * s * s * s * s;
      a4 = a[4] * s * s * s * s;
      a3 = a[3] * s * s * s;
      a2 = a[2] * s * s;
      a1 = a[1] * s;
      a0 = a[0];
      n = 4.0 / 11.0 * (a5 * a5 + a0 * a0) + 8.0 / 27.0 * (a5 * a3 + a2 * a0)
        + 4.0 / 27.0 * (a4 * a4 + a1 * a1) + 4.0 / 35.0 * (a3 * a3 + a2 * a2)
        + 8.0 / 35.0 * (a5 * a1 + a4 * a2 + a4 * a0 + a3 * a1);
      return 0.5 * log(n / (s * s * s * s * s));
    }
  else if (d == 6)
    {
      double a6, a5, a4, a3, a2, a1, a0;

      a6 = a[6] * s * s * s * s * s * s;
      a5 = a[5] * s * s * s * s * s;
      a4 = a[4] * s * s * s * s;
      a3 = a[3] * s * s * s;
      a2 = a[2] * s * s;
      a1 = a[1] * s;
      a0 = a[0];
      n = 4.0 / 13.0 * (a6 * a6 + a0 * a0) + 8.0 / 33.0 * (a6 * a4 + a2 * a0)
        + 4.0 / 33.0 * (a5 * a5 + a1 * a1) + 4.0 / 45.0 * (a4 * a4 + a2 * a2)
        + 8.0 / 45.0 * (a6 * a2 + a5 * a3 + a4 * a0 + a3 * a1)
        + 8.0 / 49.0 * (a6 * a0 + a5 * a1 + a4 * a2) + 4.0 / 49.0 * a3 * a3;
      return 0.5 * log(n / (s * s * s * s * s * s));
    }
  else
    {
      fprintf (stderr, "L2norm not yet implemented for degree %lu\n", d);
      exit (1);
    }
}

/* Returns the logarithm of the L2-norm as defined by Kleinjung, i.e.,
   log(1/2 sqrt(int(int((F(sx,y)/s^(d/2))^2, x=-1..1), y=-1..1))).
   Since we only want to compare norms, we don't consider the log(1/2) term,
   and compute only 1/2 log(int(int(...))) [here the 1/2 factor is important,
   since it is added to the alpha root property term].
*/
double
L2_lognorm (mpz_t *f, unsigned long d, double s)
{
  double a[7];
  unsigned long i;

  for (i = 0; i <= d; i++)
    a[i] = mpz_get_d (f[i]);
  return L2_lognorm_d (a, d, s);
}

/* returns the optimal skewness corresponding to L2_lognorm */
double
L2_skewness (mpz_t *f, int deg, int prec)
{
  double s, n0, n1, a, b, c, d, nc, nd, fd[7];
  int i;

  /* convert once for all to double's to avoid expensive mpz_get_d() */
  for (i = 0; i <= deg; i++)
    fd[i] = mpz_get_d (f[i]);

  /* first isolate the minimum in an interval [s, 2s] by dichotomy */
  
  s = 1.0;
  n0 = L2_lognorm_d (fd, deg, s);
  while ((n1 = L2_lognorm_d (fd, deg, 2 * s)) < n0)
    {
      s = 2.0 * s;
      n0 = n1;
    }

  /* We have L2(s/2) > L2(s) <= L2(2s)
     Assuming L2(s) is first decreasing, then increasing, the minimum is
     attained in [s/2, 2s]. */
  a = (s == 1.0) ? 1.0 : 0.5 * s;
  b = 2.0 * s;

  /* use trichotomy */
  while (prec--)
    {
      c = (2.0 * a + b) / 3.0;
      d = (a + 2.0 * b) / 3.0;
      nc = L2_lognorm_d (fd, deg, c);
      nd = L2_lognorm_d (fd, deg, d);
      if (nc < nd) /* the minimum should be in [a,d] */
        b = d;
      else /* L2(c) > L2(d): the minimum should be in [c, b] */
        a = c;
    }

  return (a + b) * 0.5;
}

/*****************************************************************************/

/* Decompose p->n in base m/b, rounding to nearest
   (linear polynomial is b*x-m).
   For b=1, it suffices to check the coefficient f[d-1] is divisible by b,
   then all successive coefficients will be.
   Assumes b^2 fits in an unsigned long.
 */
void
generate_base_mb (cado_poly p, mpz_t m, unsigned long b)
{
  unsigned long i;
  int d = p->degree;

  ASSERT (d >= 0);

  if (b == 1)
    {
      mpz_ndiv_qr (p->f[1], p->f[0], p->n, m);
      for (i = 1; i < (unsigned long) d; i++)
        mpz_ndiv_qr (p->f[i+1], p->f[i], p->f[i], m);
    }
  else /* b >= 2 */
    {
      unsigned long m_mod_b;
      long k, im;
      mpz_t powm;
      int s;

      mpz_init (powm);
      
      /* we need that m is invertible mod b */
      while (mpz_gcd_ui (powm, m, b) != 1)
        mpz_add_ui (m, m, 1);

      mpz_pow_ui (powm, m, d);

      mpz_set_ui (p->f[0], b);
      mpz_invert (p->f[0], powm, p->f[0]);
      im = mpz_get_ui (p->f[0]); /* 1/m^d mod b */

      mpz_ndiv_qr (p->f[d], p->f[d - 1], p->n, powm);
      /* N = f[d] * m^d + f[d-1]. We want f[d-1] to be divisible by b,
         which may need to consider f[d]+k instead of f[d], i.e., we want
         f[d-1]-k*m^d = 0 (mod b), i.e., k = f[d-1]/m^d (mod b). */
      k = mpz_fdiv_ui (p->f[d - 1], b); /* f[d-1] mod b */
      k = (k * im) % b;                 /* k = f[d-1]/m^d (mod b) */
      if (k != 0)
        {
          s = mpz_sgn (p->f[d - 1]);
          if ((s >= 0 && 2 * k + 1 > (long) b) ||
              (s < 0 && 2 * k - 1 > (long) b))
            k = k - (long) b;
          /* f[d] -> f[d]+l, f[d-1] -> f[d-1]-l*m^d */
          mpz_add_si (p->f[d], p->f[d], k);
          mpz_submul_si (p->f[d - 1], powm, k);
        }
      /* now N = f[d] * m^d + f[d-1], with f[d-1] divisible by b */
      m_mod_b = mpz_fdiv_ui (m, b); /* m mod b */
      mpz_divexact_ui (p->f[d - 1], p->f[d - 1], b);
      for (i = d - 1; i >= 1; i--)
        {
          /* p->f[i] is the current remainder */
          mpz_divexact (powm, powm, m); /* m^i */
          mpz_ndiv_qr (p->f[i], p->f[i - 1], p->f[i], powm);
          /* f[i] = q*m^i + r: we want to write f[i] = (q+k) * m^i + (r-k*m^i)
             with r-k*m^i divisible by b, i.e., k = r/m^i mod b */
          k = mpz_fdiv_ui (p->f[i - 1], b); /* r mod b */
          im = (im * m_mod_b) % b; /* 1/m^i mod b */
          k = (k * im) % b; /* r/m^i mod b */
          s = mpz_sgn (p->f[i - 1]);
          if ((s >= 0 && 2 * k + 1 > (long) b) ||
              (s < 0 && 2 * k - 1 > (long) b))
            k = k - (long) b;
          mpz_add_si (p->f[i], p->f[i], k);
          mpz_submul_si (p->f[i - 1], powm, k);
          mpz_divexact_ui (p->f[i - 1], p->f[i - 1], b);
        }
      mpz_clear (powm);
    }
  mpz_neg (p->g[0], m);
  mpz_set_ui (p->g[1], b);
  mpz_set (p->m, m);
  p->skew = SKEWNESS (p->f, d, SKEWNESS_DEFAULT_PREC);
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
}

/* Return the smallest value of alpha(f + k*(b*x-m)) for |k| <= K,
   and modify f[] accordingly. */
double
rotate (mpz_t *f, int d, unsigned long alim, long K, mpz_t m, unsigned long b,
        int verbose)
{
  long k, k0, kmin = 0, K0, K1;
  double alpha, alpha0 = 0.0, best_alpha = 9.0;

  K0 = -K;
  K1 = K;

  ASSERT_ALWAYS(K0 <= 0);
  ASSERT_ALWAYS(K1 >= 0);

  k0 = 0; /* the coefficients f[] correspond to f+k*g */

  /* search for the smallest alpha[k] */
  for (k = K0; k <= K1; k++)
    {
      /* we now have f + k0*(bx-m) and we want f + k*(bx-m), thus we have
         to add (k-k0)*(bx-m) */
      mpz_add_si (f[1], f[1], (long) b * (k - k0));
      mpz_submul_si (f[0], m, k - k0);
      k0 = k;
      alpha = get_alpha (f, d, alim);
      if (k == 0)
        alpha0 = alpha;
      if (alpha < best_alpha)
        {
          best_alpha = alpha;
          kmin = k;
        }
    }

  /* we now have f + k0*(bx-m) and we want f + kmin*(bx-m), thus we have
     to add (kmin-k0)*(bx-m) */
  mpz_add_si (f[1], f[1], (long) b * (kmin - k0));
  mpz_submul_si (f[0], m, kmin - k0);

  if (verbose)
    printf ("# Rotate by %ld: alpha improved from %f to %f\n", kmin,
            alpha0, best_alpha);

  return best_alpha;
}

/* Returns k such that f(x+k) has the smallest 1-norm, with the corresponding
   skewness.
   The coefficients f[] are modified to those of f(x+k).

   The linear polynomial b*x-m is changed into b*(x+k)-m, thus m is
   changed into m-k*b.
*/
long
translate (mpz_t *f, int d, mpz_t m, unsigned long b)
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
  mpz_sub_ui (m, m, b);
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
        mpz_sub_ui (m, m, b);
      else
        mpz_add_ui (m, m, b);
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
    mpz_add_ui (m, m, b);
  else
    mpz_sub_ui (m, m, b);
  k = k - dir;

  return k;
}

/* Method described in reference [2], slide 8
   (Bernstein considered degree d=5, here we generalize to arbitrary d):
   Take m \approx B n^{1/(d+1)}
   Then f_d \approx B^{-d} n^{1/(d+1)}.
   Decompose n in base m, and control f_{d-1} as in Method 2.

   The number of tries T is related to the decrease B by:
   For d=5, T = B^4.5; for d=4, T = B^3.3.

   b is the leading coefficient of the linear polynomial.
*/
void
generate_poly (cado_poly out, double T, unsigned long b)
{
  unsigned long d = out->degree, alim;
  double B, logmu, alpha, E, mB;
  double best_E = DBL_MAX;
  mpz_t best_m, k, t, r;
  /* value of T = B^eff[d] for d <= 7 */
  static double eff[] = {0.0, 0.0, 3.0, 3.0, 3.333, 4.5, 6.6, 9.0};
  int given_m;
  m_logmu_t *M; /* stores values of m and log(mu) found in 1st phase */
  unsigned long Malloc; /* number of allocated entries in M */
  unsigned long Msize;  /* number of stored entries in M */
  unsigned long Malloc2;
  unsigned long Msize2, i;
  double st;

  ASSERT (T <= 9007199254740992.0); /* if T > 2^53, T-- will loop */

  mpz_init (best_m);

  alim = (unsigned long) cbrt ((double) default_alim (out->n));

  given_m = mpz_cmp_ui (out->m, 0) != 0;
  if (given_m) /* use given value of m */
    {
      mpz_set (best_m, out->m);
      goto end;
    }
  
  Malloc = (unsigned long) sqrt (T);
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
  while (T-- > 0)
    {
      mpz_pow_ui (t, out->m, d - 1); /* m^(d-1) */
      mpz_ndiv_qr (t, r, out->n, t); /* t = round (N / m^(d-1)) */
      /* we use mpz_fdiv_qr here (round towards -infinity) to ensure
         f[d-1] >= 0, which in turn ensures f[d-1]/(d f[d]) >= 0 */
      mpz_fdiv_qr (out->f[d], out->f[d - 1], t, out->m); /* f[d]*m+f[d-1]*/
      if (mpz_cmp_ui (out->f[d], 0) == 0)
        {
          gmp_fprintf (stderr, "# Stopping selection at m=%Zd since leading coefficient is zero\n", out->m);
          break;
        }
      mpz_mul_ui (t, out->f[d], d);
      mpz_fdiv_q (k, out->f[d - 1], t);
      mpz_add (out->m, out->m, k);
      generate_base_mb (out, out->m, b);
      logmu = LOGNORM (out->f, d, out->skew);
      Msize = m_logmu_insert (M, Malloc, Msize, out->m, logmu);
      if (T > 0 && mpz_sgn (out->f[d - 1]) > 0)
        { /* try another small f[d-1], avoiding a mpz_pow_ui call */
          mpz_add_ui (out->m, out->m, 1);
          generate_base_mb (out, out->m, b);
          T --;
          logmu = LOGNORM (out->f, d, out->skew);
          Msize = m_logmu_insert (M, Malloc, Msize, out->m, logmu);
        }
      mpz_add_ui (out->m, out->m, 1);
    }
  fprintf (stderr, "# First phase took %.2fs and kept %lu candidates\n",
           (seconds () - st), Msize);

  /* Second phase: loop over entries in M database. We do this in two
     subphases: first compute an estimate of alpha(F) up to alim^(1/3),
     and keep only the Malloc^(1/2) best values, then compute a more
     accurate estimate of alpha(F) on the remaining candidates. */
  Malloc2 = (unsigned long) sqrt ((double) Malloc);
  Msize2 = 0;
  st = seconds ();
  for (i = 0; i < Msize; i++)
    {
      logmu = M[i].logmu;
      generate_base_mb (out, M[i].m, b);
      alpha = get_alpha (out->f, d, alim);
      E = logmu + alpha;
      /* the following overwrites the M database */
      Msize2 = m_logmu_insert (M, Malloc2, Msize2, out->m, E);
    }
  /* In principle we should compute the contribution to alpha(F) for primes
     up to the algebraic factor base bound alim, but since it is too expensive,
     we use the heuristic of using sqrt(alim). */
  alim = (unsigned long) sqrt ((double) default_alim (out->n));
  alim = 100;
  for (i = 0; i < Msize2; i++)
    {
      generate_base_mb (out, M[i].m, b);
      translate (out->f, d, out->m, b);
      logmu = LOGNORM (out->f, d, out->skew);
      alpha = rotate (out->f, d, alim, MAX_ROTATE, out->m, b, 0);
      /* if there was some translation or rotation, we need to recompute
         skewness and norm */
      out->skew = SKEWNESS (out->f, out->degree, SKEWNESS_DEFAULT_PREC);
      logmu = LOGNORM (out->f, d, out->skew);
      E = logmu + alpha;
      if (E < best_E)
        {
          best_E = E;
          gmp_fprintf (stderr, "# m=%Zd skew=%1.2f logmu=%1.2f alpha~%1.2f E=%1.2f\n",
                       out->m, out->skew, logmu, alpha, E);
          mpz_set (best_m, out->m);
        }
    }
  fprintf (stderr, "# Second phase took %.2fs for %lu candidates\n",
           (seconds () - st), Msize);

  mpz_clear (k);
  mpz_clear (t);
  mpz_clear (r);
  m_logmu_clear (M, Malloc);

 end:
  generate_base_mb (out, best_m, b);
  /* rotate */
  rotate (out->f, d, alim, MAX_ROTATE, out->m, b, 1);
  /* translate */
  translate (out->f, d, out->m, b);
  /* recompute skewness, since it might differ from the original
     polynomial with no rotation */
  out->skew = SKEWNESS (out->f, out->degree, 20);
  generate_sieving_parameters (out);
  mpz_clear (best_m);
}

/* st is the initial value cputime() at the start of the program */
void
print_poly (FILE *fp, cado_poly p, int argc, char *argv[], double st, int raw)
{
  int i;
  double alpha, logmu;

  if (strlen (p->name) != 0)
    fprintf (fp, "name: %s\n", p->name);
  fprintf (fp, "n: ");
  mpz_out_str (fp, 10, p->n);
  fprintf (fp, "\n");
  fprintf (fp, "skew: %1.3f\n", p->skew);
  logmu = LOGNORM (p->f, p->degree, p->skew);
  alpha = get_alpha (p->f, p->degree, (p->alim > 1000) ? 1000 : p->alim);
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
  fprintf (fp, "m: ");
  /* if g[1]<>1, then m = -g[0]/g[1] mod n */
  if (mpz_cmp_ui (p->g[1], 1) != 0)
    {
      mpz_invert (p->m, p->g[1], p->n);
      mpz_neg (p->m, p->m);
      mpz_mul (p->m, p->m, p->g[0]);
      mpz_mod (p->m, p->m, p->n);
    }
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
    }
  fprintf (fp, "# generated by polyselect.r%s", REV);
  for (i = 1; i < argc; i++)
    fprintf (fp, " %s", argv[i]);
  fprintf (fp, " in %.2fs\n", (seconds () - st));
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
  cado_input in;
  cado_poly out;
  int verbose = 0, raw = 1, i;
  double effort = 1.0;
  char **argv0 = argv;
  int argc0 = argc;
  double st = seconds ();
  mpz_t m;
  int degree = 0;      /* if 0, use default values */
  unsigned long b = 1; /* leading coefficient of linear polynomial */

  fprintf (stderr, "# %s.r%s", argv[0], REV);
  for (i = 1; i < argc; i++)
    fprintf (stderr, " %s", argv[i]);
  fprintf (stderr, "\n");

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
      else if (strcmp (argv[1], "-full") == 0)
	{
	  raw = 0;
	  argv ++;
	  argc --;
	}
      else if (argc > 2 && strcmp (argv[1], "-e") == 0)
	{
	  effort = atof (argv[2]);
	  argv += 2;
	  argc -= 2;
	}
      else if (argc > 2 && strcmp (argv[1], "-m") == 0)
	{
          mpz_set_str (m, argv[2], 0);
	  argv += 2;
	  argc -= 2;
	}
      else if (argc > 2 && strcmp (argv[1], "-degree") == 0)
	{
          degree = atoi (argv[2]);
	  argv += 2;
	  argc -= 2;
	}
      else if (argc > 2 && strcmp (argv[1], "-b") == 0)
	{
          b = atoi (argv[2]);
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

  if (mpz_cmp_ui (m, 0) != 0 && effort != 1.0)
    {
      fprintf (stderr, "Warning: option -m xxx implies -e 1\n");
      effort = 1.0;
    }

  /* checks b^2 fits in an unsigned long */
  if (b != 1)
    {
      mpz_t bb;
      
      mpz_init_set_ui (bb, b);
      mpz_mul_ui (bb, bb, b);
      if (mpz_fits_ulong_p (bb) == 0)
        {
          fprintf (stderr, "Error, b^2 must fit in an unsigned long\n");
          exit (1);
        }
      mpz_clear (bb);
    }

  init (in);
  parse_input (in);
  if (verbose)
    fprintf (stderr, "L[1/3,(64/9)^(1/3)](n)=%1.2e\n", L (in->n));
  if (degree != 0)
    in->degree = degree;
  init_poly (out, in);
  mpz_set (out->m, m); /* if non zero, use given value of m */
  clear (in);

  generate_poly (out, effort, b);
  print_poly (stdout, out, argc0, argv0, st, raw);
  clear_poly (out);

  mpz_clear (m);

  return 0;
}
