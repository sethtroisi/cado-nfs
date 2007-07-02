/* arithmetic on polynomials over Z/pZ, with coefficients represented by
   the 'long' type */

#include <stdio.h>
#include <stdlib.h>
#include "utils.h"

/* allocate an array of d coefficients, and initialize it */
static long*
alloc_long_array (int d)
{
  long *f;
  int i;

  f = (long*) malloc (d * sizeof (long));
  return f;
}

/* reallocate an array to d coefficients */
static long*
realloc_long_array (long *f, int d)
{
  return (long*) realloc (f, d * sizeof (long));
}

/* initialize a polynomial with maximal degree d */
void
long_poly_init (long_poly_t f, int d)
{
  f->alloc = d + 1;
  f->degree = -1; /* initialize to 0 */
  f->coeff = alloc_long_array (d + 1);
}

/* clear a polynomial */
void
long_poly_clear (long_poly_t f)
{
  free (f->coeff);
}

/* realloc f to (at least) n coefficients */
void
long_poly_realloc (long_poly_t f, int n)
{
  if (f->alloc < n)
    {
      f->coeff = realloc_long_array (f->coeff, n);
      f->alloc = n;
    }
}

/* return 1/s mod t */
long
invert_si (long s, long t)
{
  long u1, v1, q;
  long u2, v2;

  assert (t > 0);
 
  u1 = 1;
  v1 = 0;
  u2 = (s > 0) ? s : -s;
  v2 = t;

  while (v2 != 0)
    {
      /* unroll twice and swap u/v */
      q = u2 / v2;
      u1 = u1 - q * v1;
      u2 = u2 - q * v2;
      
      if (u2 == 0)
	{
	  u1 = v1;
	  break;
	}

      q = v2 / u2;
      v1 = v1 - q * u1;
      v2 = v2 - q * u2;
    }
 
  if (u1 < 0)
    u1 = u1 - t * (-u1 / t - 1);
  
  return (s > 0) ? u1 : t-u1;
}

/* f <- g */
void
long_poly_set (long_poly_t f, const long_poly_t g)
{
  int i, d = g->degree;

  long_poly_realloc (f, d + 1);
  f->degree = d;
  for (i = 0; i <= d; i++)
    f->coeff[i] = g->coeff[i];
}

/* fp <- f/lc(f) mod p. Return degree of fp (-1 if fp=0). */
int
long_poly_set_mod (long_poly_t fp, mpz_t *f, int d, long p)
{
  int i;
  long ifd;

  while (d >= 0 && mpz_divisible_ui_p (f[d], p))
    d --;
  ASSERT (d >= 0); /* f is 0 mod p: should not happen since otherwise p would
		      divide N, because f(m)=N */
  ifd = invert_si (mpz_fdiv_ui (f[d], p), p);
  for (i = 0; i < d; i++)
    fp->coeff[i] = (mpz_fdiv_ui (f[i], p) * ifd) % p;
  fp->coeff[d] = 1;
  fp->degree = d;

  return d;
}

/* f <- a*x+b, a <> 0 */
void
long_poly_set_linear (long_poly_t f, long a, long b)
{
  long_poly_realloc (f, 2);
  f->degree = 1;
  f->coeff[1] = a;
  f->coeff[0] = b;
}

/* swap f and g */
void
long_poly_swap (long_poly_t f, long_poly_t g)
{
  int i;
  long *t;

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
long_poly_sqr (long_poly_t h, const long_poly_t g)
{
  int i, j, dg = g->degree;
  long *gc = g->coeff, *hc = h->coeff;

  long_poly_realloc (h, 2 * g->degree + 1);
  for (i = 0; i <= dg; i++)
    for (j = i + 1; j <= dg; j++)
      if (i == 0 || j == dg)
	hc[i + j] = gc[i] * gc[j];
      else
	hc[i + j] += gc[i] * gc[j];
  for (i = 1; i < 2 * dg; i++)
    hc[i] <<= 1;
  hc[0] = gc[0] * gc[0];
  hc[2 * dg] = gc[dg] * gc[dg];
  for (i = 1; i < dg; i++)
    hc[2 * i] += gc[i] * gc[i];
  h->degree = 2 * dg;
}

/* normalize h so that h->coeff[deg(h)] <> 0 */
static void
long_poly_normalize (long_poly_t h)
{
  int dh = h->degree;

  while (dh >= 0 && h->coeff[dh] == 0)
    dh --;
  h->degree = dh;
}

/* g <- h mod p */
void
long_poly_mod_ui (long_poly_t g, const long_poly_t h, long p)
{
  int i;

  long_poly_realloc (g, h->degree + 1);
  for (i = 0; i <= h->degree; i++)
    g->coeff[i] = h->coeff[i] % p;
  g->degree = h->degree;
  long_poly_normalize (g);
}

/* h <- (x+a)*h mod p */
void
long_poly_mul_x (long_poly_t h, long a, long p)
{
  int i, d = h->degree;

  long_poly_realloc (h, d + 2);
  h->coeff[d + 1] = h->coeff[d];
  for (i = d - 1; i >= 0; i--)
    h->coeff[i + 1] = (h->coeff[i] + a * h->coeff[i + 1]) % p;
  h->coeff[0] = (a * h->coeff[0]) % p;
  h->degree ++;
}

/* h <- g - x mod p */
void
long_poly_sub_x (long_poly_t h, const long_poly_t g, long p)
{
  int i, d = g->degree;

  long_poly_realloc (h, (d < 2) ? 2 : d + 1);
  for (i = 0; i <= d; i++)
    h->coeff[i] = g->coeff[i];
  for (i = d + 1; i <= 1; i++)
    h->coeff[i] = 0;
  h->coeff[1] = (h->coeff[1] - 1) % p;
  h->degree = (d < 2) ? 2 : d;
  long_poly_normalize (h);
}

/* h <- g - 1 mod p */
void
long_poly_sub_1 (long_poly_t h, const long_poly_t g, long p)
{
  int i, d = g->degree;

  long_poly_realloc (h, (d < 1) ? 1 : d + 1);
  for (i = 0; i <= d; i++)
    h->coeff[i] = g->coeff[i];
  for (i = d + 1; i <= 0; i++)
    h->coeff[i] = 0;
  h->coeff[0] = (h->coeff[0] - 1) % p;
  h->degree = (d < 1) ? 1 : d;
  long_poly_normalize (h);
}

/* h <- rem(h, f) mod p, f not necessarily monic */
void
long_poly_div_r (long_poly_t h, long_poly_t f, long p)
{
  int i, d = f->degree, dh = h->degree, monic;
  long *hc = h->coeff, t;

  monic = f->coeff[d] == 1;
  if (!monic)
    t = invert_si (f->coeff[d], p); /* 1/f[d] mod p */
  while (dh >= d)
    {
      /* subtract h[dh]/f[d]*x^(dh-d)*f from h */
      if (!monic)
	hc[dh] = (hc[dh] * t) % p;
      for (i = 0; i < d; i++)
	hc[dh - d + i] -= hc[dh] * f->coeff[i];
      h->degree = dh - 1;
      /* it is not necessary to reduce h[j] mod p here for j < dh */
      do
	{
	  dh --;
	  if (dh >= 0)
	    hc[dh] = hc[dh] % p;
	}
      while (dh >= 0 && hc[dh] == 0);
      h->degree = dh;
    }
  h->degree = dh;
  /* reduce mod p */
  long_poly_mod_ui (h, h, p);
}

/* q <- divexact(h, f) mod p, f not necessarily monic. Clobbers h. */
void
long_poly_divexact (long_poly_t q, long_poly_t h, const long_poly_t f, long p)
{
  int i, d = f->degree, dh = h->degree, monic;
  long *hc = h->coeff, t;

  assert (d >= 0);
  assert (dh >= 0);
  assert (dh >= d);
  long_poly_realloc (q, dh + 1 - d);
  q->degree = dh - d;
  monic = f->coeff[d] == 1;
  if (!monic)
    t = invert_si (f->coeff[d], p); /* 1/f[d] mod p */
  while (dh >= d)
    {
      /* subtract h[dh]/f[d]*x^(dh-d)*f from h */
      if (!monic)
	hc[dh] = (hc[dh] * t) % p;
      q->coeff[dh - d] = hc[dh];
      /* we only need to update the coefficients of degree >= d of h,
	 i.e., we want i >= 2d - dh */
      for (i = (2 * d > dh) ? 2 * d - dh : 0; i < d; i++)
	hc[dh - d + i] -= hc[dh] * f->coeff[i];
      h->degree = dh - 1;
      /* it is not necessary to reduce h[j] mod p here for j < dh */
      dh --;
      if (dh >= 0)
	hc[dh] = hc[dh] % p;
    }
  h->degree = dh;
  long_poly_normalize (q);
  long_poly_normalize (h);
}

/* fp <- gcd (fp, g), clobbers g */
void
long_poly_gcd (long_poly_t fp, long_poly_t g, unsigned long p)
{
  while (g->degree >= 0)
    {
      long_poly_div_r (fp, g, p);
      long_poly_mod_ui (fp, fp, p);
      /* now deg(fp) < deg(g): swap f and g */
      long_poly_swap (fp, g);
    }
}

void
long_poly_out (FILE *fp, long_poly_t f)
{
  if (f->degree < 0)
    fprintf (fp, "0\n");
  else
    {
      int i;
      for (i = 0; i <= f->degree; i++)
	if (f->coeff[i] > 0)
	  fprintf (fp, "+%d*x^%d", f->coeff[i], i);
	else if (f->coeff[i] < 0)
	  fprintf (fp, "%d*x^%d", f->coeff[i], i);
      fprintf (fp, ";\n");
    }
}

/* returns f(x) mod p */
long
long_poly_eval (long_poly_t f, long x, long p)
{
  int i, d = f->degree;

  long v = f->coeff[d];
  for (i = d - 1; i >= 0; i--)
    v = (v * x + f->coeff[i]) % p;
  return v;
}

/* Return the number n of roots of f mod p using a naive algorithm.
   If r is not NULL, put the roots in r[0], ..., r[n-1]. */
int
long_poly_roots_mod (long *r, long_poly_t f, const long p)
{
  int n = 0;
  long x;
  
  for (x = 0; x < p; x++)
    if (long_poly_eval (f, x, p) == 0)
      {
	if (r != NULL)
	  r[n] = x;
	n ++;
      }
  return n;
}

/* g <- (x+a)^e mod (fp, e), using auxiliary polynomial h */
void
long_poly_powmod_ui (long_poly_t g, long_poly_t fp, long_poly_t h, long a,
		     long e, long p)
{
  int k = nbits (e);

  /* initialize g to x */
  long_poly_set_linear (g, 1, a);

  assert (e > 0);
  for (k -= 2; k >= 0; k--)
    {
      long_poly_sqr (h, g);             /* h <- g^2 */
      if (e & (1 << k))
	long_poly_mul_x (h, a, p);            /* h <- x*h */

      long_poly_div_r (h, fp, p);       /* h -> rem(h, fp) */
      long_poly_mod_ui (g, h, p);       /* g <- h mod p */
    }
}

/* Assuming f is a (squarefree) product of linear factors mod p, splits it
   and put the corresponding roots mod p in r[]. 
   Return number of roots found (should be degree of f).
   Assumes p is odd, and deg(f) >= 1.
*/
int
long_poly_cantor_zassenhaus (long *r, long_poly_t f, long p, int depth)
{
  long a;
  long_poly_t q, h, ff;
  int d = f->degree, dq, n, m;
  
  assert (p & 1);
  assert (d >= 1);

  if (d == 1) /* easy case: linear factor x + a ==> root -a mod p */
    {
      a = invert_si (-f->coeff[1], p); /* -1/f[1] mod p */
      r[0] = (a * f->coeff[0]) % p;
      return 1;
    }

  /* if f has degree d, then q,h may have up to degree 2d-1 in the powering
     algorithm */
  long_poly_init (q, 2 * d - 1);
  long_poly_init (h, 2 * d - 1);
  long_poly_init (ff, d);
  for (a = 0; a < p; a++)
    {
      long_poly_set_linear (q, 1, a);
      long_poly_powmod_ui (q, f, h, a, (p - 1) / 2, p);
      long_poly_sub_1 (q, q, p);
      long_poly_set (h, f);
      long_poly_gcd (q, h, p);
      dq = q->degree;
      if (0 < dq && dq < d)
	{
	  n = long_poly_cantor_zassenhaus (r, q, p, depth + 1);
	  assert (n == dq);
	  long_poly_set (ff, f); /* long_poly_divexact clobbers its 2nd arg */
	  long_poly_divexact (h, ff, q, p);
	  m = long_poly_cantor_zassenhaus (r + n, h, p, depth + 1);
	  assert (m == h->degree);
	  n += m;
	  break;
	}
    }
  assert (a != p);
  long_poly_clear (q);
  long_poly_clear (h);
  long_poly_clear (ff);
  return n;
}
