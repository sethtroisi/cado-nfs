/* Arithmetic on polynomials over Z/pZ, with coefficients represented by
   the 'long' type.
   Since products are performed with the '*' operator, and several products
   are accumulated before a reduction by '% p', p should be at most half the
   size of the 'long' type.
   FIXME: get a precise bound on p.
*/

#include <stdio.h>
#include <stdlib.h>
#include "utils.h"

/* allocate an array of d coefficients, and initialize it */
static LONG*
alloc_long_array (int d)
{
  LONG *f;

  f = (LONG*) malloc (d * sizeof (LONG));
  if (f == NULL)
    {
      fprintf (stderr, "Error, not enough memory\n");
      exit (1);
    }
  return f;
}

/* reallocate an array to d coefficients */
static LONG*
realloc_long_array (LONG *f, int d)
{
  f = (LONG*) realloc (f, d * sizeof (LONG));
  if (f == NULL)
    {
      fprintf (stderr, "Error, not enough memory\n");
      exit (1);
    }
  return f;
}

/* initialize a polynomial with maximal degree d */
void
long_poly_init (long_poly_t f, int d)
{
  assert (d >= 0);
  f->degree = -1; /* initialize to 0 */
  f->coeff = alloc_long_array (d + 1);
  f->alloc = d + 1;
}

/* clear a polynomial */
void
long_poly_clear (long_poly_t f)
{
  free (f->coeff);
  f->alloc = 0;
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
LONG
invert_si (LONG s, LONG t)
{
  LONG u1, v1, q;
  LONG u2, v2;

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

/* f <- f/lc(f) mod p */
void
long_poly_make_monic (long_poly_t f, LONG p)
{
  int d = f->degree, i;
  LONG ilc;

  if (f->coeff[d] == 1)
    return;

  ilc = invert_si (f->coeff[d], p);
  for (i = 0; i < d; i++)
    f->coeff[i] = (f->coeff[i] * ilc) % p;
  f->coeff[d] = 1;
}

/* fp <- f/lc(f) mod p. Return degree of fp (-1 if fp=0). */
int
long_poly_set_mod (long_poly_t fp, mpz_t *f, int d, LONG p)
{
  int i;

  while (d >= 0 && mpz_divisible_ui_p (f[d], p))
    d --;
  ASSERT (d >= 0); /* f is 0 mod p: should not happen in the CADO-NFS context
                      since otherwise p would divide N, indeed f(m)=N */
  long_poly_realloc (fp, d + 1);
  fp->degree = d;
  for (i = 0; i <= d; i++)
    fp->coeff[i] = mpz_fdiv_ui (f[i], p);
  long_poly_make_monic (fp, p);

  return d;
}

/* f <- a*x+b, a <> 0 */
void
long_poly_set_linear (long_poly_t f, LONG a, LONG b)
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
  LONG *t;

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

/* h <- f*g mod p */
void
long_poly_mul(long_poly_t h, const long_poly_t f, const long_poly_t g, LONG p) {
  int df = f->degree;
  int dg = g->degree;
  int i, j;
  long_poly_t res;

  long_poly_init(res, df+dg);
  for (i = 0; i <= df+dg; ++i)
    res->coeff[i] = 0;
  for (i = 0; i <= df; ++i)
    for (j = 0; j <= dg; ++j)
      res->coeff[i+j] += f->coeff[i]*g->coeff[j];
  /* reduce mod p */
  for (i = 0; i <= df+dg; i++)
    res->coeff[i] %= p;
  res->degree=df+dg;
  long_poly_set(h, res);
  long_poly_clear(res);
}


/* h <- g^2 mod p, g and h must differ */
void
long_poly_sqr (long_poly_t h, const long_poly_t g, LONG p)
{
  int i, j, dg = g->degree;
  LONG *gc, *hc;

  assert (dg >= -1);
  if (dg == -1) /* g is zero */
    {
      h->degree = -1;
      return;
    }
  long_poly_realloc (h, 2 * dg + 1);
  gc = g->coeff;
  hc = h->coeff;
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
  /* reduce mod p */
  for (i = 0; i <= 2 * dg; i++)
    hc[i] %= p;
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
long_poly_mod_ui (long_poly_t g, const long_poly_t h, LONG p)
{
  int i;
  int dh = h->degree;

  long_poly_realloc (g, dh + 1);
  for (i = 0; i <= dh; i++)
    g->coeff[i] = h->coeff[i] % p;
  g->degree = dh;
  long_poly_normalize (g);
}

/* h <- (x+a)*h mod p */
void
long_poly_mul_x (long_poly_t h, LONG a, LONG p)
{
  int i, d = h->degree;
  LONG *hc;

  long_poly_realloc (h, d + 2); /* (x+a)*h has degree d+1, thus d+2 coeffs */
  hc = h->coeff;
  hc[d + 1] = hc[d];
  for (i = d - 1; i >= 0; i--)
    hc[i + 1] = (hc[i] + a * hc[i + 1]) % p;
  hc[0] = (a * hc[0]) % p;
  h->degree = d + 1;
}

/* h <- g - x mod p */
void
long_poly_sub_x (long_poly_t h, const long_poly_t g, LONG p)
{
  int i, d = g->degree;

  /* if g has degree d >= 2, then g-x has degree d too;
     otherwise g-x has degree <= 1 */
  long_poly_realloc (h, ((d < 2) ? 1 : d) + 1);
  for (i = 0; i <= d; i++)
    h->coeff[i] = g->coeff[i];
  for (i = d + 1; i <= 1; i++)
    h->coeff[i] = 0;
  h->coeff[1] = (h->coeff[1] - 1) % p;
  h->degree = (d < 1) ? 1 : d;
  long_poly_normalize (h);
}

/* h <- g - 1 mod p */
void
long_poly_sub_1 (long_poly_t h, const long_poly_t g, LONG p)
{
  int i, d = g->degree;

  /* g-1 has degree d if d >= 1, and degree 0 otherwise */
  long_poly_realloc (h, ((d < 1) ? 0 : d) + 1);
  for (i = 0; i <= d; i++)
    h->coeff[i] = g->coeff[i];
  for (i = d + 1; i <= 0; i++)
    h->coeff[i] = 0;
  h->coeff[0] = (h->coeff[0] - 1) % p;
  h->degree = (d < 0) ? 0 : d;
  long_poly_normalize (h);
}

/* h <- rem(h, f) mod p, f not necessarily monic */
void
long_poly_div_r (long_poly_t h, const long_poly_t f, LONG p)
{
  int i, d = f->degree, dh = h->degree, monic;
  LONG *hc = h->coeff, t = 1;

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
  /* reduce mod p */
  long_poly_mod_ui (h, h, p);
}

/* q <- divexact(h, f) mod p, f not necessarily monic. Clobbers h. */
void
long_poly_divexact (long_poly_t q, long_poly_t h, const long_poly_t f, LONG p)
{
  int i, d = f->degree, dh = h->degree, monic;
  LONG *hc = h->coeff, t = 1;

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
  long_poly_normalize (q);
}

/* fp <- gcd (fp, g), clobbers g */
void
long_poly_gcd (long_poly_t fp, long_poly_t g, ULONG p)
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
	  fprintf (fp, "+%lld*x^%d", (long long int) f->coeff[i], i);
	else if (f->coeff[i] < 0)
	  fprintf (fp, "%lld*x^%d", (long long int) f->coeff[i], i);
      fprintf (fp, ";\n");
    }
}

/* returns f(x) mod p */
LONG
long_poly_eval (long_poly_t f, LONG x, LONG p)
{
  int i, d = f->degree;
  LONG v;
  
  assert (d >= 0);
  v = f->coeff[d];
  for (i = d - 1; i >= 0; i--)
    v = (v * x + f->coeff[i]) % p;
  return v;
}

/* Return the number n of roots of f mod p using a naive algorithm.
   If r is not NULL, put the roots in r[0], ..., r[n-1]. */
int
long_poly_roots_mod (LONG *r, long_poly_t f, const LONG p)
{
  int n = 0;
  LONG x;
  
  for (x = 0; x < p; x++)
    if (long_poly_eval (f, x, p) == 0)
      {
	if (r != NULL)
	  r[n] = x;
	n ++;
      }
  return n;
}

/* return the number of bits of p */
int
nbits (ULONG p)
{
  int k;

  for (k = 0; p != 0; p >>= 1, k ++);
  return k;
}

/* g <- (x+a)^e mod (fp, p), using auxiliary polynomial h */
void
long_poly_powmod_ui (long_poly_t g, long_poly_t fp, long_poly_t h, LONG a,
		     LONG e, LONG p)
{
  int k = nbits (e);

  long_poly_make_monic (fp, p);

  /* initialize g to x */
  long_poly_set_linear (g, 1, a);

  assert (e > 0);
  for (k -= 2; k >= 0; k--)
    {
      long_poly_sqr (h, g, p);             /* h <- g^2 */
      if (e & (1L << k))
        long_poly_mul_x (h, a, p);            /* h <- x*h */

      long_poly_div_r (h, fp, p);       /* h -> rem(h, fp) */
      long_poly_mod_ui (g, h, p);       /* g <- h mod p */
    }
}

/* g <- g^e mod (fp, p), using auxiliary polynomial h */
void
long_poly_general_powmod_ui (long_poly_t g, long_poly_t fp, long_poly_t h,
		     LONG e, LONG p)
{
  int k = nbits (e);
  long_poly_t g_sav;

  long_poly_init(g_sav, g->degree);
  long_poly_set(g_sav, g);

  long_poly_make_monic (fp, p);

  assert (e > 0);
  for (k -= 2; k >= 0; k--)
    {
      long_poly_sqr (h, g, p);             /* h <- g^2 */
      long_poly_div_r (h, fp, p);       /* h -> rem(h, fp) */
      long_poly_mod_ui (g, h, p);       /* g <- h mod p */
      if (e & (1L << k))
        long_poly_mul (h, h, g_sav, p);            /* h <- g_sav*h */
      long_poly_div_r (h, fp, p);       /* h -> rem(h, fp) */
      long_poly_mod_ui (g, h, p);       /* g <- h mod p */
    }
  long_poly_clear(g_sav);
}
/* Assuming f is a (squarefree) product of linear factors mod p, splits it
   and put the corresponding roots mod p in r[]. 
   Return number of roots found (should be degree of f).
   Assumes p is odd, and deg(f) >= 1.
*/
int
long_poly_cantor_zassenhaus (LONG *r, long_poly_t f, LONG p, int depth)
{
  LONG a;
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
  a = lrand48 () % p;
  for (;;)
    {
      long_poly_powmod_ui (q, f, h, a, (p - 1) / 2, p);
      long_poly_sub_1 (q, q, p);
      long_poly_set (h, f);
      long_poly_gcd (q, h, p);
      dq = q->degree;
      assert (dq >= 0);
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
      if (++a >= p)
        a -= p;
    }
  long_poly_clear (q);
  long_poly_clear (h);
  long_poly_clear (ff);
  return n;
}

#define ROOTS_MOD_THRESHOLD  43 /* if only the number of roots is needed */
#define ROOTS_MOD_THRESHOLD2 67 /* if the roots are needed too */

/* The following function returns the number of roots of f(x) mod p,
   where f has degree d.
   If r is not NULL, put the n roots of f(x) mod p in r[0]...r[n-1].

   For p <= ROOTS_MOD_THRESHOLD, determine the number of roots of f mod p
   by evaluating f on 0, 1, ..., p-1. In theory, this threshold should 
   depend on the degree of f. The above thresholds seem to be optimal
   on a Pentium M. We must have ROOTS_MOD_THRESHOLD2 >= 2.

   The following ideas are from Emmanuel Thome:
   FIXME1: instead of first computing f1 = gcd(x^p-x, f) which is the
   product of linear factors, and then splitting f1, we can directly
   compute f1 = gcd(x^((p-1)/2)-1, f) and f2 = gcd(x^((p-1)/2)+1, f),
   assuming the roots 0 has been taken out.
   FIXME2: in the EDF, instead of dividing out f by gcd((x+a)^((p-1)/2)-1, f),
   we can simply compute gcd((x+a)^((p-1)/2)+1, f).
   FIXME3: instead of computing (x+a)^((p-1)/2), we can translate f by -a,
   and keep x^((p-1)/2) unchanged. [But then why not translate x^((p-1)/2),
   which has lower degree? Warning: if f has root -a, we might miss it.]
*/
int
roots_mod_long (LONG *r, mpz_t *f, int d, const LONG p)
{
  long_poly_t fp, g, h;
  int df, n;

  /* the number of roots is the degree of gcd(x^p-x,f) */

  long_poly_init (fp, d);
  d = long_poly_set_mod (fp, f, d, p);
  /* d is the degree of fp (-1 if fp=0) */

  if (d == 0)
    {
      df = 0;
      goto clear_fp;
    }

  assert (d > 0);

  if ((r == NULL && p <= ROOTS_MOD_THRESHOLD) ||
      (r != NULL && p <= ROOTS_MOD_THRESHOLD2))
    {
      df = long_poly_roots_mod (r, fp, p);
      goto clear_fp;
    }

  long_poly_init (g, 2 * d - 1);
  long_poly_init (h, 2 * d - 1);

  /* we first compute x^p mod fp; since fp has degree d, all operations can
     be done with polynomials of degree < 2d: a square give at most degree
     2d-2, and multiplication by x gives 2d-1. */

  /* g <- x^p mod fp */
  long_poly_powmod_ui (g, fp, h, 0, p, p);

  /* subtract x */
  long_poly_sub_x (g, g, p);

  /* fp <- gcd (fp, g) */
  long_poly_gcd (fp, g, p);

  /* now fp contains gcd(x^p-x, f) */

  df = fp->degree;
  assert (df >= 0);
  
  if (r != NULL && df > 0)
    {
      n = long_poly_cantor_zassenhaus (r, fp, p, 0);
      assert (n == df);
    }

  long_poly_clear (g);

  long_poly_clear (h);

 clear_fp:
  long_poly_clear (fp);

  return df;
}

int isirreducible_mod_long(long_poly_t fp, const LONG p) {
  long_poly_t g, gmx, h;
  int d, i;

  long_poly_make_monic (fp, p);
  d = fp->degree;

  if (d == 0)
    return 1;

  assert (d > 0);

  long_poly_init (g, 2 * d - 1);
  long_poly_init (gmx, 2 * d - 1);
  long_poly_init (h, 2 * d - 1);

  for (i = 1; 2*i <= d; ++i) {
    /* we first compute x^(p^i) mod fp; since fp has degree d, all operations can
       be done with polynomials of degree < 2d: a square give at most degree
       2d-2, and multiplication by x gives 2d-1. */

    if (i == 1) {
      /* g <- x^p mod fp */
      long_poly_powmod_ui (g, fp, h, 0, p, p);
    } else {
      /* g <- g^p mod fp */
      long_poly_general_powmod_ui (g, fp, h, p, p);
    }

    /* subtract x */
    long_poly_sub_x (gmx, g, p);

    /* h <- gcd (fp, x^(p^i)-x) */
    long_poly_set(h, fp);
    long_poly_gcd (h, gmx, p);

    if (h->degree > 0)
      return 0;
  }

  long_poly_clear (g);
  long_poly_clear (gmx);
  long_poly_clear (h);
  return 1;
}
