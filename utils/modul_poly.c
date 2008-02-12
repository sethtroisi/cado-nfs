/* arithmetic on polynomials over Z/pZ, with coefficients represented by
   the 'long' type */

#include <stdio.h>
#include <stdlib.h>
#include "utils.h"
#include "mod_ul.h"



#ifdef LONG
#undef LONG
#endif
#define LONG unsigned long

typedef struct {
  int alloc;    /* number of allocated coefficients */
  int degree;   /* degree < alloc */
  LONG *coeff; /* coefficient list */
} __modul_poly_struct;
typedef __modul_poly_struct modul_poly_t[1];


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
modul_poly_init (modul_poly_t f, int d)
{
  assert (d >= 0);
  f->degree = -1; /* initialize to 0 */
  f->coeff = alloc_long_array (d + 1);
  f->alloc = d + 1;
}

/* clear a polynomial */
void
modul_poly_clear (modul_poly_t f)
{
  free (f->coeff);
  f->alloc = 0;
}

/* realloc f to (at least) n coefficients */
void
modul_poly_realloc (modul_poly_t f, int n)
{
  if (f->alloc < n)
    {
      f->coeff = realloc_long_array (f->coeff, n);
      f->alloc = n;
    }
}

#if 0
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
#endif

/* f <- g */
void
modul_poly_set (modul_poly_t f, const modul_poly_t g)
{
  int i, d = g->degree;

  modul_poly_realloc (f, d + 1);
  f->degree = d;
  for (i = 0; i <= d; i++)
    f->coeff[i] = g->coeff[i];
}

/* f <- f/lc(f) mod p */
void
modul_poly_make_monic (modul_poly_t f, modulus_t p)
{
  int d = f->degree, i;
  residueul_t ilc;

  if (f->coeff[d] == 1)
    return;

  modul_inv(ilc, &(f->coeff[d]), p);
  for (i = 0; i < d; i++)
      modul_mul(&(f->coeff[i]), &(f->coeff[i]), ilc, p);
  f->coeff[d] = 1;
}

/* fp <- f/lc(f) mod p. Return degree of fp (-1 if fp=0). */
int
modul_poly_set_mod (modul_poly_t fp, mpz_t *f, int d, modulus_t p)
{
  int i;

  while (d >= 0 && mpz_divisible_ui_p (f[d], p[0]))
    d --;
  ASSERT (d >= 0); /* f is 0 mod p: should not happen in the CADO-NFS context
                      since otherwise p would divide N, indeed f(m)=N */
  modul_poly_realloc (fp, d + 1);
  fp->degree = d;
  for (i = 0; i <= d; i++)
    fp->coeff[i] = mpz_fdiv_ui (f[i], p[0]);
  modul_poly_make_monic (fp, p);

  return d;
}

/* f <- a*x+b, a <> 0 */
void
modul_poly_set_linear (modul_poly_t f, LONG a, LONG b)
{
  modul_poly_realloc (f, 2);
  f->degree = 1;
  f->coeff[1] = a;
  f->coeff[0] = b;
}

/* swap f and g */
void
modul_poly_swap (modul_poly_t f, modul_poly_t g)
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
modul_poly_mul(modul_poly_t h, const modul_poly_t f, const modul_poly_t g, modulus_t p) {
  int df = f->degree;
  int dg = g->degree;
  int i, j;
  modul_poly_t res;
  residueul_t aux;

  modul_poly_init(res, df+dg);
  for (i = 0; i <= df+dg; ++i)
    res->coeff[i] = 0;
  for (i = 0; i <= df; ++i)
    for (j = 0; j <= dg; ++j) {
        modul_mul(aux, &(f->coeff[i]), &(g->coeff[j]), p);
        modul_add(&(res->coeff[i+j]), &(res->coeff[i+j]), aux, p);
    }
  res->degree=df+dg;
  modul_poly_set(h, res);
  modul_poly_clear(res);
}


/* h <- g^2 mod p, g and h must differ */
void
modul_poly_sqr (modul_poly_t h, const modul_poly_t g, modulus_t p)
{
  int i, j, dg = g->degree;
  LONG *gc, *hc;
  residueul_t aux;

  assert (dg >= -1);
  if (dg == -1) /* g is zero */
    {
      h->degree = -1;
      return;
    }
  modul_poly_realloc (h, 2 * dg + 1);
  gc = g->coeff;
  hc = h->coeff;
  for (i = 0; i <= dg; i++)
    for (j = i + 1; j <= dg; j++)
      if (i == 0 || j == dg) 
          modul_mul(&(hc[i + j]), &(gc[i]), &(gc[j]), p);
      else {
          modul_mul(aux, &(gc[i]), &(gc[j]), p);
          modul_add(&(hc[i + j]), &(hc[i + j]), aux, p);
      }
  for (i = 1; i < 2 * dg; i++)
      modul_add(&(hc[i]), &(hc[i]), &(hc[i]), p);
  modul_mul(&(hc[0]), &(gc[0]), &(gc[0]), p);
  modul_mul(&(hc[2*dg]), &(gc[dg]), &(gc[dg]), p);
  for (i = 1; i < dg; i++) {
      modul_mul(aux, &(gc[i]), &(gc[i]), p);
      modul_add(&(hc[2*i]), &(hc[2*i]), aux, p);
  }
  h->degree = 2 * dg;
}

/* normalize h so that h->coeff[deg(h)] <> 0 */
static void
modul_poly_normalize (modul_poly_t h)
{
  int dh = h->degree;

  while (dh >= 0 && h->coeff[dh] == 0)
    dh --;
  h->degree = dh;
}

#if 0
/* g <- h mod p */
void
modul_poly_mod_ui (modul_poly_t g, const modul_poly_t h, LONG p)
{
  int i;
  int dh = h->degree;

  modul_poly_realloc (g, dh + 1);
  for (i = 0; i <= dh; i++)
    g->coeff[i] = h->coeff[i] % p;
  g->degree = dh;
  modul_poly_normalize (g);
}
#endif

/* h <- (x+a)*h mod p */
void
modul_poly_mul_x (modul_poly_t h, LONG a, modulus_t p)
{
  int i, d = h->degree;
  LONG *hc;
  residueul_t aux;

  modul_poly_realloc (h, d + 2); /* (x+a)*h has degree d+1, thus d+2 coeffs */
  hc = h->coeff;
  hc[d + 1] = hc[d];
  for (i = d - 1; i >= 0; i--) {
      modul_mul(aux, &a, &(hc[i+1]), p);
      modul_add(&(hc[i+1]), aux, &(hc[i]), p);
  }
  modul_mul(&(hc[0]), &(hc[0]), &a, p);
  h->degree = d + 1;
}

/* h <- g - x mod p */
void
modul_poly_sub_x (modul_poly_t h, const modul_poly_t g, modulus_t p)
{
  int i, d = g->degree;

  /* if g has degree d >= 2, then g-x has degree d too;
     otherwise g-x has degree <= 1 */
  modul_poly_realloc (h, ((d < 2) ? 1 : d) + 1);
  for (i = 0; i <= d; i++)
    h->coeff[i] = g->coeff[i];
  for (i = d + 1; i <= 1; i++)
    h->coeff[i] = 0;
  modul_sub_ul(&(h->coeff[1]), &(h->coeff[1]), 1, p);
  h->degree = (d < 1) ? 1 : d;
  modul_poly_normalize (h);
}

/* h <- g - 1 mod p */
void
modul_poly_sub_1 (modul_poly_t h, const modul_poly_t g, modulus_t p)
{
  int i, d = g->degree;

  /* g-1 has degree d if d >= 1, and degree 0 otherwise */
  modul_poly_realloc (h, ((d < 1) ? 0 : d) + 1);
  for (i = 0; i <= d; i++)
    h->coeff[i] = g->coeff[i];
  for (i = d + 1; i <= 0; i++)
    h->coeff[i] = 0;
  modul_sub_ul(&(h->coeff[0]), &(h->coeff[0]), 1, p);
  h->degree = (d < 0) ? 0 : d;
  modul_poly_normalize (h);
}

/* h <- rem(h, f) mod p, f not necessarily monic */
void
modul_poly_div_r (modul_poly_t h, const modul_poly_t f, modulus_t p)
{
  int i, d = f->degree, dh = h->degree, monic;
  LONG *hc = h->coeff, t = 1;
  residueul_t aux;

  monic = f->coeff[d] == 1;
  if (!monic)
      modul_inv(&t, &(f->coeff[d]), p); /* 1/f[d] mod p */
  while (dh >= d)
    {
      /* subtract h[dh]/f[d]*x^(dh-d)*f from h */
      if (!monic)
          modul_mul(&(hc[dh]), &(hc[dh]), &t, p);
      for (i = 0; i < d; i++) {
          modul_mul(aux, &(hc[dh]), &(f->coeff[i]), p);
          modul_sub(&(hc[dh - d + i]), &(hc[dh - d + i]), aux, p);
      }
      do {
	  dh --;
      } while (dh >= 0 && hc[dh] == 0);
      h->degree = dh;
    }
}

/* q <- divexact(h, f) mod p, f not necessarily monic. Clobbers h. */
void
modul_poly_divexact (modul_poly_t q, modul_poly_t h, const modul_poly_t f, modulus_t p)
{
  int i, d = f->degree, dh = h->degree, monic;
  LONG *hc = h->coeff, t = 1;
  residueul_t aux;

  assert (d >= 0);
  assert (dh >= 0);
  assert (dh >= d);
  modul_poly_realloc (q, dh + 1 - d);
  q->degree = dh - d;
  monic = f->coeff[d] == 1;
  if (!monic)
      modul_inv(&t, &(f->coeff[d]), p);
  while (dh >= d)
    {
      /* subtract h[dh]/f[d]*x^(dh-d)*f from h */
      if (!monic)
          modul_mul(&(hc[dh]), &(hc[dh]), &t, p);
      q->coeff[dh - d] = hc[dh];
      /* we only need to update the coefficients of degree >= d of h,
	 i.e., we want i >= 2d - dh */
      for (i = (2 * d > dh) ? 2 * d - dh : 0; i < d; i++) {
          modul_mul(aux, &(hc[dh]), &(f->coeff[i]), p);
          modul_sub(&(hc[dh - d + i]), &(hc[dh - d + i]), aux, p);
      }
      dh --;
    }
  modul_poly_normalize (q);
}

/* fp <- gcd (fp, g), clobbers g */
void
modul_poly_gcd (modul_poly_t fp, modul_poly_t g, modulus_t p)
{
  while (g->degree >= 0)
    {
      modul_poly_div_r (fp, g, p);
//      modul_poly_mod_ui (fp, fp, p);
      /* now deg(fp) < deg(g): swap f and g */
      modul_poly_swap (fp, g);
    }
}

void
modul_poly_out (FILE *fp, modul_poly_t f)
{
  if (f->degree < 0)
    fprintf (fp, "0\n");
  else
    {
      int i;
      for (i = 0; i <= f->degree; i++)
	if (f->coeff[i] > 0)
	  fprintf (fp, "+%lld*x^%d", (long long int) f->coeff[i], i);
      fprintf (fp, ";\n");
    }
}

/* returns f(x) mod p */
LONG
modul_poly_eval (modul_poly_t f, LONG x, modulus_t p)
{
  int i, d = f->degree;
  LONG v;
  residueul_t aux;
  
  assert (d >= 0);
  v = f->coeff[d];
  for (i = d - 1; i >= 0; i--) {
      modul_mul(aux, &v, &x, p);
      modul_add(&v, aux, &(f->coeff[i]), p);
  }
  return v;
}

/* Return the number n of roots of f mod p using a naive algorithm.
   If r is not NULL, put the roots in r[0], ..., r[n-1]. */
int
modul_poly_roots_mod (LONG *r, modul_poly_t f, modulus_t p)
{
  int n = 0;
  LONG x;
  
  for (x = 0; x < p[0]; x++)
    if (modul_poly_eval (f, x, p) == 0)
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
modul_poly_powmod_ui (modul_poly_t g, modul_poly_t fp, modul_poly_t h, LONG a,
		     LONG e, modulus_t p)
{
  int k = nbits (e);

  modul_poly_make_monic (fp, p);

  /* initialize g to x */
  modul_poly_set_linear (g, 1, a);

  assert (e > 0);
  for (k -= 2; k >= 0; k--)
    {
      modul_poly_sqr (h, g, p);             /* h <- g^2 */
      if (e & (1 << k))
        modul_poly_mul_x (h, a, p);            /* h <- x*h */

      modul_poly_div_r (h, fp, p);       /* h -> rem(h, fp) */
      modul_poly_set (g, h);       /* g <- h  */
    }
}

/* g <- g^e mod (fp, p), using auxiliary polynomial h */
void
modul_poly_general_powmod_ui (modul_poly_t g, modul_poly_t fp, modul_poly_t h,
		     LONG e, modulus_t p)
{
  int k = nbits (e);
  modul_poly_t g_sav;

  modul_poly_init(g_sav, g->degree);
  modul_poly_set(g_sav, g);

  modul_poly_make_monic (fp, p);

  assert (e > 0);
  for (k -= 2; k >= 0; k--)
    {
      modul_poly_sqr (h, g, p);             /* h <- g^2 */
      modul_poly_div_r (h, fp, p);       /* h -> rem(h, fp) */
      modul_poly_set (g, h);       /* g <- h */
      if (e & (1 << k))
        modul_poly_mul (h, h, g_sav, p);            /* h <- g_sav*h */
      modul_poly_div_r (h, fp, p);       /* h -> rem(h, fp) */
      modul_poly_set (g, h);       /* g <- h */
    }
  modul_poly_clear(g_sav);
}
/* Assuming f is a (squarefree) product of linear factors mod p, splits it
   and put the corresponding roots mod p in r[]. 
   Return number of roots found (should be degree of f).
   Assumes p is odd, and deg(f) >= 1.
*/
int
modul_poly_cantor_zassenhaus (LONG *r, modul_poly_t f, modulus_t p, int depth)
{
  LONG a;
  modul_poly_t q, h, ff;
  int d = f->degree, dq, n, m;

  assert (p[0] & 1);
  assert (d >= 1);

  if (d == 1) /* easy case: linear factor x + a ==> root -a mod p */
    {
        residueul_t aux;
        modul_neg(aux, &(f->coeff[1]), p);
        modul_inv(&a, aux, p);
        modul_mul(&(r[0]), &a, &(f->coeff[0]), p);
        return 1;
    }

  /* if f has degree d, then q,h may have up to degree 2d-1 in the powering
     algorithm */
  modul_poly_init (q, 2 * d - 1);
  modul_poly_init (h, 2 * d - 1);
  modul_poly_init (ff, d);
  a = lrand48 () % p[0];
  for (;;)
    {
      modul_poly_powmod_ui (q, f, h, a, (p[0] - 1) / 2, p);
      modul_poly_sub_1 (q, q, p);
      modul_poly_set (h, f);
      modul_poly_gcd (q, h, p);
      dq = q->degree;
      assert (dq >= 0);
      if (0 < dq && dq < d)
	{
	  n = modul_poly_cantor_zassenhaus (r, q, p, depth + 1);
	  assert (n == dq);
	  modul_poly_set (ff, f); /* modul_poly_divexact clobbers its 2nd arg */
	  modul_poly_divexact (h, ff, q, p);
	  m = modul_poly_cantor_zassenhaus (r + n, h, p, depth + 1);
	  assert (m == h->degree);
	  n += m;
          break;
	}
      if (++a >= p[0])
        a -= p[0];
    }
  modul_poly_clear (q);
  modul_poly_clear (h);
  modul_poly_clear (ff);
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
modul_roots_mod_long (LONG *r, mpz_t *f, int d, modulus_t p)
{
  modul_poly_t fp, g, h;
  int df, n;

  /* the number of roots is the degree of gcd(x^p-x,f) */

  modul_poly_init (fp, d);
  d = modul_poly_set_mod (fp, f, d, p);
  /* d is the degree of fp (-1 if fp=0) */

  if (d == 0)
    {
      df = 0;
      goto clear_fp;
    }

  assert (d > 0);

  if ((r == NULL && p[0] <= ROOTS_MOD_THRESHOLD) ||
      (r != NULL && p[0] <= ROOTS_MOD_THRESHOLD2))
    {
      df = modul_poly_roots_mod (r, fp, p);
      goto clear_fp;
    }

  modul_poly_init (g, 2 * d - 1);
  modul_poly_init (h, 2 * d - 1);

  /* we first compute x^p mod fp; since fp has degree d, all operations can
     be done with polynomials of degree < 2d: a square give at most degree
     2d-2, and multiplication by x gives 2d-1. */

  /* g <- x^p mod fp */
  modul_poly_powmod_ui (g, fp, h, 0, p[0], p);

  /* subtract x */
  modul_poly_sub_x (g, g, p);

  /* fp <- gcd (fp, g) */
  modul_poly_gcd (fp, g, p);

  /* now fp contains gcd(x^p-x, f) */

  df = fp->degree;
  assert (df >= 0);
  
  if (r != NULL && df > 0)
    {
      n = modul_poly_cantor_zassenhaus (r, fp, p, 0);
      assert (n == df);
    }

  modul_poly_clear (g);

  modul_poly_clear (h);

 clear_fp:
  modul_poly_clear (fp);

  return df;
}

int modul_isirreducible_mod_long(modul_poly_t fp, modulus_t p)
{
  modul_poly_t g, gmx, h;
  int d, i;

  modul_poly_make_monic (fp, p);
  d = fp->degree;

  if (d == 0)
    return 1;

  assert (d > 0);

  modul_poly_init (g, 2 * d - 1);
  modul_poly_init (gmx, 2 * d - 1);
  modul_poly_init (h, 2 * d - 1);

  for (i = 1; 2*i <= d; ++i) {
    /* we first compute x^(p^i) mod fp; since fp has degree d, all operations can
       be done with polynomials of degree < 2d: a square give at most degree
       2d-2, and multiplication by x gives 2d-1. */

    if (i == 1) {
      /* g <- x^p mod fp */
      modul_poly_powmod_ui (g, fp, h, 0, p[0], p);
    } else {
      /* g <- g^p mod fp */
      modul_poly_general_powmod_ui (g, fp, h, p[0], p);
    }

    /* subtract x */
    modul_poly_sub_x (gmx, g, p);

    /* h <- gcd (fp, x^(p^i)-x) */
    modul_poly_set(h, fp);
    modul_poly_gcd (h, gmx, p);

    if (h->degree > 0)
      return 0;
  }

  modul_poly_clear (g);
  modul_poly_clear (gmx);
  modul_poly_clear (h);
  return 1;
}
