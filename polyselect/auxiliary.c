/* common routines to polyselect and kleinjung */

#include <stdio.h>
#include <stdlib.h>
#include <float.h> /* for DBL_MAX */
#include <math.h>
#include "gmp.h"
#include "macros.h" /* for ASSERT_ALWAYS */
#include "cado.h"
#include "utils.h"
#include "auxiliary.h"
#include "rootsieve.h"
#include "murphyE.h"

// #define NEW_ROOTSIEVE

/* for the rotation, we try (j*x+k) for |k| <= 2^MAX_k */
int MAX_k = 16;

/************************* norm and skewness *********************************/

/* Same as L2_lognorm, but takes 'double' instead of 'mpz_t' as coefficients.
   Returns 1/2*log(int(int(F^2(xs,y)/s^d, x=-1..1), y=-1..1)) for the
   'rectangular' method, and
   1/2*log(int(int(F^2(r cos(t) s,r sin(t))*r/s^d, r=0..1), t=0..2*Pi))
   with the 'circular' method.
 */
static double
L2_lognorm_d (double *a, unsigned long d, double s, int method)
{
  double n;

  /* coefficients for degree 2 to 10:
    sage: [[4/(2*i+1)/(2*(d-i)+1) for i in [0..d]] for d in [2..10]]

    [[4/5, 4/9, 4/5],
     [4/7, 4/15, 4/15, 4/7],
     [4/9, 4/21, 4/25, 4/21, 4/9],
     [4/11, 4/27, 4/35, 4/35, 4/27, 4/11],
     [4/13, 4/33, 4/45, 4/49, 4/45, 4/33, 4/13],
     [4/15, 4/39, 4/55, 4/63, 4/63, 4/55, 4/39, 4/15],
     [4/17, 4/45, 4/65, 4/77, 4/81, 4/77, 4/65, 4/45, 4/17],
     [4/19, 4/51, 4/75, 4/91, 4/99, 4/99, 4/91, 4/75, 4/51, 4/19],
     [4/21, 4/57, 4/85, 4/105, 4/117, 4/121, 4/117, 4/105, 4/85, 4/57, 4/21]]

     (to be multiplied by the coefficients of the even part of the square
     of the de-skewed polynomial)
   */

#define USE_CIRCULAR
  if (d == 3)
    {
      double a3, a2, a1, a0;
      a3 = a[3] * s * s * s;
      a2 = a[2] * s * s;
      a1 = a[1] * s;
      a0 = a[0];
      if (method == RECTANGULAR)
        n = 4.0/7.0*(a3*a3+a0*a0)+8.0/15.0*(a1*a3+a0*a2)+4.0/15.0*(a2*a2+a1*a1);
      else /* circular method */
        {
          /* f := a3*x^3+a2*x^2+a1*x+a0:
             d := degree(f,x):
             F := expand(y^d*subs(x=x/y,f));
             subs(x=r*cos(t), y=r*sin(t), %);
             int(int(%^2/s^d*r, r=0..1),t=0..2*Pi); */
          n = 5.0 * (a3 * a3 + a0 * a0) + 2.0 * (a3 * a1 + a0 * a2)
            + a1 * a1 + a2 * a2;
          n = n * 0.049087385212340519352; /* Pi/64 */
        }
      return 0.5 * log(n / (s * s * s));
    }
  else if (d == 4)
    {
      double a4, a3, a2, a1, a0;

      a4 = a[4] * s * s * s * s;
      a3 = a[3] * s * s * s;
      a2 = a[2] * s * s;
      a1 = a[1] * s;
      a0 = a[0];
      if (method == RECTANGULAR)
        n = 4.0 / 9.0 * (a4 * a4 + a0 * a0) + 8.0 / 21.0 * (a4 * a2 + a2 * a0)
          + 4.0 / 21.0 * (a3 * a3 + a1 * a1) + 8.0 / 25.0 * (a4 * a0 + a3 * a1)
          + 4.0 / 25.0 * a2 * a2;
      else /* circular method */
        {
          /* f := a4*x^4+a3*x^3+a2*x^2+a1*x+a0:
             d := degree(f,x):
             F := expand(y^d*subs(x=x/y,f));
             subs(x=r*cos(t), y=r*sin(t), %);
             int(int(%^2/s^d*r, r=0..1),t=0..2*Pi); */
          n = 35.0 * (a4 * a4 + a0 * a0) + 10.0 * (a4 * a2 + a2 * a0)
            + 5.0 * (a3 * a3 + a1 * a1) + 6.0 * (a4 * a0 + a3 * a1)
            + 3.0 * a2 * a2;
          n = n * 0.0049087385212340519352; /* Pi/640 */
        }
      return 0.5 * log(n / (s * s * s * s));
    }
  else if (d == 5)
    {
      double a5, a4, a3, a2, a1, a0, s2 = s * s, s3 = s2 * s, s4 = s3 * s,
        s5 = s4 * s;

      /*
        f := a5*x^5+a4*x^4+a3*x^3+a2*x^2+a1*x+a0:
        F := expand(y^5*subs(x=x/y,f));
        int(int(subs(x=x*s,F)^2/s^5, x=-1..1), y=-1..1);
       */
      a0 = a[0];
      a1 = a[1] * s;
      a2 = a[2] * s2;
      a3 = a[3] * s3;
      a4 = a[4] * s4;
      a5 = a[5] * s5;
      if (method == RECTANGULAR)
        n = 4.0 / 11.0 * (a5 * a5 + a0 * a0) + 8.0 / 27.0 * (a5 * a3 + a2 * a0)
          + 4.0 / 27.0 * (a4 * a4 + a1 * a1) + 4.0 / 35.0 * (a3 * a3 + a2 * a2)
          + 8.0 / 35.0 * (a5 * a1 + a4 * a2 + a4 * a0 + a3 * a1);
      else
        {
          /* use circular integral:
             f := a5*x^5+a4*x^4+a3*x^3+a2*x^2+a1*x+a0:
             F := expand(y^5*subs(x=x/y,f));
             subs(x=r*cos(t), y=r*sin(t), %);
             int(int(%^2/s^5*r, r=0..1),t=0..2*Pi); */
          n = 6.0 * (a3 * a1 + a1 * a5 + a4 * a2 + a0 * a4)
            + 14.0 * (a0 * a2 + a3 * a5) + 63.0 * (a0 * a0 + a5 * a5)
            + 7.0 * (a4 * a4 + a1 * a1) + 3.0 * (a3 * a3 + a2 * a2);
          n = n * 0.0020453077171808549730; /* Pi/1536 */
        }
      return 0.5 * log(n / s5);
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
      if (method == RECTANGULAR)
        n = 4.0 / 13.0 * (a6 * a6 + a0 * a0) + 8.0 / 33.0 * (a6 * a4 + a2 * a0)
          + 4.0 / 33.0 * (a5 * a5 + a1 * a1) + 4.0 / 45.0 * (a4 * a4 + a2 * a2)
          + 8.0 / 45.0 * (a6 * a2 + a5 * a3 + a4 * a0 + a3 * a1)
          + 8.0 / 49.0 * (a6 * a0 + a5 * a1 + a4 * a2) + 4.0 / 49.0 * a3 * a3;
      else
        { /* use circular integral:
             f := a6*x^6+a5*x^5+a4*x^4+a3*x^3+a2*x^2+a1*x+a0:
             F := expand(y^6*subs(x=x/y,f));
             subs(x=r*cos(t), y=r*sin(t), %);
             int(int(%^2/s^6*r, r=0..1),t=0..2*Pi); */
          n = 231.0 * (a6 * a6 + a0 * a0) + 42.0 * (a6 * a4 + a2 * a0)
            + 21.0 * (a5 * a5 + a1 * a1) + 7.0 * (a4 * a4 + a2 * a2)
            + 14.0 * (a6 * a2 + a5 * a3 + a4 * a0 + a3 * a1)
            + 10.0 * (a6 * a0 + a5 * a1 + a4 * a2) + 5.0 * a3 * a3;
          n = n * 0.00043828022511018320850; /* Pi/7168 */
        }
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

   If method=RECTANGULAR: integrate over [-1,1] x [-1,1].
   If method=CIRCULAR: integrate over the unit circle.
*/
double
L2_lognorm (mpz_t *f, unsigned long d, double s, int method)
{
  double a[MAX_DEGREE + 1];
  unsigned long i;

  for (i = 0; i <= d; i++)
    a[i] = mpz_get_d (f[i]);
  return L2_lognorm_d (a, d, s, method);
}

/* returns the optimal skewness corresponding to L2_lognorm */
double
L2_skewness (mpz_t *f, int deg, int prec, int method)
{
  double s, n0, n1, a, b, c, d, nc, nd, fd[MAX_DEGREE + 1];
  int i;

  /* convert once for all to double's to avoid expensive mpz_get_d() */
  for (i = 0; i <= deg; i++)
    fd[i] = mpz_get_d (f[i]);

  /* first isolate the minimum in an interval [s, 2s] by dichotomy */
  
  s = 1.0;
  n0 = L2_lognorm_d (fd, deg, s, method);
  while ((n1 = L2_lognorm_d (fd, deg, 2 * s, method)) < n0)
    {
      s = 2.0 * s;
      n0 = n1;
    }

  /* We have L2(s/2) > L2(s) <= L2(2s)
     Assuming L2(s) is first decreasing, then increasing, the minimum is
     attained in [s/2, 2s]. */
  a = (s == 1.0) ? 1.0 : 0.5 * s;
  b = 2.0 * s;

  /* since we use trichotomy, the intervals shrink by 3/2 instead of 2 at each
     iteration, thus we multiply the precision (in bits) by log(2)/log(3/2) */
  prec = (int) (1.70951129135145 * (double) prec);

  /* use trichotomy */
  while (prec--)
    {
      c = (2.0 * a + b) / 3.0;
      d = (a + 2.0 * b) / 3.0;
      nc = L2_lognorm_d (fd, deg, c, method);
      nd = L2_lognorm_d (fd, deg, d, method);
      if (nc < nd) /* the minimum should be in [a,d] */
        b = d;
      else /* L2(c) > L2(d): the minimum should be in [c, b] */
        a = c;
    }

  s = (a + b) * 0.5;

  return s;
}

/********************* data structures for first phase ***********************/

m_logmu_t*
m_logmu_init (unsigned long Malloc)
{
  unsigned long i;
  m_logmu_t* M;

  M = (m_logmu_t*) malloc (Malloc * sizeof (m_logmu_t));
  for (i = 0; i < Malloc; i++)
    {
      mpz_init (M[i].b);
      mpz_init (M[i].m);
    }
  return M;
}

void
m_logmu_clear (m_logmu_t* M, unsigned long Malloc)
{
  unsigned long i;

  for (i = 0; i < Malloc; i++)
    {
      mpz_clear (M[i].b);
      mpz_clear (M[i].m);
    }
  free (M);
}

/* Insert (b, m, logmu) in the M database. Current implementation of M is
   a sorted list, with increasing values of logmu.
   alloc: number of allocated entries in M.
   size:  number of stored entries in M (size <= alloc).
   Changes size into min (size + 1, alloc).
   Returns non-zero if new entry is kept in the database.
*/
int
m_logmu_insert (m_logmu_t* M, unsigned long alloc, unsigned long *psize,
                mpz_t b, mpz_t m, double logmu, char *str)
{
  unsigned long size = *psize;

  ASSERT(size <= alloc);

  if (size < alloc)
    {
      mpz_set (M[size].b, b);
      mpz_set (M[size].m, m);
      M[size].logmu = logmu;
      *psize = size + 1;
      return 1;
    }
  else /* size=alloc: database is full, remove entry with smallest logmu */
    {
      if (logmu < M[alloc - 1].logmu)
        {
          size --;
          while (size > 0 && logmu < M[size - 1].logmu)
            {
              mpz_swap (M[size - 1].b, M[size].b);
              mpz_swap (M[size - 1].m, M[size].m);
              M[size].logmu = M[size - 1].logmu;
              size --;
            }
          /* now either size = 0 or M[size - 1].logmu <= logmu */
          mpz_set (M[size].b, b);
          mpz_set (M[size].m, m);
          M[size].logmu = logmu;
          if (size == 0)
            gmp_fprintf (stderr, "# p=%Zd m=%Zd %s%1.2f\n", b, m, str, logmu);
          return 1;
        }
      else
        return 0;
    }
}

/******************** base (p,m) decomposition *******************************/

/* round to nearest, assume d > 0 */
void
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

#if 0
/* Decompose p->n in base m/b, rounding to nearest
   (linear polynomial is b*x-m).
   For b=1, it suffices to check the coefficient f[d-1] is divisible by b,
   then all successive coefficients will be.
   Assumes b^2 fits in an unsigned long.
 */
void
generate_base_mb (cado_poly p, mpz_t m, mpz_t b)
{
  unsigned long i;
  int d = p->degree;

  ASSERT (d >= 0);

  if (mpz_cmp_ui (b, 1) == 0)
    {
      mpz_ndiv_qr (p->f[1], p->f[0], p->n, m);
      for (i = 1; i < (unsigned long) d; i++)
        mpz_ndiv_qr (p->f[i+1], p->f[i], p->f[i], m);
    }
  else /* b >= 2 */
    {
      mpz_t powm, im, k, m_mod_b;
      int s, tst;

      mpz_init (powm);
      mpz_init (im);
      mpz_init (k);
      mpz_init (m_mod_b);
      
      /* we need that m is invertible mod b */
      do
        {
          mpz_gcd (powm, m, b);
          if (mpz_cmp_ui (powm, 1) == 0)
            break;
          mpz_add_ui (m, m, 1);
        }
      while (1);

      mpz_pow_ui (powm, m, d);
      mpz_invert (im, powm, b); /* 1/m^d mod b */

      mpz_ndiv_qr (p->f[d], p->f[d - 1], p->n, powm);
      /* N = f[d] * m^d + f[d-1]. We want f[d-1] to be divisible by b,
         which may need to consider f[d]+k instead of f[d], i.e., we want
         f[d-1]-k*m^d = 0 (mod b), i.e., k = f[d-1]/m^d (mod b). */
      mpz_fdiv_r (k, p->f[d - 1], b);   /* f[d-1] mod b */
      mpz_mul (k, k, im);
      mpz_mod (k, k, b);                /* k = f[d-1]/m^d (mod b) */
      if (mpz_cmp_ui (k, 0) != 0)
        {
          int tst;

          s = mpz_sgn (p->f[d - 1]);
          /* if (s >= 0 and 2k >= b) or (s < 0 and 2k-2 >= b) then k <- k-b */
          mpz_mul_2exp (k, k, 1);
          if (s >= 0)
            tst = mpz_cmp (k, b) >= 0;
          else
            {
              mpz_sub_ui (k, k, 2);
              tst = mpz_cmp (k, b) >= 0;
              mpz_add_ui (k, k, 2);
            }
          mpz_div_2exp (k, k, 1);
          if (tst)
            mpz_sub (k, k, b);

          /* f[d] -> f[d]+l, f[d-1] -> f[d-1]-l*m^d */
          mpz_add (p->f[d], p->f[d], k);
          mpz_submul (p->f[d - 1], powm, k);
        }
      /* now N = f[d] * m^d + f[d-1], with f[d-1] divisible by b */
      mpz_fdiv_r (m_mod_b, m, b); /* m mod b */
      mpz_divexact (p->f[d - 1], p->f[d - 1], b);
      for (i = d - 1; i >= 1; i--)
        {
          /* p->f[i] is the current remainder */
          mpz_divexact (powm, powm, m); /* m^i */
          mpz_ndiv_qr (p->f[i], p->f[i - 1], p->f[i], powm);
          /* f[i] = q*m^i + r: we want to write f[i] = (q+k) * m^i + (r-k*m^i)
             with r-k*m^i divisible by b, i.e., k = r/m^i mod b */
          mpz_fdiv_r (k, p->f[i - 1], b); /* r mod b */
          mpz_mul (im, im, m_mod_b);
          mpz_mod (im, im, b); /* 1/m^i mod b */
          mpz_mul (k, k, im);
          mpz_mod (k, k, b); /* r/m^i mod b */

          s = mpz_sgn (p->f[i - 1]);
          mpz_mul_2exp (k, k, 1);
          if (s >= 0)
            tst = mpz_cmp (k, b) >= 0;
          else
            {
              mpz_sub_ui (k, k, 2);
              tst = mpz_cmp (k, b) >= 0;
              mpz_add_ui (k, k, 2);
            }
          mpz_div_2exp (k, k, 1);
          if (tst)
            mpz_sub (k, k, b);
          mpz_add (p->f[i], p->f[i], k);
          mpz_submul (p->f[i - 1], powm, k);
          mpz_divexact (p->f[i - 1], p->f[i - 1], b);
        }
      mpz_clear (powm);
      mpz_clear (im);
      mpz_clear (k);
      mpz_clear (m_mod_b);
    }
  mpz_neg (p->g[0], m);
  mpz_set (p->g[1], b);
  mpz_set (p->m, m);
  p->skew = L2_skewness (p->f, d, SKEWNESS_DEFAULT_PREC, 0);
}
#endif

/************************** polynomial arithmetic ****************************/

/* g <- content(f) where deg(f)=d */
static void
content_poly (mpz_t g, mpz_t *f, int d)
{
  int i;

  ASSERT(d >= 1);

  mpz_gcd (g, f[0], f[1]);
  for (i = 2; i <= d; i++)
    mpz_gcd (g, g, f[i]);
}

#if 0
/* v <- f(r), where f is of degree d */
static void
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
#endif

/* v <- f'(r), where f is of degree d */
static void
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
static void
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

/********************* computation of alpha **********************************/

/* Auxiliary routine for special_valuation(), see below. It returns the
   average p-valuation of the polynomial f. Works recursively. */
double
special_val0 (mpz_t *f, int d, unsigned long p)
{
  double v;
  mpz_t c, *g, *h;
  int i, r0, nroots;
  unsigned long *roots, r;
  mpz_array_t *G = NULL, *H;

  mpz_init (c);
  content_poly (c, f, d);
  for (v = 0.0; mpz_divisible_ui_p (c, p); v++, mpz_divexact_ui (c, c, p));

  /* g <- f/p^v */
  if (v != 0.0)
    {
      G = alloc_mpz_array (d + 1);
      g = G->data;
      mpz_ui_pow_ui (c, p, (unsigned long) v); /* p^v */
      for (i = 0; i <= d; i++)
        mpz_divexact (g[i], f[i], c);
    }
  else
    g = f;

  H = alloc_mpz_array (d + 1);
  h = H->data;
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
  clear_mpz_array (H);

  if (G != NULL)
    clear_mpz_array (G);
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
static double
special_valuation (mpz_t * f, int d, unsigned long p, mpz_t disc)
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
    } else if (pvaluation_disc == 1) {
      /* special case where p^2 does not divide disc */
	int e;
	e = poly_roots_ulong(NULL, f, d, p);
        if (p_divides_lc)
          e ++;
	/* something special here. */
	return (pd * e - 1) / (pd * pd - 1);
    } else {
	v = special_val0(f, d, p) * pd;
	if (p_divides_lc) {
	    /* compute g(x) = f(1/(px))*(px)^d, i.e., g[i] = f[d-i]*p^i */
	    /* IOW, the reciprocal polynomial evaluated at px */
            mpz_array_t *G;
	    mpz_t *g;
	    mpz_t t;
	    int i;

	    G = alloc_mpz_array (d + 1);
            g = G->data;
	    mpz_init_set_ui(t, 1);	/* will contains p^i */
	    for (i = 0; i <= d; i++) {
		mpz_mul(g[i], f[d - i], t);
		mpz_mul_ui(t, t, p);
	    }
	    v += special_val0(g, d, p);
	    clear_mpz_array (G);
            mpz_clear(t);
	}
	v /= pd + 1.0;
	return v;
    }
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

/**************************** rotation ***************************************/

/* D <- discriminant (f+k*g), which has degree d */
static void
discriminant_k (mpz_t *D, mpz_t *f, int d, mpz_t m, mpz_t b)
{
  int i, j, k;
  uint32_t **M, pivot;

  ASSERT_ALWAYS(d <= 9);

  /* we first put in D[i] the value of disc(f + i*g) for 0 <= i <= d,
     thus if disc(f + k*g) = a[d]*k^d + ... + a[0], then
          D[0] = a[0]
          D[1] = a[0] + a[1] + ... + a[d]
          ...
          D[d] = a[0] + a[1]*d + ... + a[d]*d^d */

  discriminant (D[0], f, d); /* a[0] */
  for (i = 1; i <= d; i++)
    {
      /* add b*x - m */
      mpz_add (f[1], f[1], b);
      mpz_sub (f[0], f[0], m);
      discriminant (D[i], f, d);
    }

  /* initialize matrix coefficients */
  M = (uint32_t**) malloc ((d + 1) * sizeof(uint32_t*));
  for (i = 0; i <= d; i++)
    {
      M[i] = (uint32_t*) malloc ((d + 1) * sizeof(uint32_t));
      M[i][0] = 1;
      for (j = 1; j <= d; j++)
        M[i][j] = i * M[i][j-1];
    }

  for (j = 0; j < d; j++)
    {
      /* invariant: D[i] = M[i][0] * a[0] + ... + M[i][d] * a[d]
         with M[i][k] = 0 for k < j and k < i */
      for (i = j + 1; i <= d; i++)
        {
          /* eliminate M[i][j] */
          pivot = M[i][j] / M[j][j];
          mpz_submul_ui (D[i], D[j], pivot);
          for (k = j; k <= d; k++)
            M[i][k] -= pivot * M[j][k];
        }
    }

  /* now we have an upper triangular matrix */
  for (j = d; j > 0; j--)
    {
      for (k = j + 1; k <= d; k++)
        mpz_submul_ui (D[j], D[k], M[j][k]);
      ASSERT_ALWAYS(mpz_divisible_ui_p (D[j], M[j][j]));
      mpz_divexact_ui (D[j], D[j], M[j][j]);
    }

  /* restore the original f[] */
  mpz_submul_ui (f[1], b, d);
  mpz_addmul_ui (f[0], m, d);

  for (i = 0; i <= d; i++)
    free (M[i]);
  free (M);
}

/* replace f + k0 * x^t * (b*x - m) by f + k * x^t * (b*x - m), and return k */
long
rotate_aux (mpz_t *f, mpz_t b, mpz_t m, long k0, long k, unsigned int t)
{
  mpz_addmul_si (f[t + 1], b, k - k0);
  mpz_submul_si (f[t], m, k - k0);
  return k;
}

/* replace f by f + k * x^t * (b*x + g0) */
static void
rotate_auxg_z (mpz_t *f, mpz_t b, mpz_t g0, mpz_t k, unsigned int t)
{
  mpz_addmul (f[t + 1], b, k);
  mpz_addmul (f[t], g0, k);
}

unsigned long
rotate_area (long K0, long K1, long J0, long J1)
{
  return (unsigned long) (J1 - J0 + 1) * (unsigned long) (K1 - K0 + 1);
}

/* Use Emmanuel Thome's idea: assuming alpha(f + k*g) admits a Gaussian
   distribution with mean 'm' and standard deviation 's', then the probability
   that one polynomial has a value of alpha >= A is
   1/2*(1 - erf((A-m)/s/sqrt(2))), thus the probability that K polynomials
   have all their alpha values >= A is [1/2*(1 - erf((A-m)/s/sqrt(2)))]^K.
   For 'm' and 's' fixed, and a given K, we define the value of A for which
   this probability is 1/2 to be the expected value of A for K polynomials.

   We assume lognorm(f + k*g) + E(alpha(f + k*g)) is first decreasing, then
   increasing, then the optimal K corresponds to the minimum of that function.
*/
void
rotate_bounds (mpz_t *f, int d, mpz_t b, mpz_t m, long *K0, long *K1,
               long *J0, long *J1, int verbose, int method)
{
  /* exp_alpha[i] corresponds to K=2^i polynomials for m=0, s=1
     f := x -> 1/2*(1 - erf(x/sqrt(2)));
     seq(fsolve(f(x)^(2^k) = 1/2, x=-log[2](1.0+k)), k=0..30);
   */
  double exp_alpha[] = {0.0 /* K=1 */, -0.5449521356 /* 2 */,
                        -0.9981488825 /* 4 */, -1.385198061 /* 8 */,
                        -1.723526050 /* 16 */, -2.025111894 /* 32 */,
                        -2.298313131 /* 64 */, -2.549067054 /* 128 */,
                        -2.781676726 /* 256 */, -2.999326227 /* 512 */,
                        -3.204420841 /* 1024 */, -3.398814100 /* 2048 */,
                        -3.583961388 /* 4096 */, -3.761025425 /* 8192 */,
                        -3.930949902 /* 16384 */, -4.094511733 /* 32768 */,
                        -4.252358774 /* 65536 */, -4.405037486 /* 131072 */,
                        -4.553013560 /* 262144 */, -4.696687518 /* 524288 */,
                        -4.836406692 /* 2^20 */, -4.972474538 /* 2^21 */,
                        -5.105157963 /* 2^22 */, -5.234693169 /* 2^23 */,
                        -5.361290351 /* 2^24 */, -5.485137511 /* 2^25 */,
                        -5.606403590 /* 2^26 */, -5.725241052 /* 2^27 */,
                        -5.841788041 /* 2^27 */, -5.956170181 /* 2^29 */,
                        DBL_MAX};
  int i;
  long k0 = 0, j0 = 0;
  double lognorm, alpha, E0, E, best_E;
  double skewness = L2_skewness (f, d, SKEWNESS_DEFAULT_PREC, method);
  long jmax = (long) ((double) (1L << MAX_k) / skewness);
  unsigned long max_area = 1UL << MAX_k;

#define MARGIN 0.12 /* we allow a small error margin in the expected lognorm
                       + alpha values, to get a larger search range */

  E0 = L2_lognorm (f, d, skewness, method);

  *K0 = -2;
  *J0 = -16;
  *J1 = 16;

  /* look for negative k: -2, -4, -8, ... */
  best_E = E0;
  for (i = 1; rotate_area (*K0, -*K0, *J0, *J1) < max_area; i++, *K0 *= 2)
    {
      k0 = rotate_aux (f, b, m, k0, *K0, 0);
      lognorm = L2_lognorm (f, d, L2_skewness (f, d, SKEWNESS_DEFAULT_PREC, method),
                         method);
      alpha = exp_alpha[i];
      E = lognorm + alpha;
      if (E < best_E + MARGIN)
        {
          if (E < best_E)
            best_E = E;
        }
      else
        break;
    }
  /* go back to k=0 */
  k0 = rotate_aux (f, b, m, k0, 0, 0);
  *K1 = -*K0;

  /* now try negative j: -1, -3, -7, ... */
  for (i++; exp_alpha[i] != DBL_MAX && rotate_area (*K0, *K1, *J0, -*J0)
	 < max_area; i++, *J0 = 2 * *J0 - 1)
    {
      j0 = rotate_aux (f, b, m, j0, *J0, 1);
      lognorm = L2_lognorm (f, d, L2_skewness (f, d, SKEWNESS_DEFAULT_PREC, method),
                         method);
      alpha = exp_alpha[i];
      E = lognorm + alpha;
      if (E < best_E + MARGIN)
        {
          if (E < best_E)
            best_E = E;
        }
      else
        break;
      if (1 - 2 * *J0 > jmax)
        break;
    }
  *J1 = -*J0;

  if (verbose)
    fprintf (stderr, "# Rotate bounds: %ld <= j <= %ld, %ld <= k <= %ld\n",
             *J0, *J1, *K0, *K1);

  /* rotate back to j=k=0 */
  rotate_aux (f, b, m, k0, 0, 0);
  rotate_aux (f, b, m, j0, 0, 1);
}

/* res <- f(k) */
static void
mpz_poly_eval_si (mpz_t res, mpz_t *f, int d, long k)
{
  int i;
  
  mpz_set (res, f[d]);
  for (i = d - 1; i >= 0; i--)
    {
      mpz_mul_si (res, res, k);
      mpz_add (res, res, f[i]);
    }
}

long
old_rotate (unsigned long p, double *A, long K0, long K1, long k0, mpz_t *f,
	    int d, mpz_t b, mpz_t m, mpz_array_t *D)
{
  unsigned long pp;
  long k, i;
  mpz_t v;
  double e, alpha;
  long l;
  double one_over_pm1 = 1.0 / (double) (p - 1);
  double logp = log ((double) p);

  mpz_init (v);
  for (pp = p * p, k = K0; (k < K0 + (long) p) && (k <= K1); k++)
    {
      /* translate from k0 to k */
      k0 = rotate_aux (f, b, m, k0, k, 0);
      
      /* compute contribution for k */
      mpz_poly_eval_si (v, D->data, d, k);
      e = special_valuation (f, d, p, v);
      alpha = (one_over_pm1 - e) * logp;
      
      /* and alpha is the contribution for k */
      if (!mpz_divisible_ui_p (v, p))
	{
	  /* then any k + t*p has the same contribution */
	  for (i = k; i <= K1; i += p)
	    A[i - K0] += alpha;
	}
      else
	{
	  /* consider classes mod p^2 */
	  for (i = k;;)
	    {
	      /* invariant: v = disc (f + i*g), and alpha is the
		 contribution for i */
	      
	      A[i - K0] += alpha;
	      if (!mpz_divisible_ui_p (v, pp))
		{
		  /* then any i + t*p^2 has the same contribution */
		  for (l = i + pp; l <= K1; l += pp)
		    A[l - K0] += alpha;
		}
	      else
		{
		  for (l = i + pp; l <= K1; l += pp)
		    {
		      mpz_poly_eval_si (v, D->data, d, l);
		      /* translate from k0 to l */
		      k0 = rotate_aux (f, b, m, k0, l, 0);
		      e = special_valuation (f, d, p, v);
		      alpha = (one_over_pm1 - e) * logp;
		      A[l - K0] += alpha;
		    }
		}
	      i += p;
	      if (i >= k + (long) pp || i > K1)
		break;
	      mpz_poly_eval_si (v, D->data, d, i);
	      /* translate from k0 to i */
	      k0 = rotate_aux (f, b, m, k0, i, 0);
	      e = special_valuation (f, d, p, v);
	      alpha = (one_over_pm1 - e) * logp;
	    }
	}
    }
  mpz_clear (v);
  return k0;
}

/* Return the smallest value of lognorm + alpha(f + (j*x+k)*(b*x-m)) for
   j and k small enough such that the norm does not increase too much, and
   modify f[] accordingly.
   The parameter "multi" means that several polynomials are wanted. If
   multi=0 or 1, then only 1 polynomial is returned (classical behavior).
   Otherwise, multi polynomials are stored in jmin and kmin (that
   must be initialized arrays with at least multi elements). This option
   might be useful for Coppersmith variant (a.k.a. MNFS). 
   In the multi case, the smallest of the returned values of lognorm + alpha
   is returned (and f[] accordingly).
   Warning: the caller is responsible to update the skewness if needed.
   */
double
rotate (mpz_t *f, int d, unsigned long alim, mpz_t m, mpz_t b,
        long *jmin, long *kmin, int multi, int verbose, int method)
{
  mpz_array_t *D;
  long K0, K1, J0, J1, k0, k, i, j, j0;
  double *A, alpha, lognorm, best_alpha = DBL_MAX, best_lognorm = DBL_MAX;
  double corr;
  double alpha0 = 0.0;
  unsigned long p;
  double *best_E = NULL; /* set to NULL to avoid warning... */
  double time_alpha = 0.0, time_norm = 0.0;
#ifdef NEW_ROOTSIEVE
  unsigned long pp;
  double one_over_pm1, logp, average_alpha = 0.0;
#endif

  /* allocate best_E, to store the best (lognorm+alpha) in multi mode */
  if (multi > 1)
    {
      best_E = (double *) malloc (multi * sizeof (double));
      for (i = 0; i < multi; ++i)
	best_E[i] = DBL_MAX;
    }

  /* allocate D(k) = disc(f + (j*x+k)*g, x) */
  D = alloc_mpz_array (d + 1);

  /* compute range for k */
  rotate_bounds (f, d, b, m, &K0, &K1, &J0, &J1, verbose, method);
  ASSERT_ALWAYS(K0 <= 0 && 0 <= K1);

  /* allocate sieving zone for computing alpha */
  A = (double*) malloc ((K1 + 1 - K0) * sizeof (double));
  j0 = k0 = 0; /* the current coefficients f[] correspond to f+(j*x+k)*g */
  
  *jmin = *kmin = 0;

  ASSERT_ALWAYS(J0 < 0 && 0 < J1);
  for (j = 0;;)
    {
      /* we consider j=0, 1, ..., J1, then J0, J0+1, ..., -1 */
      j0 = rotate_aux (f, b, m, j0, j, 1);
      /* go back to k=0 for the discriminant */
      k0 = rotate_aux (f, b, m, k0, 0, 0);
      /* D(k) = disc(f + (j*x+k)*g, x) (j is now fixed) */
      discriminant_k (D->data, f, d, m, b);

      for (k = K0; k <= K1; k++)
	A[k - K0] = 0.0; /* A[k - K0] will store the value alpha(f + k*g) */

  for (p = 2; p <= alim; p += 1 + (p & 1))
    if (isprime (p))
      {
	if (k0 != 0)
	  k0 = rotate_aux (f, b, m, k0, 0, 0);

        time_alpha -= seconds ();
#ifdef NEW_ROOTSIEVE
	one_over_pm1 = 1.0 / (double) (p - 1);
	logp = log ((double) p);
	for (pp = p; pp <= alim; pp *= p)
	  {
	    /* Murphy (page 48) defines cont_p(F) = q_p*p/(p^2-1)
	       = q_p*p/(p+1)*(1/p+1/p^2+...)
	       The contribution for p^k is thus q_p*p/(p+1)/p^k. */
	    alpha = logp / (double) (p + 1) * (double) p / (double) pp;
	    /* the following is the average contribution for a prime not
	       dividing the discriminant, cf alpha.sage, function alpha_p.
	       We take it into account only for p, not for p^2, p^3, ... */
	    if (p == pp)
	      average_alpha += logp * one_over_pm1;
	    /* we do not consider the projective roots here, since the
	       corresponding correction will be considered separately for each
	       row below */
	    /* + alpha_p_projective (f, d, (D->data)[0], p); */
	    update_table (f, d, m, b, A, K0, K1, pp, alpha);
	  }
#else
	k0 = old_rotate (p, A, K0, K1, k0, f, d, b, m, D);
#endif
	time_alpha += seconds ();
      } /* end of loop on primes p */

#ifdef NEW_ROOTSIEVE
  /* Correction to apply to the current row (takes into account the projective
     roots). FIXME: we are lazy here, we should only consider the contribution
     from the projective roots. */
  corr = get_alpha (f, d, alim) - A[k0 - K0];
#else
  corr = 0.0;
#endif

  if (j == 0)
    alpha0 = A[0 - K0] + corr;

  /* determine the best alpha in each row */
  {
    long bestk = K0;
    A[K0 - K0] += corr;
    for (k = K0 + 1; k <= K1; k++)
      {
	A[k - K0] += corr;
	if (A[k - K0] < A[bestk - K0])
	  bestk = k;
      }
    if (verbose > 1)
      fprintf (stderr, "# best alpha for j=%ld: k=%ld with %f\n",
	       j, bestk, A[bestk - K0]);
  }

  /* now finds the best lognorm+alpha */
  time_norm -= seconds ();
  for (k = K0; k <= K1; k++)
    {
      alpha = A[k - K0];
      if (alpha < best_alpha + 2.0)
        {
          /* check lognorm + alpha < best_lognorm + best_alpha */

          /* translate from k0 to k */
          k0 = rotate_aux (f, b, m, k0, k, 0);
          lognorm = L2_lognorm (f, d, L2_skewness (f, d, SKEWNESS_DEFAULT_PREC,
                                             method), method);
          if (multi <= 1) {
              if (lognorm + alpha < best_lognorm + best_alpha) {
                  best_lognorm = lognorm;
                  best_alpha = alpha;
                  *kmin = k;
                  *jmin = j;
              } 
          } else { /* multi mode */
              /* Rem: best_lognorm + best_alpha is the worse of the
                 preselected */
              double newE = lognorm + alpha;
              if (newE < best_E[multi-1]) {
                  int ii;
                  /* find position; assume list of preselected is sorted */
                  for (ii = 0; ii < multi; ++ii) {
                      if (best_E[ii] > newE)
                          break;
                  }
                  ASSERT_ALWAYS(ii < multi); 
                  /* insert */
                  for (i = multi - 1; i > ii; --i) {
                      kmin[i] = kmin[i-1];
                      jmin[i] = jmin[i-1];
                      best_E[i] = best_E[i-1];
                  }
                  kmin[ii] = k;
                  jmin[ii] = j;
                  best_E[ii] = newE;
              }
          }
        }
    }
  time_norm += seconds ();

  j++;
  if (j > J1)
    j = J0;
  else if (j == 0)
    break;
  
    } /* end of loop on j */

  /* we now have f + (j0*x+k0)*(bx-m) and we want f + (jmin*x+kmin)*(bx-m),
     thus we have to add ((jmin-j0)*x+(kmin-k0)*(bx-m) */
  /* if you are in multi-mode, we use the best polynomial */
  rotate_aux (f, b, m, k0, *kmin, 0);
  rotate_aux (f, b, m, j0, *jmin, 1);

  if (verbose && (multi <= 1))
    {
      fprintf (stderr, "# Rotate by ");
      if (*jmin != 0)
        {
          if (*jmin == -1)
            fprintf (stderr, "-");
          else if (*jmin != 1)
            fprintf (stderr, "%ld*", *jmin);
          fprintf (stderr, "x");
          if (*kmin >= 0)
            fprintf (stderr, "+");
        }
      fprintf (stderr, "%ld: alpha improved from %1.2f to %1.2f (alpha %1.2fs, norm %1.2fs)\n", *kmin, alpha0, best_alpha, time_alpha, time_norm);
    }

  if (verbose && (multi > 1)) {
      fprintf(stderr, "Found the following polynomials  (j, k, E):\n");
      for (i = 0; i < multi; ++i) {
          fprintf(stderr, "  %ld\t%ld\t%1.2f\n", jmin[i], kmin[i], best_E[i]);
      }
  }

  free (A);

  clear_mpz_array (D);

  {
      double ret_val = best_lognorm + best_alpha;
      if (multi>1) {
          ret_val = best_E[0]; /* we return the smallest */
          free(best_E);
      }
      return ret_val;
  }
}

/*****************************************************************************/

/* st is the initial value cputime() at the start of the program */
void
print_poly (FILE *fp, cado_poly p, int argc, char *argv[], double st, int raw)
{
  int i;
  double alpha, logmu;

  fprintf (fp, "n: ");
  mpz_out_str (fp, 10, p->n);
  fprintf (fp, "\n");
  fprintf (fp, "skew: %1.3f\n", p->skew);
  logmu = L2_lognorm (p->f, p->degree, p->skew, DEFAULT_L2_METHOD);
  alpha = get_alpha (p->f, p->degree, ALPHA_BOUND);
  fprintf (fp, "# lognorm: %1.2f, alpha: %1.2f E=%1.2f\n", logmu, alpha,
           logmu + alpha);
  fprintf (fp, "# Murphy's E(Bf=%.0f,Bg=%.0f,area=%.2e)=%1.2e\n",
	   BOUND_F, BOUND_G, AREA,
	   MurphyE (p, BOUND_F, BOUND_G, AREA, MURPHY_K));
  for (i = p->degree; i >= 0; i--)
    {
      fprintf (fp, "c%d: ", i);
      mpz_out_str (fp, 10, p->f[i]);
      fprintf (fp, "\n");
    }
#ifdef DEBUG
  fprintf (fp, "# ");
  fprint_polynomial (fp, p->f, p->degree);
#endif
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
  else
    mpz_neg (p->m, p->g[0]);
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
  fprintf (fp, "# generated by %s.r%s", argv[0], CADO_REV);
  for (i = 1; i < argc; i++)
    fprintf (fp, " %s", argv[i]);
  fprintf (fp, " in %.2fs\n", (seconds () - st));
}

/* Returns k such that f(x+k) has the smallest 1-norm, with the corresponding
   skewness.
   The coefficients f[] are modified to those of f(x+k).

   The linear polynomial b*x-m is changed into b*(x+k)-m, thus m is
   changed into m-k*b.
*/
long
translate (mpz_t *f, int d, mpz_t *g, mpz_t m, mpz_t b, int verbose,
           int method)
{
  double logmu0, logmu;
  int i, j, dir;
  long k;
  int prec = 2 * SKEWNESS_DEFAULT_PREC;

  logmu0 = L2_lognorm (f, d, L2_skewness (f, d, prec, method), method);

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
  logmu = L2_lognorm (f, d, L2_skewness (f, d, prec, method), method);

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
      logmu = L2_lognorm (f, d, L2_skewness (f, d, prec, method), method);
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

  if (verbose > 0)
    fprintf (stderr, "# Translate by %ld\n", k);

  /* change linear polynomial */
  mpz_neg (g[0], m);

  return k;
}

/* f <- f(x+k), g <- g(x+k) */
void
do_translate (mpz_t *f, int d, mpz_t *g, long k)
{
  int i, j;

  for (i = d - 1; i >= 0; i--)
    for (j = i; j < d; j++)
      mpz_addmul_si (f[j], f[j+1], k);
  mpz_addmul_si (g[0], g[1], k);
}

/* f <- f(x+k), g <- g(x+k) */
void
do_translate_z (mpz_t *f, int d, mpz_t *g, mpz_t k)
{
  int i, j;

  for (i = d - 1; i >= 0; i--)
    for (j = i; j < d; j++)
      mpz_addmul (f[j], f[j+1], k);
  mpz_addmul (g[0], g[1], k);
}

/* Use rotation and translation to find a polynomial with smaller norm
   (local minimum). Modify f and g accordingly.
   If skew = 0.0, recompute the skewness each time,
   otherwise optimize for the given skewness.
*/
static void
optimize_aux (mpz_t *f, int d, mpz_t *g, int verbose, double skew, int method)
{
  mpz_t k, l; /* current offset */
  int changed;
  double logmu00, logmu0, logmu;
  int prec = SKEWNESS_DEFAULT_PREC;
  int recompute = skew == 0.0;

  if (recompute)
    skew = L2_skewness (f, d, prec, method);
  logmu00 = logmu0 = L2_lognorm (f, d, skew, method);
  mpz_init_set_ui (k, 1);
  mpz_init (l);
  while (1)
    {
      changed = 0;

      mpz_mul_si (l, k, -2); /* l = -2k */

      /* first try translation by k */
      do_translate_z (f, d, g, k); /* f(x+k) */
      skew = L2_skewness (f, d, prec, method);
      logmu = L2_lognorm (f, d, skew, method);
      if (logmu < logmu0)
        {
          changed = 1;
          logmu0 = logmu;
        }
      else
        {
          do_translate_z (f, d, g, l); /* f(x-k) */
          if (recompute)
            skew = L2_skewness (f, d, prec, method);
          logmu = L2_lognorm (f, d, skew, method);
          if (logmu < logmu0)
            {
              changed = 1;
              logmu0 = logmu;
            }
          else
            do_translate_z (f, d, g, k); /* original f */
        }
      
      /* then do rotation by k*x^2*g if d >= 6 */
      if (d >= 6)
        {
          rotate_auxg_z (f, g[1], g[0], k, 2);
          if (recompute)
            skew = L2_skewness (f, d, prec, method);
          logmu = L2_lognorm (f, d, skew, method);
          if (logmu < logmu0)
            {
              changed = 1;
              logmu0 = logmu;
            }
          else
            {
              rotate_auxg_z (f, g[1], g[0], l, 2); /* f - k*x^2*g */
              if (recompute)
                skew = L2_skewness (f, d, prec, method);
              logmu = L2_lognorm (f, d, skew, method);
              if (logmu < logmu0)
                {
                  changed = 1;
                  logmu0 = logmu;
                }
              else
                rotate_auxg_z (f, g[1], g[0], k, 2); /* previous f */
            }
        }

      /* then do rotation by k*x*g */
      rotate_auxg_z (f, g[1], g[0], k, 1);
      if (recompute)
        skew = L2_skewness (f, d, prec, method);
      logmu = L2_lognorm (f, d, skew, method);
      if (logmu < logmu0)
        {
          changed = 1;
          logmu0 = logmu;
        }
      else
        {
          rotate_auxg_z (f, g[1], g[0], l, 1); /* f - k*x*g */
          if (recompute)
            skew = L2_skewness (f, d, prec, method);
          logmu = L2_lognorm (f, d, skew, method);
          if (logmu < logmu0)
            {
              changed = 1;
              logmu0 = logmu;
            }
          else
            rotate_auxg_z (f, g[1], g[0], k, 1); /* previous f */
        }

      /* then do rotation by k*g */
      rotate_auxg_z (f, g[1], g[0], k, 0);
      if (recompute)
        skew = L2_skewness (f, d, prec, method);
      logmu = L2_lognorm (f, d, skew, method);
      if (logmu < logmu0)
        {
          changed = 1;
          logmu0 = logmu;
        }
      else
        {
          rotate_auxg_z (f, g[1], g[0], l, 0); /* f - k*g */
          if (recompute)
            skew = L2_skewness (f, d, prec, method);
          logmu = L2_lognorm (f, d, skew, method);
          if (logmu < logmu0)
            {
              changed = 1;
              logmu0 = logmu;
            }
          else
            rotate_auxg_z (f, g[1], g[0], k, 0); /* previous f */
        }

      if (changed == 1)
        mpz_neg (k, l);               /* k <- 2k */
      else if (mpz_cmp_ui (k, 1) > 0)
        mpz_div_2exp (k, k, 1);       /* k <- k/2 */
      else /* changed = 0 and k = 1 */
        break;
    }
  mpz_clear (k);
  mpz_clear (l);

  if (verbose > 0) {
    gmp_fprintf (stderr, "# ad=%Zd: optimized lognorm from %.2f to %.2f\n",
		 f[d], logmu00, logmu0);
    if (verbose > 1) {
      fprintf (stderr, "# "); fprint_polynomial (stderr, f, d);
      fprintf (stderr, "# "); fprint_polynomial (stderr, g, d);
    }
  }
}

void
optimize (mpz_t *f, int d, mpz_t *g, int verbose)
{
  optimize_aux (f, d, g, verbose, 0.0, CIRCULAR);
  /* if we want to optimize for the rectangular method, it seems better to
     first optimize for the circular method */
  if (DEFAULT_L2_METHOD == RECTANGULAR)
    optimize_aux (f, d, g, verbose, 0.0, RECTANGULAR);
}
