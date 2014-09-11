/* arithmetic on polynomial with double-precision coefficients */
#include "cado.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>   /* for fabs */
#include <float.h> /* for DBL_MAX */
#include "utils.h"
#include "portability.h"

/* Initialize a polynomial of degree d */
void
double_poly_init (double_poly_ptr p, unsigned int d)
{
  p->coeff = malloc ((d + 1) * sizeof (double));
  FATAL_ERROR_CHECK(p->coeff == NULL, "malloc failed");
  p->deg = d;
}

/* Clear a polynomial */
void
double_poly_clear (double_poly_ptr p)
{
  free (p->coeff);
}

/* Set r = s. Assumes r has enough memory allocated. */
void
double_poly_set (double_poly_ptr r, double_poly_srcptr s)
{
  r->deg = s->deg;
  for (unsigned int i = 0; i <= s->deg; i++) {
    r->coeff[i] = s->coeff[i];
  }
}

/* Evaluate the polynomial p at point x */
double
double_poly_eval (double_poly_srcptr p, const double x)
{
  double r;
  unsigned int k;
  const double *f = p->coeff;
  const unsigned int deg = p->deg;

  switch (deg) {
  case 0: return f[0];
  case 1: return f[0]+x*f[1];
  case 2: return f[0]+x*(f[1]+x*f[2]);
  case 3: return f[0]+x*(f[1]+x*(f[2]+x*f[3]));
  case 4: return f[0]+x*(f[1]+x*(f[2]+x*(f[3]+x*f[4])));
  case 5: return f[0]+x*(f[1]+x*(f[2]+x*(f[3]+x*(f[4]+x*f[5]))));
  case 6: return f[0]+x*(f[1]+x*(f[2]+x*(f[3]+x*(f[4]+x*(f[5]+x*f[6])))));
  case 7: return f[0]+x*(f[1]+x*(f[2]+x*(f[3]+x*(f[4]+x*(f[5]+x*(f[6]+x*f[7]))))));
  case 8: return f[0]+x*(f[1]+x*(f[2]+x*(f[3]+x*(f[4]+x*(f[5]+x*(f[6]+x*(f[7]+x*f[8])))))));
  case 9: return f[0]+x*(f[1]+x*(f[2]+x*(f[3]+x*(f[4]+x*(f[5]+x*(f[6]+x*(f[7]+x*(f[8]+x*f[9]))))))));
  default: for (r = f[deg], k = deg - 1; k != UINT_MAX; r = r * x + f[k--]); return r;
  }
}

/* assuming g(a)*g(b) < 0, and g has a single root in [a, b],
   refines that root by dichotomy with n iterations.
   Assumes sa is of same sign as g(a).
*/
double
double_poly_dichotomy (double_poly_srcptr p, double a, double b, double sa,
                       unsigned int n MAYBE_UNUSED)
{
  double s;

  for(;;) {
#if defined(__i386)
      /* See comment in utils/usp.c on this. We want to avoid comparison
       * involving an extended precision double ! */
      { volatile double ms = (a + b) * 0.5; s = ms; }
#else
      s = (a + b) * 0.5;
#endif
      if (s == a || s == b) return s;
      if (double_poly_eval (p, s) * sa > 0)
	a = s;
      else
	b = s;
  }
  return s;
}

/* assuming g(a)*g(b) < 0, and g has a single root in [a, b],
   refines that root by the weighted false position method
   Assumes sa is of same sign as g(a).
*/
double
double_poly_falseposition (double_poly_srcptr p, double a, double b, double pa)
{
  double pb;
  int side=0;

  pb = double_poly_eval(p, b);

  for(;;) {
      double s, middle;
#if defined(__i386)
      /* See above */
      { volatile double ms = (a*pb-b*pa)/(pb-pa); s = ms; }
      { volatile double ms = (a + b) * 0.5; middle = ms; }
#else
      s = (a*pb-b*pa)/(pb-pa);
      middle = (a + b) * 0.5;
#endif
      /* It may happen that because of overflow, (a*pb-b*pa)/(pb-pa)
       * reaches s==a or s==b too early. If it so happens that we're
       * doing this, while the middle cut doesn't behave this way, use
       * the middle cut instead.
       *
       * Note that almost by design, this countermeasure also cancels
       * some of the benefit of the false position method.
       */
      if (s < a || s > b || ((s == a || s == b) && !(middle == a || middle == b)))
          s = middle;
      if (s == a || s == b) return s;
      double ps = double_poly_eval (p, s);
      if (ps * pa > 0) {
          a = s; pa = ps;
          if (side==1) pb /= 2;
          side=1;
      } else {
          b = s; pb = ps;
          if (side==-1) pa /= 2;
          side=-1;
      }
      ASSERT_ALWAYS (!isnan(b));
  }
}

/* Stores the derivative of f in df. f and df may be identical.
   Assumes df has been initialized with degree at least f->deg-1. */
void
double_poly_derivative(double_poly_ptr df, double_poly_srcptr f)
{
  unsigned int n;
  if (f->deg == 0) {
    df->deg = 0; /* How do we store deg -\infty polynomials? */
    df->coeff[0] = 0.;
    return;
  }
  // at this point, f->deg >=1
  for (n = 1; n <= f->deg; n++)
    df->coeff[n - 1] = f->coeff[n] * (double) n;
  df->deg = f->deg - 1;
}

/* Stores the product of f and g in h (h = f * g).
   Assumes h != f && h != g (f *= g is not accepted for example).
   Assumes h has been initialized with degree at least f->deg + g->deg.
*/
void
double_poly_product(double_poly_ptr h, double_poly_srcptr f, double_poly_srcptr g)
{
  ASSERT(h->coeff != f->coeff && h->coeff != g->coeff);
  ASSERT(h->deg >= f->deg + g->deg);
  h->deg = f->deg + g->deg;
  memset (h->coeff, 0, sizeof(double) * (h->deg + 1));
  for (size_t i_f = f->deg + 1; i_f--; ) {
    double fcoeff = f->coeff[i_f], *hcoeff = h->coeff + i_f;
    for (size_t i_g = g->deg + 1; i_g--; )
      hcoeff[i_g] += fcoeff * g->coeff[i_g];
  }
}

/* Stores the sum of f and g in h (h = f + g).
   Assumes h has been initialized with degree at least MAX(f->deg, g->deg).
*/
void
MAYBE_UNUSED double_poly_sum(double_poly_ptr h, double_poly_srcptr f, double_poly_srcptr g)
{
  size_t i;
  if (f->deg <= g->deg) {
    ASSERT(h->deg >= g->deg);
    h->deg = g->deg;
    for (i = 0; i <= f->deg; ++i) h->coeff[i] = f->coeff[i] + g->coeff[i];
    for (     ; i <= g->deg; ++i) h->coeff[i] = g->coeff[i];
  } else {
    ASSERT(h->deg >= f->deg);
    h->deg = f->deg;
    for (i = 0; i <= g->deg; ++i) h->coeff[i] = f->coeff[i] + g->coeff[i];
    for (     ; i <= f->deg; ++i) h->coeff[i] = f->coeff[i];
  }
}

/* Stores the subtraction of g to f in h (h = f - g).
   Assumes h has been initialized with degree at least MAX(f->deg, g->deg).
*/
void
double_poly_subtract(double_poly_ptr h, double_poly_srcptr f, double_poly_srcptr g)
{
  size_t i;
  if (f->deg <= g->deg) {
    ASSERT(h->deg >= g->deg);
    h->deg = g->deg;
    for (i = 0; i <= f->deg; ++i) h->coeff[i] = f->coeff[i] - g->coeff[i];
    for (     ; i <= g->deg; ++i) h->coeff[i] =             - g->coeff[i];
  } else {
    ASSERT(h->deg >= f->deg);
    h->deg = f->deg;
    for (i = 0; i <= g->deg; ++i) h->coeff[i] = f->coeff[i] - g->coeff[i];
    for (     ; i <= f->deg; ++i) h->coeff[i] = f->coeff[i];
  }
}

/* Revert the coefficients in-place: f(x) => f(1/x) * x^degree */
void
double_poly_revert (double_poly_ptr f)
{
  const unsigned int d = f->deg;

  if (d <= 0)
    return;

  /* if d is even, nothing to do for k=d/2 */
  for (unsigned int k = 0; k <= (d - 1) / 2; k++)
    {
      double tmp = f->coeff[k];
      f->coeff[k] = f->coeff[d - k];
      f->coeff[d - k] = tmp;
    }
}

/* Change sign of variable: x -> -x */
static void
double_poly_neg_x (double_poly_ptr r, double_poly_srcptr s)
{
  unsigned int i;
  r->deg = s->deg;
  for (i = 0; i < s->deg; i += 2) {
    r->coeff[i]     = s->coeff[i]; /* Even power coeff */
    r->coeff[i + 1] = -s->coeff[i + 1]; /* Odd power */
  }
  if (i == s->deg) /* iff deg is even */
    r->coeff[i] = s->coeff[i];
}

static unsigned int
recurse_roots(double_poly_srcptr poly, double *roots,
              const unsigned int sign_changes, const double s)
{
  unsigned int new_sign_changes = 0;
  if (poly->deg <= 0) {
      /* A constant polynomial (degree 0 or -\infty) has no sign changes */
  } else if (poly->deg == 1) {
      /* A polynomial of degree 1 can have at most one sign change in (0, s),
         this happens iff poly(0) = poly[0] and poly(s) have different signs */
      if (poly->coeff[0] * double_poly_eval(poly, s) < 0) {
          new_sign_changes = 1;
          roots[0] = - poly->coeff[0] / poly->coeff[1];
      }
  } else {
      /* invariant: sign_changes is the number of sign changes of the
         (k+1)-th derivative, with corresponding roots in roots[0]...
         roots[sign_changes-1], and roots[sign_changes] = s. */
      double a = 0.0;
      double va = poly->coeff[0]; /* value of poly at x=0 */
      for (unsigned int l = 0; l <= sign_changes; l++)
        {
          /* b is a root of dg[k+1], or s, the end of the interval */
          const double b = (l < sign_changes) ? roots[l] : s;
          const double vb = double_poly_eval (poly, b);
          if (va * vb < 0) /* root in interval [va, vb] */
            roots[new_sign_changes++] = double_poly_falseposition (poly, a, b, va);
          a = b;
          va = vb;
        }
  }

  return new_sign_changes;
}

/* Return a bound on the positive roots of p.
   Assume the leading coefficient of p is positive, then for a positive root
   r we have
   p[d]*r^d + ... + p[1]*r + p[0] = 0 thus
   p[d]*r^d <= -p[d-1]*r^(d-1) - ... - p[1]*r - p[0]
            <= max(-p[d-1],0)*r^(d-1) + ... + max(-p[1],0)*r + max(-p[0],0)
   thus q(r) <= 0 where q is the degree-d polynomial formed from p as follows:
   * q[d] = p[d]
   * q[i] = p[i] if p[i] < 0, and 0 otherwise for 0 <= i < d.
   Since q has a unique positive root, say r0, and q(r) < 0 iff r < r0,
   then positive roots of p are bounded by r0.
*/
double
double_poly_bound_roots (double_poly_srcptr p)
{
  unsigned int d = p->deg, i;
  double_poly_t q;
  double s;

  s = p->coeff[d] > 0 ? 1.0 : -1.0;
  double_poly_init (q, d);
  for (i = 0; i < d; i++)
    q->coeff[i] = (s * p->coeff[i] < 0) ? s * p->coeff[i] : 0.0;
  q->coeff[d] = fabs (p->coeff[d]);
  s = 1.0;
  while (double_poly_eval (q, s) < 0)
    s = s + s;
  double_poly_clear (q);
  return s;
}

/* put in roots[0], roots[1], ..., roots[k-1] all k roots of poly in [0, s],
   and return the number k of roots */
unsigned int
double_poly_compute_roots(double *roots, double_poly_srcptr poly, double s)
{
  const unsigned int d = poly->deg;
  double_poly_t *dg; /* derivatives of poly */
  
  /* The roots of the zero polynomial are ill-defined. Bomb out */
  ASSERT_ALWAYS(d > 0 || poly->coeff[0] != 0.);

  /* Handle constant polynomials separately */
  if (d == 0)
    return 0; /* Constant non-zero poly -> no roots */

  dg = (double_poly_t *) malloc (d * sizeof (double_poly_t));
  FATAL_ERROR_CHECK(dg == NULL, "malloc failed");

  dg[0]->deg = poly->deg;
  dg[0]->coeff = poly->coeff;
  
  for (unsigned int k = 1; k < d; k++) {
    /* dg[k] is the k-th derivative, thus has degree d-k, i.e., d-k+1
       coefficients */
    double_poly_init (dg[k], d - k);
    double_poly_derivative (dg[k], dg[k - 1]);
  }
  
  unsigned int sign_changes = 0;
  for (unsigned int k = d; k > 0; k--)
    sign_changes = recurse_roots(dg[k - 1], roots, sign_changes, s);

  for (unsigned int k = 1; k < d; k++)
    double_poly_clear (dg[k]);
  free (dg);

  return sign_changes;
}

/* compute all roots whose absolute value is <= B */
unsigned int
double_poly_compute_all_roots_with_bound (double *roots,
                                          double_poly_srcptr poly,
                                          double B)
{
  /* Positive roots */
  double bound = double_poly_bound_roots (poly);
  if (B < bound)
    bound = B;
  unsigned int nr_roots_pos = double_poly_compute_roots (roots, poly, bound);
  /* Negative roots */
  double_poly_t t; /* Copy of poly which gets sign-flipped */
  double_poly_init (t, poly->deg);
  double_poly_neg_x (t, poly);
  bound = double_poly_bound_roots (t);
  if (B < bound)
    bound = B;
  unsigned int nr_roots_neg =
    double_poly_compute_roots (roots + nr_roots_pos, t, bound);
  double_poly_clear(t);
  /* Flip sign of negative roots */
  for (unsigned int i = 0; i < nr_roots_neg; i++)
    roots[nr_roots_pos + i] *= -1.;
  return nr_roots_pos + nr_roots_neg;
}

unsigned int
double_poly_compute_all_roots (double *roots, double_poly_srcptr poly)
{
  return double_poly_compute_all_roots_with_bound (roots, poly, DBL_MAX);
}

/* Print polynomial with floating point coefficients. Assumes f[deg] != 0
   if deg > 0. */
void 
double_poly_print (FILE *stream, double_poly_srcptr p, char *name)
{
  int i;
  const double *f = p->coeff;
  const unsigned int deg = p->deg;

  fprintf (stream, "%s", name);

  if (deg == 0)
    fprintf (stream, "%.16e", f[0]);

  if (deg == 1)
    fprintf (stream, "%.16e*x", f[1]);

  if (deg > 1)
    fprintf (stream, "%.16e*x^%d", f[deg], deg);

  for (i = deg - 1; i >= 0; i--)
    {
      if (f[i] == 0.)
	continue;
      if (i == 0)
	fprintf (stream, " %s %.16e", (f[i] > 0) ? "+" : "-", fabs(f[i]));
      else if (i == 1)
	fprintf (stream, " %s %.16e*x", (f[i] > 0) ? "+" : "-", fabs(f[i]));
      else 
	fprintf (stream, " %s %.16e*x^%d", (f[i] > 0) ? "+" : "-", fabs(f[i]), i);
    }

  fprintf (stream, "\n");
}

void
double_poly_set_mpz_poly (double_poly_ptr p, mpz_poly_ptr q)
{
  unsigned int d = q->deg, i;

  for (i = 0; i <= d; i++)
    p->coeff[i] = mpz_get_d (q->coeff[i]);
  p->deg = d;
}
