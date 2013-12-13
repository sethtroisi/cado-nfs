/* arithmetic on polynomial with double-precision coefficients */
#include "cado.h"

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>   /* for fabs */
#include "macros.h"
#include "double_poly.h"
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
                       unsigned int n)
{
  double s;

  do
    {
      s = (a + b) * 0.5;
      if (double_poly_eval (p, s) * sa > 0)
	a = s;
      else
	b = s;
    }
  while (n-- > 0);
  return (a + b) * 0.5;
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
    fprintf (stream, "%f", f[0]);

  if (deg == 1)
    fprintf (stream, "%f*x", f[1]);

  if (deg > 1)
    fprintf (stream, "%f*x^%d", f[deg], deg);

  for (i = deg - 1; i >= 0; i--)
    {
      if (f[i] == 0.)
	continue;
      if (i == 0)
	fprintf (stream, " %s %f", (f[i] > 0) ? "+" : "-", fabs(f[i]));
      else if (i == 1)
	fprintf (stream, " %s %f*x", (f[i] > 0) ? "+" : "-", fabs(f[i]));
      else 
	fprintf (stream, " %s %f*x^%d", (f[i] > 0) ? "+" : "-", fabs(f[i]), i);
    }

  fprintf (stream, "\n");
}
