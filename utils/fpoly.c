/* arithmetic on polynomial with double-precision coefficients */
#include "cado.h"

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>   /* for fabs */
#include "fpoly.h"

/* Evaluate the poynomial f of degree deg at point x */
double
fpoly_eval (const double *f, const int deg, const double x)
{
  double r;
  int i;

  r = f[deg];
  for (i = deg - 1; i >= 0; i--)
    r = r * x + f[i];
  
  return r;
}

/* assuming g(a)*g(b) < 0, and g has a single root in [a, b],
   refines that root by dichotomy with n iterations.
   Assumes sa is of same sign as g(a).
*/
double
fpoly_dichotomy (double *g, int d, double a, double b, double sa,
                 unsigned int n)
{
  double s;

  do
    {
      s = (a + b) * 0.5;
      if (fpoly_eval (g, d, s) * sa > 0)
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
fpoly_print (FILE *stream, const double *f, const int deg, char *name)
{
  int i;

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

