#include "cado.h"
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "double_poly.h"


/* Parse a space-separated string as polynomial coefficients.
   Exactly one space between coefficients, none at start or end of string.
   If poly is NULL, nothing is written, but number of coefficients is
   returned.
   Returns the number of parsed coefficients (not the degree!) */
static int
parse_poly_str (double_poly_ptr poly, const char *str)
{
  const char *next = str;
  int i = 0;
  double coeff;
  
  while (next != NULL) {
    sscanf(next, "%lf", &coeff);
    // printf("%lf\n", coeff);
    if (poly) {
      assert((unsigned int) i <= poly->deg);
      poly->coeff[i] = coeff;
    }
    i++;
    next = strchr(next, ' '); /* Pointer to next ' ' or NULL if not found */
    if (next != NULL)
      next++;
  }
  return i;
}

static void
parse_poly(double_poly_ptr poly, const char *str)
{
  int n = parse_poly_str(NULL, str);
  if (n == 0) {
    double_poly_init (poly, 0);
    poly->coeff[0] = 0.;
  } else {
    double_poly_init (poly, n - 1);
  }
  parse_poly_str(poly, str);
}

int cmp_double(const double d1, const double d2, const double err_margin)
{
  return fabs(d1) * (1. - err_margin) <= fabs(d2) && fabs(d2) <= fabs(d1) * (1. + err_margin);
}

void 
test_double_poly_compute_roots1(const char *poly_str, const char *roots_str, 
                               const double err_margin, const double s,
                               const int verbose)
{
  double_poly_t poly, roots1, roots2;
  int nr_roots, i;
  parse_poly(poly, poly_str);
  parse_poly(roots1, roots_str);
  double_poly_init(roots2, poly->deg);

  if (verbose)
    double_poly_print (stdout, poly, "Testing polynomial ");

  nr_roots = double_poly_compute_roots(roots2->coeff, poly, s);
  if (nr_roots != roots1->deg + 1) {
    fprintf (stderr, "double_poly_compute_roots() produced wrong number of roots %d, reference has %d\n",
             nr_roots, roots1->deg + 1);
    abort();
  }
  for (i = 0; i < nr_roots; i++) {
    if (!cmp_double(roots1->coeff[i], roots2->coeff[i], err_margin)) {
      fprintf (stderr, "double_poly_compute_roots() produced wrong roots %f, reference has %f\n",
               roots2->coeff[i], roots1->coeff[i]);
    
      abort();
    }
  }

  double_poly_clear(poly);
  double_poly_clear(roots1);
  double_poly_clear(roots2);
}

static void
test_double_poly_compute_roots(const int verbose)
{
  /* A few roots of 2 */
  test_double_poly_compute_roots1("-2 1", "2", 1e-9, 3., verbose);
  test_double_poly_compute_roots1("-2 0 1", "1.41421356237310", 1e-6, 3., verbose);
  test_double_poly_compute_roots1("-2 0 0 1", "1.25992104989487", 1e-6, 3., verbose);
  test_double_poly_compute_roots1("-2 0 0 0 1", "1.18920711500272", 1e-6, 3., verbose);
  test_double_poly_compute_roots1("-2 0 0 0 0 1", "1.14869835499704", 1e-6, 3., verbose);

  /* (x-1)*(x-2) */
  test_double_poly_compute_roots1("2 -3 1", "1 2", 1e-6, 3., verbose);
  /* (x-1)*(x-2)*(x-3) */
  test_double_poly_compute_roots1("-6 11 -6 1", "1 2 3", 1e-6, 4., verbose);
  /* (x-1)*(x-2)*(x-3)*(x-4) */
  test_double_poly_compute_roots1("24 -50 35 -10 1", "1 2 3 4", 1e-6, 5., verbose);
  /* (x-1)*(x-2)*(x-3)*(x-4)*(x-5) */
  test_double_poly_compute_roots1("-120 274 -225 85 -15 1", "1 2 3 4 5", 1e-6, 6., verbose);
  
  /* Let f(x+1/x) * x^6 == (x^13-1)/(x-1). Test the positive roots */
  test_double_poly_compute_roots1("-1 3 6 -4 -5 1 1", "0.241073360510646, 1.13612949346231, 1.77091205130642", 1e-6, 2., verbose);
  /* and the negative ones */
  test_double_poly_compute_roots1("-1 3 6 -4 -5 1 1", "-0.709209774085071 -1.49702149634220 -1.94188363485210", 1e-6, -2., verbose);
  
  /* this is f(x-2). 6 real roots, 0 rational */
  test_double_poly_compute_roots1("1, -21, 70, -84, 45, -11, 1", "0.0581163651478959 0.502978503657798 1.29079022591493 2.24107336051065 3.13612949346231 3.77091205130642", 1e-6, 4., verbose);
  
}


int main()
{
  test_double_poly_compute_roots(0);
  exit(EXIT_SUCCESS);
}
