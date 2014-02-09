#include "cado.h"
#include "utils.h"
#include "polyselect/auxiliary.h"
#include "tests_common.h"

int
check_num (double x, double y, double emax)
{
  double e = fabs (x - y) / fabs (y);

  if (e > emax)
    {
      printf ("expected %.16e, got %.16e (rel. error %e)\n", y, x, e);
      return 0;
    }
  return 1;
}

static void
test_L2_lognorm (void)
{
  mpz_poly_t p;
  double n;

  mpz_poly_init (p, MAXDEGREE);

  /* degree 1 */
  mpz_set_ui (p->coeff[0], 1);
  mpz_set_ui (p->coeff[1], 2);
  p->deg = 1;
  n = L2_lognorm (p, 1.0);
  check_num (n, 0.68393671871075337595, 2e-10);
  n = L2_lognorm (p, 2.0);
  check_num (n, 0.94925084419690070780, 2e-10);

  /* degree 2 */
  mpz_set_ui (p->coeff[2], 3);
  p->deg = 2;
  n = L2_lognorm (p, 1.0);
  check_num (n, 0.82777775480769542870, 2e-10);
  n = L2_lognorm (p, 2.0);
  check_num (n, 1.3718482492081025724, 2e-10);

  /* degree 3 */
  mpz_set_si (p->coeff[3], -4);
  p->deg = 3;
  n = L2_lognorm (p, 1.0);
  check_num (n, 0.73159180848396739500, 2e-10);
  n = L2_lognorm (p, 2.0);
  check_num (n, 1.7170713330509004847, 2e-10);

  /* degree 4 */
  mpz_set_si (p->coeff[4], 5);
  p->deg = 4;
  n = L2_lognorm (p, 1.0);
  check_num (n, 0.88625243226159843905, 2e-10);
  n = L2_lognorm (p, 2.0);
  check_num (n, 2.1476529791715452362, 2e-10);

  /* degree 5 */
  mpz_set_si (p->coeff[5], -6);
  p->deg = 5;
  n = L2_lognorm (p, 1.0);
  check_num (n, 0.90490889517894624795, 2e-10);
  n = L2_lognorm (p, 2.0);
  check_num (n, 2.5284494790416146421, 2e-10);

  /* degree 6 */
  mpz_set_si (p->coeff[6], 7);
  p->deg = 6;
  n = L2_lognorm (p, 1.0);
  check_num (n, 0.94130996783641860020, 2e-10);
  n = L2_lognorm (p, 2.0);
  check_num (n, 2.9064780318501317870, 2e-10);

  /* degree 7 */
  mpz_set_si (p->coeff[7], -8);
  p->deg = 7;
  n = L2_lognorm (p, 1.0);
  check_num (n, 0.95405307631727864405, 2e-10);
  n = L2_lognorm (p, 2.0);
  check_num (n, 3.2796366844822268800, 2e-10);

  mpz_poly_clear (p);
}

/* t=0: generate polynomial with skewness < 1
   t=1: generate polynomial with skewness > 1
   t=2: generate random polynomial */
void
test_L2_skewness (int t)
{
  mpz_poly_t p;
  int d, i;
  double s, n, sl, nl, sh, nh, eps;
  int prec = 10;

  mpz_poly_init (p, MAXDEGREE);
  if (t == 0)
    mpz_set_ui (p->coeff[0], 1);
  else if (t == 1)
    mpz_set_ui (p->coeff[0], 4294967295UL);
  else
    mpz_set_ui (p->coeff[0], lrand48 ());
  for (d = 1; d <= 7; d++)
    {
      if ((t == 0 || t == 1) && d > 1)
        mpz_set_ui (p->coeff[d-1], 1);
      if (t == 0)
        mpz_set_ui (p->coeff[d], 4294967295UL);
      else if (t == 1)
        mpz_set_ui (p->coeff[0], 1);
      else
        {
          do
            mpz_set_ui (p->coeff[d], lrand48 ());
          while (mpz_cmp_ui (p->coeff[d], 0) == 0);
        }
      p->deg = d;
      s = L2_skewness (p, prec);
      n = L2_lognorm (p, s);
      eps = ldexp (fabs (n), -prec);
      /* check that skewness to the left and to the right is worse */
      for (i = 0; i < 53; i++)
        {
          sl = s - ldexp (s, i - 53);
          nl = L2_lognorm (p, sl);
          if (nl < n - eps)
            {
              printf ("Non-optimal skewness for polynomial\n");
              mpz_poly_fprintf (stdout, p);
              printf ("For skewness %.16e, norm is %.16e\n", s, n);
              printf ("For skewness %.16e, norm is %.16e\n", sl, nl);
              abort ();
            }
          sh = s + ldexp (s, i - 53);
          nh = L2_lognorm (p, sh);
          if (nh < n - eps)
            {
              printf ("Non-optimal skewness for polynomial\n");
              mpz_poly_fprintf (stdout, p);
              printf ("For skewness %.16e, norm is %.16e\n", s, n);
              printf ("For skewness %.16e, norm is %.16e\n", sh, nh);
              abort ();
            }
        }
    }
  mpz_poly_clear (p);
}

int
main (int argc, const char *argv[])
{
  tests_common_cmdline (&argc, &argv, PARSE_SEED | PARSE_ITER);
  test_L2_lognorm ();
  test_L2_skewness (0);
  test_L2_skewness (1);
  test_L2_skewness (2);
  tests_common_clear ();
  exit (EXIT_SUCCESS);
}
