#include "cado.h"
#include "utils.h"
#include "polyselect/auxiliary.h"
#include "tests_common.h"

int
check_num (double x, double y)
{
  double e = fabs (x - y) / fabs (y);

  if (e > 2e-10)
    {
      printf ("expected %.16e, got %.16e (rel. error %e)\n", y, x, e);
      return 0;
    }
  return 1;
}

static void
test_L2_lognorm ()
{
  mpz_poly_t p;
  double n;

  mpz_poly_init (p, MAXDEGREE);

  /* degree 1 */
  mpz_set_ui (p->coeff[0], 1);
  mpz_set_ui (p->coeff[1], 2);
  p->deg = 1;
  n = L2_lognorm (p, 1.0);
  check_num (n, 0.68393671871075337595);
  n = L2_lognorm (p, 2.0);
  check_num (n, 0.94925084419690070780);

  /* degree 2 */
  mpz_set_ui (p->coeff[2], 3);
  p->deg = 2;
  n = L2_lognorm (p, 1.0);
  check_num (n, 0.82777775480769542870);
  n = L2_lognorm (p, 2.0);
  check_num (n, 1.3718482492081025724);

  /* degree 3 */
  mpz_set_si (p->coeff[3], -4);
  p->deg = 3;
  n = L2_lognorm (p, 1.0);
  check_num (n, 0.73159180848396739500);
  n = L2_lognorm (p, 2.0);
  check_num (n, 1.7170713330509004847);

  /* degree 4 */
  mpz_set_si (p->coeff[4], 5);
  p->deg = 4;
  n = L2_lognorm (p, 1.0);
  check_num (n, 0.88625243226159843905);
  n = L2_lognorm (p, 2.0);
  check_num (n, 2.1476529791715452362);

  /* degree 5 */
  mpz_set_si (p->coeff[5], -6);
  p->deg = 5;
  n = L2_lognorm (p, 1.0);
  check_num (n, 0.90490889517894624795);
  n = L2_lognorm (p, 2.0);
  check_num (n, 2.5284494790416146421);

  /* degree 6 */
  mpz_set_si (p->coeff[6], 7);
  p->deg = 6;
  n = L2_lognorm (p, 1.0);
  check_num (n, 0.94130996783641860020);
  n = L2_lognorm (p, 2.0);
  check_num (n, 2.9064780318501317870);

  /* degree 7 */
  mpz_set_si (p->coeff[7], -8);
  p->deg = 7;
  n = L2_lognorm (p, 1.0);
  check_num (n, 0.95405307631727864405);
  n = L2_lognorm (p, 2.0);
  check_num (n, 3.2796366844822268800);

  mpz_poly_clear (p);
}

int
main (int argc, const char *argv[])
{
  tests_common_cmdline (&argc, &argv, PARSE_SEED | PARSE_ITER);
  test_L2_lognorm ();
  tests_common_clear ();
  exit (EXIT_SUCCESS);
}
