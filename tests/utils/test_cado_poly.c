#include "cado.h"
#include <time.h>
#include "tests_common.h"
#include "cado_poly.h"

void
test_cado_poly_set ()
{
  cado_poly p, q;

  cado_poly_init (p);
  cado_poly_init (q);

  const double s = 3.1415;

  mpz_set_ui (q->n, 123456789);
  q->skew = s;
  q->pols[0]->deg = 1;
  mpz_set_ui (q->pols[0]->coeff[0], 17);
  mpz_set_ui (q->pols[0]->coeff[1], 42);
  q->pols[1]->deg = 0;
  mpz_set_ui (q->pols[1]->coeff[0], 59);
  mpz_set_ui (q->m, 3613);

  cado_poly_set (p, q);

  ASSERT_ALWAYS (mpz_cmp_ui (p->n, 123456789) == 0);
  ASSERT_ALWAYS (p->skew == s);
  ASSERT_ALWAYS (mpz_cmp_ui (p->pols[0]->coeff[0], 17) == 0);
  ASSERT_ALWAYS (mpz_cmp_ui (p->pols[0]->coeff[1], 42) == 0);
  ASSERT_ALWAYS (mpz_cmp_ui (p->pols[1]->coeff[0], 59) == 0);
  ASSERT_ALWAYS (mpz_cmp_ui (p->m, 3613) == 0);

  cado_poly_clear (p);
  cado_poly_clear (q);
}

int
main (int argc, const char *argv[])
{
  tests_common_cmdline(&argc, &argv, 0);
  test_cado_poly_set ();
  tests_common_clear();
  exit (EXIT_SUCCESS);
}
