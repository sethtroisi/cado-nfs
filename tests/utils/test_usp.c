#include "cado.h"
#include <stdio.h>
#include "macros.h"
#include "tests_common.h"
#include "usp.h"

#define MAX_DEGREE 100

void
test_usp ()
{
  mpz_t p[MAX_DEGREE], u;
  int n, i, d;

  for (i = 0; i < MAX_DEGREE; i++)
    mpz_init (p[i]);

  /* polynomial x */
  mpz_set_ui (p[0], 0);
  mpz_set_ui (p[1], 1);
  n = numberOfRealRoots (p, 1, 2, 1, NULL);
  assert (n == 1);

  /* polynomial x+1 */
  mpz_set_ui (p[0], 1);
  mpz_set_ui (p[1], 1);
  n = numberOfRealRoots (p, 1, 2, 0, NULL);
  assert (n == 1);

  /* polynomial x-1 */
  mpz_set_si (p[0], -1);
  mpz_set_ui (p[1], 1);
  n = numberOfRealRoots (p, 1, 2, 0, NULL);
  assert (n == 1);

  /* polynomial x^2+1 */
  mpz_set_si (p[0], 1);
  mpz_set_ui (p[1], 0);
  mpz_set_ui (p[2], 1);
  n = numberOfRealRoots (p, 2, 2, 0, NULL);
  assert (n == 0);

  /* polynomial x^2-1 */
  mpz_set_si (p[0], -1);
  mpz_set_ui (p[1], 0);
  mpz_set_ui (p[2], 1);
  n = numberOfRealRoots (p, 2, 2, 0, NULL);
  assert (n == 2);

  mpz_init (u);
  mpz_set_ui (u, 1);
  mpz_mul_2exp (u, u, 127);
  for (d = 1; d < MAX_DEGREE; d++)
    {
      for (i = 0; i <= d; i++)
        {
          mpz_urandomb (p[i], state, 128);
          mpz_sub (p[i], p[i], u);
        }
      if (mpz_cmp_ui (p[d], 0) == 0)
        mpz_set_ui (p[d], 1);
      n = numberOfRealRoots (p, d, 1000, 0, NULL);
      assert (0 <= n && n <= d);
    }

  for (i = 0; i < MAX_DEGREE; i++)
    mpz_clear (p[i]);
  mpz_clear (u);
}

int
main (int argc, const char *argv[])
{
  tests_common_cmdline(&argc, &argv, PARSE_SEED | PARSE_ITER);
  test_usp ();
  tests_common_clear();
  exit (EXIT_SUCCESS);
}

