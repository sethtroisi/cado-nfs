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
  root_struct R[1];

  for (i = 0; i < MAX_DEGREE; i++)
    mpz_init (p[i]);

  mpz_init (R[0].a);
  mpz_init (R[0].b);

  /* polynomial x */
  mpz_set_ui (p[0], 0);
  mpz_set_ui (p[1], 1);
  n = numberOfRealRoots (p, 1, 2, 1, R);
  assert (n == 1);
  /* check [a/2^ka, b/2^kb] contains root 0 */
  assert (mpz_cmp_ui (R[0].a, 0) <= 0);
  assert (0 <= mpz_cmp_ui (R[0].b, 0));

  mpz_clear (R[0].a);
  mpz_clear (R[0].b);

  /* polynomial x+1 */
  mpz_set_ui (p[0], 1);
  mpz_set_ui (p[1], 1);
  n = numberOfRealRoots (p, 1, 0, 0, NULL);
  assert (n == 1);

  /* polynomial x-1 */
  mpz_set_si (p[0], -1);
  mpz_set_ui (p[1], 1);
  n = numberOfRealRoots (p, 1, 0, 0, NULL);
  assert (n == 1);

  /* polynomial x^2+1 */
  mpz_set_si (p[0], 1);
  mpz_set_ui (p[1], 0);
  mpz_set_ui (p[2], 1);
  n = numberOfRealRoots (p, 2, 0, 0, NULL);
  assert (n == 0);

  /* polynomial x^2-1 */
  mpz_set_si (p[0], -1);
  mpz_set_ui (p[1], 0);
  mpz_set_ui (p[2], 1);
  n = numberOfRealRoots (p, 2, 0, 0, NULL);
  assert (n == 2);

  /* polynomial x^2-2*x+1 */
  mpz_set_si (p[0], 2);
  mpz_set_si (p[1], -3);
  mpz_set_ui (p[2], 1);
  n = numberOfRealRoots (p, 2, 0, 0, NULL);
  assert (n == 2);
  exit (1);

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
      n = numberOfRealRoots (p, d, 0, 0, NULL);
      assert (0 <= n && n <= d);
    }

  /* check with large coefficients */
  mpz_urandomb (p[0], state, 2048);
  mpz_urandomb (p[1], state, 2048);
  mpz_urandomb (p[2], state, 2048);
  mpz_urandomb (p[3], state, 2048);
  n = numberOfRealRoots (p, 3, 0, 0, NULL);
  assert (0 <= n && n <= 3);

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

