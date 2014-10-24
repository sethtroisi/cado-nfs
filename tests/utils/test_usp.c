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
  double root;

  for (i = 0; i < MAX_DEGREE; i++)
    mpz_init (p[i]);

  mpz_init (R[0].a);
  mpz_init (R[0].b);

  /* polynomial x */
  mpz_set_ui (p[0], 0);
  mpz_set_ui (p[1], 1);
  n = numberOfRealRoots (p, 1, 2, 1, R);
  ASSERT_ALWAYS (n == 1);
  /* check [a/2^ka, b/2^kb] contains root 0 */
  ASSERT_ALWAYS (mpz_cmp_ui (R[0].a, 0) <= 0);
  ASSERT_ALWAYS (0 <= mpz_cmp_ui (R[0].b, 0));

  /* polynomial 2*x */
  mpz_set_ui (p[0], 0);
  mpz_set_ui (p[1], 2);
  n = numberOfRealRoots (p, 1, 2, 1, R);
  ASSERT_ALWAYS (n == 1);
  root = rootRefine (R, p, 1, 1e-9);
  ASSERT_ALWAYS(-0.000001 <= root && root <= 0.000001);

  /* polynomial x+1 */
  mpz_set_ui (p[0], 1);
  mpz_set_ui (p[1], 1);
  n = numberOfRealRoots (p, 1, 0, 0, NULL);
  ASSERT_ALWAYS (n == 1);

  /* polynomial x-1 */
  mpz_set_si (p[0], -1);
  mpz_set_ui (p[1], 1);
  n = numberOfRealRoots (p, 1, 0, 0, NULL);
  ASSERT_ALWAYS (n == 1);

  /* polynomial x^2+1 */
  mpz_set_si (p[0], 1);
  mpz_set_ui (p[1], 0);
  mpz_set_ui (p[2], 1);
  n = numberOfRealRoots (p, 2, 0, 0, NULL);
  ASSERT_ALWAYS (n == 0);

  /* polynomial x^2+2 */
  mpz_set_si (p[0], 2);
  mpz_set_ui (p[1], 0);
  mpz_set_ui (p[2], 1);
  n = numberOfRealRoots (p, 2, 0, 0, NULL);
  ASSERT_ALWAYS (n == 0);

  /* polynomial x^2-1 */
  mpz_set_si (p[0], -1);
  mpz_set_ui (p[1], 0);
  mpz_set_ui (p[2], 1);
  n = numberOfRealRoots (p, 2, 0, 0, NULL);
  ASSERT_ALWAYS (n == 2);

  /* polynomial x^2-3*x+2 */
  mpz_set_si (p[0], 2);
  mpz_set_si (p[1], -3);
  mpz_set_ui (p[2], 1);
  n = numberOfRealRoots (p, 2, 0, 0, NULL);
  ASSERT_ALWAYS (n == 2);

  /* polynomial x^2-2*x+2 */
  mpz_set_si (p[0], 2);
  mpz_set_si (p[1], -2);
  mpz_set_ui (p[2], 1);
  n = numberOfRealRoots (p, 2, 0, 0, NULL);
  ASSERT_ALWAYS (n == 0);

  /* polynomial x^2+1000*x+2 */
  mpz_set_si (p[0], 2);
  mpz_set_si (p[1], 1000);
  mpz_set_ui (p[2], 1);
  n = numberOfRealRoots (p, 2, 0, 0, NULL);
  ASSERT_ALWAYS (n == 2);

  /* polynomial (x-1)*(x-2)*(x-3) */
  mpz_set_si (p[0], -6);
  mpz_set_si (p[1], 11);
  mpz_set_si (p[2], -6);
  mpz_set_si (p[3], 1);
  n = numberOfRealRoots (p, 3, 0, 0, NULL);
  ASSERT_ALWAYS (n == 3);

  /* polynomial (x-1)*(x-2)*(x-3)*(x-4) */
  mpz_set_si (p[0], 24);
  mpz_set_si (p[1], -50);
  mpz_set_si (p[2], 35);
  mpz_set_si (p[3], -10);
  mpz_set_si (p[4], 1);
  n = numberOfRealRoots (p, 4, 0, 0, NULL);
  ASSERT_ALWAYS (n == 4);

  /* polynomial (x-1)*(x-2)*(x-3)*(x-4)*(x-5) */
  mpz_set_si (p[0], -120);
  mpz_set_si (p[1], 274);
  mpz_set_si (p[2], -225);
  mpz_set_si (p[3], 85);
  mpz_set_si (p[4], -15);
  mpz_set_si (p[5], 1);
  n = numberOfRealRoots (p, 5, 0, 0, NULL);
  ASSERT_ALWAYS (n == 5);

  /* polynomial (x-1)*(x-2)*(x-3)*(x-4)*(x-5)*(x-6) */
  mpz_set_si (p[0], 720);
  mpz_set_si (p[1], -1764);
  mpz_set_si (p[2], 1624);
  mpz_set_si (p[3], -735);
  mpz_set_si (p[4], 175);
  mpz_set_si (p[5], -21);
  mpz_set_si (p[6], 1);
  n = numberOfRealRoots (p, 6, 0, 0, NULL);
  ASSERT_ALWAYS (n == 6);

  /* polynomial 2*(x-1)*(x-2)*(x-3/2) */
  mpz_set_si (p[0], -6);
  mpz_set_si (p[1], 13);
  mpz_set_si (p[2], -9);
  mpz_set_si (p[3], 2);
  n = numberOfRealRoots (p, 3, 0, 0, NULL);
  ASSERT_ALWAYS (n == 3);

  /* polynomial 8*(x-1)*(x-2)*(x-3/2)*(x-5/4) */
  mpz_set_si (p[0], 30);
  mpz_set_si (p[1], -89);
  mpz_set_si (p[2], 97);
  mpz_set_si (p[3], -46);
  mpz_set_si (p[4], 8);
  n = numberOfRealRoots (p, 4, 0, 0, NULL);
  ASSERT_ALWAYS (n == 4);

  /* polynomial 2*(x-1)*(x-3/2)*(x-2)*(x-3)*(x-4) */
  mpz_set_si (p[0], -72);
  mpz_set_si (p[1], 198);
  mpz_set_si (p[2], -205);
  mpz_set_si (p[3], 100);
  mpz_set_si (p[4], -23);
  mpz_set_si (p[5], 2);
  n = numberOfRealRoots (p, 5, 0, 0, NULL);
  ASSERT_ALWAYS (n == 5);

  mpz_set_si (p[0], -3);
  mpz_set_si (p[1], 2);
  n = numberOfRealRoots (p, 1, 0, 0, NULL);
  ASSERT_ALWAYS (n == 1);

  /* 2*(x-17)*(x-17/5) */
  mpz_set_si (p[0], 595);
  mpz_set_si (p[1], -69);
  mpz_set_si (p[2], 2);
  n = numberOfRealRoots (p, 2, 0, 0, NULL);
  ASSERT_ALWAYS (n == 2);

  /* 2000*x^3 - 35000*x^2 + 2*x - 35 */
  mpz_set_si (p[0], -35);
  mpz_set_si (p[1], 2);
  mpz_set_si (p[2], -35000);
  mpz_set_si (p[3], 2000);
  n = numberOfRealRoots (p, 3, 0, 0, NULL);
  ASSERT_ALWAYS (n == 1);

  /* 33*x^4-124*x^3+137*x^2-32*x-11 */
  mpz_set_si (p[0], -11);
  mpz_set_si (p[1], -32);
  mpz_set_si (p[2], 137);
  mpz_set_si (p[3], -124);
  mpz_set_si (p[4], 33);
  n = numberOfRealRoots (p, 4, 0, 0, NULL);
  ASSERT_ALWAYS (n == 4);

  /* infinite loop that appeared on 32-bit computers, 2nd April 2014 */
  mpz_set_str (p[0], "202156081496031128767910903808", 10);
  mpz_set_str (p[1], "-1014774808369763543554925264896", 10);
  mpz_set_str (p[2], "305592973565950609559904059392", 10);
  mpz_set_str (p[3], "532741928198739321928042938368", 10);
  mpz_set_str (p[4], "-265275048860303928171627020288", 10);
  mpz_set_str (p[5], "152821965764546794524860481536", 10);
  mpz_set_str (p[6], "-40951570592501524318484692992", 10);
  mpz_set_ui (R[0].a, 2);
  mpz_set_ui (R[0].b, 4);
  R[0].ka = R[0].kb = 0;
  root = rootRefine (R, p, 6, 1e-9); /* root is near 3.00763029864372 */
  ASSERT_ALWAYS(3.00763 <= root && root <= 3.00764);

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
      ASSERT_ALWAYS (0 <= n && n <= d);
    }

  /* check with large coefficients */
  mpz_urandomb (p[0], state, 2048);
  mpz_urandomb (p[1], state, 2048);
  mpz_urandomb (p[2], state, 2048);
  mpz_urandomb (p[3], state, 2048);
  n = numberOfRealRoots (p, 3, 0, 0, NULL);
  ASSERT_ALWAYS (0 <= n && n <= 3);

  for (i = 0; i < MAX_DEGREE; i++)
    mpz_clear (p[i]);
  mpz_clear (u);
  mpz_clear (R[0].a);
  mpz_clear (R[0].b);
}

int
main (int argc, const char *argv[])
{
  tests_common_cmdline(&argc, &argv, PARSE_SEED | PARSE_ITER);
  test_usp ();
  tests_common_clear();
  exit (EXIT_SUCCESS);
}

