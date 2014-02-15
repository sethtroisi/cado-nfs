#include "cado.h"
#include "discriminant.h"
#include "macros.h"
#include "test_iter.h"
#include "tests_common.h"

#define N 10

void
test_discriminant (unsigned long iter)
{
  mpz_t f[N], D;
  int i, d;

  for (i = 0; i < N; i++)
    mpz_init (f[i]);
  mpz_init (D);

  mpz_set_ui (f[3], 1);
  mpz_set_ui (f[2], 0);
  mpz_set_ui (f[1], 3);
  mpz_set_ui (f[0], 4);
  discriminant (D, f, 3);
  ASSERT_ALWAYS(mpz_cmp_ui (D, 540) == 0);

  mpz_set_ui (f[3], 2);
  mpz_set_ui (f[2], 0);
  mpz_set_ui (f[1], 0);
  mpz_set_ui (f[0], 0);
  discriminant (D, f, 3);
  ASSERT_ALWAYS(mpz_cmp_ui (D, 0) == 0);

  while (iter--)
    {
      d = 1 + (lrand48 () % (N-1));
      for (i = 0; i <= d; i++)
        mpz_set_si (f[i], (lrand48 () % 5) - 2);
      while (mpz_cmp_ui (f[d], 0) == 0)
        mpz_set_si (f[d], (lrand48 () % 5) - 2);
      discriminant (D, f, d);
    }

  for (i = 0; i < N; i++)
    mpz_clear (f[i]);
  mpz_clear (D);
}

int
main (int argc, const char *argv[])
{
  unsigned long iter = 20000;
  tests_common_cmdline (&argc, &argv, PARSE_SEED | PARSE_ITER);
  tests_common_get_iter (&iter);
  test_discriminant (iter);
  tests_common_clear();
  exit (EXIT_SUCCESS);
}

