#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include "tests_common.h"
#include "sieve/ecm/mpqs_doit.h"

int
main (int argc, const char *argv[])
{
  mpz_t N, f;
  unsigned long bits, iter = 100, i, found = 0;

  tests_common_cmdline (&argc, &argv, PARSE_SEED | PARSE_ITER);
  tests_common_get_iter (&iter);

  if (argc > 1)
    bits = atoi (argv[1]);
  else
    /* generate a random size between 64 and 128 bits */
    bits = 64 + (rand () % 65);
  printf ("bits=%lu iter=%lu\n", bits, iter);

  mpz_init (N);
  mpz_init (f);
  for (i = 0; i < iter; i++)
    {
      do
        {
          mpz_urandomb (N, state, bits / 2);
          mpz_nextprime (N, N);
          mpz_urandomb (f, state, (bits + 1) / 2);
          mpz_nextprime (f, f);
          mpz_mul (N, N, f);
        }
      while (mpz_sizeinbase (N, 2) != bits);
      gmp_printf ("N=%Zd\n", N);
      mpqs_doit (f, N, 1);
      if (mpz_cmp_ui (f, 1) > 0 && mpz_cmp (f, N) < 0)
        found ++;
    }
  mpz_clear (N);
  mpz_clear (f);

  printf ("Found %lu factors out of %lu composites\n", found, iter);

  tests_common_clear ();

  return 0;
}
