#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <sys/time.h>

#include "sieve/ecm/prac_bc.h"
#include "sieve/ecm/addchain_bc.h"
#include "utils.h"
#include "test_iter.h"
#include "tests_common.h"

void
mpz_prod_primes_below_B1 (mpz_t E, unsigned int B1)
{
  mpz_set_ui (E, 1);
  prime_info pi;
  prime_info_init (pi);
  unsigned long p = getprime_mt (pi);
  ASSERT (p == 3);

  for ( ; p <= B1; p = (unsigned int) getprime_mt (pi))
  {
    for (unsigned long q = 1; q <= B1 / p; q *= p)
      mpz_mul_ui (E, E, p);
  }
  prime_info_clear (pi);
}

unsigned int
test_prac (unsigned int B1min, unsigned int B1max, int compress,
           const prac_cost_t *opcost, int verbose)
{
  unsigned int n = 0;
  /* Test all values of B1 in [B1min..B1max] */
  for (unsigned int B1 = B1min; B1 <= B1max; B1++)
  {
    if (verbose)
      printf ("##### %s with B1 = %u\n", __func__, B1);
    char *bc;
    mpz_t E;
    unsigned int len = prac_bytecode (&bc, B1, 0, 0, opcost, compress, verbose);

    mpz_init (E);
    mpz_prod_primes_below_B1 (E, B1);

    if (prac_bytecode_check (bc, len, E, verbose) != 0)
    {
      printf ("##### Test with B1 = %u failed with PRAC\n", B1);
      n++;
    }
    free (bc);
    mpz_clear (E);
  }
  return n;
}

unsigned int
test_addchain (unsigned int B1min, unsigned int B1max, int compress,
               const addchain_cost_t *opcost, int verbose)
{
  unsigned int n = 0;
  /* Test all values of B1 in [B1min..B1max] */
  for (unsigned int B1 = B1min; B1 <= B1max; B1++)
  {
    if (verbose)
      printf ("##### %s with B1 = %u\n", __func__, B1);
    char *bc;
    mpz_t E;
    unsigned int len = addchain_bytecode (&bc, B1, 0, 0, opcost, compress,
                                                                       verbose);

    mpz_init (E);
    mpz_prod_primes_below_B1 (E, B1);

    if (addchain_bytecode_check (bc, len, E, verbose) != 0)
    {
      printf ("##### Test with B1 = %u failed with addchain\n", B1);
      n++;
    }
    free (bc);
    mpz_clear (E);
  }
  return n;
}

int main (int argc, const char **argv)
{
  unsigned int nerrors = 0;

  tests_common_cmdline (&argc, &argv, PARSE_SEED | PARSE_VERBOSE);

  int verbose = tests_common_get_verbose ();

  { /* PRAC test */
    /* Random opcost (between 0 and 16) [ 0x1p60 == double 2^60 ] */
    prac_cost_t cost;
    cost.dbl = (double) random_uint64 () / 0x1p60;
    cost.dadd = (double) random_uint64 () / 0x1p60;
    printf ("PRAC cost: dbl = %f; dadd = %f\n", cost.dbl, cost.dadd);

    /* compress prac chains */
    nerrors += test_prac (100, 1000, 1, &cost, verbose);
    /* uncompress prac chains */
    nerrors += test_prac (100, 1000, 0, &cost, verbose);
  }

  { /* addchain test */
    /* Random opcost (between 0 and 16) [ 0x1p60 == double 2^60 ] */
    addchain_cost_t cost;
    cost.dbl = (double) random_uint64 () / 0x1p60;
    cost.add = (double) random_uint64 () / 0x1p60;
    cost.dbladd = (double) random_uint64 () / 0x1p60;
    cost.dbl_precomp = (double) random_uint64 () / 0x1p60;
    cost.add_precomp = (double) random_uint64 () / 0x1p60;
    printf ("Addchain cost: dbl = %f; add = %f; dbladd = %f; dbl_precomp = %f; "
            "add_precomp = %f\n", cost.dbl, cost.add, cost.dbladd,
            cost.dbl_precomp, cost.add_precomp);
    

    /* compress additions chains */
    nerrors += test_addchain (100, 1000, 1, &cost, verbose);
    /* uncompress additions chains */
    nerrors += test_addchain (100, 1000, 0, &cost, verbose);
  }

  tests_common_clear ();
  return (nerrors == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
