#include "cado.h"
#include <stdint.h>
#include <time.h>
#include <gmp.h>
#include "gcd.h"
#include "macros.h"
#include "gmp_aux.h"
#include "test_iter.h"
#include "tests_common.h"

static void
cmp_mpz_gcd_i64(const int64_t a, const int64_t b, const int64_t g)
{
  mpz_t ma, mb, mg;
  mpz_init (ma);
  mpz_init (mb);
  mpz_init (mg);

  mpz_set_int64 (ma, a);
  mpz_set_int64 (mb, b);
  mpz_gcd (mg, ma, mb);
  
  if (mpz_get_int64 (mg) != g) {
    fprintf (stderr, "GCD(%lld, %lld) = %lld is incorrect, GMP has %lld",
             (long long) a, (long long) b, (long long) g,
             (long long) mpz_get_int64 (mg));
    abort();
  }

  mpz_clear (ma);
  mpz_clear (mb);
  mpz_clear (mg);
}

static void
cmp_mpz_gcd_ui64(const uint64_t a, const uint64_t b, const uint64_t g)
{
  mpz_t ma, mb, mg;
  mpz_init (ma);
  mpz_init (mb);
  mpz_init (mg);

  mpz_set_uint64 (ma, a);
  mpz_set_uint64 (mb, b);
  mpz_gcd (mg, ma, mb);
  
  if (mpz_get_uint64 (mg) != g) {
    fprintf (stderr, "GCD(%llu, %llu) = %llu is incorrect, GMP has %llu",
             (unsigned long long) a, (unsigned long long int) b,
             (unsigned long long int) g,
             (unsigned long long) mpz_get_uint64 (mg));
    abort();
  }

  mpz_clear (ma);
  mpz_clear (mb);
  mpz_clear (mg);
}

void
test_gcd_int64 (void)
{
  int64_t a, b, g;
  int i;
  
  for (i = 0; i < 200000; i++)
    {
      a = (i == 0 || i == 1) ? 0 : random_int64 ();
      b = (i == 0 || i == 2) ? 0 : random_int64 ();
      g = gcd_int64 (a, b);
      cmp_mpz_gcd_i64(a, b, g);
    }
}

void
test_gcd_uint64 (void)
{
  uint64_t a, b, g;
  int i;
  
  for (i = 0; i < 200000; i++)
    {
      a = (i == 0 || i == 1) ? 0 : random_uint64 ();
      b = (i == 0 || i == 2) ? 0 : random_uint64 ();
      g = gcd_uint64 (a, b);
      cmp_mpz_gcd_ui64(a, b, g);
    }
}

void
test_gcd_ul (void)
{
  unsigned long a, b, g;
  int i;

  assert (sizeof(unsigned long) <= sizeof(uint64_t));
  for (i = 0; i < 200000; i++)
    {
      a = (unsigned long) (i == 0 || i == 1) ? 0 : random_uint64 ();
      b = (unsigned long) (i == 0 || i == 2) ? 0 : random_uint64 ();
      g = gcd_ul (a, b);
      ASSERT(sizeof(unsigned long) <= sizeof(uint64_t));
      cmp_mpz_gcd_ui64(a, b, g);
    }
}

void
test_bin_gcd_int64_safe_ab (const int64_t a, const int64_t b)
{
  const int64_t g = bin_gcd_int64_safe (a, b);
  cmp_mpz_gcd_i64 (a, b, g);
}

void
test_bin_gcd_int64_safe (void)
{
  test_iter_t iter_a, iter_b;
  int i;

  test_iter_init(iter_a, 100, test_iter_int64_next);
  while (!test_iter_isdone(iter_a)) {
    int64_t a;
    test_iter_next(&a, iter_a);
    test_iter_init(iter_b, 100, test_iter_int64_next);
    while (!test_iter_isdone(iter_b)) {
      int64_t b;
      test_iter_next(&b, iter_b);
      test_bin_gcd_int64_safe_ab(a, b);
    }
  }

  for (i = 0; i < 200000; i++)
    {
      int64_t a = (i == 0 || i == 1) ? 0 : random_int64 ();
      int64_t b = (i == 0 || i == 2) ? 0 : random_int64 ();
      test_bin_gcd_int64_safe_ab(a,b);
    }
}

int
main (int argc, const char *argv[])
{
  tests_common_cmdline(&argc, &argv, PARSE_SEED);
  test_gcd_int64 ();
  test_gcd_uint64 ();
  test_gcd_ul ();
  test_bin_gcd_int64_safe ();
  tests_common_clear();
  exit (EXIT_SUCCESS);
}
