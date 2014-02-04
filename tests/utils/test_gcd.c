#include "cado.h"
#include <stdint.h>
#include <time.h>
#include <gmp.h>
#include "gcd.h"
#include "macros.h"
#include "gmp_aux.h"
#include "test_iter.h"

int64_t
random_int64 (void)
{
  return (mrand48 () << 32) + mrand48 ();
}

int64_t
random_uint64 (void)
{
  return (lrand48 () << 33) + (lrand48 () << 2) + (lrand48 () & 3);
}

void
test_gcd_int64 (void)
{
  int64_t a, b, g, c, d, h;
  int i;
  
  for (i = 0; i < 200000; i++)
    {
      a = (i == 0 || i == 1) ? 0 : random_int64 ();
      b = (i == 0 || i == 2) ? 0 : random_int64 ();
      g = gcd_int64 (a, b);
      if (g == 0)
        {
          assert (a == 0 && b == 0);
          continue;
        }
      assert ((a % g) == 0);
      assert ((b % g) == 0);
      c = a / g;
      d = b / g;
      h = gcd_int64 (c, d);
      assert (h == 1);
    }
}

void
test_gcd_uint64 (void)
{
  uint64_t a, b, g, c, d, h;
  int i;
  
  for (i = 0; i < 200000; i++)
    {
      a = (i == 0 || i == 1) ? 0 : random_uint64 ();
      b = (i == 0 || i == 2) ? 0 : random_uint64 ();
      g = gcd_uint64 (a, b);
      if (g == 0)
        {
          assert (a == 0 && b == 0);
          continue;
        }
      assert ((a % g) == 0);
      assert ((b % g) == 0);
      c = a / g;
      d = b / g;
      h = gcd_uint64 (c, d);
      assert (h == 1);
    }
}

void
test_gcd_ul (void)
{
  unsigned long a, b, g, c, d, h;
  int i;

  assert (sizeof(unsigned long) <= sizeof(uint64_t));
  for (i = 0; i < 200000; i++)
    {
      a = (unsigned long) (i == 0 || i == 1) ? 0 : random_uint64 ();
      b = (unsigned long) (i == 0 || i == 2) ? 0 : random_uint64 ();
      g = gcd_ul (a, b);
      if (g == 0)
        {
          assert (a == 0 && b == 0);
          continue;
        }
      assert ((a % g) == 0);
      assert ((b % g) == 0);
      c = a / g;
      d = b / g;
      h = gcd_ul (c, d);
      assert (h == 1);
    }
}

void

cmp_mpz_gcd_i64(int64_t a, int64_t b, int64_t g)
{
  mpz_t ma, mb, mg;
  mpz_init (ma);
  mpz_init (mb);
  mpz_init (mg);

  mpz_set_int64 (ma, a);
  mpz_set_int64 (mb, b);
  mpz_gcd (mg, ma, mb);
  mpz_set_int64 (ma, a);
  
  if (mpz_get_int64 (mg) != g) {
    fprintf (stderr, "GCD(%lld, %lld) = %lld is incorrect, GMP has %lld",
             (long long int) a, (long long int) b, (long long int) g,
             (long long int) mpz_get_int64 (mg));
    abort();
  }

  mpz_clear (ma);
  mpz_clear (mb);
  mpz_clear (mg);
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
  int64_t a, b;
  test_iter_t iter_a, iter_b;
  int i;

  test_iter_init(iter_a, 100);
  while (!test_iter_isdone(iter_a)) {
    int64_t a = test_iter_int64_next(iter_a);
    test_iter_init(iter_b, 100);
    while (!test_iter_isdone(iter_b)) {
      int64_t b = test_iter_int64_next(iter_b);
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
main (int argc, char *argv[])
{
  long int seed;

  seed = (argc >= 2) ? atoi (argv[1]) : time (NULL);
  fprintf (stderr, "Using random seed=%ld\n", seed);
  srand48 (seed);
  test_gcd_int64 ();
  test_gcd_uint64 ();
  test_gcd_ul ();
  test_bin_gcd_int64_safe ();
  exit (EXIT_SUCCESS);
}
