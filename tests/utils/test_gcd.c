#include "cado.h"
#include <stdint.h>
#include <time.h>
#include "gcd.h"
#include "macros.h"

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
test_bin_gcd_int64_safe (void)
{
  int64_t a, b, g, c, d, h;
  int i;

  for (i = 0; i < 200000; i++)
    {
      a = (i == 0 || i == 1) ? 0 : random_int64 ();
      b = (i == 0 || i == 2) ? 0 : random_int64 ();
      /* the current bin_gcd_int64_safe code might overflow when a, b are
         larger than 2^62 in absolute value */
      a >>= 1;
      b >>= 1;
      g = bin_gcd_int64_safe (a, b);
      if (g == 0)
        {
          assert (a == 0 && b == 0);
          continue;
        }
      assert ((a % g) == 0);
      if (b % g)
        printf ("a=%ld b=%ld g=%ld\n", a, b, g);
      assert ((b % g) == 0);
      c = a / g;
      d = b / g;
      h = gcd_int64 (c, d);
      assert (h == 1);
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
