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
      a = (i == 0) ? 0 : random_int64 ();
      do b = random_int64 (); while (b == 0);
      g = gcd_int64 (a, b);
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
      a = (i == 0) ? 0 : random_uint64 ();
      do b = random_uint64 (); while (b == 0);
      g = gcd_uint64 (a, b);
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
  exit (EXIT_SUCCESS);
}
