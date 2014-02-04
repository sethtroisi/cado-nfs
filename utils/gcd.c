#include "cado.h"
#include <stdint.h>
#include "gcd.h"
#include "macros.h"
#include "misc.h" /* for ctzl */

int64_t
gcd_int64 (int64_t a, int64_t b)
{
  int64_t t;

  if (b == 0)
    return a;

  if (a < 0)
    a = -a;
  if (b < 0)
    b = -b;
  
  if (a >= b)
    a %= b;

  while (a > 0)
    {
      /* Here 0 < a < b */
      t = b % a; /* 0 <= t < a */
      b = a;
      a = t; /* 0 <= a < b */
    }

  return b;
}

uint64_t
gcd_uint64 (uint64_t a, uint64_t b)
{
  uint64_t t;

  if (b == 0)
    return a;

  if (a >= b)
    a %= b;

  while (a > 0)
    {
      /* Here 0 < a < b */
      t = b % a; /* 0 <= t < a */
      b = a;
      a = t; /* 0 <= a < b */
    }

  return b;
}

unsigned long
gcd_ul (unsigned long a, unsigned long b)
{
  unsigned long t;

  if (b == 0)
    return a;

  if (a >= b)
    a %= b;

  while (a > 0)
    {
      /* Here 0 < a < b */
      t = b % a;
      b = a;
      a = t;
    }

  return b;
}

/* Binary gcd; any input allowed. */
int64_t
bin_gcd_int64_safe (int64_t a, int64_t b)
{
  int s, t;

  if (a == 0)
    return b;

  if (b == 0)
    return a;

  s = ctzl (a);
  t = ctzl (b);
  a >>= s;
  b >>= t;
  if (t < s)
    s = t;
  
  while (b != 0)
    {
      /* if a is odd, reduce a wrt b, i.e., cancel the two low bits of a,
         so that a + q*b = 0 (mod 4) */
      b >>= ctzl (b);
      a = ((a ^ b) & 2) ? a + b : a - b;
      /* if a was even, then since b is now odd, the new a is odd */
      if (a == 0)
        return (b > 0) ? (b << s) : ((-b) << s);
      a >>= ctzl (a);
      /* from here on, a and b are odd (or zero) */
      ASSERT(a & 1);
      /* reduce b wrt a */
      b = ((b ^ a) & 2) ? b + a : b - a;
    }
  return (a > 0) ? (a << s) : ((-a) << s);
}
