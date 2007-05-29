/*****************************************************************
 *       Some basic number theory functions for inlining         *
 *****************************************************************/

#include "config.h"
#include "cado.h"

#if 0
/* Buggy ? */
static inline unsigned long
gcd (unsigned long a, unsigned long b)
{
  unsigned long t;
  uint32_t a32, b32, t32;

  ASSERT (b > 0);
  
  if (a >= b)
    a %= b;

  while (a > UINT32_MAX)
    {
      if (b - a < a)
	t = b - a;
      else
	t = b % a;
      b = a;
      a = t;
    }

  if (a == 0)
    return b;

  a32 = a;
  b32 = b;
  t32 = t;

  while (a32 > 0)
    {
      t32 = b32 % a32; /* The 32 bit DIV is much faster */
      b32 = a32;
      a32 = t32;
    }

  return b32;
}
#else
static inline unsigned long
gcd (unsigned long a, unsigned long b)
{
  unsigned long t;

  ASSERT (b > 0);
  
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
#endif

/* Returns 0 if n is prime, otherwise the smallest prime factor of n */
static inline unsigned long
iscomposite (const unsigned long n)
{
  unsigned long i, i2;

  if (n % 2 == 0)
    return (n == 2) ? 0 : 2;

  /* (i + 2)^2 = i^2 + 4*i + 4 */
  for (i = 3, i2 = 9; i2 <= n; i2 += (i+1) * 4, i += 2)
    if (n % i == 0)
	return i;

  return 0;
}

static inline uint32_t
signed_mod_longto32 (long a, uint32_t p)
{
  uint32_t amodp;
  if (a < 0)
    {
      amodp = ((unsigned long)(-a)) % p;
      if (amodp > 0)
	amodp = p - amodp;
    }
  else
    amodp = ((unsigned long) a) % p;

  return amodp;
}
