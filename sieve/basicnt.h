/*****************************************************************
 *       Some basic number theory functions for inlining         *
 *****************************************************************/

#include <stdint.h>

static inline unsigned long 
eulerphi_ul (unsigned long n)
{
  unsigned long p, r = 1UL;
  
  if (n == 0UL) /* Undefined, we return 0 */
    return 0UL;

  if (n % 2UL == 0UL)
    {
      n /= 2UL;
      while (n % 2UL == 0UL)
        {
          n /= 2UL;
          r *= 2UL;
        }
    }

  for (p = 3UL; p*p <= n; p += 2UL)
    {
      if (n % p == 0UL)
        {
          n /= p;
          r *= p - 1UL;
          while (n % p == 0UL)
            {
              n /= p;
              r *= p;
            }
        }
    }
  /* Now n is either 1 or a prime */
  if (n > 1UL)
    r *= n - 1UL;
  
  return r;
}

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

