#include "mod_ul.h"

#ifndef ASSERT
 #ifdef WANT_ASSERT
  #include <assert.h>
  #define ASSERT(x) assert(x)
 #else
  #define ASSERT(x)
 #endif
#define MOD_UL_ASSERT
#endif

void
modul_div3 (residueul r, residueul a, modulusul m)
{
  const unsigned long a3 = a[0] % 3;
  const unsigned long m3 = m[0] % 3;
#ifdef WANT_ASSERT
  residueul t;
#endif

  ASSERT(m3 != 0UL);
  if (a3 == 0UL)
    r[0] = a[0] / 3UL;
  else if (a3 + m3 == 3UL) /* Hence a3 == 1, m3 == 2 or a3 == 3, m3 == 1 */
    r[0] = a[0] / 3UL + m[0] / 3UL + 1UL;
  else /* a3 == 1, m3 == 1 or a3 == 2, m3 == 2 */
    r[0] = m[0] - (m[0] - a[0]) / 3UL;

#ifdef WANT_ASSERT
  modul_add (t, r, r, m);
  modul_add (t, t, r, m);
  assert (t[0] == a[0]);
#endif
}

unsigned long
modul_gcd (residueul r, modulusul m)
{
  unsigned long a = r[0], b = m[0], t;

  ASSERT (b > 0UL);
  
  if (a >= b)
    a %= b;

  while (a > 0UL)
    {
      /* Here 0 < a < b */
      t = b % a;
      b = a;
      a = t;
    }

  return b;
}


/* Compute 1/t (mod 2^wordsize) */
/* FIXME: not optimised yet */
unsigned long
modul_invmodlong (modulusul m)
{
  unsigned long p, r; 
  int i;

  ASSERT (m[0] % 2 != 0);
  
  r = 1UL; p = m[0];
  for (i = 1; p != 1UL; i++) /* Invariant: r * t == p (mod 2^wordsize) */
    if ((1UL << i) & p)
      {
        r += 1UL << i;
        p += m[0] << i;
      }

  ASSERT (r * m[0] == 1UL);
  return r;
}

/* Compute r = b^e, with r and a in Montgomery representation */
void
modul_powredc_ul (residueul r, residueul b, unsigned long e, 
                  unsigned long invm, modulusul m)
{
  unsigned long mask = ~0UL - (~0UL >> 1); /* Only MSB set */

  if (e == 0UL)
    {
      r[0] = 1UL;
      return;
    }

  /* Assume r = 1 here for the invariant.
     r^mask * b^e is invariant, and is the result we want */

  /* Find highest set bit in e. */
  while ((e & mask) == 0UL)
    mask >>= 1; /* r = 1, so r^(mask/2) * b^e = r^mask * b^e  */

  /* Exponentiate */

  modul_set (r, b, m);       /* e -= mask; */
#ifdef WANT_ASSERT
  e -= mask;                 /* (r*b)^mask * b^(e-mask) = r^mask * b^e */
#endif

  while (mask > 1UL)
    {
      modul_mulredc (r, r, r, invm, m);
      mask >>= 1;            /* (r^2)^(mask/2) * b^e = r^mask * b^e */
      if (e & mask)
        {
	  modul_mulredc (r, r, b, invm, m);
#ifdef WANT_ASSERT
          e -= mask;         /* e -= mask; As above */
#endif
        }
    }
  ASSERT (e == 0 && mask == 1);
  /* Now e = 0, mask = 1, and r^mask * b^0 = r^mask is the result we want */
}


/* Put 1/s (mod t) in r and return 1 if s is invertible, 
   or return 0 if s is not invertible */

int
modul_inv (residueul r, residueul s, modulusul t)
{
  long u1, v1;
  unsigned long q, u2, v2;

  ASSERT (t[0] > 0UL);
  ASSERT (s[0] < t[0]);

  if (s[0] == 0UL)
    {
      r[0] = 0UL; /* Not invertible */
      return 0;
    }

  if (s[0] == 1UL)
    {
      r[0] = 1UL;
      return 1;
    }

#if 0
  u1 = 1L;
  v1 = 0L;
  u2 = s;
  v2 = t;

  q = s / t; /* == 0, due to s < t */
  u1 = u1 - (long) q * v1; /* == u1 - 0 == u1 == 1*/
  u2 = u2 - q * v2; /* == u2 - 0 == u2 == s */

  q = v2 / u2; /* == t / s <= t / 2*/
  v1 = v1 - (long) q * u1; /* == v1 - q * 1 == 0 - q == -q */
  v2 = v2 - q * u2; /* == v2 - q * s == t - q * s == t % s */
#endif

  u1 = 1L;
  u2 = s[0];
  v1 = - (t[0] / s[0]); /* No overflow, since s >= 2 */
  v2 = t[0] % s[0];

  while (v2 != 0UL)
    {
      /* unroll twice and swap u/v */
      q = u2 / v2;
      u1 = u1 - (long) q * v1;
      u2 = u2 - q * v2;

      if (u2 == 0UL)
        {
          u1 = v1;
          u2 = v2;
          break;
        }

      q = v2 / u2;
      v1 = v1 - (long) q * u1;
      v2 = v2 - q * u2;
    }

  if (u2 != 1UL)
    {
      /* printf ("s=%lu t=%lu found %lu\n", s[0], t[0], u2); */
      r[0] = 0UL; /* non-trivial gcd */
      return 0;
    }

  if (u1 < 0L)
    u1 = u1 - t[0] * (-u1 / t[0] - 1L);

#ifdef WANT_ASSERT
  modul_mul (&u2, (unsigned long *)&u1, s, t);
  ASSERT(u2 == 1UL);
#endif

  r[0] = (unsigned long) u1;
  return 1;
}

int
modul_jacobi (residueul a_par, modulusul m_par)
{
  unsigned long a, m, s;
  int t = 1;

  ASSERT (a_par[0] < m_par[0]);
  a = a_par[0];
  m = m_par[0];
  
  while (a != 0UL)
  {
    while (a % 2 == 0)
    {
      a /= 2UL;
      if (m % 8 == 3 || m % 8 == 5)
        t = -t;
    }
    s = a; a = m; m = s; /* swap */
    if (a % 4 == 3 && m % 4 == 3)
      t = -t;
    a %= m;
  }
  if (m != 1UL)
    t = 0;
  
#ifdef MODTRACE
  printf ("kronecker(%lu, %lu) == %d\n", a_par[0], m_par[0], t);
#endif
  return t;
}

#ifdef MOD_UL_ASSERT
#undef ASSERT
#undef MOD_UL_ASSERT
#endif
