#include "mod_ul.h"

void
modul_div3 (residueul_t r, residueul_t a, modulusul_t m)
{
  const unsigned long a3 = a[0] % 3UL;
  const unsigned long m3 = m[0] % 3UL;
#ifndef NDEBUG
  residueul_t t;
#endif

  if (a3 == 0UL)
    r[0] = a[0] / 3UL;
  else 
    {
      ASSERT(m3 != 0UL);
      if (a3 + m3 == 3UL) /* Hence a3 == 1, m3 == 2 or a3 == 3, m3 == 1 */
	r[0] = a[0] / 3UL + m[0] / 3UL + 1UL;
      else /* a3 == 1, m3 == 1 or a3 == 2, m3 == 2 */
	r[0] = m[0] - (m[0] - a[0]) / 3UL;
    }

#ifndef NDEBUG
  modul_add (t, r, r, m);
  modul_add (t, t, r, m);
  ASSERT (t[0] == a[0]);
#endif
}


/* Compute 1/t (mod 2^wordsize) */
unsigned long
modul_invmodlong (modulusul_t m)
{
  unsigned long r;

  ASSERT (m[0] % 2UL != 0UL);
  
  /* The square of an odd integer is always 1 (mod 8). So by
     initing r = m, the low three bits in the approximate inverse
     are correct. 
     When r = 1/m (mod 16), the 4th bit of r happens to be the
     XOR of bits 2, 3 and 4 of m. This gives us an approximate 
     inverse with the 4 lowest bits correct, so 3 (for 32 bit) or
     4 (for 64 bit) Newton iterations are enough. */
  r = m[0] ^ ((m[0] & 4UL) << 1) ^ ((m[0] & 2UL) << 2);
  r = 2UL * r - r * r * m[0]; /* Newton iteration */
  r = 2UL * r - r * r * m[0];
  r = 2UL * r - r * r * m[0];
  if (sizeof (unsigned long) > 4)
    r = 2UL * r - r * r * m[0];

  ASSERT (r * m[0] == 1UL);

  return r;
}


unsigned long
modul_gcd (residueul_t r, modulusul_t m)
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


/* Compute r = b^e, with r and b in Montgomery representation */
void
modul_powredc_ul (residueul_t r, const residueul_t b, const unsigned long e, 
                  const unsigned long invm, const modulusul_t m)
{
  unsigned long mask = ~0UL - (~0UL >> 1); /* Only MSB set */
#ifndef NDEBUG
  unsigned long e1 = e;
#endif
  if (e == 0UL)
    {
      const residueul_t one = {1UL};
      modul_tomontgomery (r, one, m);
      return;
    }

  /* Assume r = 1 here for the invariant.
     r^mask * b^e is invariant, and is the result we want */

  /* Find highest set bit in e. */
  while ((e & mask) == 0UL)
    mask >>= 1; /* r = 1, so r^(mask/2) * b^e = r^mask * b^e  */

  /* Exponentiate */

  modul_set (r, b, m);       /* (r*b)^mask * b^(e-mask) = r^mask * b^e */
#ifndef NDEBUG
  e1 -= mask;
#endif

  while (mask > 1UL)
    {
      modul_mulredc (r, r, r, invm, m);
      mask >>= 1;            /* (r^2)^(mask/2) * b^e = r^mask * b^e */
      if (e & mask)
        {
	  modul_mulredc (r, r, b, invm, m);
#ifndef NDEBUG
          e1 -= mask;
#endif
        }
    }
  ASSERT (e1 == 0UL && mask == 1UL);
  /* Now e = 0, mask = 1, and r^mask * b^0 = r^mask is the result we want */
}


/* Compute r = b^e, with r and b in Montgomery representation. 
   Here e is a multiple precision integer 
   sum_{i=0}^{e_nrwords} e[i] * (machine word base)^i */
void
modul_powredc_mp (residueul_t r, const residueul_t b, const unsigned long *e, 
                  const int e_nrwords, const unsigned long invm, 
                  const modulusul_t m)
{
  unsigned long mask = ~0UL - (~0UL >> 1); /* Only MSB set */
  int i = e_nrwords - 1;

  if (e_nrwords == 0 || e[i] == 0UL)
    {
      const residueul_t one = {1UL};
      modul_tomontgomery (r, one, m);
      return;
    }

  /* Find highest set bit in e. */
  while ((e[i] & mask) == 0UL)
    mask >>= 1; /* r = 1, so r^(mask/2) * b^e = r^mask * b^e  */

  /* Exponentiate */

  modul_set (r, b, m);       /* (r*b)^mask * b^(e-mask) = r^mask * b^e */
  mask >>= 1;

  for ( ; i >= 0; i--)
    {
      while (mask > 0UL)
        {
          modul_mulredc (r, r, r, invm, m);
          if (e[i] & mask)
            modul_mulredc (r, r, b, invm, m);
          mask >>= 1;            /* (r^2)^(mask/2) * b^e = r^mask * b^e */
        }
      mask = ~0UL - (~0UL >> 1);
    }
}


void
modul_2powredc_mp (residueul_t r, const unsigned long *e, 
                   const int e_nrwords, const unsigned long invm, 
                   const modulusul_t m)
{
  unsigned long mask = ~0UL - (~0UL >> 1); /* Only MSB set */
  int i = e_nrwords - 1;

  if (e_nrwords == 0 || e[i] == 0UL)
    {
      const residueul_t one = {1UL};
      modul_tomontgomery (r, one, m);
      return;
    }

  /* Find highest set bit in e. */
  while ((e[i] & mask) == 0UL)
    mask >>= 1; /* r = 1, so r^(mask/2) * b^e = r^mask * b^e  */

  /* Set r to 2 in Montgomery representation. We use the fact that
     (m*invm + 1)/w == 1/w (mod m) */
  {
    unsigned long dummy;
    mul_ul_ul_2ul (r , &dummy, m[0], invm);
    r[0]++;
  }

  /* Exponentiate */
  {
    const residueul_t two = {2UL};
    modul_tomontgomery (r, two, m);  
  }
  mask >>= 1;

  for ( ; i >= 0; i--)
    {
      while (mask > 0UL)
        {
          modul_mulredc (r, r, r, invm, m);
          if (e[i] & mask)
            modul_add (r, r, r, m);
          mask >>= 1;            /* (r^2)^(mask/2) * b^e = r^mask * b^e */
        }
      mask = ~0UL - (~0UL >> 1);
    }
}


/* Compute r = V_e(b), with r and b in Montgomery representation. 
   Here e is a multiple precision integer 
   sum_{i=0}^{e_nrwords} e[i] * (machine word base)^i */

void
modul_Vredc_mp (residueul_t r, const residueul_t b, const unsigned long *e, 
                const int e_nrwords, const unsigned long invm, 
                const modulusul_t m)
{
  unsigned long mask = ~0UL - (~0UL >> 1); /* Only MSB set */
  int i = e_nrwords - 1;
  residueul_t r1, two;

  modul_init_noset0 (two, m);
  modul_set_ul_reduced (two, 2UL, m);
  modul_tomontgomery (two, two, m);

  if (e_nrwords == 0 || e[i] == 0UL)
    {
      modul_set (r, two, m);
      return;
    }

  /* Find highest set bit in e. */
  while ((e[i] & mask) == 0UL)
    mask >>= 1; /* r = 1, so r^(mask/2) * b^e = r^mask * b^e  */

  /* Exponentiate */

  modul_init (r1, m);
  modul_set (r, b, m);         /* r = b = V_1 (b) */
  modul_mulredc (r1, b, b, invm, m);
  modul_sub (r1, r1, two, m);  /* r1 = b^2 - 2 = V_2 (b) */
  mask >>= 1;

  /* Here r = V_j (b) and r1 = V_{j+1} (b) for j = 1 */

  for ( ; i >= 0; i--)
    {
      while (mask > 0UL)
        {
          if (e[i] & mask)
	    {
	      /* j -> 2*j+1. Compute V_{2j+1} and V_{2j+2} */
	      modul_mulredc (r, r, r1, invm, m);
	      modul_sub (r, r, b, m); /* V_j * V_{j+1} - V_1 = V_{2j+1} */
	      modul_mulredc (r1, r1, r1, invm, m);
	      modul_sub (r1, r1, two, m); /* (V_{j+1})^2 - 2 = V_{2j+2} */
	    }
	  else
	    {
	      /* j -> 2*j. Compute V_{2j} and V_{2j+1} */
	      modul_mulredc (r1, r1, r, invm, m);
	      modul_sub (r1, r1, b, m); /* V_j * V_{j+1} - V_1 = V_{2j+1}*/
	      modul_mulredc (r, r, r, invm, m);
	      modul_sub (r, r, two, m);
	    }
          mask >>= 1;
        }
      mask = ~0UL - (~0UL >> 1);
    }
}


/* Returns 1 if m is a strong probable prime wrt base b, 0 otherwise.
   b must be < m.
   b is NOT passed in Montgomery form, even though the exponentiation
   in this function uses REDC - this function converts b by itself.
   invn must be passed so that n * invn  == -1 (mod 2^wordsize).  */
int
modul_sprp (const residueul_t b, const unsigned long invn, 
	    const modulusul_t m)
{
  residueul_t b1, r1, t;
  int i = 0, po2 = 1;
  unsigned long mm1 = (m[0] - 1UL) >> 1;

  ASSERT (b[0] < m[0]);
  ASSERT (invn * m[0] == (unsigned long) -1);

  if (m[0] <= 3UL)
    return (m[0] >= 2UL);

  if (m[0] % 2UL == 0UL)
    return 0;

  /* Set mm1 to the odd part of m-1 */
  while (mm1 % 2UL == 0UL)
    {
      po2++;
      mm1 >>= 1;
    }
  /* Hence, m-1 = mm1 * 2^po2 */

  modul_init_noset0 (b1, m);
  modul_init_noset0 (r1, m);
  modul_init_noset0 (t, m);
  mod_tomontgomery (b1, b, m);

  /* Exponentiate */
  modul_powredc_ul (r1, b1, mm1, invn, m);
  modul_clear (b1, m);

  modul_frommontgomery (t, r1, invn, m);
  /* Now t == b^mm1 (mod m) */
#if defined(PARI)
  printf ("(Mod(%lu,%lu) ^ %lu) == %lu /* PARI */\n", 
	  b, n, mm1, modul_get_ul (t, m));
#endif
  
  /* If m is prime, then b^mm1 might be == 1 or == -1 (mod m) here */
  if (modul_get_ul (t, m) == 1UL || modul_get_ul (t, m) == m[0] - 1UL)
    i = 1;
  else
    {
      /* If m is a prime, then one of b^(2*mm1), b^(2^2*mm1), ..., 
	 b^(2^(po2 - 1)*mm1)  must be == -1 (mod m) */
      for ( ; po2 > 1; po2--)
	{
	  modul_mulredc (r1, r1, r1, invn, m);
	  modul_frommontgomery (t, r1, invn, m);
	  if (modul_get_ul (t, m) == m[0] - 1UL)
	    {
	      i = 1;
	      break;
	    }
	}
    }

  modul_clear (r1, m);
  modul_clear (t, m);
  return i;
}


/* Put 1/s (mod t) in r and return 1 if s is invertible, 
   or return 0 if s is not invertible */

int
modul_inv (residueul_t r, residueul_t s, modulusul_t t)
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

#ifndef NDEBUG
  modul_mul (&u2, (unsigned long *)&u1, s, t);
  ASSERT(u2 == 1UL);
#endif

  r[0] = (unsigned long) u1;
  return 1;
}

int
modul_jacobi (residueul_t a_par, modulusul_t m_par)
{
  unsigned long a, m, s;
  int t = 1;

  ASSERT (a_par[0] < m_par[0]);
  a = a_par[0];
  m = m_par[0];
  
  while (a != 0UL)
  {
    while (a % 2UL == 0UL)
    {
      a /= 2UL;
      if (m % 8UL == 3UL || m % 8UL == 5UL)
        t = -t;
    }
    s = a; a = m; m = s; /* swap */
    if (a % 4UL == 3UL && m % 4UL == 3UL)
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
