
#if defined(WANT_ASSERT) || defined (WANT_ASSERT_EXPENSIVE)
#include <limits.h>
#endif

/* This file defines some functions that work more or less the same 
   with mod_ul.h and modredc_ul.h. I.e. mod_div3() and mod_gcd() work
   unchanged with plain and Montgomery representation (so we can work on 
   the stored residue directly, whatever its representation is); 
   mod_inv() and mod_jacobi() convert to plain "unsigned long" first
   and to residue_t again at the end (this could be speeded up, e.g. by
   precomputing 2^(3*wordsize) % m for mod_inv); the others use only
   mod_*() inline functions.
   Speed-critical functions need to be rewritten in assembly for REDC,
   but this is a start.
*/

void
mod_div3 (residue_t r, const residue_t a, const modulus_t m)
{
  const unsigned long a3 = a[0] % 3UL;
  unsigned long ml, m3;
#ifndef WANT_ASSERT_EXPENSIVE
  residue_t t;
#endif

#ifndef WANT_ASSERT_EXPENSIVE
  mod_init_noset0 (t, m);
  mod_set (t, a, m); /* r might overwrite a */
#endif

  if (a3 == 0UL)
    r[0] = a[0] / 3UL;
  else 
    {
      ml = mod_getmod_ul (m);
      m3 = ml % 3UL;
      ASSERT(m3 != 0UL);
      
      if (a3 + m3 == 3UL) /* Hence a3 == 1, m3 == 2 or a3 == 2, m3 == 1 */
	r[0] = a[0] / 3UL + ml / 3UL + 1UL;
      else /* a3 == 1, m3 == 1 or a3 == 2, m3 == 2 */
	r[0] = ml - (ml - a[0]) / 3UL;
    }

#ifndef WANT_ASSERT_EXPENSIVE
  mod_sub (t, t, r, m);
  mod_sub (t, t, r, m);
  mod_sub (t, t, r, m);
  ASSERT_EXPENSIVE (mod_is0 (t, m));
  mod_clear (t, m);
#endif
}


void 
mod_gcd (unsigned long *g, const residue_t r, const modulus_t m)
{
  unsigned long a, b, t;

  a = r[0]; /* This works the same for "a" in plain or Montgomery 
               representation */
  b = mod_getmod_ul (m);
  /* ASSERT (a < b); Should we require this? */
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

  g[0] = b;
}


/* Compute r = b^e. Here, e is an unsigned long */
void
mod_pow_ul (residue_t r, const residue_t b, const unsigned long e, 
            const modulus_t m)
{
  unsigned long mask = ~0UL - (~0UL >> 1); /* Only MSB set */
#ifndef NDEBUG
  unsigned long e1 = e;
#endif
  if (e == 0UL)
    {
      mod_set_ul (r, 1UL, m);
      return;
    }

  /* Assume r = 1 here for the invariant.
     r^mask * b^e is invariant, and is the result we want */

  /* Find highest set bit in e. */
  while ((e & mask) == 0UL)
    mask >>= 1; /* r = 1, so r^(mask/2) * b^e = r^mask * b^e  */

  /* Exponentiate */

  mod_set (r, b, m);       /* (r*b)^mask * b^(e-mask) = r^mask * b^e */
#ifndef NDEBUG
  e1 -= mask;
#endif

  while (mask > 1UL)
    {
      mod_mul (r, r, r, m);
      mask >>= 1;            /* (r^2)^(mask/2) * b^e = r^mask * b^e */
      if (e & mask)
        {
	  mod_mul (r, r, b, m);
#ifndef NDEBUG
          e1 -= mask;
#endif
        }
    }
  ASSERT (e1 == 0UL && mask == 1UL);
  /* Now e = 0, mask = 1, and r^mask * b^0 = r^mask is the result we want */
}


/* Compute r = b^e. Here e is a multiple precision integer 
   sum_{i=0}^{e_nrwords} e[i] * (machine word base)^i */
void
mod_pow_mp (residue_t r, const residue_t b, const unsigned long *e, 
            const int e_nrwords, const modulus_t m)
{
  unsigned long mask = ~0UL - (~0UL >> 1); /* Only MSB set */
  int i = e_nrwords - 1;

  if (e_nrwords == 0 || e[i] == 0UL)
    {
      mod_set_ul (r, 1UL, m);
      return;
    }

  /* Find highest set bit in e. */
  while ((e[i] & mask) == 0UL)
    mask >>= 1; /* r = 1, so r^(mask/2) * b^e = r^mask * b^e  */

  /* Exponentiate */

  mod_set (r, b, m);       /* (r*b)^mask * b^(e-mask) = r^mask * b^e */
  mask >>= 1;

  for ( ; i >= 0; i--)
    {
      while (mask > 0UL)
        {
          mod_mul (r, r, r, m);
          if (e[i] & mask)
            mod_mul (r, r, b, m);
          mask >>= 1;            /* (r^2)^(mask/2) * b^e = r^mask * b^e */
        }
      mask = ~0UL - (~0UL >> 1);
    }
}


/* Computes 2^e (mod m), where e is a multiple precision integer.
   Requires e != 0. The value of 2 in Montgomery representation 
   (i.e. 2*2^w (mod m) must be passed. */

void
mod_2pow_mp (residue_t r, const residue_t two, const unsigned long *e, 
             const int e_nrwords, const unsigned long e_mask, 
             const modulus_t m)
{
  residue_t t;
  unsigned long mask = e_mask;
  int i = e_nrwords - 1;

  ASSERT (e_nrwords != 0 && e[e_nrwords - 1] != 0);
  ASSERT ((e[e_nrwords - 1] & e_mask) == e_mask);

  mod_init (t, m);
  mod_set (t, two, m);
  mask >>= 1;

  for ( ; i >= 0; i--)
    {
      while (mask > 0UL)
        {
          mod_mul (t, t, t, m);
          if (e[i] & mask)
            mod_add (t, t, t, m);
          mask >>= 1;            /* (r^2)^(mask/2) * b^e = r^mask * b^e */
        }
      mask = ~0UL - (~0UL >> 1);
    }
    
  mod_set (r, t, m);
  mod_clear (t, m);
}


/* Compute r = V_e(b), where V_e(x) is the Chebyshev polynomial defined by
   V_e(x + 1/x) = x^e + 1/x^e. Here e is an unsigned long. */

void
mod_V_ul (residue_t r, const residue_t b, const residue_t two, 
          const unsigned long e, const modulus_t m)
{
  unsigned long mask = ~0UL - (~0UL >> 1); /* Only MSB set */
  residue_t r1;

  if (e == 0UL)
    {
      mod_set (r, two, m);
      return;
    }

  /* Find highest set bit in e. */
  while ((e & mask) == 0UL)
    mask >>= 1; /* r = 1, so r^(mask/2) * b^e = r^mask * b^e  */

  /* Exponentiate */

  mod_init_noset0 (r1, m);
  mod_set (r, b, m);         /* r = b = V_1 (b) */
  mod_mul (r1, b, b, m);
  mod_sub (r1, r1, two, m);  /* r1 = b^2 - 2 = V_2 (b) */
  mask >>= 1;

  /* Here r = V_j (b) and r1 = V_{j+1} (b) for j = 1 */

  while (mask > 0UL)
    {
      if (e & mask)
        {
          /* j -> 2*j+1. Compute V_{2j+1} and V_{2j+2} */
          mod_mul (r, r, r1, m);
          mod_sub (r, r, b, m); /* V_j * V_{j+1} - V_1 = V_{2j+1} */
          mod_mul (r1, r1, r1, m);
          mod_sub (r1, r1, two, m); /* (V_{j+1})^2 - 2 = V_{2j+2} */
        }
      else
        {
          /* j -> 2*j. Compute V_{2j} and V_{2j+1} */
          mod_mul (r1, r1, r, m);
          mod_sub (r1, r1, b, m); /* V_j * V_{j+1} - V_1 = V_{2j+1}*/
          mod_mul (r, r, r, m);
          mod_sub (r, r, two, m);
        }
      mask >>= 1;
    }

  mod_clear (r1, m);
}


/* Compute r = V_e(b), where V_e(x) is the Chebyshev polynomial defined by
   V_e(x + 1/x) = x^e + 1/x^e. Here e is a multiple precision integer 
   sum_{i=0}^{e_nrwords} e[i] * (machine word base)^i */

void
mod_V_mp (residue_t r, const residue_t b, const unsigned long *e, 
          const int e_nrwords, const modulus_t m)
{
  unsigned long mask = ~0UL - (~0UL >> 1); /* Only MSB set */
  int i = e_nrwords - 1;
  residue_t r1, two;

  mod_init_noset0 (two, m);
  mod_set_ul (two, 2UL, m);

  if (e_nrwords == 0 || e[i] == 0UL)
    {
      mod_set (r, two, m);
      mod_clear (two, m);
      return;
    }

  /* Find highest set bit in e. */
  while ((e[i] & mask) == 0UL)
    mask >>= 1; /* r = 1, so r^(mask/2) * b^e = r^mask * b^e  */

  /* Exponentiate */

  mod_init_noset0 (r1, m);
  mod_set (r, b, m);         /* r = b = V_1 (b) */
  mod_mul (r1, b, b, m);
  mod_sub (r1, r1, two, m);  /* r1 = b^2 - 2 = V_2 (b) */
  mask >>= 1;

  /* Here r = V_j (b) and r1 = V_{j+1} (b) for j = 1 */

  for ( ; i >= 0; i--)
    {
      while (mask > 0UL)
        {
          if (e[i] & mask)
	    {
	      /* j -> 2*j+1. Compute V_{2j+1} and V_{2j+2} */
	      mod_mul (r, r, r1, m);
	      mod_sub (r, r, b, m); /* V_j * V_{j+1} - V_1 = V_{2j+1} */
	      mod_mul (r1, r1, r1, m);
	      mod_sub (r1, r1, two, m); /* (V_{j+1})^2 - 2 = V_{2j+2} */
	    }
	  else
	    {
	      /* j -> 2*j. Compute V_{2j} and V_{2j+1} */
	      mod_mul (r1, r1, r, m);
	      mod_sub (r1, r1, b, m); /* V_j * V_{j+1} - V_1 = V_{2j+1}*/
	      mod_mul (r, r, r, m);
	      mod_sub (r, r, two, m);
	    }
          mask >>= 1;
        }
      mask = ~0UL - (~0UL >> 1);
    }
  mod_clear (two, m);
  mod_clear (r1, m);
}


/* Returns 1 if m is a strong probable prime wrt base b, 0 otherwise.
   b must be < m. */
int
mod_sprp (const residue_t b, const modulus_t m)
{
  residue_t r1;
  int i = 0, po2 = 1;
  unsigned long mm1, t;

  mm1 = mod_getmod_ul (m);
  ASSERT (b[0] < mm1);

  if (mm1 <= 3UL)
    return (mm1 >= 2UL);

  if (mm1 % 2UL == 0UL)
    return 0;

  /* Set mm1 to the odd part of m-1 */
  mm1 = (mm1 - 1) >> 1;
  while (mm1 % 2UL == 0UL)
    {
      po2++;
      mm1 >>= 1;
    }
  /* Hence, m-1 = mm1 * 2^po2 */

  mod_init_noset0 (r1, m);

  /* Exponentiate */
  mod_pow_ul (r1, b, mm1, m);

  t = mod_get_ul (r1, m);
  /* Now t == b^mm1 (mod m) */
#if defined(PARI)
  printf ("(Mod(%lu,%lu) ^ %lu) == %lu /* PARI */\n", 
	  b, mod_getmod_ul (m), mm1, t);
#endif
  
  /* If m is prime, then b^mm1 might be == 1 or == -1 (mod m) here */
  if (t == 1UL || t == mod_getmod_ul (m) - 1UL)
    i = 1;
  else
    {
      /* If m is a prime, then one of b^(2*mm1), b^(2^2*mm1), ..., 
	 b^(2^(po2 - 1)*mm1)  must be == -1 (mod m) */
      for ( ; po2 > 1; po2--)
	{
	  mod_mul (r1, r1, r1, m);
	  t = mod_get_ul (r1, m);
	  if (t == mod_getmod_ul (m) - 1UL)
	    {
	      i = 1;
	      break;
	    }
	}
    }

  mod_clear (r1, m);
  return i;
}


/* Put 1/s (mod t) in r and return 1 if s is invertible, 
   or set r to 0 and return 0 if s is not invertible */

/* We apply REDC twice on the input.
   Let sp = a*w (mod m), i.e. sp is $a$ in Montgomery representation.
   Then REDC(REDC(sp)) = a/w, and the inverse of that is w/a, which
   is the inverse we want in Montgomery representation. The extra REDC
   is faster than computing 1/a (mod m) first, then w/a (mod m) which 
   takes an integer divide instruction to get the remainder w*(1/a) (mod m). 
*/

int
mod_inv (residue_t r, const residue_t sp, const modulus_t m)
{
  long u1, v1;
  unsigned long u2, v2, s, t;
#ifndef NDEBUG
  residue_t tmp;
#endif

#ifndef NDEBUG
  /* Remember input in case r overwrites it */
  mod_init_noset0 (tmp, m);
  mod_set (tmp, sp, m);
#endif

  s = mod_get_ul (sp, m);
  s = mod_get_ul (&s, m); /* If we use REDC, do REDC again */
  t = mod_getmod_ul (m);

  ASSERT (t > 0UL);
  ASSERT (s < t);

  if (s == 0UL)
    {
      r[0] = 0UL; /* Not invertible */
#ifndef NDEBUG
      mod_clear (tmp, m);
#endif
      return 0;
    }

  if (s == 1UL)
    {
      r[0] = 1UL;
#ifndef NDEBUG
      mod_clear (tmp, m);
#endif
      return 1;
    }

  u1 = 1L;
  u2 = s;
  v1 = - (long) (t / s); /* No overflow, since s >= 2 */
  v2 = t % s;

  if (v2 == 1UL)
    {
       u1 = v1 + t;
    }
  else 
    {
      while (v2 != 0UL)
	{
	  unsigned long q;
	  /* unroll twice and swap u/v */
	  q = u2 / v2;
	  ASSERT_EXPENSIVE (q <= (unsigned long) LONG_MAX);
	  u1 = u1 - (long) q * v1;
	  u2 = u2 - q * v2;
	  
	  if (u2 == 0UL)
	    {
	      u1 = v1;
	      u2 = v2;
	      break;
	    }
	  
	  q = v2 / u2;
	  ASSERT_EXPENSIVE (q <= (unsigned long) LONG_MAX);
	  v1 = v1 - (long) q * u1;
	  v2 = v2 - q * u2;
	}
  
      if (u2 != 1UL)
	{
	  /* printf ("s=%lu t=%lu found %lu\n", s[0], t[0], u2); */
	  r[0] = 0UL; /* non-trivial gcd */
#ifndef NDEBUG
          mod_clear (tmp, m);
#endif
	  return 0;
	}

      if (u1 < 0L)
        u1 = u1 + t;
    }

    ASSERT ((unsigned long) u1 < t);
    
  r[0] = (unsigned long) u1;

#ifndef NDEBUG
  mod_mul (tmp, tmp, r, m);
  ASSERT(mod_get_ul (tmp, m) == 1UL);
  mod_clear (tmp, m);
#endif

  return 1;
}

int
mod_jacobi (const residue_t a_par, const modulus_t m_par)
{
  unsigned long a, m, s;
  int t = 1;

  a = mod_get_ul (a_par, m_par);
  m = mod_getmod_ul (m_par);
  ASSERT (a < m);
  
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
  printf ("kronecker(%lu, %lu) == %d\n", 
          mod_get_ul (a_par, m_par), mod_getmod_ul (m_par), t);
#endif
  return t;
}
