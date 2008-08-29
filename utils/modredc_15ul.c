#include <stdio.h>
#include "modredc_15ul.h"

#if defined(__GNUC__) && (__GNUC__ >= 4 || __GNUC__ >= 3 && __GNUC_MINOR__ >= 4)
/* Opteron prefers LOOKUP_TRAILING_ZEROS 1, 
   Core2 prefers LOOKUP_TRAILING_ZEROS 0 */
#ifndef LOOKUP_TRAILING_ZEROS
#define LOOKUP_TRAILING_ZEROS 1
#endif
#define ctzl(x) __builtin_ctzl(x)
#define clzl(x) __builtin_clzl(x)
#else
/* If we have no ctzl(), we always use the table lookup */
#define LOOKUP_TRAILING_ZEROS 1
#endif

#if LOOKUP_TRAILING_ZEROS
static const unsigned char trailing_zeros[256] = 
  {8,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
   5,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
   6,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
   5,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
   7,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
   5,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
   6,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
   5,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0};
#endif


void
modredc15ul_div3 (residueredc15ul_t r, const residueredc15ul_t a, 
		  const modulusredc15ul_t m)
{
  residueredc15ul_t t;
  unsigned long a3 = (a[1] % 256UL + a[1] / 256UL + 
		      a[0] % 256UL + a[0] / 256UL) % 3UL;
  const unsigned long m3 = (m[0].m[0] % 256UL + m[0].m[0] / 256UL +
			    m[0].m[1] % 256UL + m[0].m[1] / 256UL) % 3UL;

  ASSERT(m3 != 0UL);

  modredc15ul_init_noset0 (t, m);
  t[1] = a[1];
  t[0] = a[0];

  /* Make t[1]:t[0] divisible by 3 */
  if (a3 != 0UL)
    {
      if (a3 + m3 == 3UL) /* Hence a3 == 1, m3 == 2 or a3 == 2, m3 == 1 */
	{
	  ularith_add_2ul_2ul (&(t[0]), &(t[1]), m[0].m[0], m[0].m[1]);
	}
      else /* a3 == 1, m3 == 1 or a3 == 2, m3 == 2 */
	{
	  ularith_add_2ul_2ul (&(t[0]), &(t[1]), m[0].m[0], m[0].m[1]);
	  ularith_add_2ul_2ul (&(t[0]), &(t[1]), m[0].m[0], m[0].m[1]);
	}
      
      /* Now t[1]:t[0] is divisible by 3 */
      ASSERT_EXPENSIVE ((t[0] % 3UL + t[1] % 3UL) % 3UL == 0UL);
    }
  
  /* a = a1 * 2^w + a0, 3|a
     Let a = a' * 3 * 2^w + a'', a'' < 3 * 2^w. 
     3 | a'', a'' / 3 < 2^w
     So a / 3 = a' * w + a'' / 3
     a' = trunc(a1 / 3)
     a'' = a0 * 3^{-1} (mod 2^w)
     Hence we get the correct result with one one-word multiplication
     and one one-word truncating division by a small constant.
  */
  
  t[1] = t[1] / 3UL;
  if (sizeof (unsigned long) == 4)
    t[0] *= 0xaaaaaaabUL; /* 1/3 (mod 2^32) */
  else
    t[0] *= 0xaaaaaaaaaaaaaaabUL; /* 1/3 (mod 2^64) */

#ifdef WANT_ASSERT_EXPENSIVE
  modredc15ul_sub (r, a, t, m);
  modredc15ul_sub (r, r, t, m);
  modredc15ul_sub (r, r, t, m);
  ASSERT_EXPENSIVE (modredc15ul_is0 (r, m));
#endif
  r[0] = t[0];
  r[1] = t[1];
  modredc15ul_clear (t, m);
}


void
modredc15ul_div7 (residueredc15ul_t r, const residueredc15ul_t a, 
		  const modulusredc15ul_t m)
{
  const unsigned long w_mod_7 = (sizeof (unsigned long) == 4) ? 4UL : 2UL;
  const unsigned long a7 = ((a[1] % 7UL) * w_mod_7 + a[0] % 7UL) % 7UL;
  const unsigned long inv7[7] = {0,6,3,2,5,4,1}; /* inv7[i] = -1/i (mod 7) */
  unsigned long m7 = ((m[0].m[0] % 7UL) * w_mod_7 + m[0].m[1] % 7UL) % 7UL;
  residueredc15ul_t t;
  
  ASSERT(m7 != 0UL);
  
  modredc15ul_init_noset0 (t, m);
  t[1] = a[1];
  t[0] = a[0];
  
  /* Make t[1]:t[0] divisible by 7 */
  if (a7 != 0UL)
    {
      /* We want a+km == 0 (mod 7), so k = -a*m^{-1} (mod 7) */
      m7 = (inv7[m7] * a7) % 7UL;
      
      switch (m7) 
	{
	case 6: ularith_add_2ul_2ul (&(t[0]), &(t[1]), m[0].m[0], m[0].m[1]);
	case 5: ularith_add_2ul_2ul (&(t[0]), &(t[1]), m[0].m[0], m[0].m[1]);
	case 4: ularith_add_2ul_2ul (&(t[0]), &(t[1]), m[0].m[0], m[0].m[1]);
	case 3: ularith_add_2ul_2ul (&(t[0]), &(t[1]), m[0].m[0], m[0].m[1]);
	case 2: ularith_add_2ul_2ul (&(t[0]), &(t[1]), m[0].m[0], m[0].m[1]);
	case 1: ularith_add_2ul_2ul (&(t[0]), &(t[1]), m[0].m[0], m[0].m[1]);
	case 0: ;
	}
      
      /* Now t[1]:t[0] is divisible by 7 */
      ASSERT_EXPENSIVE (((t[1] % 7UL) * w_mod_7 + t[0] % 7UL) % 7UL == 0UL);
    }
  
  t[1] = t[1] / 7UL;
  if (sizeof (unsigned long) == 4)
    t[0] *= 0xb6db6db7UL; /* 1/7 (mod 2^32) */
  else
    t[0] *= 0x6db6db6db6db6db7UL; /* 1/7 (mod 2^64) */

#ifdef WANT_ASSERT_EXPENSIVE
  modredc15ul_sub (r, a, t, m);
  modredc15ul_sub (r, r, t, m);
  modredc15ul_sub (r, r, t, m);
  modredc15ul_sub (r, r, t, m);
  modredc15ul_sub (r, r, t, m);
  modredc15ul_sub (r, r, t, m);
  modredc15ul_sub (r, r, t, m);
  ASSERT_EXPENSIVE (modredc15ul_is0 (r, m));
#endif
  r[0] = t[0];
  r[1] = t[1];
  modredc15ul_clear (t, m);
}


void
modredc15ul_gcd (unsigned long *r, const residueredc15ul_t A, 
		 const modulusredc15ul_t m)
{
  unsigned long a[2], b[2];
#ifdef LOOKUP_TRAILING_ZEROS
  int sh;
#endif

  /* Since we do REDC arithmetic, we must have m odd */
  ASSERT_EXPENSIVE (m[0].m[0] & 1UL);

  if (A[0] == 0UL && A[1] == 0UL)
    {
      r[0] = m[0].m[0];
      r[1] = m[0].m[1];
      return;
    }

  a[0] = A[0];
  a[1] = A[1];
  b[0] = m[0].m[0];
  b[1] = m[0].m[1];

  while (b[1] != 0UL || b[0] != 0UL)
    {
      /* Make b odd */
#ifdef LOOKUP_TRAILING_ZEROS
      do {
	sh = trailing_zeros [(unsigned char) b[0]];
	ularith_shrd (&(b[0]), b[1], sh);
	*(long *) &(b[1]) >>= sh;
      } while (sh == 8);
#else
      if (b[0] == 0UL) /* ctzl does not like zero input */
	{
	  b[0] = b[1];
	  b[1] = 0UL;
	}
      sh = ctzl (b[0]);
      ularith_shrd (&(b[0]), b[1], sh);
      *(long *) &(b[1]) >>= sh;
#endif
      
      /* Try to make the low two bits of a[0] zero */
      ASSERT_EXPENSIVE (b[0] % 2UL == 1UL);
      if ((a[0] ^ b[0]) & 2UL)
	ularith_add_2ul_2ul (&(a[0]), &(a[1]), b[0], b[1]);
      else
	ularith_sub_2ul_2ul (&(a[0]), &(a[1]), b[0], b[1]);

      if (a[0] == 0UL && a[1] == 0UL)
	{
	  if ((long) b[1] < 0)
	    {
	      b[1] = -b[1];
	      if (b[0] != 0UL)
		b[1]--;
	      b[0] = -b[0];
	    }
	  r[0] = b[0];
	  r[1] = b[1];
	  return;
	}

#ifdef LOOKUP_TRAILING_ZEROS
      do {
	sh = trailing_zeros [(unsigned char) a[0]];
	ularith_shrd (&(a[0]), a[1], sh);
	*(long *) &(a[1]) >>= sh;
      } while (sh == 8);
#else
      if (a[0] == 0UL) /* ctzl does not like zero input */
	{
	  a[0] = a[1];
	  a[1] = 0UL;
	}
      sh = ctzl (a[0]);
      ularith_shrd (&(a[0]), a[1], sh);
      *(long *) &(a[1]) >>= sh;
#endif
      ASSERT_EXPENSIVE (a[0] & 1);

      if ((a[0] ^ b[0]) & 2)
	ularith_add_2ul_2ul (&(b[0]), &(b[1]), a[0], a[1]);
      else
	ularith_sub_2ul_2ul (&(b[0]), &(b[1]), a[0], a[1]);
    }

  if ((long) a[1] < 0)
    {
      a[1] = -a[1];
      if (a[0] != 0UL)
	a[1]--;
      a[0] = -a[0];
    }
  r[0] = a[0];
  r[1] = a[1];

  return;
}


/* Compute r = b^e. Here, e is an unsigned long */
void
modredc15ul_pow_ul (residueredc15ul_t r, const residueredc15ul_t b, 
		    const unsigned long e, const modulusredc15ul_t m)
{
  unsigned long mask = ~0UL - (~0UL >> 1); /* Only MSB set */
#ifndef NDEBUG
  unsigned long e1 = e;
#endif
  if (e == 0UL)
    {
      modredc15ul_set_ul (r, 1UL, m);
      return;
    }

  /* Assume r = 1 here for the invariant.
     r^mask * b^e is invariant, and is the result we want */

  /* Find highest set bit in e. */
  while ((e & mask) == 0UL)
    mask >>= 1; /* r = 1, so r^(mask/2) * b^e = r^mask * b^e  */

  /* Exponentiate */

  modredc15ul_set (r, b, m);       /* (r*b)^mask * b^(e-mask) = r^mask * b^e */
#ifndef NDEBUG
  e1 -= mask;
#endif

  while (mask > 1UL)
    {
      modredc15ul_mul (r, r, r, m);
      mask >>= 1;            /* (r^2)^(mask/2) * b^e = r^mask * b^e */
      if (e & mask)
        {
	  modredc15ul_mul (r, r, b, m);
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
modredc15ul_pow_mp (residueredc15ul_t r, const residueredc15ul_t b, 
		    const unsigned long *e, const int e_nrwords, 
		    const modulusredc15ul_t m)
{
  unsigned long mask = ~0UL - (~0UL >> 1); /* Only MSB set */
  int i = e_nrwords - 1;

  if (e_nrwords == 0 || e[i] == 0UL)
    {
      modredc15ul_set_ul (r, 1UL, m);
      return;
    }

  /* Find highest set bit in e. */
  while ((e[i] & mask) == 0UL)
    mask >>= 1; /* r = 1, so r^(mask/2) * b^e = r^mask * b^e  */

  /* Exponentiate */

  modredc15ul_set (r, b, m);       /* (r*b)^mask * b^(e-mask) = r^mask * b^e */
  mask >>= 1;

  for ( ; i >= 0; i--)
    {
      while (mask > 0UL)
        {
          modredc15ul_mul (r, r, r, m);
          if (e[i] & mask)
            modredc15ul_mul (r, r, b, m);
          mask >>= 1;            /* (r^2)^(mask/2) * b^e = r^mask * b^e */
        }
      mask = ~0UL - (~0UL >> 1);
    }
}


/* Computes 2^e (mod m), where e is a multiple precision integer.
   Requires e != 0. The value of 2 in Montgomery representation 
   (i.e. 2*2^w (mod m) must be passed. */

void
modredc15ul_2pow_mp (residueredc15ul_t r, const residueredc15ul_t two, 
		     const unsigned long *e, const int e_nrwords, 
		     const unsigned long e_mask, const modulusredc15ul_t m)
{
  residueredc15ul_t t;
  unsigned long mask = e_mask;
  int i = e_nrwords - 1;

  ASSERT (e_nrwords != 0 && e[e_nrwords - 1] != 0);
  ASSERT ((e[e_nrwords - 1] & e_mask) == e_mask);

  modredc15ul_init (t, m);
  modredc15ul_set (t, two, m);
  mask >>= 1;

  for ( ; i >= 0; i--)
    {
      while (mask > 0UL)
        {
          modredc15ul_mul (t, t, t, m);
          if (e[i] & mask)
            modredc15ul_add (t, t, t, m);
          mask >>= 1;            /* (r^2)^(mask/2) * b^e = r^mask * b^e */
        }
      mask = ~0UL - (~0UL >> 1);
    }
    
  modredc15ul_set (r, t, m);
  modredc15ul_clear (t, m);
}


/* Compute r = V_e(b), where V_e(x) is the Chebyshev polynomial defined by
   V_e(x + 1/x) = x^e + 1/x^e. Here e is an unsigned long. */

void
modredc15ul_V_ul (residueredc15ul_t r, const residueredc15ul_t b, 
		  const residueredc15ul_t two, const unsigned long e, 
		  const modulusredc15ul_t m)
{
  unsigned long mask = ~0UL - (~0UL >> 1); /* Only MSB set */
  residueredc15ul_t r1;

  if (e == 0UL)
    {
      modredc15ul_set (r, two, m);
      return;
    }

  /* Find highest set bit in e. */
  while ((e & mask) == 0UL)
    mask >>= 1; /* r = 1, so r^(mask/2) * b^e = r^mask * b^e  */

  /* Exponentiate */

  modredc15ul_init_noset0 (r1, m);
  modredc15ul_set (r, b, m);         /* r = b = V_1 (b) */
  modredc15ul_mul (r1, b, b, m);
  modredc15ul_sub (r1, r1, two, m);  /* r1 = b^2 - 2 = V_2 (b) */
  mask >>= 1;

  /* Here r = V_j (b) and r1 = V_{j+1} (b) for j = 1 */

  while (mask > 0UL)
    {
      if (e & mask)
        {
          /* j -> 2*j+1. Compute V_{2j+1} and V_{2j+2} */
          modredc15ul_mul (r, r, r1, m);
          modredc15ul_sub (r, r, b, m); /* V_j * V_{j+1} - V_1 = V_{2j+1} */
          modredc15ul_mul (r1, r1, r1, m);
          modredc15ul_sub (r1, r1, two, m); /* (V_{j+1})^2 - 2 = V_{2j+2} */
        }
      else
        {
          /* j -> 2*j. Compute V_{2j} and V_{2j+1} */
          modredc15ul_mul (r1, r1, r, m);
          modredc15ul_sub (r1, r1, b, m); /* V_j * V_{j+1} - V_1 = V_{2j+1}*/
          modredc15ul_mul (r, r, r, m);
          modredc15ul_sub (r, r, two, m);
        }
      mask >>= 1;
    }

  modredc15ul_clear (r1, m);
}


/* Compute r = V_e(b), where V_e(x) is the Chebyshev polynomial defined by
   V_e(x + 1/x) = x^e + 1/x^e. Here e is a multiple precision integer 
   sum_{i=0}^{e_nrwords} e[i] * (machine word base)^i */

void
modredc15ul_V_mp (residueredc15ul_t r, const residueredc15ul_t b, 
		  const unsigned long *e, const int e_nrwords, 
		  const modulusredc15ul_t m)
{
  unsigned long mask = ~0UL - (~0UL >> 1); /* Only MSB set */
  int i = e_nrwords - 1;
  residueredc15ul_t r1, two;

  modredc15ul_init_noset0 (two, m);
  modredc15ul_set_ul (two, 2UL, m);

  if (e_nrwords == 0 || e[i] == 0UL)
    {
      modredc15ul_set (r, two, m);
      modredc15ul_clear (two, m);
      return;
    }

  /* Find highest set bit in e. */
  while ((e[i] & mask) == 0UL)
    mask >>= 1; /* r = 1, so r^(mask/2) * b^e = r^mask * b^e  */

  /* Exponentiate */

  modredc15ul_init_noset0 (r1, m);
  modredc15ul_set (r, b, m);         /* r = b = V_1 (b) */
  modredc15ul_mul (r1, b, b, m);
  modredc15ul_sub (r1, r1, two, m);  /* r1 = b^2 - 2 = V_2 (b) */
  mask >>= 1;

  /* Here r = V_j (b) and r1 = V_{j+1} (b) for j = 1 */

  for ( ; i >= 0; i--)
    {
      while (mask > 0UL)
        {
          if (e[i] & mask)
	    {
	      /* j -> 2*j+1. Compute V_{2j+1} and V_{2j+2} */
	      modredc15ul_mul (r, r, r1, m);
	      modredc15ul_sub (r, r, b, m); /* V_j * V_{j+1} - V_1 = V_{2j+1} */
	      modredc15ul_mul (r1, r1, r1, m);
	      modredc15ul_sub (r1, r1, two, m); /* (V_{j+1})^2 - 2 = V_{2j+2} */
	    }
	  else
	    {
	      /* j -> 2*j. Compute V_{2j} and V_{2j+1} */
	      modredc15ul_mul (r1, r1, r, m);
	      modredc15ul_sub (r1, r1, b, m); /* V_j * V_{j+1} - V_1 = V_{2j+1}*/
	      modredc15ul_mul (r, r, r, m);
	      modredc15ul_sub (r, r, two, m);
	    }
          mask >>= 1;
        }
      mask = ~0UL - (~0UL >> 1);
    }
  modredc15ul_clear (two, m);
  modredc15ul_clear (r1, m);
}


/* Returns 1 if m is a strong probable prime wrt base b, 0 otherwise.
   b must be < m. */
int
modredc15ul_sprp (const residueredc15ul_t b, const modulusredc15ul_t m)
{
  residueredc15ul_t r1;
  int i = 0, po2 = 0;
  unsigned long mm1[2];

  modredc15ul_getmod_uls (mm1, m);

  if (mm1[1] == 0UL && mm1[0] <= 3UL)
    return (mm1[0] >= 2UL);

  if (mm1[0] % 2UL == 0UL)
    return 0;

  /* Let m-1 = 2^l * k, k odd. Set mm1 = k, po2 = l */
  mm1[0]--; /* No borrow since m is odd */
  while (mm1[0] % 2UL == 0UL) /* TODO: use ctzl */
    {
      po2++;
      ularith_shrd (&(mm1[0]), mm1[1], 1);
      mm1[1] >>= 1;
    }

  modredc15ul_init_noset0 (r1, m);

  /* Exponentiate */
  if (mm1[1] > 0UL)
    modredc15ul_pow_mp (r1, b, mm1, 2UL, m);
  else
    modredc15ul_pow_ul (r1, b, mm1[0], m);

  modredc15ul_get_uls (mm1, r1, m);
  
  /* If m is prime, then b^mm1 might be == 1 or == -1 (mod m) here */
  if ((mm1[1] == 0UL && mm1[0] == 1UL) || 
      (mm1[1] == m[0].m[1] && mm1[0] == m[0].m[0] - 1UL))
    i = 1;
  else
    {
      /* If m is a prime, then one of b^(2*mm1), b^(2^2*mm1), ..., 
	 b^(2^(po2 - 1)*mm1)  must be == -1 (mod m) */
      for ( ; po2 > 1; po2--)
	{
	  modredc15ul_mul (r1, r1, r1, m);
	  modredc15ul_get_uls (mm1, r1, m);
	  if (mm1[1] == m[0].m[1] && mm1[0] == m[0].m[0] - 1UL)
	    {
	      i = 1;
	      break;
	    }
	}
    }

  modredc15ul_clear (r1, m);
  return i;
}


int
modredc15ul_inv (residueredc15ul_t r, const residueredc15ul_t A, 
		 const modulusredc15ul_t m) 
{
  unsigned long a[2], b[2], u[2], v[2];
  int t, lsh;

  ASSERT_EXPENSIVE (modredc15ul_lt_2ul (A, m[0].m));
  ASSERT_EXPENSIVE (m[0].m[0] & 1UL);

  if (A[0] == 0UL && A[1] == 0UL)
    return 0;

  b[0] = m[0].m[0];
  b[1] = m[0].m[1];

  /* Let A = x*2^w, so we want the Montgomery representation of 1/x, 
     which is 2^w/x. We start by getting a = x */ 
  modredc15ul_get_uls (a, A, m);

  /* We simply set a = x/2^w and t=0. The result before 
     correction will be 2^(w+t)/x so we have to divide by t, which
     may be >64, so we may have to do a full and a variable width REDC. */
  /* Since t >= log_2(a+b) - 1, we can skip one REDC here if b[1] > 1 */
  if (b[1] > 1UL)
    t = -LONG_BIT;
  else
    {
      modredc15ul_frommontgomery (a, a, m);
      /* Now a = x/2^w */
      t = 0;
    }

  u[0] = 1UL; u[1] = 0UL; v[0] = 0UL; v[1] = 0UL;

  /* make a odd */
#if LOOKUP_TRAILING_ZEROS
  do {
    lsh = trailing_zeros [(unsigned char) a[0]];
    ularith_shrd (&(a[0]), a[1], lsh);
    a[1] >>= lsh;
    t += lsh;
  } while (lsh == 8);
#else
  if (a[0] == 0UL)
    {
      /* x86 bsf gives undefined result for zero input */
      a[0] = a[1];
      a[1] = 0UL;
      t += LONG_BIT;
    }
  lsh = ctzl(a[0]);
  ularith_shrd (&(a[0]), a[1], lsh);
  a[1] >>= lsh;
  t += lsh;
#endif

  // Here a and b are odd, and a < b
  do {
    /* Here, a and b are odd, 0 < a < b, u is odd and v is even */
    ASSERT_EXPENSIVE ((a[0] & 1UL) == 1UL);
    ASSERT_EXPENSIVE ((b[0] & 1UL) == 1UL);
    ASSERT_EXPENSIVE ((u[0] & 1UL) == 1UL);
    ASSERT_EXPENSIVE ((v[0] & 1UL) == 0UL);
    
    do {
      ularith_sub_2ul_2ul (&(b[0]), &(b[1]), a[0], a[1]);
      ularith_add_2ul_2ul (&(v[0]), &(v[1]), u[0], u[1]);
      if (b[0] == 0UL)
	{
	  b[0] = b[1];
	  b[1] = 0UL;
	  u[1] = u[0]; /* Shift left u by LONG_BIT */
	  u[0] = 0UL;
	  t += LONG_BIT; /* b[0] can be odd now, so lsh might be 0 below! */
	}
#if LOOKUP_TRAILING_ZEROS
      do {
	lsh = trailing_zeros [(unsigned char) b[0]];
#else
	lsh = ctzl(b[0]);
#endif
	ularith_shrd (&(b[0]), b[1], lsh);
	b[1] >>= lsh;
	t += lsh;
	ularith_shld (&(u[1]), u[0], lsh);
	u[0] <<= lsh;
#if LOOKUP_TRAILING_ZEROS
      } while (lsh == 8);
#endif
    } while (modredc15ul_lt_2ul(a, b)); /* ~50% branch taken :( */
    
    /* Here, a and b are odd, 0 < b =< a, u is even and v is odd */
    ASSERT_EXPENSIVE ((a[0] & 1UL) == 1UL);
    ASSERT_EXPENSIVE ((b[0] & 1UL) == 1UL);
    ASSERT_EXPENSIVE ((u[0] & 1UL) == 0UL);
    ASSERT_EXPENSIVE ((v[0] & 1UL) == 1UL);
    
    if (a[0] == b[0] && a[1] == b[1])
      break;
    
    /* Here, a and b are odd, 0 < b < a, u is even and v is odd */
    do {
      ularith_sub_2ul_2ul (&(a[0]), &(a[1]), b[0], b[1]);
      ularith_add_2ul_2ul (&(u[0]), &(u[1]), v[0], v[1]);
      
      if (a[0] == 0UL)
	{
	  a[0] = a[1];
	  a[1] = 0UL;
	  v[1] = v[0]; /* Shift left u by LONG_BIT */
	  v[0] = 0UL;
	  t += LONG_BIT;
	}
#if LOOKUP_TRAILING_ZEROS
      do {
	lsh = trailing_zeros [(unsigned char) a[0]];
#else
	lsh = ctzl(a[0]);
#endif
	ularith_shrd (&(a[0]), a[1], lsh);
	a[1] >>= lsh;
	t += lsh;
	ularith_shld (&(v[1]), v[0], lsh);
	v[0] <<= lsh;
#if LOOKUP_TRAILING_ZEROS
      } while (lsh == 8);
#endif

    } while (modredc15ul_lt_2ul(b, a)); /* about 50% branch taken :( */
    /* Here, a and b are odd, 0 < a =< b, u is odd and v is even */
  } while (a[0] != b[0] || a[1] != b[1]);
  
  if (a[1] != 0UL || a[0] != 1UL) /* Non-trivial GCD */
    return 0;

  ASSERT (t >= 0);

  /* Here, u = 2^w * 2^t / x. We want 2^w / x. */

  /* Here, the inverse of a is u/2^t mod b. To do the division by t,
     we use a variable-width REDC. We want to add a multiple of m to u
     so that the low t bits of the sum are 0 and we can right-shift by t
     with impunity. */
  while (t >= LONG_BIT)
    {
      /* Do full REDC */
      unsigned long s[4], k;
      k = u[0] * m[0].invm;
      ularith_mul_ul_ul_2ul (&(s[0]), &(s[1]), k, m[0].m[0]);
      if (s[0] != 0UL)
	s[1]++;
      s[2] = 0UL;
      ularith_add_ul_2ul (&(s[1]), &(s[2]), u[1]); /* s[2] <= 1 */
      ularith_mul_ul_ul_2ul (&(s[0]), &(s[3]), k, m[0].m[1]); /* s[3]<2^(w/2) */
      ularith_add_2ul_2ul (&(s[1]), &(s[2]), s[0], s[3]);
      u[0] = s[1];
      u[1] = s[2];
      t -= LONG_BIT;
    }

  if (t > 0)
    {
      unsigned long s[5], k;
      k = ((u[0] * m[0].invm) & ((1UL << t) - 1UL)); /* tlow <= 2^t-1 */
      ularith_mul_ul_ul_2ul (&(s[0]), &(s[1]), k, m[0].m[0]); 
      /* s[1]:s[0] <= (2^w-1)*(2^t-1) <= (2^w-1)*(2^(w-1)-1) */
      ularith_add_2ul_2ul (&(s[0]), &(s[1]), u[0], u[1]);
      /* s[1]:s[0] <= (2^w-1)*(2^(w-1)-1) + (m-1) < 2^(2w) */
      /* s[0] == 0 (mod 2^t) */
      ASSERT_EXPENSIVE ((s[0] & ((1UL << t) - 1UL)) == 0UL);
      s[2] = 0;
      ularith_mul_ul_ul_2ul (&(s[3]), &(s[4]), k, m[0].m[1]);
      ularith_add_2ul_2ul (&(s[1]), &(s[2]), s[3], s[4]);

      /* Now shift s[2]:s[1]:s[0] right by t */
      ularith_shrd (&(s[0]), s[1], t);
      ularith_shrd (&(s[1]), s[2], t);

      u[0] = s[0];
      u[1] = s[1];
    }

  r[0] = u[0];
  r[1] = u[1];
  return 1;
}


int
modredc15ul_jacobi (const residueredc15ul_t a_par, 
		    const modulusredc15ul_t m_par)
{
  unsigned long a[2], m[2], s;
  int t = 1;

  modredc15ul_get_uls (a, a_par, m_par);
  modredc15ul_getmod_uls (m, m_par);
  
  while (a[1] != 0UL || a[0] != 0UL)
  {
    while (a[0] % 2UL == 0UL) /* TODO speedup */
    {
      ularith_shrd (&(a[0]), a[1], 1);
      a[1] >>= 1;
      if (m[0] % 8UL == 3UL || m[0] % 8UL == 5UL)
        t = -t;
    }
    s = a[0]; a[0] = m[0]; m[0] = s; /* swap */
    s = a[1]; a[1] = m[1]; m[1] = s;
    if (a[0] % 4UL == 3UL && m[0] % 4UL == 3UL)
      t = -t;
    
    /* m is odd here */
    while (!modredc15ul_lt_2ul (a, m)) {
      int sh;
      ularith_sub_2ul_2ul (&(a[0]), &(a[1]), m[0], m[1]);
#if LOOKUP_TRAILING_ZEROS
      do {
	sh = trailing_zeros [(unsigned char) a[0]];
	ularith_shrd (&(a[0]), a[1], sh);
	a[1] >>= sh;
      } while (sh == 8);
#else
      if (a[0] == 0UL) /* ctzl does not like zero input */
	{
	  a[0] = a[1];
	  a[1] = 0UL;
	}
      sh = ctzl(a[0]);
      ularith_shrd (&(a[0]), a[1], lsh);
      a[1] >>= sh;
#endif
    }

  }
  if (m[1] != 0UL || m[0] != 1UL)
    t = 0;
  
#ifdef MODTRACE
  printf ("kronecker(%lu, %lu) == %d\n", 
          modredc15ul_get_ul (a_par, m_par), modredc15ul_getmod_ul (m_par), t);
#endif
  return t;
}
