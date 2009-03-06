
/* This file defines some functions that work more or less the same 
   with mod_ul.h and modredc_ul.h. I.e. mod_div3() and mod_gcd() work
   unchanged with plain and Montgomery representation (so we can work on 
   the stored residue directly, whatever its representation is); 
   mod_jacobi() converts to plain "unsigned long" first, the others use 
   only mod_*() inline functions.
   Speed-critical functions need to be rewritten in assembly for REDC,
   but this is a start.
*/

void
mod_div3 (residue_t r, const residue_t a, const modulus_t m)
{
  const unsigned long a3 = a[0] % 3UL;
  unsigned long ml, m3;
  residue_t t;

  mod_init_noset0 (t, m);
  
  mod_set (t, a, m);
  if (a3 != 0UL)
    {
      ml = mod_getmod_ul (m);
      m3 = ml % 3UL;
      ASSERT(m3 != 0UL);
      
      if (a3 + m3 == 3UL) /* Hence a3 == 1, m3 == 2 or a3 == 2, m3 == 1 */
	t[0] = t[0] + ml;
      else /* a3 == 1, m3 == 1 or a3 == 2, m3 == 2 */
	t[0] = t[0] + 2UL * ml;
    }

  /* Now t[0] == a+k*m (mod 2^w) so that a+k*m is divisible by 3.
     (a+k*m)/3 < 2^w, so doing a division (mod 2^w) produces the 
     correct result. */
  
  if (sizeof (unsigned long) == 4)
    t[0] *= 0xaaaaaaabUL; /* 1/3 (mod 2^32) */
  else
    t[0] *= 0xaaaaaaaaaaaaaaabUL; /* 1/3 (mod 2^64) */
  
#ifdef WANT_ASSERT_EXPENSIVE
  mod_sub (r, a, t, m);
  mod_sub (r, r, t, m);
  mod_sub (r, r, t, m);
  ASSERT_EXPENSIVE (mod_is0 (r, m));
#endif

  mod_set (r, t, m);
  mod_clear (t, m);
}


void
mod_div7 (residue_t r, const residue_t a, const modulus_t m)
{
  unsigned long ml, m7, k;
  residue_t t;
  const unsigned long a7 = a[0] % 7UL;
  const unsigned long inv7[7] = {0,6,3,2,5,4,1}; /* inv7[i] = -1/i (mod 7) */

  ASSERT_EXPENSIVE (a[0] < mod_getmod_ul (m));

  mod_init_noset0 (t, m);
  mod_set (t, a, m);
  if (a7 != 0UL)
    {
      ml = mod_getmod_ul (m);
      m7 = ml % 7UL;
      ASSERT (m7 != 0UL);
      
      /* We want a+km == 0 (mod 7), so k = -a*m^{-1} (mod 7) */
      k = (a7 * inv7[m7]) % 7UL;
      ASSERT_EXPENSIVE ((k*m7 + a7) % 7UL == 0UL);
      t[0] = a[0] + k * ml;
    }

  /* Now t[0] == a+k*m (mod 2^w) so that a+k*m is divisible by 7.
     (a+k*m)/7 < 2^w, so doing a division (mod 2^w) produces the 
     correct result. */
  
  if (sizeof (unsigned long) == 4)
    t[0] *= 0xb6db6db7UL; /* 1/7 (mod 2^32) */
  else
    t[0] *= 0x6db6db6db6db6db7UL; /* 1/7 (mod 2^64) */
  
#ifdef WANT_ASSERT_EXPENSIVE
  ASSERT_EXPENSIVE (t[0] < mod_getmod_ul (m));
  mod_sub (r, a, t, m);
  mod_sub (r, r, t, m);
  mod_sub (r, r, t, m);
  mod_sub (r, r, t, m);
  mod_sub (r, r, t, m);
  mod_sub (r, r, t, m);
  mod_sub (r, r, t, m);
  ASSERT_EXPENSIVE (mod_is0 (r, m));
#endif

  mod_set (r, t, m);
  mod_clear (t, m);
}


void 
mod_gcd (modint_t g, const residue_t r, const modulus_t m)
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
      mod_set1 (r, m);
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


/* Compute r = 2^e. Here, e is an unsigned long */
void
mod_2pow_ul (residue_t r, const unsigned long e, const modulus_t m)
{
  unsigned long mask = ~0UL - (~0UL >> 1); /* Only MSB set */
#ifndef NDEBUG
  unsigned long e1 = e;
#endif
  if (e == 0UL)
    {
      mod_set1 (r, m);
      return;
    }

  /* Assume r = 1 here for the invariant.
     r^mask * b^e is invariant, and is the result we want */

  /* Find highest set bit in e. */
  while ((e & mask) == 0UL)
    mask >>= 1; /* r = 1, so r^(mask/2) * b^e = r^mask * b^e  */

  /* Exponentiate */

  mod_set1 (r, m);       
  mod_add (r, r, r, m);  /* (r*b)^mask * b^(e-mask) = r^mask * b^e */
#ifndef NDEBUG
  e1 -= mask;
#endif

  while (mask > 1UL)
    {
      mod_mul (r, r, r, m);
      mask >>= 1;            /* (r^2)^(mask/2) * b^e = r^mask * b^e */
      if (e & mask)
        {
	  mod_add (r, r, r, m);
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
      mod_set1 (r, m);
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
  mod_set1 (two, m);
  mod_add (two, two, two, m);

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


/* Returns 1 if m is a strong probable prime wrt base b, 0 otherwise. */
int
mod_sprp (const residue_t b, const modulus_t m)
{
  residue_t r1, minusone;
  int i = 0, po2 = 1;
  unsigned long mm1;

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
  mod_init_noset0 (minusone, m);
  mod_set1 (minusone, m);
  mod_neg (minusone, minusone, m);

  /* Exponentiate */
  mod_pow_ul (r1, b, mm1, m);

  /* Now r1 == b^mm1 (mod m) */
#if defined(PARI)
  printf ("(Mod(%lu,%lu) ^ %lu) == %lu /* PARI */\n", 
	  mod_get_ul (b, m), mod_getmod_ul (m), mm1, mod_get_ul (r1, m));
#endif
  
  /* If m is prime, then b^mm1 might be == 1 or == -1 (mod m) here */
  if (mod_is1 (r1, m) || mod_equal (r1, minusone, m))
    i = 1;
  else
    {
      /* If m is a prime, then one of b^(2*mm1), b^(2^2*mm1), ..., 
	 b^(2^(po2 - 1)*mm1)  must be == -1 (mod m) */
      for ( ; po2 > 1; po2--)
	{
	  mod_mul (r1, r1, r1, m);
	  if (mod_equal (r1, minusone, m))
	    {
	      i = 1;
	      break;
	    }
	}
    }

  mod_clear (r1, m);
  mod_clear (minusone, m);
  return i;
}


/* Returns 1 if m is a strong probable prime wrt base 2, 0 otherwise. */
int
mod_sprp2 (const modulus_t m)
{
  residue_t r1, minusone;
  int i = 0, po2 = 1;
  unsigned long mm1;

  mm1 = mod_getmod_ul (m);

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
  mod_init_noset0 (minusone, m);
  mod_set1 (minusone, m);
  mod_neg (minusone, minusone, m);

  /* Exponentiate */
  mod_2pow_ul (r1, mm1, m);

  /* Now r1 == b^mm1 (mod m) */
#if defined(PARI)
  printf ("(Mod(2,%lu) ^ %lu) == %lu /* PARI */\n", 
	  mod_getmod_ul (m), mm1, mod_get_ul (r1, m));
#endif
  
  /* If m is prime, then b^mm1 might be == 1 or == -1 (mod m) here */
  if (mod_is1 (r1, m) || mod_equal (r1, minusone, m))
    i = 1;
  else
    {
      /* If m is a prime, then one of b^(2*mm1), b^(2^2*mm1), ..., 
	 b^(2^(po2 - 1)*mm1)  must be == -1 (mod m) */
      for ( ; po2 > 1; po2--)
	{
	  mod_mul (r1, r1, r1, m);
	  if (mod_equal (r1, minusone, m))
	    {
	      i = 1;
	      break;
	    }
	}
    }

  mod_clear (r1, m);
  mod_clear (minusone, m);
  return i;
}


int
mod_isprime (const modulus_t m)
{
  residue_t b;
  const unsigned long n = mod_getmod_ul (m);
  int r = 0;
  
  if (n % 2UL == 0UL)
    {
      r = (n == 2UL);
#if defined(PARI)
      printf ("isprime(%lu) == %d /* PARI */\n", n, r);
#endif
      return r;
    }

  mod_init_noset0 (b, m);

  if (!mod_sprp2 (m))
    goto end;

  if (n < 2047UL)
    {
      r = 1;
      goto end;
    }

  if (n % 3UL == 1UL)
    {
      mod_set_ul_reduced (b, 7UL, m);
      if (mod_sprp (b, m))
        {
          mod_set_ul_reduced (b, 61UL, m);
          if (mod_sprp (b, m))
	    {
#if (ULONG_MAX > 4294967295UL)
	      if (n < 4759123141UL)
		{
		  return 1;
		  goto end;
		}

              mod_set_ul_reduced (b, 5UL, m);
              if (mod_sprp (b, m))
                {
                  /* These are the base 5,7,61 SPSP */
                  r = (n != 30926647201UL && n != 45821738881UL && 
                       n != 74359744201UL && n != 90528271681UL && 
                       n != 110330267041UL && n != 373303331521UL && 
                       n != 440478111067UL && n != 1436309367751UL && 
                       n != 1437328758421UL && n != 1858903385041UL && 
                       n != 4897239482521UL && n != 5026103290981UL && 
                       n != 5219055617887UL && n != 5660137043641UL && 
                       n != 6385803726241UL);
                }
#else
	      r = 1;
#endif
	    }
        }
    }
  else
    {
      /* Case n % 3 == 0, 2 */
      
      mod_set_ul_reduced (b, 3UL, m);
      if (mod_sprp (b, m))
        {
          if (n < 1373653UL)
            {
              r = 1;
            }
          else
            {
              mod_set_ul_reduced (b, 5UL, m);
              if (mod_sprp (b, m))
                {
#if (ULONG_MAX > 4294967295UL)
	          r = (n != 244970876021UL && n != 405439595861UL && 
		       n != 1566655993781UL && n != 3857382025841UL &&
		       n != 4074652846961UL && n != 5783688565841UL);
#else
	          r = 1;
#endif
                }
	    }
        }
    }
  
 end:
#if defined(PARI)
  printf ("isprime(%lu) == %d /* PARI */\n", n, r);
#endif
  mod_clear (b, m);
  return r;
}

int
mod_jacobi (const residue_t a_par, const modulus_t m_par)
{
  unsigned long a, m, s;
  int t = 1;

  /* We probably could stay in Montgomery representation here,
     a_par = a * 2^w, w even, so 2^w is a square and won't change
     the quadratic character */
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
