/* Routines that are the same for modredc_15ul.c and modredc_2ul2.c */


/* Divide residue by 3. Returns 1 if division is possible, 0 otherwise.
   Assumes that a+3m does not overflow */

int
mod_div3 (residue_t r, const residue_t a, const modulus_t m)
{
  residue_t t;
  unsigned long a3 = (a[1] % 256UL + a[1] / 256UL + 
		      a[0] % 256UL + a[0] / 256UL) % 3UL;
  const unsigned long m3 = (m[0].m[0] % 256UL + m[0].m[0] / 256UL +
			    m[0].m[1] % 256UL + m[0].m[1] / 256UL) % 3UL;

  if (m3 == 0)
    return 0;

  mod_init_noset0 (t, m);
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
  
  r[1] = t[1] / 3UL;
#if LONG_BIT == 32
    r[0] = t[0] * 0xaaaaaaabUL; /* 1/3 (mod 2^32) */
#elif LONG_BIT == 64
    r[0] = t[0] * 0xaaaaaaaaaaaaaaabUL; /* 1/3 (mod 2^64) */
#else
#error Unknown word size
#endif

#ifdef WANT_ASSERT_EXPENSIVE
  mod_add (t, r, r, m);
  mod_add (t, t, r, m);
  ASSERT_EXPENSIVE (mod_equal (a, t, m));
#endif
  mod_clear (t, m);

  return 1;
}


/* Divide residue by 5. Returns 1 if division is possible, 0 otherwise */

int
mod_div5 (residue_t r, const residue_t a, const modulus_t m)
{
  /* inv_5[i] = -1/i (mod 5) */
  const unsigned long inv_5[5] = {0,4,2,3,1};
#if LONG_BIT == 32
  const unsigned long c = 0xcccccccdUL; /* 1/5 (mod 2^32) */
#elif LONG_BIT == 64
  const unsigned long c = 0xcccccccccccccccdUL; /* 1/5 (mod 2^64) */
#else
#error Unknown word size
#endif
  
  return mod_divn (r, a, 5UL, 1UL, inv_5, c, m);
}


/* Divide residue by 7. Returns 1 if division is possible, 0 otherwise */

int
mod_div7 (residue_t r, const residue_t a, const modulus_t m)
{
  const unsigned long w_mod_7 = (sizeof (unsigned long) == 4) ? 4UL : 2UL;
  /* inv_7[i] = -1/i (mod 7) */
  const unsigned long inv_7[7] = {0,6,3,2,5,4,1};
#if LONG_BIT == 32
  const unsigned long c = 0xb6db6db7UL; /* 1/7 (mod 2^32) */
#elif LONG_BIT == 64
  const unsigned long c = 0x6db6db6db6db6db7UL; /* 1/7 (mod 2^64) */
#else
#error Unknown word size
#endif

  return mod_divn (r, a, 7UL, w_mod_7, inv_7, c, m);
}


/* Divide residue by 11. Returns 1 if division is possible, 0 otherwise */

int
mod_div11 (residue_t r, const residue_t a, const modulus_t m)
{
  const unsigned long w_mod_11 = (sizeof (unsigned long) == 4) ? 4UL : 5UL;
  /* inv_11[i] = -1/i (mod 11) */
  const unsigned long inv_11[11] = {0, 10, 5, 7, 8, 2, 9, 3, 4, 6, 1}; 
#if LONG_BIT == 32
  const unsigned long c = 0xba2e8ba3UL; /* 1/11 (mod 2^32) */
#elif LONG_BIT == 64
  const unsigned long c = 0x2e8ba2e8ba2e8ba3UL; /* 1/11 (mod 2^64) */
#else
#error Unknown word size
#endif

  return mod_divn (r, a, 11UL, w_mod_11, inv_11, c, m);
}


/* Divide residue by 13. Returns 1 if division is possible, 0 otherwise */

int
mod_div13 (residue_t r, const residue_t a, const modulus_t m)
{
  const unsigned long w_mod_13 = (sizeof (unsigned long) == 4) ? 9UL : 3UL;
  /* inv_13[i] = -1/i (mod 13) */
  const unsigned long inv_13[13] = {0, 12, 6, 4, 3, 5, 2, 11, 8, 10, 9, 7, 1}; 
#if LONG_BIT == 32
  const unsigned long c = 0xc4ec4ec5UL; /* 1/13 (mod 2^32) */
#elif LONG_BIT == 64
  const unsigned long c = 0x4ec4ec4ec4ec4ec5UL; /* 1/13 (mod 2^64) */
#else
#error Unknown word size
#endif

  return mod_divn (r, a, 13UL, w_mod_13, inv_13, c, m);
}


void
mod_gcd (unsigned long *r, const residue_t A, const modulus_t m)
{
  unsigned long a[2], b[2];
  int sh;

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

  while (a[1] != 0UL || a[0] != 0UL)
    {
      /* Make a odd */
#if LOOKUP_TRAILING_ZEROS
      do {
	sh = trailing_zeros [(unsigned char) a[0]];
	ularith_shrd (&(a[0]), a[1], sh);
	*(long *) &(a[1]) >>= sh;
      } while (sh == 8);
#else
      if (a[0] == 0UL) /* ctzl does not like zero input */
	{
	  a[0] = a[1];
	  a[1] = ((long)a[1] < 0L) ? (unsigned long) (-1L) : 0UL;
	}
      sh = ularith_ctz (a[0]);
      ularith_shrd (&(a[0]), a[1], sh);
      *(long *) &(a[1]) >>= sh;
#endif
      
      /* Try to make the low two bits of b[0] zero */
      ASSERT_EXPENSIVE (a[0] % 2UL == 1UL);
      ASSERT_EXPENSIVE (b[0] % 2UL == 1UL);
      if ((a[0] ^ b[0]) & 2UL)
	ularith_add_2ul_2ul (&(b[0]), &(b[1]), a[0], a[1]);
      else
	ularith_sub_2ul_2ul (&(b[0]), &(b[1]), a[0], a[1]);

      if (b[0] == 0UL && b[1] == 0UL)
	{
	  if ((long) a[1] < 0L)
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

      /* Make b odd */
#if LOOKUP_TRAILING_ZEROS
      do {
	sh = trailing_zeros [(unsigned char) b[0]];
	ularith_shrd (&(b[0]), b[1], sh);
	*(long *) &(b[1]) >>= sh;
      } while (sh == 8);
#else
      if (b[0] == 0UL) /* ctzl does not like zero input */
	{
	  b[0] = b[1];
	  b[1] = ((long)b[1] < 0) ? (unsigned long) (-1L) : 0UL;
	}
      sh = ularith_ctz (b[0]);
      ularith_shrd (&(b[0]), b[1], sh);
      *(long *) &(b[1]) >>= sh;
#endif
      ASSERT_EXPENSIVE (a[0] % 2UL == 1UL);
      ASSERT_EXPENSIVE (b[0] % 2UL == 1UL);

      if ((a[0] ^ b[0]) & 2)
	ularith_add_2ul_2ul (&(a[0]), &(a[1]), b[0], b[1]);
      else
	ularith_sub_2ul_2ul (&(a[0]), &(a[1]), b[0], b[1]);
    }

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


/* Compute r = b^e. Here, e is an unsigned long */
void
mod_pow_ul (residue_t r, const residue_t b, const unsigned long e, 
	    const modulus_t m)
{
  unsigned long mask;
  residue_t t;
#ifndef NDEBUG
  unsigned long e1 = e;
#endif
  
  if (e == 0UL)
    {
      mod_set1 (r, m);
      return;
    }

  /* Assume t = 1 here for the invariant.
     r^mask * b^e is invariant, and is the result we want */

  /* Find highest set bit in e. */
  mask = (1UL << (LONG_BIT - 1)) >> ularith_clz (e);
  /* r = 1, so r^(mask/2) * b^e = r^mask * b^e  */

  /* Exponentiate */

  mod_init (t, m);
  mod_set (t, b, m);       /* (r*b)^mask * b^(e-mask) = r^mask * b^e */
#ifndef NDEBUG
  e1 -= mask;
#endif

  while (mask > 1UL)
    {
      mod_mul (t, t, t, m);
      mask >>= 1;            /* (r^2)^(mask/2) * b^e = r^mask * b^e */
      if (e & mask)
        {
	  mod_mul (t, t, b, m);
#ifndef NDEBUG
          e1 -= mask;
#endif
        }
    }
  ASSERT (e1 == 0UL && mask == 1UL);
  /* Now e = 0, mask = 1, and r^mask * b^0 = r^mask is the result we want */
  mod_set (r, t, m);
  mod_clear (t, m);
}


/* Compute r = 2^e. Here, e is an unsigned long */
void
mod_2pow_ul (residue_t r, const unsigned long e, const modulus_t m)
{
  unsigned long mask;
  residue_t t;

  if (e == 0UL)
    {
      mod_set1 (r, m);
      return;
    }

  mask = (1UL << (LONG_BIT - 1)) >> ularith_clz (e);
  ASSERT (e & mask);

  mod_init (t, m);
  mod_set1 (t, m);
  mod_add (t, t, t, m); /* t = 2 */
  mask >>= 1;

  while (mask > 0UL)
    {
      mod_mul (t, t, t, m);
      mod_intshl (t, t, (e & mask) ? 1 : 0);
      ularith_sub_2ul_2ul_ge (&(t[0]), &(t[1]), m[0].m[0], m[0].m[1]);
      mask >>= 1;
    }
  mod_set (r, t, m);
  mod_clear (t, m);
}


/* Compute r = b^e. Here e is a multiple precision integer 
   sum_{i=0}^{e_nrwords-1} e[i] * (machine word base)^i */
void
mod_pow_mp (residue_t r, const residue_t b, const unsigned long *e, 
	    const int e_nrwords, const modulus_t m)
{
  unsigned long mask;
  residue_t t;
  int i = e_nrwords - 1;

  if (e_nrwords == 0 || e[i] == 0UL)
    {
      mod_set1 (r, m);
      return;
    }

  /* Find highest set bit in e[i]. */
  mask = (1UL << (LONG_BIT - 1)) >> ularith_clz (e[i]);
  /* t = 1, so t^(mask/2) * b^e = t^mask * b^e  */

  /* Exponentiate */

  mod_init (t, m);
  mod_set (t, b, m);       /* (r*b)^mask * b^(e-mask) = r^mask * b^e */
  mask >>= 1;

  for ( ; i >= 0; i--)
    {
      while (mask > 0UL)
        {
          mod_mul (t, t, t, m);
          if (e[i] & mask)
            mod_mul (t, t, b, m);
          mask >>= 1;            /* (r^2)^(mask/2) * b^e = r^mask * b^e */
        }
      mask = ~0UL - (~0UL >> 1);
    }
  mod_set (r, t, m);
  mod_clear (t, m);
}


/* Compute r = 2^e. Here e is a multiple precision integer 
   sum_{i=0}^{e_nrwords-1} e[i] * (machine word base)^i */
void
mod_2pow_mp (residue_t r, const unsigned long *e, const int e_nrwords, 
	     const modulus_t m)
{
  unsigned long mask;
  residue_t t;
  int i = e_nrwords - 1;

  if (e_nrwords == 0 || e[i] == 0UL)
    {
      mod_set1 (r, m);
      return;
    }

  ASSERT (e[e_nrwords - 1] != 0);

  /* Find highest set bit in e[i]. */
  mask = (1UL << (LONG_BIT - 1)) >> ularith_clz (e[i]);
  /* t = 1, so t^(mask/2) * b^e = t^mask * b^e  */
  
  /* Exponentiate */

  mod_init (t, m);
  mod_set1 (t, m);
  mod_add (t, t, t, m); /* t = 2 */
  /* (t*2)^mask * b^(e-mask) = t^mask * b^e */
  mask >>= 1;

  for ( ; i >= 0; i--)
    {
      while (mask > 0UL)
        {
          mod_mul (t, t, t, m);
          mod_intshl (t, t, (e[i] & mask) ? 1 : 0);
          ularith_sub_2ul_2ul_ge (&(t[0]), &(t[1]), m[0].m[0], m[0].m[1]);
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
mod_V_ul (residue_t r, const residue_t b, 
		  const unsigned long e, const modulus_t m)
{
  unsigned long mask;
  residue_t t, t1, two;

  if (e == 0UL)
    {
      mod_set1 (r, m);
      mod_add (r, r, r, m);
      return;
    }
  
  /* Find highest set bit in e. */
  mask = (1UL << (LONG_BIT - 1)) >> ularith_clz (e);
  /* r = 1, so r^(mask/2) * b^e = r^mask * b^e  */

  /* Exponentiate */
  
  mod_init_noset0 (t, m);
  mod_init_noset0 (t1, m);
  mod_init_noset0 (two, m);
  mod_set1 (two, m);
  mod_add (two, two, two, m);
  mod_set (t, b, m);        /* t = b = V_1 (b) */
  mod_mul (t1, b, b, m);
  mod_sub (t1, t1, two, m);  /* r1 = b^2 - 2 = V_2 (b) */
  mask >>= 1;

  /* Here t = V_j (b) and t1 = V_{j+1} (b) for j = 1 */

  while (mask > 0UL)
    {
      if (e & mask)
        {
          /* j -> 2*j+1. Compute V_{2j+1} and V_{2j+2} */
          mod_mul (t, t, t1, m);
          mod_sub (t, t, b, m); /* V_j * V_{j+1} - V_1 = V_{2j+1} */
          mod_mul (t1, t1, t1, m);
          mod_sub (t1, t1, two, m); /* (V_{j+1})^2 - 2 = V_{2j+2} */
        }
      else
        {
          /* j -> 2*j. Compute V_{2j} and V_{2j+1} */
          mod_mul (t1, t1, t, m);
          mod_sub (t1, t1, b, m); /* V_j * V_{j+1} - V_1 = V_{2j+1}*/
          mod_mul (t, t, t, m);
          mod_sub (t, t, two, m);
        }
      mask >>= 1;
    }

  mod_set (r, t, m);
  mod_clear (t, m);
  mod_clear (t1, m);
  mod_clear (two, m);
}


/* Compute r = V_e(b), where V_e(x) is the Chebyshev polynomial defined by
   V_e(x + 1/x) = x^e + 1/x^e. Here e is a multiple precision integer 
   sum_{i=0}^{e_nrwords} e[i] * (machine word base)^i */

void
mod_V_mp (residue_t r, const residue_t b, 
		  const unsigned long *e, const int e_nrwords, 
		  const modulus_t m)
{
  unsigned long mask;
  int i = e_nrwords - 1;
  residue_t t, t1, two;

  if (e_nrwords == 0 || e[i] == 0UL)
    {
      mod_set1 (r, m);
      mod_add (r, r, r, m);
      return;
    }

  /* Find highest set bit in e. */
  mask = (1UL << (LONG_BIT - 1)) >> ularith_clz (e[i]);
  /* t = 1, so r^(mask/2) * b^e = t^mask * b^e  */

  /* Exponentiate */

  mod_init_noset0 (t, m);
  mod_init_noset0 (t1, m);
  mod_init_noset0 (two, m);
  mod_set1 (two, m);
  mod_add (two, two, two, m);
  mod_set (t, b, m);         /* t = b = V_1 (b) */
  mod_mul (t1, b, b, m);
  mod_sub (t1, t1, two, m);  /* t1 = b^2 - 2 = V_2 (b) */
  mask >>= 1;

  /* Here t = V_j (b) and t1 = V_{j+1} (b) for j = 1 */

  for ( ; i >= 0; i--)
    {
      while (mask > 0UL)
        {
          if (e[i] & mask)
	    {
	      /* j -> 2*j+1. Compute V_{2j+1} and V_{2j+2} */
	      mod_mul (t, t, t1, m);
	      mod_sub (t, t, b, m); /* V_j * V_{j+1} - V_1 = V_{2j+1} */
	      mod_mul (t1, t1, t1, m);
	      mod_sub (t1, t1, two, m); /* (V_{j+1})^2 - 2 = V_{2j+2} */
	    }
	  else
	    {
	      /* j -> 2*j. Compute V_{2j} and V_{2j+1} */
	      mod_mul (t1, t1, t, m);
	      mod_sub (t1, t1, b, m); /* V_j * V_{j+1} - V_1 = V_{2j+1}*/
	      mod_mul (t, t, t, m);
	      mod_sub (t, t, two, m);
	    }
          mask >>= 1;
        }
      mask = ~0UL - (~0UL >> 1);
    }
  
  mod_set (r, t, m);
  mod_clear (t, m);
  mod_clear (t1, m);
  mod_clear (two, m);
}


/* Returns 1 if m is a strong probable prime wrt base b, 0 otherwise.
   We assume m is odd. */
int
mod_sprp (const residue_t b, const modulus_t m)
{
  residue_t r, minusone;
  int i = 0, po2 = 0;
  modint_t mm1;

  mod_getmod_uls (mm1, m);

  if (mod_intequal_ul (mm1, 1UL))
    return 0;

  /* Let m-1 = 2^l * k, k odd. Set mm1 = k, po2 = l */
  mm1[0]--; /* No borrow since m is odd */
  if (mm1[0] == 0UL)
    {
      mm1[0] = mm1[1];
      mm1[1] = 0UL;
      po2 += LONG_BIT;
    }
  ASSERT (mm1[0] != 0UL);
  i = ularith_clz (mm1[0]);
  mod_intshr (mm1, mm1, i);
  po2 += i;

  mod_init_noset0 (r, m);
  mod_init_noset0 (minusone, m);
  mod_set1 (minusone, m);
  mod_neg (minusone, minusone, m);

  /* Exponentiate */
  if (mm1[1] > 0UL)
    mod_pow_mp (r, b, mm1, 2UL, m);
  else
    mod_pow_ul (r, b, mm1[0], m);

  /* If m is prime, then b^mm1 might be == 1 or == -1 (mod m) here */
  if (mod_is1 (r, m) || mod_equal (r, minusone, m))
    i = 1;
  else
    {
      /* If m is a prime, then one of b^(2*mm1), b^(2^2*mm1), ..., 
	 b^(2^(po2 - 1)*mm1)  must be == -1 (mod m) */
      for ( ; po2 > 1; po2--)
	{
	  mod_mul (r, r, r, m);
	  if (mod_equal (r, minusone, m))
	    {
	      i = 1;
	      break;
	    }
	}
    }

  mod_clear (r, m);
  mod_clear (minusone, m);
  return i;
}


/* Returns 1 if m is a strong probable prime wrt base 2, 0 otherwise. 
   We assume m is odd. */
int
mod_sprp2 (const modulus_t m)
{
  residue_t r, minusone;
  int i = 0, po2 = 0;
  modint_t mm1;

  mod_getmod_uls (mm1, m);

  /* Let m-1 = 2^l * k, k odd. Set mm1 = k, po2 = l */
  mm1[0]--; /* No borrow since m is odd */

  if (mm1[0] == 0UL)
    {
      mm1[0] = mm1[1];
      mm1[1] = 0UL;
      po2 += LONG_BIT;
    }
  ASSERT (mm1[0] != 0UL);
  i = ularith_clz (mm1[0]);
  mod_intshr (mm1, mm1, i);
  po2 += i;

  mod_init_noset0 (r, m);
  mod_init_noset0 (minusone, m);
  mod_set1 (minusone, m);
  mod_neg (minusone, minusone, m);

  /* Exponentiate */
  if (mm1[1] > 0UL)
    mod_2pow_mp (r, mm1, 2UL, m);
  else
    mod_2pow_ul (r, mm1[0], m);

  /* If m is prime, then b^mm1 might be == 1 or == -1 (mod m) here */
  if (mod_is1 (r, m) || mod_equal (r, minusone, m))
    i = 1;
  else
    {
      /* If m is a prime, then one of b^(2*mm1), b^(2^2*mm1), ..., 
	 b^(2^(po2 - 1)*mm1)  must be == -1 (mod m) */
      for ( ; po2 > 1; po2--)
	{
	  mod_mul (r, r, r, m);
	  if (mod_equal (r, minusone, m))
	    {
	      i = 1;
	      break;
	    }
	}
    }

  mod_clear (r, m);
  mod_clear (minusone, m);
  return i;
}


int
mod_jacobi (const residue_t a_par, const modulus_t m_par)
{
  modint_t a, m, s;
  int t = 1;

  mod_get_uls (a, a_par, m_par);
  mod_getmod_uls (m, m_par);
  
  while (!mod_intequal_ul (a, 0UL))
  {
    while (a[0] % 2UL == 0UL) /* TODO speedup */
    {
      mod_intshr (a, a, 1);
      if (m[0] % 8UL == 3UL || m[0] % 8UL == 5UL)
        t = -t;
    }
    mod_intset (s, a); /* swap a and m */
    mod_intset (a, m);
    mod_intset (m, s);
    if (a[0] % 4UL == 3UL && m[0] % 4UL == 3UL)
      t = -t;
    
    /* m is odd here */
    if (mod_intcmp (a, m) >= 0)
      {
        if (a[1] == 0UL)
          {
            a[0] %= m[0];
          }
        else
          {
            /* FIXME, slow and stupid */
            modint_t t;
            mod_intset (t, m);
            while (mod_intcmp (t, a) < 0)
              mod_intshl (t, t, 1);
            while (mod_intcmp (a, m) >= 0)
              {
		ularith_sub_2ul_2ul_ge (&(a[0]), &(a[1]), t[0], t[1]);
                mod_intshr (t, t, 1);
              }
          }
      }
  }
  if (m[1] != 0UL || m[0] != 1UL)
    t = 0;
  
#ifdef MODTRACE
  printf ("kronecker(%lu, %lu) == %d\n", 
          mod_get_ul (a_par, m_par), mod_getmod_ul (m_par), t);
#endif
  return t;
}
