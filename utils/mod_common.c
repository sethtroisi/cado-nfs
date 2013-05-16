#include "ularith.h"

/* Some functions that are implemented entirely on mod_*() functions, and which 
   thus can share source code for all arithmetic implentations */

#ifndef MOD_NO_SHARED_MOD_POW_UL
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
      mod_sqr (t, t, m);
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
#endif /* MOD_NO_SHARED_MOD_POW_UL */


#ifndef MOD_NO_SHARED_MOD_POW_MP
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
          mod_sqr (t, t, m);
          if (e[i] & mask)
            mod_mul (t, t, b, m);
          mask >>= 1;            /* (r^2)^(mask/2) * b^e = r^mask * b^e */
        }
      mask = ~0UL - (~0UL >> 1);
    }
  mod_set (r, t, m);
  mod_clear (t, m);
}
#endif /* MOD_NO_SHARED_MOD_POW_MP */


/* Returns 1 if r1 == 1 (mod m) or if r1 == -1 (mod m) or if
   one of r1^(2^1), r1^(2^2), ..., r1^(2^(po2-1)) == -1 (mod m),
   zero otherwise. Requires -1 (mod m) in minusone. */

static inline int
find_minus1 (residue_t r1, const residue_t minusone, const int po2, 
             const modulus_t m)
{
  int i;

  if (mod_is1 (r1, m) || mod_equal (r1, minusone, m))
    return 1;

  for (i = 1 ; i < po2; i++)
    {
      mod_sqr (r1, r1, m);
      if (mod_equal (r1, minusone, m))
        break;
    }

  return (i < po2) ? 1 : 0;
}


/* Compute r = V_e(b), where V_e(x) is the Chebyshev polynomial defined by
   V_e(x + 1/x) = x^e + 1/x^e. Here e is an unsigned long. */

void
mod_V_ul (residue_t r, const residue_t b, const unsigned long e, 
          const modulus_t m)
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
  mod_sqr (t1, b, m);
  mod_sub (t1, t1, two, m);  /* t1 = b^2 - 2 = V_2 (b) */
  mask >>= 1;

  /* Here t = V_j (b) and t1 = V_{j+1} (b) for j = 1 */

  while (mask > 0UL)
    {
      if (e & mask)
        {
          /* j -> 2*j+1. Compute V_{2j+1} and V_{2j+2} */
          mod_mul (t, t, t1, m);
          mod_sub (t, t, b, m); /* V_j * V_{j+1} - V_1 = V_{2j+1} */
          mod_sqr (t1, t1, m);
          mod_sub (t1, t1, two, m); /* (V_{j+1})^2 - 2 = V_{2j+2} */
        }
      else
        {
          /* j -> 2*j. Compute V_{2j} and V_{2j+1} */
          mod_mul (t1, t1, t, m);
          mod_sub (t1, t1, b, m); /* V_j * V_{j+1} - V_1 = V_{2j+1}*/
          mod_sqr (t, t, m);
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
mod_V_mp (residue_t r, const residue_t b, const unsigned long *e, 
          const int e_nrwords, const modulus_t m)
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
  mod_sqr (t1, b, m);
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
	      mod_sqr (t1, t1, m);
	      mod_sub (t1, t1, two, m); /* (V_{j+1})^2 - 2 = V_{2j+2} */
	    }
	  else
	    {
	      /* j -> 2*j. Compute V_{2j} and V_{2j+1} */
	      mod_mul (t1, t1, t, m);
	      mod_sub (t1, t1, b, m); /* V_j * V_{j+1} - V_1 = V_{2j+1}*/
	      mod_sqr (t, t, m);
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
