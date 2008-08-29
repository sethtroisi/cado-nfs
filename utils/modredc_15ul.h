/* Some functions for modular arithmetic with residues and modulus in
   unsigned long variables. The modulus can be up to 1.5 unsigned longs 
   in size (meaning 48 bits if unsigned long has 32 bits, or 96 bits if 
   unsigned long has 64 bits). Residues are stored in Montgomery form,
   reduction after multiplication is done with REDC. Due to inlining, 
   this file must be included in the caller's source code with #include */

/* Naming convention: all function start with modredc15ul, for 
   MODulus REDC 1.5 Unsigned Longs, followed by underscore, functionality of 
   function (add, mul, etc), and possibly underscore and specification of 
   what argument types the function takes (_ul, etc). */

#ifndef __MODREDC_15UL_H

#define __MODREDC_15UL_H

/**********************************************************************/
#include <assert.h>
#if defined(MODTRACE)
#include <stdio.h>
#endif
#include <limits.h>
#include "ularith.h"

#ifndef ASSERT
#define ASSERT(x)	assert(x)
#endif

/* Even simple assertions are relatively expensive in very simple functions.
   If we want them anyway to hunt a bug, define WANT_ASSERT_EXPENSIVE */
#ifdef WANT_ASSERT_EXPENSIVE
#define ASSERT_EXPENSIVE(x) ASSERT(x)
#else
#define ASSERT_EXPENSIVE(x)
#endif

/*********************************************************************/
/* Helper macros */
#ifndef	MAYBE_UNUSED
#if defined(__GNUC__)
#define MAYBE_UNUSED __attribute__ ((unused))
#else
#define MAYBE_UNUSED
#endif
#endif

#define MODREDC15UL_SIZE 2
#define MODREDC15UL_MAXBITS (LONG_BIT + LONG_BIT/2)

typedef unsigned long residueredc15ul_t[MODREDC15UL_SIZE];
typedef unsigned long modintredc15ul_t[MODREDC15UL_SIZE];
typedef struct { 
  unsigned long m[MODREDC15UL_SIZE];
  unsigned long invm;
} __modulusredc15ul_t;
typedef __modulusredc15ul_t modulusredc15ul_t[1];


/* ==================== Functions used internally ==================== */

MAYBE_UNUSED
static int
modredc15ul_lt_2ul (const unsigned long *a, const unsigned long *b)
{
  return (a[1] < b[1] || (a[1] == b[1] && a[0] < b[0]));
}

MAYBE_UNUSED
static inline void
modredc15ul_tomontgomery (residueredc15ul_t r, const residueredc15ul_t s,
			  const modulusredc15ul_t m)
{
  int i;

  r[0] = s[0];
  r[1] = s[1];
  /* TODO FIXME: ridiculously slow */
  for (i = 0; i < 2 * LONG_BIT; i++)
    {
      ularith_add_2ul_2ul (&(r[0]), &(r[1]), r[0], r[1]);
      if (! modredc15ul_lt_2ul (r, m[0].m))
	ularith_sub_2ul_2ul (&(r[0]), &(r[1]), m[0].m[0], m[0].m[1]);
    }
}

MAYBE_UNUSED
static inline void
modredc15ul_frommontgomery (residueredc15ul_t r, const residueredc15ul_t s,
			    const modulusredc15ul_t m)
{
  unsigned long t[5], k;
  
  /* Do two REDC steps */
  k = s[0] * m[0].invm;
  ularith_mul_ul_ul_2ul (&(t[0]), &(t[1]), k, m[0].m[0]);
  if (t[0] != 0UL)
    t[1]++;
  t[2] = 0;
  ularith_add_ul_2ul (&(t[1]), &(t[2]), s[1]); /* t[2] <= 1 */
  ularith_mul_ul_ul_2ul (&(t[0]), &(t[4]), k, m[0].m[1]); /* t[4] < 2^(w/2) */
  ularith_add_2ul_2ul (&(t[1]), &(t[2]), t[0], t[4]);
  /* Here, t2:t1:0 = s + k*m */
  k = t[1] * m[0].invm;
  ularith_mul_ul_ul_2ul (&(t[0]), &(t[4]), k, m[0].m[0]);
  if (t[0] != 0UL)
    t[4]++;
  t[3] = 0UL;
  ularith_add_ul_2ul (&(t[2]), &(t[3]), t[4]); /* t[3] <= 1 */
  ularith_mul_ul_ul_2ul (&(t[0]), &(t[4]), k, m[0].m[1]);
  ularith_add_2ul_2ul (&(t[2]), &(t[3]), t[0], t[4]);
  r[0] = t[2];
  r[1] = t[3];
}


/* Compute 1/n (mod 2^wordsize) */
MAYBE_UNUSED
static inline unsigned long
modredc15ul_invmodul (unsigned long n)
{
  unsigned long r;

  ASSERT (n % 2UL != 0UL);
  
  /* The square of an odd integer is always 1 (mod 8). So by
     initing r = m, the low three bits in the approximate inverse
     are correct. 
     When r = 1/m (mod 16), the 4th bit of r happens to be the
     XOR of bits 2, 3 and 4 of m. This gives us an approximate 
     inverse with the 4 lowest bits correct, so 3 (for 32 bit) or
     4 (for 64 bit) Newton iterations are enough. */
  r = n ^ ((n & 4UL) << 1) ^ ((n & 2UL) << 2);
  r = 2UL * r - r * r * n; /* Newton iteration */
  r = 2UL * r - r * r * n;
  r = 2UL * r - r * r * n;
  if (sizeof (unsigned long) > 4)
    r = 2UL * r - r * r * n;

  ASSERT_EXPENSIVE (r * n == 1UL);

  return r;
}


/* ================= Functions that are part of the API ================= */

/* Some functions for integers of the same width as the modulus */

MAYBE_UNUSED
static inline void
modredc15ul_intset (modintredc15ul_t r, const modintredc15ul_t s)
{
  r[0] = s[0];
  r[1] = s[1];
}

MAYBE_UNUSED
static inline void
modredc15ul_intset_ul (modintredc15ul_t r, const unsigned long s)
{
  r[0] = s;
  r[1] = 0UL;
}

MAYBE_UNUSED
static inline int
modredc15ul_intequal (const modintredc15ul_t a, const modintredc15ul_t b)
{
  return (a[0] == b[0] && a[1] == b[1]);
}

MAYBE_UNUSED
static inline int
modredc15ul_intequal_ul (const modintredc15ul_t a, const unsigned long b)
{
  return (a[0] == b && a[1] == 0UL);
}

MAYBE_UNUSED
static inline int
modredc15ul_intcmp (const modintredc15ul_t a, const modintredc15ul_t b)
{
  if (a[1] < b[1])
    return -1;
  if (a[1] > b[1])
    return 1;
  return (a[0] < b[0]) ? -1 : (a[0] == b[0]) ? 0 : 1;
}

MAYBE_UNUSED
static inline int
modredc15ul_intcmp_ul (const modintredc15ul_t a, const unsigned long b)
{
  if (a[1] > 0UL)
    return 1;
  return (a[0] < b) ? -1 : (a[0] == b) ? 0 : 1;
}

MAYBE_UNUSED
static inline int
modredc15ul_intfits_ul (const modintredc15ul_t a)
{
  return (a[1] == 0UL);
}

/* Returns the number of bits in a, that is, floor(log_2(n))+1. 
   For n==0 returns 0. */
MAYBE_UNUSED
static inline int
modredc15ul_intbits (const modintredc15ul_t a)
{
  int bits = 0;
  unsigned long n = a[0];
  if (a[1] > 0UL)
    {
      bits = LONG_BIT;
      n = a[1];
    }
  while (n > 0UL)
    {
      bits++;
      n >>= 1;
    }
  return bits;
}

/* r = n/d. We require d|n */
MAYBE_UNUSED
static inline void
modredc15ul_intdivexact (modintredc15ul_t r, const modintredc15ul_t n,
                         const modintredc15ul_t d)
{
  modintredc15ul_t n1, d1;
  unsigned long invf, r0, k0, k1;
#ifdef WANT_ASSERT_EXPENSIVE
  unsigned long s0 = n[0], s1 = n[1];
#endif
  
  n1[0] = n[0];
  n1[1] = n[1];
  d1[0] = d[0];
  d1[1] = d[1];

  /* Make d odd. TODO: use ctzl */
  while (d1[0] % 2 == 0UL)
    {
      ASSERT_EXPENSIVE (n1[0] % 2 == 0UL);
      ularith_shrd (&(d1[0]), d1[1], 1);
      d1[1] >>= 1;
      ularith_shrd (&(n1[0]), n1[1], 1);
      n1[1] >>= 1;
    }
  
  invf = modredc15ul_invmodul (d1[0]);
  r0 = invf * n1[0];
  ularith_mul_ul_ul_2ul (&k0, &k1, r0, d1[0]);
  ularith_sub_2ul_2ul (&(n1[0]), &(n1[1]), k0, k1);
  n1[1] -= r0 * d1[1];
  ASSERT (n1[0] == 0UL);
  r[0] = r0;
  r[1] = invf * n1[1];
#ifdef WANT_ASSERT_EXPENSIVE
  ASSERT_EXPENSIVE (d[1] == 0UL || r[1] == 0UL);
  ularith_mul_ul_ul_2ul (&k0, &k1, r[0], d[0]);
  ularith_sub_2ul_2ul (&s0, &s1, k0, k1);
  ASSERT_EXPENSIVE (s0 == 0UL);
  ularith_mul_ul_ul_2ul (&k0, &k1, r[0], d[1]);
  ASSERT_EXPENSIVE (k1 == 0UL);
  s1 -= k0;
  ularith_mul_ul_ul_2ul (&k0, &k1, r[1], d[0]);
  ASSERT_EXPENSIVE (k1 == 0UL);
  s1 -= k0;
  ASSERT_EXPENSIVE (s1 == 0UL);
#endif
}


/* Functions for the modulus */

MAYBE_UNUSED
static inline void
modredc15ul_initmod_ul (modulusredc15ul_t m, const unsigned long s)
{
  m[0].m[0] = s;
  m[0].m[1] = 0UL;
  m[0].invm = -modredc15ul_invmodul (s);
}


/* Init the modulus from a multi-word integer. s must point to an array of
   at least two unsigned longs, where s[0] is the low word of the modulus, 
   and s[1] is the high word. */
MAYBE_UNUSED
static inline void
modredc15ul_initmod_uls (modulusredc15ul_t m, const modintredc15ul_t s)
{
  m[0].m[0] = s[0];
  m[0].m[1] = s[1];
  m[0].invm = -modredc15ul_invmodul (s[0]);
}


/* Returns the modulus as an unsigned long if the modulus fits in an 
   unsigned long. Return 0 if it doesn't fit. */
MAYBE_UNUSED
static inline unsigned long
modredc15ul_getmod_ul (const modulusredc15ul_t m)
{
  if (m[0].m[1] > 0UL)
    return 0UL;
  return m[0].m[0];
}


/* Returns the modulus to an array of unsigned longs. */
MAYBE_UNUSED
static inline void
modredc15ul_getmod_uls (modintredc15ul_t r, const modulusredc15ul_t m)
{
  r[0] = m[0].m[0];
  r[1] = m[0].m[1];
}


MAYBE_UNUSED
static inline void
modredc15ul_clearmod (modulusredc15ul_t m MAYBE_UNUSED)
{
  return;
}


/* Functions for residues */

/* Initialises a residueredc15ul_t and sets it to zero */
MAYBE_UNUSED
static inline void
modredc15ul_init (residueredc15ul_t r, const modulusredc15ul_t m MAYBE_UNUSED)
{
  r[0] = 0UL;
  r[1] = 0UL;
}


/* Initialises a residueredc15ul_t, but does not set it to zero. For fixed 
   length residueredc15ul_t, that leaves nothing to do at all. */
MAYBE_UNUSED
static inline void
modredc15ul_init_noset0 (residueredc15ul_t r MAYBE_UNUSED, 
			 const modulusredc15ul_t m MAYBE_UNUSED)
{
  return;
}


MAYBE_UNUSED
static inline void
modredc15ul_clear (residueredc15ul_t r MAYBE_UNUSED, 
		   const modulusredc15ul_t m MAYBE_UNUSED)
{
  return;
}


MAYBE_UNUSED
static inline void
modredc15ul_set (residueredc15ul_t r, const residueredc15ul_t s, 
		 const modulusredc15ul_t m MAYBE_UNUSED)
{
  ASSERT_EXPENSIVE (modredc15ul_lt_2ul (s, m[0].m));
  r[0] = s[0];
  r[1] = s[1];
}


MAYBE_UNUSED
static inline void
modredc15ul_set_ul (residueredc15ul_t r, const unsigned long s, 
		    const modulusredc15ul_t m)
{
  r[0] = s;
  if (m[0].m[1] == 0UL)
    r[0] %= m[0].m[0];
  r[1] = 0UL;
  modredc15ul_tomontgomery (r, r, m);
}


/* Sets the residueredc15ul_t to the class represented by the integer s. 
   Assumes that s is reduced (mod m), i.e. 0 <= s < m */

MAYBE_UNUSED
static inline void
modredc15ul_set_ul_reduced (residueredc15ul_t r, const unsigned long s, 
			    const modulusredc15ul_t m MAYBE_UNUSED)
{
  ASSERT (m[0].m[1] > 0UL || s < m[0].m[0]);
  r[0] = s;
  r[1] = 0UL;
  modredc15ul_tomontgomery (r, r, m);
}


MAYBE_UNUSED
static inline void
modredc15ul_set_uls (residueredc15ul_t r, const modintredc15ul_t s, 
		     const modulusredc15ul_t m)
{
  r[0] = s[0];
  r[1] = s[1];
  if (! modredc15ul_lt_2ul (s, m[0].m))
    {
      /* Do reduction */
      if (m[0].m[1] == 0UL)
	{
	  r[1] %= m[0].m[0]; /* To avoid quotient overflow in next division */
	  ularith_div_2ul_ul_ul_r (&(r[0]), r[0], r[1], m[0].m[0]);
	  r[1] = 0UL;
	}
      else
	{
	  /* FIXME, slow and stupid */
	  unsigned long t[2];
	  t[0] = m[0].m[0];
	  t[1] = m[0].m[1];
	  while (t[1] < r[1] / 2)
	    {
	      ularith_shld (&(t[1]), t[0], 1);
	      t[0] <<= 1;
	    }
	  while (! modredc15ul_lt_2ul (r, m[0].m))
	    {
	      if (! modredc15ul_lt_2ul (r, t))
		ularith_sub_2ul_2ul (&(r[0]), &(r[1]), t[0], t[1]);
	      ularith_shrd (&(t[0]), t[1], 1);
	      t[1] >>= 1;
	    }
	}
    }
  ASSERT (modredc15ul_lt_2ul (r, m[0].m));
  modredc15ul_tomontgomery (r, r, m);
}


MAYBE_UNUSED
static inline void
modredc15ul_set_uls_reduced (residueredc15ul_t r, const modintredc15ul_t s, 
			     const modulusredc15ul_t m)
{
  ASSERT (modredc15ul_lt_2ul (s, m[0].m));
  r[0] = s[0];
  r[1] = s[1];
  modredc15ul_tomontgomery (r, r, m);
}


/* This one is so trivial that we don't really require m in the
 * interface. For interface homogeneity we make it take the m parameter 
 * anyway.
 */
MAYBE_UNUSED 
static inline void 
modredc15ul_set0 (residueredc15ul_t r, const modulusredc15ul_t m MAYBE_UNUSED) 
{ 
  r[0] = 0UL; 
  r[1] = 0UL;
}


MAYBE_UNUSED 
static inline void 
modredc15ul_set1 (residueredc15ul_t r, const modulusredc15ul_t m) 
{ 
  r[0] = 1UL;
  r[1] = 0UL;
  modredc15ul_tomontgomery (r, r, m);
}


/* Exchanges the values of the two arguments */

MAYBE_UNUSED
static inline void
modredc15ul_swap (residueredc15ul_t a, residueredc15ul_t b, 
		  const modulusredc15ul_t m MAYBE_UNUSED)
{
  unsigned long t0, t1;
  ASSERT_EXPENSIVE (modredc15ul_lt_2ul (a, m[0].m));
  ASSERT_EXPENSIVE (modredc15ul_lt_2ul (b, m[0].m));
  t0 = a[0];
  t1 = a[1];
  a[0] = b[0];
  a[1] = b[1];
  b[0] = t0;
  b[1] = t1;
}
                          

/* Returns the least significant unsigned long of the residue. How to signal
   if the residue does not fit in one unsigned long? */

MAYBE_UNUSED
static inline unsigned long
modredc15ul_get_ul (const residueredc15ul_t s, 
		    const modulusredc15ul_t m MAYBE_UNUSED)
{
  unsigned long t[2];
  ASSERT_EXPENSIVE (modredc15ul_lt_2ul (s, m[0].m));
  modredc15ul_frommontgomery (t, s, m);
  ASSERT (t[1] == 0UL);
  return t[0];
}


/* Returns the residue into an array of unsigned longs */

MAYBE_UNUSED
static inline void
modredc15ul_get_uls (modintredc15ul_t r, const residueredc15ul_t s, 
		     const modulusredc15ul_t m MAYBE_UNUSED)
{
  ASSERT_EXPENSIVE (modredc15ul_lt_2ul (s, m[0].m));
  modredc15ul_frommontgomery (r, s, m);
}


MAYBE_UNUSED
static inline int
modredc15ul_equal (const residueredc15ul_t a, const residueredc15ul_t b, 
		   const modulusredc15ul_t m MAYBE_UNUSED)
{
  ASSERT_EXPENSIVE (modredc15ul_lt_2ul (a, m[0].m));
  ASSERT_EXPENSIVE (modredc15ul_lt_2ul (b, m[0].m));
  return (a[1] == b[1] && a[0] == b[0]);
}


MAYBE_UNUSED
static inline int
modredc15ul_is0 (const residueredc15ul_t a, const modulusredc15ul_t m MAYBE_UNUSED)
{
  ASSERT_EXPENSIVE (modredc15ul_lt_2ul (a, m[0].m));
  return (a[1] == 0UL && a[0] == 0UL);
}


MAYBE_UNUSED
static inline int
modredc15ul_is1 (const residueredc15ul_t a, const modulusredc15ul_t m MAYBE_UNUSED)
{
  unsigned long t[2];
  ASSERT_EXPENSIVE (modredc15ul_lt_2ul (a, m[0].m));
  modredc15ul_frommontgomery (t, a, m);
  return (t[1] == 0UL && t[0] == 1UL);
}


MAYBE_UNUSED
static inline void
modredc15ul_add (residueredc15ul_t r, const residueredc15ul_t a, 
		 const residueredc15ul_t b, const modulusredc15ul_t m)
{
  unsigned long t1 = a[0], t2 = a[1];
  
  ASSERT_EXPENSIVE (modredc15ul_lt_2ul (a, m[0].m));
  ASSERT_EXPENSIVE (modredc15ul_lt_2ul (b, m[0].m));

  ularith_add_2ul_2ul (&t1, &t2, b[0], b[1]);
  
  r[0] = t1;
  r[1] = t2;
  if (! modredc15ul_lt_2ul (r, m[0].m))
    ularith_sub_2ul_2ul (&(r[0]), &(r[1]), m[0].m[0], m[0].m[1]);
}


MAYBE_UNUSED
static inline void
modredc15ul_sub (residueredc15ul_t r, const residueredc15ul_t a, 
		 const residueredc15ul_t b, const modulusredc15ul_t m)
{
  unsigned long t1 = a[0], t2 = a[1];
  
  ASSERT_EXPENSIVE (modredc15ul_lt_2ul (a, m[0].m));
  ASSERT_EXPENSIVE (modredc15ul_lt_2ul (b, m[0].m));

  ularith_sub_2ul_2ul (&t1, &t2, b[0], b[1]);
  
  if (t2 > a[1] || (t2 == a[1] && t1 > a[0]))
    ularith_add_2ul_2ul (&t1, &t2, m[0].m[0], m[0].m[1]);
  
  r[0] = t1;
  r[1] = t2;
}


MAYBE_UNUSED
static inline void
modredc15ul_add_ul (residueredc15ul_t r, const residueredc15ul_t a,
		    const unsigned long b, const modulusredc15ul_t m)
{
  residueredc15ul_t t;

  /* TODO: speed up */
  ASSERT_EXPENSIVE (modredc15ul_lt_2ul (a, m[0].m));
  modredc15ul_init_noset0 (t, m);
  modredc15ul_set_ul (t, b, m);
  modredc15ul_add (r, a, t, m);
  modredc15ul_clear (t, m);
}


MAYBE_UNUSED
static inline void
modredc15ul_sub_ul (residueredc15ul_t r, const residueredc15ul_t a,
		    const unsigned long b, const modulusredc15ul_t m)
{
  residueredc15ul_t t;

  /* TODO: speed up */
  ASSERT_EXPENSIVE (modredc15ul_lt_2ul (a, m[0].m));
  modredc15ul_init_noset0 (t, m);
  modredc15ul_set_ul (t, b, m);
  modredc15ul_sub (r, a, t, m);
  modredc15ul_clear (t, m);
}


MAYBE_UNUSED
static inline void
modredc15ul_neg (residueredc15ul_t r, const residueredc15ul_t a, 
		 const modulusredc15ul_t m)
{
  ASSERT_EXPENSIVE (modredc15ul_lt_2ul (a, m[0].m));
  if (a[0] == 0UL && a[1] == 0UL)
    modredc15ul_set (r, a, m);
  else
    {
      unsigned long t1 = m[0].m[0], t2 = m[0].m[1];
      ularith_sub_2ul_2ul (&t1, &t2, a[0], a[1]);
      r[0] = t1;
      r[1] = t2;
    }
}


MAYBE_UNUSED
static inline void
modredc15ul_div2 (residueredc15ul_t r, const residueredc15ul_t a, 
		  const modulusredc15ul_t m)
{
  ASSERT_EXPENSIVE (m[0].m[0] % 2UL != 0UL);
  r[0] = a[0];
  r[1] = a[1];
  if (r[0] % 2UL == 1UL)
    ularith_add_2ul_2ul (&(r[0]), &(r[1]), m[0].m[0], m[0].m[1]);
  ularith_shrd (&(r[0]), r[1], 1);
  r[1] >>= 1;
}


MAYBE_UNUSED
static inline void
modredc15ul_mul (residueredc15ul_t r, const residueredc15ul_t a, 
               const residueredc15ul_t b, const modulusredc15ul_t m)
{
  unsigned long pl, ph, t[4], k;
  
  ASSERT_EXPENSIVE (modredc15ul_lt_2ul (a, m[0].m));
  ASSERT_EXPENSIVE (modredc15ul_lt_2ul (b, m[0].m));
#if defined(MODTRACE)
  printf ("((%lu * 2^%d + %lu) * (%lu * 2^%d + %lu) / 2^%d) %% "
	  "(%lu * 2^%d + %lu)", 
          a[1], LONG_BIT, a[0], b[1], LONG_BIT, b[0], 2 * LONG_BIT, 
	  m[0].m[1], LONG_BIT, m[0].m[0]);
#endif
  
  ularith_mul_ul_ul_2ul (&(t[0]), &(t[1]), a[0], b[0]);
  k = t[0] * m[0].invm;
  ularith_mul_ul_ul_2ul (&pl, &ph, k, m[0].m[0]);
  if (pl != 0UL)
    ph++;
  t[2] = 0UL;
  ularith_add_ul_2ul (&(t[1]), &(t[2]), ph); /* t[2] <= 1 */
  ularith_mul_ul_ul_2ul (&pl, &ph, k, m[0].m[1]);
  ularith_add_2ul_2ul (&(t[1]), &(t[2]), pl, ph); 
  ularith_mul_ul_ul_2ul (&pl, &ph, a[1], b[0]);
  ularith_add_2ul_2ul (&(t[1]), &(t[2]), pl, ph);
  ularith_mul_ul_ul_2ul (&pl, &ph, a[0], b[1]);
  ularith_add_2ul_2ul (&(t[1]), &(t[2]), pl, ph);
  t[3] = 0UL;
  pl = a[1] * b[1];
  ularith_add_ul_2ul (&(t[2]), &(t[3]), pl);
  k = t[1] * m[0].invm;
  ularith_mul_ul_ul_2ul (&pl, &ph, k, m[0].m[0]);
  if (pl != 0UL)
    ph++;
  ularith_add_ul_2ul (&(t[2]), &(t[3]), ph);
  ularith_mul_ul_ul_2ul (&pl, &ph, k, m[0].m[1]);
  ularith_add_2ul_2ul (&(t[2]), &(t[3]), pl, ph);

  /* Result may be larger than m, but is < 2*m */

  if (!modredc15ul_lt_2ul (t + 2, m[0].m))
    ularith_sub_2ul_2ul (&(t[2]), &(t[3]), m[0].m[0], m[0].m[1]);

  r[0] = t[2];
  r[1] = t[3];
#if defined(MODTRACE)
  printf (" == (%lu * 2^%d + %lu) /* PARI */ \n", r[1], LONG_BIT, r[0]);
#endif
  ASSERT_EXPENSIVE (modredc15ul_lt_2ul (r, m[0].m));
}


MAYBE_UNUSED
static inline int
modredc15ul_next (residueredc15ul_t r, const modulusredc15ul_t m)
{
  ularith_add_ul_2ul (&(r[0]), &(r[1]), 1UL);
  return (r[1] == m[0].m[1] && r[0] == m[0].m[0]);
}


MAYBE_UNUSED
static inline int
modredc15ul_finished (const residueredc15ul_t r, const modulusredc15ul_t m)
{
  return (r[1] == m[0].m[1] && r[0] == m[0].m[0]);
}

/* prototypes of non-inline functions */
void modredc15ul_div3 (residueredc15ul_t, const residueredc15ul_t, 
		       const modulusredc15ul_t);
void modredc15ul_div7 (residueredc15ul_t, const residueredc15ul_t, 
		       const modulusredc15ul_t);
void modredc15ul_gcd (unsigned long *, const residueredc15ul_t, 
		      const modulusredc15ul_t);
int modredc15ul_inv (residueredc15ul_t, const residueredc15ul_t, 
		     const modulusredc15ul_t);
void modredc15ul_pow_ul (residueredc15ul_t, const residueredc15ul_t, 
			 const unsigned long, const modulusredc15ul_t);
void modredc15ul_pow_mp (residueredc15ul_t, const residueredc15ul_t, 
			 const unsigned long *, const int, 
			 const modulusredc15ul_t);
int modredc15ul_sprp (const residueredc15ul_t, const modulusredc15ul_t);
void modredc15ul_2pow_mp (residueredc15ul_t, const residueredc15ul_t, 
			  const unsigned long *, const int, 
			  const unsigned long, const modulusredc15ul_t);
void modredc15ul_V_ul (residueredc15ul_t, const residueredc15ul_t, 
		       const residueredc15ul_t, const unsigned long, 
		       const modulusredc15ul_t);
int modredc15ul_jacobi (const residueredc15ul_t, const modulusredc15ul_t);
#endif
