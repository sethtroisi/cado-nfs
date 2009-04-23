/* Some functions for modular arithmetic with residues and modulus in
   unsigned long variables. The modulus can be up to 2 unsigned longs 
   in size with the two most significant bits zert (meaning 62 bits if 
   unsigned long has 32 bits, or 126 bits if unsigned long has 64 bits).
   Moduli must be odd and have the upper word non-zero. Residues are stored 
   in Montgomery form, reduction after multiplication is done with REDC. 
   Due to inlining, this file must be included in the caller's source code with 
   #include */

/* Naming convention: all function start with modredc2ul2, for 
   MODulus REDC 2 Unsigned Longs minus 2 bits, followed by underscore, 
   functionality of function (add, mul, etc), and possibly underscore and 
   specification of what argument types the function takes (_ul, etc). */

#ifndef __MODREDC_2UL2_H

#define __MODREDC_2UL2_H

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

#define MODREDC2UL2_SIZE 2
#define MODREDC2UL2_MINBITS LONG_BIT
#define MODREDC2UL2_MAXBITS (2 * LONG_BIT - 2)

typedef unsigned long residueredc2ul2_t[MODREDC2UL2_SIZE];
typedef unsigned long modintredc2ul2_t[MODREDC2UL2_SIZE];
typedef struct { 
  modintredc2ul2_t m;
  residueredc2ul2_t one;
  unsigned long invm;
} __modulusredc2ul2_t;
typedef __modulusredc2ul2_t modulusredc2ul2_t[1];


/* ==================== Functions used internally ==================== */

static inline int
modredc2ul2_intlt (const modintredc2ul2_t a, const modintredc2ul2_t b);

static inline void
modredc2ul2_add (residueredc2ul2_t r, const residueredc2ul2_t a, 
		 const residueredc2ul2_t b, const modulusredc2ul2_t m);
static inline void
modredc2ul2_get_uls (modintredc2ul2_t r, const residueredc2ul2_t s, 
		     const modulusredc2ul2_t m MAYBE_UNUSED);

MAYBE_UNUSED
static inline void
modredc2ul2_tomontgomery (residueredc2ul2_t r, const residueredc2ul2_t s,
			  const modulusredc2ul2_t m)
{
  int i;

  r[0] = s[0];
  r[1] = s[1];
  /* TODO FIXME: ridiculously slow */
  for (i = 0; i < 2 * LONG_BIT; i++)
    modredc2ul2_add (r, r, r, m);
}


/* Do a one-word REDC, i.e., r == s / 2^LONG_BIT (mod m), r < m. 
   If m > 2^w, r < 2m. If s<m, then r<m */
MAYBE_UNUSED
static inline void
modredc2ul2_redc1 (residueredc2ul2_t r, const residueredc2ul2_t s,
		   const modulusredc2ul2_t m)
{
  unsigned long t[4], k;
  
  k = s[0] * m[0].invm;
  ularith_mul_ul_ul_2ul (&(t[0]), &(t[1]), k, m[0].m[0]);
  if (s[0] != 0UL)
    t[1]++;
  t[2] = 0;
  ularith_add_ul_2ul (&(t[1]), &(t[2]), s[1]); /* t[2] <= 1 */
  ularith_mul_ul_ul_2ul (&(t[0]), &(t[3]), k, m[0].m[1]); /* t[3] < 2^w-1 */
  ularith_add_2ul_2ul (&(t[1]), &(t[2]), t[0], t[3]);     /* t[2] < 2^w */

  /* r = (k*m + s) / wb, k <= wb-1. If s < m, then r < m */
  r[0] = t[1];
  r[1] = t[2];
}

/* Converts s out of Montgomery form by dividing by 2^(2*LONG_BIT).
   Requires s < m. */
MAYBE_UNUSED
static inline void
modredc2ul2_frommontgomery (residueredc2ul2_t r, const residueredc2ul2_t s,
			    const modulusredc2ul2_t m)
{
  unsigned long t[2];
  
  /* Do two REDC steps */
  modredc2ul2_redc1 (t, s, m);
  modredc2ul2_redc1 (r, t, m);
}


/* ================= Functions that are part of the API ================= */

/* Some functions for integers of the same width as the modulus */

MAYBE_UNUSED
static inline void
modredc2ul2_intset (modintredc2ul2_t r, const modintredc2ul2_t s)
{
  r[0] = s[0];
  r[1] = s[1];
}

MAYBE_UNUSED
static inline void
modredc2ul2_intset_ul (modintredc2ul2_t r, const unsigned long s)
{
  r[0] = s;
  r[1] = 0UL;
}

MAYBE_UNUSED
static inline int
modredc2ul2_intequal (const modintredc2ul2_t a, const modintredc2ul2_t b)
{
  return (a[0] == b[0] && a[1] == b[1]);
}

MAYBE_UNUSED
static inline int
modredc2ul2_intequal_ul (const modintredc2ul2_t a, const unsigned long b)
{
  return (a[0] == b && a[1] == 0UL);
}

/* Returns 1 if a < b, 0 otherwise */
MAYBE_UNUSED
static inline int
modredc2ul2_intlt (const modintredc2ul2_t a, const modintredc2ul2_t b)
{
    modintredc2ul2_t t;

    modredc2ul2_intset (t, a);
    return ularith_sub_2ul_2ul_cy (&(t[0]), &(t[1]), b[0], b[1]);
}

MAYBE_UNUSED
static inline int
modredc2ul2_intcmp (const modintredc2ul2_t a, const modintredc2ul2_t b)
{
  if (a[1] < b[1])
    return -1;
  if (a[1] > b[1])
    return 1;
  return (a[0] < b[0]) ? -1 : (a[0] == b[0]) ? 0 : 1;
}

MAYBE_UNUSED
static inline int
modredc2ul2_intcmp_ul (const modintredc2ul2_t a, const unsigned long b)
{
  if (a[1] > 0UL)
    return 1;
  return (a[0] < b) ? -1 : (a[0] == b) ? 0 : 1;
}

MAYBE_UNUSED
static inline int
modredc2ul2_intfits_ul (const modintredc2ul2_t a)
{
  return (a[1] == 0UL);
}

MAYBE_UNUSED
static inline void
modredc2ul2_intadd (modintredc2ul2_t r, const modintredc2ul2_t a,
		    const modintredc2ul2_t b)
{
  modintredc2ul2_t t;
  modredc2ul2_intset (t, a);
  ularith_add_2ul_2ul (&(t[0]), &(t[1]), b[0], b[1]);
  modredc2ul2_intset (r, t);
}

MAYBE_UNUSED
static inline void
modredc2ul2_intsub (modintredc2ul2_t r, const modintredc2ul2_t a,
		    const modintredc2ul2_t b)
{
  modintredc2ul2_t t;
  modredc2ul2_intset (t, a);
  ularith_sub_2ul_2ul (&(t[0]), &(t[1]), b[0], b[1]);
  modredc2ul2_intset (r, t);
}

/* Returns the number of bits in a, that is, floor(log_2(n))+1. 
   For n==0 returns 0. */
MAYBE_UNUSED
static inline int
modredc2ul2_intbits (const modintredc2ul2_t a)
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


MAYBE_UNUSED
static inline void
modredc2ul2_intshr (modintredc2ul2_t r, const modintredc2ul2_t s, const int i)
{
  r[0] = s[0];
  ularith_shrd (&(r[0]), s[1], i);
  r[1] = s[1] >> i;
}


MAYBE_UNUSED
static inline void
modredc2ul2_intshl (modintredc2ul2_t r, const modintredc2ul2_t s, const int i)
{
  r[1] = s[1];
  ularith_shld (&(r[1]), s[0], i);
  r[0] = s[0] << i;
}


/* r = n/d. We require d|n */
MAYBE_UNUSED
static inline void
modredc2ul2_intdivexact (modintredc2ul2_t r, const modintredc2ul2_t n,
                         const modintredc2ul2_t d)
{
  modintredc2ul2_t n1, d1;
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
  
  invf = ularith_invmod (d1[0]);
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

/* Init the modulus from a multi-word integer. s must point to an array of
   at least two unsigned longs, where s[0] is the low word of the modulus, 
   and s[1] is the high word. */
MAYBE_UNUSED
static inline void
modredc2ul2_initmod_uls (modulusredc2ul2_t m, const modintredc2ul2_t s)
{
  int i;
  ASSERT (s[1] > 0UL);
  ASSERT (s[1] < (1UL << (LONG_BIT - 2)));
  modredc2ul2_intset (m[0].m, s);
  m[0].invm = -ularith_invmod (s[0]);
  m[0].one[0] = 0UL;
  m[0].one[1] = 1UL;
  for (i = 0; i < LONG_BIT; i++)
    modredc2ul2_add (m[0].one, m[0].one, m[0].one, m);
#ifdef WANT_ASSERT_EXPENSIVE
  {
    modintredc2ul2_t t;
    modredc2ul2_get_uls (t, m[0].one, m);
    ASSERT_EXPENSIVE (t[0] == 1UL && t[1] == 0UL);
  }
#endif
}


/* Returns the modulus to an array of unsigned longs. */
MAYBE_UNUSED
static inline void
modredc2ul2_getmod_uls (modintredc2ul2_t r, const modulusredc2ul2_t m)
{
  r[0] = m[0].m[0];
  r[1] = m[0].m[1];
}


MAYBE_UNUSED
static inline void
modredc2ul2_clearmod (modulusredc2ul2_t m MAYBE_UNUSED)
{
  return;
}


/* Functions for residues */

/* Initialises a residueredc2ul2_t and sets it to zero */
MAYBE_UNUSED
static inline void
modredc2ul2_init (residueredc2ul2_t r, const modulusredc2ul2_t m MAYBE_UNUSED)
{
  r[0] = 0UL;
  r[1] = 0UL;
}


/* Initialises a residueredc2ul2_t, but does not set it to zero. For fixed 
   length residueredc2ul2_t, that leaves nothing to do at all. */
MAYBE_UNUSED
static inline void
modredc2ul2_init_noset0 (residueredc2ul2_t r MAYBE_UNUSED, 
			 const modulusredc2ul2_t m MAYBE_UNUSED)
{
  return;
}


MAYBE_UNUSED
static inline void
modredc2ul2_clear (residueredc2ul2_t r MAYBE_UNUSED, 
		   const modulusredc2ul2_t m MAYBE_UNUSED)
{
  return;
}


MAYBE_UNUSED
static inline void
modredc2ul2_set (residueredc2ul2_t r, const residueredc2ul2_t s, 
		 const modulusredc2ul2_t m MAYBE_UNUSED)
{
  ASSERT_EXPENSIVE (modredc2ul2_intcmp (s, m[0].m) < 0);
  r[0] = s[0];
  r[1] = s[1];
}


MAYBE_UNUSED
static inline void
modredc2ul2_set_ul (residueredc2ul2_t r, const unsigned long s, 
		    const modulusredc2ul2_t m)
{
  r[0] = s;
  r[1] = 0UL;
  modredc2ul2_tomontgomery (r, r, m);
#ifdef WANT_ASSERT_EXPENSIVE
  {
    modintredc2ul2_t t;
    modredc2ul2_get_uls (t, r, m);
    ASSERT_EXPENSIVE (t[0] == s && t[1] == 0UL);
  }
#endif
}


/* Sets the residueredc2ul2_t to the class represented by the integer s. 
   Assumes that s is reduced (mod m), i.e. 0 <= s < m */

MAYBE_UNUSED
static inline void
modredc2ul2_set_ul_reduced (residueredc2ul2_t r, const unsigned long s, 
			    const modulusredc2ul2_t m MAYBE_UNUSED)
{
  modredc2ul2_set_ul (r, s, m);
}


MAYBE_UNUSED
static inline void
modredc2ul2_set_uls (residueredc2ul2_t r, const modintredc2ul2_t s, 
		     const modulusredc2ul2_t m)
{
  r[0] = s[0];
  r[1] = s[1];
  if (modredc2ul2_intcmp (r, m[0].m) >= 0)
    {
      /* Do reduction */
      /* FIXME, slow and stupid */
      modintredc2ul2_t t;
      modredc2ul2_intset (t, m[0].m);
      while ((t[1] & (1UL << (LONG_BIT - 1))) == 0UL &&
             modredc2ul2_intlt (t, r))
	modredc2ul2_intshl (t, t, 1);
      while (!modredc2ul2_intlt (r, m[0].m))
	{
	  ularith_sub_2ul_2ul_ge (&(r[0]), &(r[1]), t[0], t[1]);
	  modredc2ul2_intshr (t, t, 1);
	}
    }
  ASSERT (modredc2ul2_intcmp (r, m[0].m) < 0);
  modredc2ul2_tomontgomery (r, r, m);
}


MAYBE_UNUSED
static inline void
modredc2ul2_set_uls_reduced (residueredc2ul2_t r, const modintredc2ul2_t s, 
			     const modulusredc2ul2_t m)
{
  ASSERT (modredc2ul2_intcmp (s, m[0].m) < 0);
  r[0] = s[0];
  r[1] = s[1];
  modredc2ul2_tomontgomery (r, r, m);
}


MAYBE_UNUSED 
static inline void 
modredc2ul2_set0 (residueredc2ul2_t r, const modulusredc2ul2_t m MAYBE_UNUSED) 
{ 
  r[0] = 0UL; 
  r[1] = 0UL;
}


MAYBE_UNUSED 
static inline void 
modredc2ul2_set1 (residueredc2ul2_t r, const modulusredc2ul2_t m) 
{ 
  r[0] = m[0].one[0];
  r[1] = m[0].one[1];
}


/* Exchanges the values of the two arguments */

MAYBE_UNUSED
static inline void
modredc2ul2_swap (residueredc2ul2_t a, residueredc2ul2_t b, 
		  const modulusredc2ul2_t m MAYBE_UNUSED)
{
  unsigned long t0, t1;
  ASSERT_EXPENSIVE (modredc2ul2_intcmp (a, m[0].m) < 0);
  ASSERT_EXPENSIVE (modredc2ul2_intcmp (b, m[0].m) < 0);
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
modredc2ul2_get_ul (const residueredc2ul2_t s, 
		    const modulusredc2ul2_t m MAYBE_UNUSED)
{
  unsigned long t[2];
  ASSERT_EXPENSIVE (modredc2ul2_intcmp (s, m[0].m) < 0);
  modredc2ul2_frommontgomery (t, s, m);
  ASSERT (t[1] == 0UL);
  return t[0];
}


/* Returns the residue into an array of unsigned longs */

MAYBE_UNUSED
static inline void
modredc2ul2_get_uls (modintredc2ul2_t r, const residueredc2ul2_t s, 
		     const modulusredc2ul2_t m MAYBE_UNUSED)
{
  ASSERT_EXPENSIVE (modredc2ul2_intcmp (s, m[0].m) < 0);
  modredc2ul2_frommontgomery (r, s, m);
}


MAYBE_UNUSED
static inline int
modredc2ul2_equal (const residueredc2ul2_t a, const residueredc2ul2_t b, 
		   const modulusredc2ul2_t m MAYBE_UNUSED)
{
  ASSERT_EXPENSIVE (modredc2ul2_intcmp (a, m[0].m) < 0);
  ASSERT_EXPENSIVE (modredc2ul2_intcmp (b, m[0].m) < 0);
  return (a[1] == b[1] && a[0] == b[0]);
}


MAYBE_UNUSED
static inline int
modredc2ul2_is0 (const residueredc2ul2_t a, const modulusredc2ul2_t m MAYBE_UNUSED)
{
  ASSERT_EXPENSIVE (modredc2ul2_intcmp (a, m[0].m) < 0);
  return (a[1] == 0UL && a[0] == 0UL);
}


MAYBE_UNUSED
static inline int
modredc2ul2_is1 (const residueredc2ul2_t a, const modulusredc2ul2_t m MAYBE_UNUSED)
{
  ASSERT_EXPENSIVE (modredc2ul2_intcmp (a, m[0].m) < 0);
  return (a[0] == m[0].one[0] && a[1] == m[0].one[1]);
}


MAYBE_UNUSED
static inline void
modredc2ul2_add (residueredc2ul2_t r, const residueredc2ul2_t a, 
		 const residueredc2ul2_t b, const modulusredc2ul2_t m)
{
  const unsigned long t0 = b[0], t1 = b[1]; /* r, a, and/or b may overlap */
  ASSERT_EXPENSIVE (modredc2ul2_intcmp (a, m[0].m) < 0);
  ASSERT_EXPENSIVE (modredc2ul2_intcmp (b, m[0].m) < 0);

  r[0] = a[0];
  r[1] = a[1];
  ularith_add_2ul_2ul (&(r[0]), &(r[1]), t0, t1);
  ularith_sub_2ul_2ul_ge (&(r[0]), &(r[1]), m[0].m[0], m[0].m[1]);
  ASSERT_EXPENSIVE (modredc2ul2_intcmp (r, m[0].m) < 0);
}


MAYBE_UNUSED
static inline void
modredc2ul2_sub (residueredc2ul2_t r, const residueredc2ul2_t a, 
		 const residueredc2ul2_t b, const modulusredc2ul2_t m)
{
  ASSERT_EXPENSIVE (modredc2ul2_intcmp (a, m[0].m) < 0);
  ASSERT_EXPENSIVE (modredc2ul2_intcmp (b, m[0].m) < 0);

#if defined(__x86_64__) && defined(__GNUC__)
  {
    unsigned long s1 = m[0].m[0], s2 = m[0].m[1], t1 = a[0], t2 = a[1];
    
    __asm__ (
	     "subq %4, %0\n\t"
	     "sbbq %5, %1\n\t"    /* r -= b */
	     "cmovncq %6, %2\n\t" /* If !carry, s = 0 */
	     "cmovncq %6, %3\n"
	     : "+&r" (t1), "+&r" (t2), "+&r" (s1), "+r" (s2)
	     : "g" (b[0]), "g" (b[1]), "rm" (0UL)
	     : "cc"
	     );
    ularith_add_2ul_2ul (&t1, &t2, s1, s2);
    r[0] = t1;
    r[1] = t2;
  }
#else
  {
    unsigned long t1 = a[0], t2 = a[1];
    ularith_sub_2ul_2ul (&t1, &t2, b[0], b[1]);
    
    if (t2 > a[1] || (t2 == a[1] && t1 > a[0]))
      ularith_add_2ul_2ul (&t1, &t2, m[0].m[0], m[0].m[1]);

    r[0] = t1;
    r[1] = t2;
  }
#endif
}


MAYBE_UNUSED
static inline void
modredc2ul2_add_ul (residueredc2ul2_t r, const residueredc2ul2_t a,
		    const unsigned long b, const modulusredc2ul2_t m)
{
  residueredc2ul2_t t;

  /* TODO: speed up */
  ASSERT_EXPENSIVE (modredc2ul2_intcmp (a, m[0].m) < 0);
  modredc2ul2_init_noset0 (t, m);
  modredc2ul2_set_ul (t, b, m);
  modredc2ul2_add (r, a, t, m);
  modredc2ul2_clear (t, m);
}


MAYBE_UNUSED
static inline void
modredc2ul2_sub_ul (residueredc2ul2_t r, const residueredc2ul2_t a,
		    const unsigned long b, const modulusredc2ul2_t m)
{
  residueredc2ul2_t t;

  /* TODO: speed up */
  ASSERT_EXPENSIVE (modredc2ul2_intcmp (a, m[0].m) < 0);
  modredc2ul2_init_noset0 (t, m);
  modredc2ul2_set_ul (t, b, m);
  modredc2ul2_sub (r, a, t, m);
  modredc2ul2_clear (t, m);
}


MAYBE_UNUSED
static inline void
modredc2ul2_neg (residueredc2ul2_t r, const residueredc2ul2_t a, 
		 const modulusredc2ul2_t m)
{
  ASSERT_EXPENSIVE (modredc2ul2_intcmp (a, m[0].m) < 0);
  if (a[0] == 0UL && a[1] == 0UL)
    modredc2ul2_set (r, a, m);
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
modredc2ul2_div2 (residueredc2ul2_t r, const residueredc2ul2_t a, 
		  const modulusredc2ul2_t m)
{
  ASSERT_EXPENSIVE (m[0].m[0] % 2UL != 0UL);
  ASSERT_EXPENSIVE (modredc2ul2_intcmp (a, m[0].m) < 0);
  r[0] = a[0];
  r[1] = a[1];
  if (r[0] % 2UL == 1UL)
    ularith_add_2ul_2ul (&(r[0]), &(r[1]), m[0].m[0], m[0].m[1]);
  ularith_shrd (&(r[0]), r[1], 1);
  r[1] >>= 1;
}


MAYBE_UNUSED
static inline void
modredc2ul2_mul (residueredc2ul2_t r, const residueredc2ul2_t a, 
               const residueredc2ul2_t b, const modulusredc2ul2_t m)
{
  unsigned long pl, ph, t[4], k;
  
  ASSERT_EXPENSIVE (modredc2ul2_intcmp (a, m[0].m) < 0);
  ASSERT_EXPENSIVE (modredc2ul2_intcmp (b, m[0].m) < 0);
#if defined(MODTRACE)
  printf ("((%lu * 2^%d + %lu) * (%lu * 2^%d + %lu) / 2^%d) %% "
	  "(%lu * 2^%d + %lu)", 
          a[1], LONG_BIT, a[0], b[1], LONG_BIT, b[0], 2 * LONG_BIT, 
	  m[0].m[1], LONG_BIT, m[0].m[0]);
#endif

  /* m < 1/4 W^2,  a,b < m */
  
  /* Product of the two low words */
  ularith_mul_ul_ul_2ul (&(t[0]), &(t[1]), a[0], b[0]); 

  /* One REDC step */
  modredc2ul2_redc1 (t, t, m); /* t < 2m < 1/2 W^2 */

  /* Products of one low and one high word  */
  ularith_mul_ul_ul_2ul (&pl, &ph, a[1], b[0]);   /* ph:pl < 1/4 W^2 */
  ularith_add_2ul_2ul (&(t[0]), &(t[1]), pl, ph); /* t1:t0 < 3/4 W^2 */
  ularith_mul_ul_ul_2ul (&pl, &ph, a[0], b[1]);   /* ph:pl < 1/4 W^2 */
  ularith_add_2ul_2ul (&(t[0]), &(t[1]), pl, ph); /* t1:t0 < W^2 */

  /* Product of the two high words */
  ularith_mul_ul_ul_2ul (&pl, &(t[2]), a[1], b[1]); /* t2:pl < 1/16 W^2 */
  ularith_add_ul_2ul (&(t[1]), &(t[2]), pl); /* t2:t1:t0 < 1/16 W^3 + W^2 */

  /* Compute t2:t1:t0 := t2:t1:t0 + km, km < Wm < 1/4 W^3 */
  k = t[0] * m[0].invm;
  ularith_mul_ul_ul_2ul (&pl, &ph, k, m[0].m[0]);
  if (t[0] != 0UL)
    ph++; /* t[0] = 0 */
  ularith_add_ul_2ul (&(t[1]), &(t[2]), ph);
  ularith_mul_ul_ul_2ul (&pl, &ph, k, m[0].m[1]); /* ph:pl < 1/4 W^2 */
  ularith_add_2ul_2ul (&(t[1]), &(t[2]), pl, ph);
  /* t2:t1:0 < 1/16 W^3 + W^2 + 1/4 W^3 < 5/16 W^3 + W^2 */

  /* Result may be larger than m, but is < 2*m */

  ularith_sub_2ul_2ul_ge (&(t[1]), &(t[2]), m[0].m[0], m[0].m[1]);

  r[0] = t[1];
  r[1] = t[2];
#if defined(MODTRACE)
  printf (" == (%lu * 2^%d + %lu) /* PARI */ \n", r[1], LONG_BIT, r[0]);
#endif
  ASSERT_EXPENSIVE (modredc2ul2_intlt (r, m[0].m));
}


MAYBE_UNUSED
static inline int
modredc2ul2_next (residueredc2ul2_t r, const modulusredc2ul2_t m)
{
  ularith_add_ul_2ul (&(r[0]), &(r[1]), 1UL);
  return (r[1] == m[0].m[1] && r[0] == m[0].m[0]);
}


MAYBE_UNUSED
static inline int
modredc2ul2_finished (const residueredc2ul2_t r, const modulusredc2ul2_t m)
{
  return (r[1] == m[0].m[1] && r[0] == m[0].m[0]);
}

/* prototypes of non-inline functions */
void modredc2ul2_div3 (residueredc2ul2_t, const residueredc2ul2_t, 
		       const modulusredc2ul2_t);
void modredc2ul2_div5 (residueredc2ul2_t, const residueredc2ul2_t, 
		       const modulusredc2ul2_t);
void modredc2ul2_div7 (residueredc2ul2_t, const residueredc2ul2_t, 
		       const modulusredc2ul2_t);
void modredc2ul2_div13 (residueredc2ul2_t, const residueredc2ul2_t, 
		        const modulusredc2ul2_t);
void modredc2ul2_gcd (modintredc2ul2_t, const residueredc2ul2_t, 
		      const modulusredc2ul2_t);
void modredc2ul2_pow_ul (residueredc2ul2_t, const residueredc2ul2_t, 
			 const unsigned long, const modulusredc2ul2_t);
void modredc2ul2_2pow_ul (residueredc2ul2_t, const unsigned long, 
                          const modulusredc2ul2_t);
void modredc2ul2_pow_mp (residueredc2ul2_t, const residueredc2ul2_t, 
			 const unsigned long *, const int, 
			 const modulusredc2ul2_t);
void modredc2ul2_2pow_mp (residueredc2ul2_t, const unsigned long *, const int, 
			  const modulusredc2ul2_t);
void modredc2ul2_V_ul (residueredc2ul2_t, const residueredc2ul2_t, 
		       const unsigned long, const modulusredc2ul2_t);
void modredc2ul2_V_mp (residueredc2ul2_t, const residueredc2ul2_t, 
		       const unsigned long *, const int, 
		       const modulusredc2ul2_t);
int modredc2ul2_sprp (const residueredc2ul2_t, const modulusredc2ul2_t);
int modredc2ul2_sprp2 (const modulusredc2ul2_t);
int modredc2ul2_isprime (const modulusredc2ul2_t);
int modredc2ul2_inv (residueredc2ul2_t, const residueredc2ul2_t, 
		     const modulusredc2ul2_t);
int modredc2ul2_jacobi (const residueredc2ul2_t, const modulusredc2ul2_t);
#endif
