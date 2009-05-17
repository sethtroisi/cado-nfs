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

/* A macro for function renaming. All functions here start with 
   modredc15ul_ */
#define MODREDC15UL_RENAME(x) modredc15ul_##x

#define MODREDC15UL_SIZE 2
#define MODREDC15UL_MINBITS LONG_BIT
#define MODREDC15UL_MAXBITS (LONG_BIT + LONG_BIT/2)

typedef unsigned long residueredc15ul_t[MODREDC15UL_SIZE];
typedef unsigned long modintredc15ul_t[MODREDC15UL_SIZE];
typedef struct { 
  modintredc15ul_t m;
  residueredc15ul_t one;
  unsigned long invm;
} __modulusredc15ul_t;
typedef __modulusredc15ul_t modulusredc15ul_t[1];


/* ==================== Functions used internally ==================== */

static inline void
modredc15ul_add (residueredc15ul_t r, const residueredc15ul_t a, 
		 const residueredc15ul_t b, const modulusredc15ul_t m);
static inline void
modredc15ul_get_uls (modintredc15ul_t r, const residueredc15ul_t s, 
		     const modulusredc15ul_t m MAYBE_UNUSED);

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
    modredc15ul_add (r, r, r, m);
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

/* Do a one-word REDC, i.e., divide by 2^LONG_BIT */
MAYBE_UNUSED
static inline void
modredc15ul_redc1 (residueredc15ul_t r, const residueredc15ul_t s,
		   const modulusredc15ul_t m)
{
  unsigned long t[4], k;
  
  k = s[0] * m[0].invm;
  ularith_mul_ul_ul_2ul (&(t[0]), &(t[1]), k, m[0].m[0]);
  if (t[0] != 0UL)
    t[1]++;
  t[2] = 0;
  ularith_add_ul_2ul (&(t[1]), &(t[2]), s[1]); /* t[2] <= 1 */
  ularith_mul_ul_ul_2ul (&(t[0]), &(t[3]), k, m[0].m[1]); /* t[3] < 2^(w/2) */
  ularith_add_2ul_2ul (&(t[1]), &(t[2]), t[0], t[3]);

  /* r = (k*m + s) div wb. k <= wb-1, s<m, thus r < m */
  r[0] = t[1];
  r[1] = t[2];
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

/* Returns 1 if a < b, 0 otherwise */
MAYBE_UNUSED
static inline int
modredc15ul_intlt (const modintredc15ul_t a, const modintredc15ul_t b)
{
    modintredc15ul_t t;

    modredc15ul_intset (t, a);
    return ularith_sub_2ul_2ul_cy (&(t[0]), &(t[1]), b[0], b[1]);
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

MAYBE_UNUSED
static inline void
modredc15ul_intadd (modintredc15ul_t r, const modintredc15ul_t a,
		    const modintredc15ul_t b)
{
  modintredc15ul_t t;
  modredc15ul_intset (t, a);
  ularith_add_2ul_2ul (&(t[0]), &(t[1]), b[0], b[1]);
  modredc15ul_intset (r, t);
}

MAYBE_UNUSED
static inline void
modredc15ul_intsub (modintredc15ul_t r, const modintredc15ul_t a,
		    const modintredc15ul_t b)
{
  modintredc15ul_t t;
  modredc15ul_intset (t, a);
  ularith_sub_2ul_2ul (&(t[0]), &(t[1]), b[0], b[1]);
  modredc15ul_intset (r, t);
}

/* Returns the number of bits in a, that is, floor(log_2(a))+1. 
   For a == 0 returns 0. */
MAYBE_UNUSED
static inline int
modredc15ul_intbits (const modintredc15ul_t a)
{
  if (a[1] > 0UL)
    return 2*LONG_BIT - ularith_clz (a[1]);

  if (a[0] > 0UL)
    return LONG_BIT - ularith_clz (a[0]);
  
  return 0;
}


MAYBE_UNUSED
static inline void
modredc15ul_intshr (modintredc15ul_t r, const modintredc15ul_t s, const int i)
{
  r[0] = s[0];
  ularith_shrd (&(r[0]), s[1], i);
  r[1] = s[1] >> i;
}


MAYBE_UNUSED
static inline void
modredc15ul_intshl (modintredc15ul_t r, const modintredc15ul_t s, const int i)
{
  r[1] = s[1];
  ularith_shld (&(r[1]), s[0], i);
  r[0] = s[0] << i;
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
modredc15ul_initmod_uls (modulusredc15ul_t m, const modintredc15ul_t s)
{
  int i;
  ASSERT (s[1] > 0UL);
  ASSERT (s[1] < (1UL << (LONG_BIT / 2)));
  modredc15ul_intset (m[0].m, s);
  m[0].invm = -ularith_invmod (s[0]);
  m[0].one[0] = 0UL;
  m[0].one[1] = 1UL;
  for (i = 0; i < LONG_BIT; i++)
    modredc15ul_add (m[0].one, m[0].one, m[0].one, m);
#ifdef WANT_ASSERT_EXPENSIVE
  {
    modintredc15ul_t t;
    modredc15ul_get_uls (t, m[0].one, m);
    ASSERT_EXPENSIVE (t[0] == 1UL && t[1] == 0UL);
  }
#endif
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
  ASSERT_EXPENSIVE (modredc15ul_intcmp (s, m[0].m) < 0);
  r[0] = s[0];
  r[1] = s[1];
}


MAYBE_UNUSED
static inline void
modredc15ul_set_ul (residueredc15ul_t r, const unsigned long s, 
		    const modulusredc15ul_t m)
{
  r[0] = s;
  r[1] = 0UL;
  modredc15ul_tomontgomery (r, r, m);
#ifdef WANT_ASSERT_EXPENSIVE
  {
    modintredc15ul_t t;
    modredc15ul_get_uls (t, r, m);
    ASSERT_EXPENSIVE (t[0] == s && t[1] == 0UL);
  }
#endif
}


/* Sets the residueredc15ul_t to the class represented by the integer s. 
   Assumes that s is reduced (mod m), i.e. 0 <= s < m */

MAYBE_UNUSED
static inline void
modredc15ul_set_ul_reduced (residueredc15ul_t r, const unsigned long s, 
			    const modulusredc15ul_t m MAYBE_UNUSED)
{
  modredc15ul_set_ul (r, s, m);
}


MAYBE_UNUSED
static inline void
modredc15ul_set_uls (residueredc15ul_t r, const modintredc15ul_t s, 
		     const modulusredc15ul_t m)
{
  r[0] = s[0];
  r[1] = s[1];
  if (modredc15ul_intcmp (r, m[0].m) >= 0)
    {
      /* Do reduction */
      /* FIXME, slow and stupid */
      modintredc15ul_t t;
      modredc15ul_intset (t, m[0].m);
      while ((t[1] & (1UL << (LONG_BIT - 1))) == 0UL &&
             modredc15ul_intlt (t, r))
	modredc15ul_intshl (t, t, 1);
      while (!modredc15ul_intlt (r, m[0].m))
	{
	  ularith_sub_2ul_2ul_ge (&(r[0]), &(r[1]), t[0], t[1]);
	  modredc15ul_intshr (t, t, 1);
	}
    }
  ASSERT (modredc15ul_intcmp (r, m[0].m) < 0);
  modredc15ul_tomontgomery (r, r, m);
}


MAYBE_UNUSED
static inline void
modredc15ul_set_uls_reduced (residueredc15ul_t r, const modintredc15ul_t s, 
			     const modulusredc15ul_t m)
{
  ASSERT (modredc15ul_intcmp (s, m[0].m) < 0);
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
  r[0] = m[0].one[0];
  r[1] = m[0].one[1];
}


/* Exchanges the values of the two arguments */

MAYBE_UNUSED
static inline void
modredc15ul_swap (residueredc15ul_t a, residueredc15ul_t b, 
		  const modulusredc15ul_t m MAYBE_UNUSED)
{
  unsigned long t0, t1;
  ASSERT_EXPENSIVE (modredc15ul_intcmp (a, m[0].m) < 0);
  ASSERT_EXPENSIVE (modredc15ul_intcmp (b, m[0].m) < 0);
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
  ASSERT_EXPENSIVE (modredc15ul_intcmp (s, m[0].m) < 0);
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
  ASSERT_EXPENSIVE (modredc15ul_intcmp (s, m[0].m) < 0);
  modredc15ul_frommontgomery (r, s, m);
}


MAYBE_UNUSED
static inline int
modredc15ul_equal (const residueredc15ul_t a, const residueredc15ul_t b, 
		   const modulusredc15ul_t m MAYBE_UNUSED)
{
  ASSERT_EXPENSIVE (modredc15ul_intcmp (a, m[0].m) < 0);
  ASSERT_EXPENSIVE (modredc15ul_intcmp (b, m[0].m) < 0);
  return (a[1] == b[1] && a[0] == b[0]);
}


MAYBE_UNUSED
static inline int
modredc15ul_is0 (const residueredc15ul_t a, const modulusredc15ul_t m MAYBE_UNUSED)
{
  ASSERT_EXPENSIVE (modredc15ul_intcmp (a, m[0].m) < 0);
  return (a[1] == 0UL && a[0] == 0UL);
}


MAYBE_UNUSED
static inline int
modredc15ul_is1 (const residueredc15ul_t a, const modulusredc15ul_t m MAYBE_UNUSED)
{
  ASSERT_EXPENSIVE (modredc15ul_intcmp (a, m[0].m) < 0);
  return (a[0] == m[0].one[0] && a[1] == m[0].one[1]);
}


MAYBE_UNUSED
static inline void
modredc15ul_add (residueredc15ul_t r, const residueredc15ul_t a, 
		 const residueredc15ul_t b, const modulusredc15ul_t m)
{
  const unsigned long t0 = b[0], t1 = b[1]; /* r, a, and/or b may overlap */
  ASSERT_EXPENSIVE (modredc15ul_intcmp (a, m[0].m) < 0);
  ASSERT_EXPENSIVE (modredc15ul_intcmp (b, m[0].m) < 0);

  r[0] = a[0];
  r[1] = a[1];
  ularith_add_2ul_2ul (&(r[0]), &(r[1]), t0, t1);
  ularith_sub_2ul_2ul_ge (&(r[0]), &(r[1]), m[0].m[0], m[0].m[1]);
  ASSERT_EXPENSIVE (modredc15ul_intcmp (r, m[0].m) < 0);
}


MAYBE_UNUSED
static inline void
modredc15ul_sub (residueredc15ul_t r, const residueredc15ul_t a, 
		 const residueredc15ul_t b, const modulusredc15ul_t m)
{
  ASSERT_EXPENSIVE (modredc15ul_intcmp (a, m[0].m) < 0);
  ASSERT_EXPENSIVE (modredc15ul_intcmp (b, m[0].m) < 0);

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
modredc15ul_add_ul (residueredc15ul_t r, const residueredc15ul_t a,
		    const unsigned long b, const modulusredc15ul_t m)
{
  residueredc15ul_t t;

  /* TODO: speed up */
  ASSERT_EXPENSIVE (modredc15ul_intcmp (a, m[0].m) < 0);
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
  ASSERT_EXPENSIVE (modredc15ul_intcmp (a, m[0].m) < 0);
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
  ASSERT_EXPENSIVE (modredc15ul_intcmp (a, m[0].m) < 0);
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
  ASSERT_EXPENSIVE (modredc15ul_intcmp (a, m[0].m) < 0);
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
  
  ASSERT_EXPENSIVE (modredc15ul_intcmp (a, m[0].m) < 0);
  ASSERT_EXPENSIVE (modredc15ul_intcmp (b, m[0].m) < 0);
#if defined(MODTRACE)
  printf ("((%lu * 2^%d + %lu) * (%lu * 2^%d + %lu) / 2^%d) %% "
	  "(%lu * 2^%d + %lu)", 
          a[1], LONG_BIT, a[0], b[1], LONG_BIT, b[0], 2 * LONG_BIT, 
	  m[0].m[1], LONG_BIT, m[0].m[0]);
#endif
  
  /* Product of the two low words */
  ularith_mul_ul_ul_2ul (&(t[0]), &(t[1]), a[0], b[0]); /* t1:t0 = a[0]*b[0] <= W^2 - 2W + 1 */

  /* One REDC step */
  k = t[0] * m[0].invm; /* k <= W-1 */
  ularith_mul_ul_ul_2ul (&pl, &ph, k, m[0].m[0]); /* ph:pl = k*m[0] <= W^2 - 2W + 1 */
  /* t[0] + pl == 0 (mod W) */
  if (pl != 0UL)
    ph++; /* ph <= W-1 */
  t[2] = 0UL;
  ularith_add_ul_2ul (&(t[1]), &(t[2]), ph); /* t2:t1:0 = a[0]*b[0] + k*m[0] <= 2*W^2 - 4W + 2, so
                                                t2:t1 = (a[0]*b[0] + k*m[0]) / W <= 2*W - 4 */
  ularith_mul_ul_ul_2ul (&pl, &ph, k, m[0].m[1]); /* ph:pl <= (W^(1/2)-1)*(W-1) = W^(3/2) - W - W^(1/2) + 1 */
  ularith_add_2ul_2ul (&(t[1]), &(t[2]), pl, ph); /* t2:t1 <= W^(3/2) + W - W^(1/2) - 3 */

  /* Products of one low and one high word  */
  ularith_mul_ul_ul_2ul (&pl, &ph, a[1], b[0]);   /* ph:pl <= W^(3/2) - W - W^(1/2) + 1 */
  ularith_add_2ul_2ul (&(t[1]), &(t[2]), pl, ph); /* t2:t1 <= 2W^(3/2) - 2W^(1/2) - 2 */
  ularith_mul_ul_ul_2ul (&pl, &ph, a[0], b[1]);   /* ph:pl <= W^(3/2) - W - W^(1/2) + 1 */
  ularith_add_2ul_2ul (&(t[1]), &(t[2]), pl, ph); /* t2:t1 <= 3W^(3/2) - 3W^(1/2) - W - 1 */
  t[3] = 0UL;
  pl = a[1] * b[1];                               /* pl <= (W^(1/2)-1)^2 = W - 2W^(1/2) + 1 */
  ularith_add_ul_2ul (&(t[2]), &(t[3]), pl);      /* t3:t2:t1 <= W^2 + W^(3/2) - 3W^(1/2) - 1 */
  k = t[1] * m[0].invm;
  ularith_mul_ul_ul_2ul (&pl, &ph, k, m[0].m[0]); /* ph:pl <= W^2 - 2W + 1 */
  if (pl != 0UL)
    ph++;
  ularith_add_ul_2ul (&(t[2]), &(t[3]), ph);      /* t3:t2:t1 <= 2W^2 + W^(3/2) - 2W - 3W^(1/2), t1 = 0 so
						     t3:t2 <= 2W + W^(1/2) - 3 */
  ularith_mul_ul_ul_2ul (&pl, &ph, k, m[0].m[1]); /* ph:pl <= W^(3/2) - W - W^(1/2) + 1 */
  ularith_add_2ul_2ul (&(t[2]), &(t[3]), pl, ph); /* t3:t2 <= W^(3/2) + W - 2 */

  /* Result may be larger than m, but is < 2*m */

  ularith_sub_2ul_2ul_ge (&(t[2]), &(t[3]), m[0].m[0], m[0].m[1]);

  r[0] = t[2];
  r[1] = t[3];
#if defined(MODTRACE)
  printf (" == (%lu * 2^%d + %lu) /* PARI */ \n", r[1], LONG_BIT, r[0]);
#endif
  ASSERT_EXPENSIVE (modredc15ul_intcmp (r, m[0].m) < 0);
}


MAYBE_UNUSED
static inline void
modredc15ul_sqr (residueredc15ul_t r, const residueredc15ul_t a, 
                 const modulusredc15ul_t m)
{
  unsigned long pl, ph, t[4], k;
  
  ASSERT_EXPENSIVE (modredc15ul_intcmp (a, m[0].m) < 0);
#if defined(MODTRACE)
  printf ("((%lu * 2^%d + %lu)^2 / 2^%d) %% (%lu * 2^%d + %lu)", 
          a[1], LONG_BIT, a[0], 2 * LONG_BIT, m[0].m[1], LONG_BIT, m[0].m[0]);
#endif
  
  /* Square of low word */
  ularith_mul_ul_ul_2ul (&(t[0]), &(t[1]), a[0], a[0]); /* t1:t0 = a[0]*a[0] <= W^2 - 2W + 1 */

  /* One REDC step */
  k = t[0] * m[0].invm; /* k <= W-1 */
  ularith_mul_ul_ul_2ul (&pl, &ph, k, m[0].m[0]); /* ph:pl = k*m[0] <= W^2 - 2W + 1 */
  /* t[0] + pl == 0 (mod W) */
  if (pl != 0UL)
    ph++; /* ph <= W-1 */
  t[2] = 0UL;
  ularith_add_ul_2ul (&(t[1]), &(t[2]), ph); /* t2:t1:0 = a[0]*a[0] + k*m[0] <= 2*W^2 - 4W + 2, so
                                                t2:t1 = (a[0]*a[0] + k*m[0]) / W <= 2*W - 4 */
  ularith_mul_ul_ul_2ul (&pl, &ph, k, m[0].m[1]); /* ph:pl <= (W^(1/2)-1)*(W-1) = W^(3/2) - W - W^(1/2) + 1 */
  ularith_add_2ul_2ul (&(t[1]), &(t[2]), pl, ph); /* t2:t1 <= W^(3/2) + W - W^(1/2) - 3 */

  /* Product of low and high word  */
  ularith_mul_ul_ul_2ul (&pl, &ph, a[1], a[0]);   /* ph:pl <= W^(3/2) - W - W^(1/2) + 1 */
  ularith_add_2ul_2ul (&(t[1]), &(t[2]), pl, ph); /* t2:t1 <= 2W^(3/2) - 2W^(1/2) - 2 */
  ularith_add_2ul_2ul (&(t[1]), &(t[2]), pl, ph); /* t2:t1 <= 3W^(3/2) - 3W^(1/2) - W - 1 */
  t[3] = 0UL;

  /* Square of high word */
  pl = a[1] * a[1];                               /* pl <= (W^(1/2)-1)^2 = W - 2W^(1/2) + 1 */
  ularith_add_ul_2ul (&(t[2]), &(t[3]), pl);      /* t3:t2:t1 <= W^2 + W^(3/2) - 3W^(1/2) - 1 */
  k = t[1] * m[0].invm;
  ularith_mul_ul_ul_2ul (&pl, &ph, k, m[0].m[0]); /* ph:pl <= W^2 - 2W + 1 */
  if (pl != 0UL)
    ph++;
  ularith_add_ul_2ul (&(t[2]), &(t[3]), ph);      /* t3:t2:t1 <= 2W^2 + W^(3/2) - 2W - 3W^(1/2), t1 = 0 so
						     t3:t2 <= 2W + W^(1/2) - 3 */
  ularith_mul_ul_ul_2ul (&pl, &ph, k, m[0].m[1]); /* ph:pl <= W^(3/2) - W - W^(1/2) + 1 */
  ularith_add_2ul_2ul (&(t[2]), &(t[3]), pl, ph); /* t3:t2 <= W^(3/2) + W - 2 */

  /* Result may be larger than m, but is < 2*m */

  ularith_sub_2ul_2ul_ge (&(t[2]), &(t[3]), m[0].m[0], m[0].m[1]);

  r[0] = t[2];
  r[1] = t[3];
#if defined(MODTRACE)
  printf (" == (%lu * 2^%d + %lu) /* PARI */ \n", r[1], LONG_BIT, r[0]);
#endif
  ASSERT_EXPENSIVE (modredc15ul_intcmp (r, m[0].m) < 0);
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


/* Division by small integer n, where (n-1)*m may NOT overflow the most 
   significant word. Returns 1 if n is invertible modulo m, 0 if not. 
   
   w_mod_n is word base (e.g., 2^32 or  2^64) mod n
   inv_n contains -1/i (mod n) if i is coprime to n, or 0 if i is not coprime 
   to n, for 0 <= i < n
   c = n^(-1) (mod word base)
*/

static inline int
modredc15ul_divn (residueredc15ul_t r, const residueredc15ul_t a, 
		  const unsigned long n, const unsigned long w_mod_n, 
		  const unsigned long *inv_n, const unsigned long c,
		  const modulusredc15ul_t m)
{
  const unsigned long an = ((a[1] % n)*w_mod_n + a[0] % n) % n;
  const unsigned long mn = ((m[0].m[1] % n)*w_mod_n + m[0].m[0] % n) % n;
  unsigned long k;
  residueredc15ul_t t;
  
  modredc15ul_init_noset0 (t, m);
  t[1] = a[1];
  t[0] = a[0];
  
  if (inv_n[mn] == 0)
    {
      modredc15ul_clear (t, m);
      return 0;
    }

  /* Make t[1]:t[0] divisible by n */
  /* We want a+km == 0 (mod n), so k = -a*m^{-1} (mod n) */
  k = (inv_n[mn] * an) % n;
  ularith_mul_ul_ul_2ul (&(t[0]), &(t[1]), m[0].m[0], k);
  t[1] += m[0].m[1] * k;
  ularith_add_2ul_2ul (&(t[0]), &(t[1]), a[0], a[1]);
  
  /* Now t[1]:t[0] is divisible by n */
  ASSERT_EXPENSIVE (((t[1] % n)*w_mod_n + t[0] % n) % n == 0UL);
  
  r[1] = t[1] / n;
  r[0] = t[0] * c;

#ifdef WANT_ASSERT_EXPENSIVE
  {
    unsigned long i;
    modredc15ul_set (t, r, m);
    for (i = 1; i < n; i++)
      modredc15ul_add (t, t, r, m);
    ASSERT_EXPENSIVE (modredc15ul_equal (t, a, m));
  }
#endif

  modredc15ul_clear (t, m);
  return 1;
}



/* prototypes of non-inline functions */
int modredc15ul_div3 (residueredc15ul_t, const residueredc15ul_t, 
		      const modulusredc15ul_t);
int modredc15ul_div5 (residueredc15ul_t, const residueredc15ul_t, 
		      const modulusredc15ul_t);
int modredc15ul_div7 (residueredc15ul_t, const residueredc15ul_t, 
		      const modulusredc15ul_t);
int modredc15ul_div11 (residueredc15ul_t, const residueredc15ul_t, 
		       const modulusredc15ul_t);
int modredc15ul_div13 (residueredc15ul_t, const residueredc15ul_t, 
		       const modulusredc15ul_t);
void modredc15ul_gcd (modintredc15ul_t, const residueredc15ul_t, 
		      const modulusredc15ul_t);
void modredc15ul_pow_ul (residueredc15ul_t, const residueredc15ul_t, 
			 const unsigned long, const modulusredc15ul_t);
void modredc15ul_2pow_ul (residueredc15ul_t, const unsigned long, 
                          const modulusredc15ul_t);
void modredc15ul_pow_mp (residueredc15ul_t, const residueredc15ul_t, 
			 const unsigned long *, const int, 
			 const modulusredc15ul_t);
void modredc15ul_2pow_mp (residueredc15ul_t, const unsigned long *, const int, 
			  const modulusredc15ul_t);
void modredc15ul_V_ul (residueredc15ul_t, const residueredc15ul_t, 
		       const unsigned long, const modulusredc15ul_t);
void modredc15ul_V_mp (residueredc15ul_t, const residueredc15ul_t, 
		       const unsigned long *, const int, 
		       const modulusredc15ul_t);
int modredc15ul_sprp (const residueredc15ul_t, const modulusredc15ul_t);
int modredc15ul_sprp2 (const modulusredc15ul_t);
int modredc15ul_isprime (const modulusredc15ul_t);
int modredc15ul_inv (residueredc15ul_t, const residueredc15ul_t, 
		     const modulusredc15ul_t);
int modredc15ul_jacobi (const residueredc15ul_t, const modulusredc15ul_t);
#endif
