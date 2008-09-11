/* Some functions for modular arithmetic with residues and modulus in
   unsigned long variables. Residues are stored in Montgomery form,
   reduction after multiplication is done with REDC. Due to inlining, 
   this file must be included in the caller's source code with #include */

/* Naming convention: all function start with modredcul, for 
   MODulus REDC Unsigned Long, followed by underscore, functionality of 
   function (add, mul, etc), and possibly underscore and specification of 
   what argument types the function takes (_ul, etc). */

#ifndef __MODREDC_UL_H

#define __MODREDC_UL_H

/**********************************************************************/
#include <assert.h>
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

#define MODREDCUL_SIZE 1
#define MODREDCUL_MAXBITS LONG_BIT

typedef unsigned long residueredcul_t[MODREDCUL_SIZE];
typedef unsigned long modintredcul_t[MODREDCUL_SIZE];
typedef struct { 
  unsigned long m;
  unsigned long invm;
  residueredcul_t one;
} __modulusredcul_t;
typedef __modulusredcul_t modulusredcul_t[1];


/* ==================== Functions used internally ==================== */

/* Computes (a * 2^wordsize) % m */
MAYBE_UNUSED
static inline void
modredcul_tomontgomery (residueredcul_t r, const residueredcul_t a, 
                        const modulusredcul_t m)
{
  ASSERT_EXPENSIVE (a[0] < m[0].m);
  ularith_div_2ul_ul_ul_r (r, 0UL, a[0], m[0].m);
}


/* Computes (a / 2^wordsize) % m */
MAYBE_UNUSED
static inline void
modredcul_frommontgomery (residueredcul_t r, const residueredcul_t a, 
                          const modulusredcul_t m)
{
  unsigned long tlow, thigh;
  tlow = a[0] * m[0].invm;
  ularith_mul_ul_ul_2ul (&tlow, &thigh, tlow, m[0].m);
  r[0] = thigh + ((a[0] != 0UL) ? 1UL : 0UL);
}


/* Compute 1/n (mod 2^wordsize) */
MAYBE_UNUSED
static inline unsigned long
modredcul_invmodul (unsigned long n)
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
modredcul_intset (modintredcul_t r, const modintredcul_t s)
{
  r[0] = s[0];
}


MAYBE_UNUSED
static inline void
modredcul_intset_ul (modintredcul_t r, const unsigned long s)
{
  r[0] = s;
}


MAYBE_UNUSED
static inline int
modredcul_intequal (const modintredcul_t a, const modintredcul_t b)
{
  return (a[0] == b[0]);
}


MAYBE_UNUSED
static inline int
modredcul_intequal_ul (const modintredcul_t a, const unsigned long b)
{
  return (a[0] == b);
}


MAYBE_UNUSED
static inline int
modredcul_intcmp (const modintredcul_t a, const modintredcul_t b)
{
  return (a[0] < b[0]) ? -1 : (a[0] == b[0]) ? 0 : 1;
}


MAYBE_UNUSED
static inline int
modredcul_intcmp_ul (const modintredcul_t a, const unsigned long b)
{
  return (a[0] < b) ? -1 : (a[0] == b) ? 0 : 1;
}

MAYBE_UNUSED
static inline int
modredcul_intfits_ul (const modintredcul_t a MAYBE_UNUSED)
{
  return 1;
}

/* Returns the number of bits in a, that is, floor(log_2(n))+1. 
   For n==0 returns 0. */
MAYBE_UNUSED
static inline int
modredcul_intbits (const modintredcul_t a)
{
  int bits = 0;
  unsigned long n = a[0];
  while (n > 0UL) /* TODO: use clzl */
    {
      bits++;
      n >>= 1;
    }
  return bits;
}

MAYBE_UNUSED
static inline void
modredcul_intshr (modintredcul_t r, const modintredcul_t s, const int i)
{
  r[0] = s[0] >> i;
}


MAYBE_UNUSED
static inline void
modredcul_intshl (modintredcul_t r, const modintredcul_t s, const int i)
{
  r[0] = s[0] << i;
}


/* r = n/d. We require d|n */
MAYBE_UNUSED
static inline void
modredcul_intdivexact (modintredcul_t r, const modintredcul_t n, 
                       const modintredcul_t d)
{
  r[0] = n[0] / d[0]; 
}


/* Functions for the modulus */

MAYBE_UNUSED
static inline void
modredcul_initmod_ul (modulusredcul_t m, const unsigned long s)
{
  m[0].m = s;
  m[0].invm = -modredcul_invmodul (s);
  if (m[0].m == 1UL)
    m[0].one[0] = 0UL;
  else
    {
      m[0].one[0] = 1UL;
      modredcul_tomontgomery (m[0].one, m[0].one, m);
    }
}


MAYBE_UNUSED
static inline void
modredcul_initmod_uls (modulusredcul_t m, const modintredcul_t s)
{
  m[0].m = s[0];
  m[0].invm = -modredcul_invmodul (s[0]);
  if (m[0].m == 1UL)
    m[0].one[0] = 0UL;
  else
    {
      m[0].one[0] = 1UL;
      modredcul_tomontgomery (m[0].one, m[0].one, m);
    }
}


MAYBE_UNUSED
static inline unsigned long
modredcul_getmod_ul (const modulusredcul_t m)
{
  return m[0].m;
}


MAYBE_UNUSED
static inline void
modredcul_getmod_uls (modintredcul_t r, const modulusredcul_t m)
{
  r[0] = m[0].m;
}


MAYBE_UNUSED
static inline void
modredcul_clearmod (modulusredcul_t m MAYBE_UNUSED)
{
  return;
}


/* Functions for residues */

/* Initialises a residue_t type and sets it to zero */
MAYBE_UNUSED
static inline void
modredcul_init (residueredcul_t r, const modulusredcul_t m MAYBE_UNUSED)
{
  r[0] = 0UL;
}


/* Initialises a residue_t type, but does not set it to zero. For fixed length
   residue_t types, that leaves nothing to do at all. */
MAYBE_UNUSED
static inline void
modredcul_init_noset0 (residueredcul_t r MAYBE_UNUSED, 
                       const modulusredcul_t m MAYBE_UNUSED)
{
  return;
}


MAYBE_UNUSED
static inline void
modredcul_clear (residueredcul_t r MAYBE_UNUSED, 
                 const modulusredcul_t m MAYBE_UNUSED)
{
  return;
}


MAYBE_UNUSED
static inline void
modredcul_set (residueredcul_t r, const residueredcul_t s, 
               const modulusredcul_t m MAYBE_UNUSED)
{
  ASSERT_EXPENSIVE (s[0] < m[0].m);
  r[0] = s[0];
}


MAYBE_UNUSED
static inline void
modredcul_set_ul (residueredcul_t r, const unsigned long s, 
                  const modulusredcul_t m)
{
  r[0] = s % m[0].m;
  modredcul_tomontgomery (r, r, m);
}


/* Sets the residue_t to the class represented by the integer s. Assumes that
   s is reduced (mod m), i.e. 0 <= s < m */

MAYBE_UNUSED
static inline void
modredcul_set_ul_reduced (residueredcul_t r, const unsigned long s, 
                          const modulusredcul_t m MAYBE_UNUSED)
{
  ASSERT (s < m[0].m);
  r[0] = s;
  modredcul_tomontgomery (r, r, m);
}


MAYBE_UNUSED
static inline void
modredcul_set_uls (residueredcul_t r, const modintredcul_t s, 
		   const modulusredcul_t m)
{
  r[0] = s[0] % m[0].m;
  modredcul_tomontgomery (r, r, m);
}


MAYBE_UNUSED
static inline void
modredcul_set_uls_reduced (residueredcul_t r, const modintredcul_t s, 
			   const modulusredcul_t m)
{
  ASSERT (s[0] < m[0].m);
  r[0] = s[0];
  modredcul_tomontgomery (r, r, m);
}


/* This one is so trivial that we don't really require m in the
 * interface. For interface homogeneity we make it take the m parameter 
 * anyway.
 */
MAYBE_UNUSED 
static inline void 
modredcul_set0 (residueredcul_t r, const modulusredcul_t m MAYBE_UNUSED) 
{ 
  r[0] = 0UL;
}


MAYBE_UNUSED 
static inline void 
modredcul_set1 (residueredcul_t r, const modulusredcul_t m) 
{
  r[0] = m[0].one[0];
}


/* Exchanges the values of the two arguments */

MAYBE_UNUSED
static inline void
modredcul_swap (residueredcul_t a, residueredcul_t b, 
                const modulusredcul_t m MAYBE_UNUSED)
{
  unsigned long t;
  ASSERT_EXPENSIVE (a[0] < m[0].m && b[0] < m[0].m);
  t = a[0];
  a[0] = b[0];
  b[0] = t;
}


MAYBE_UNUSED
static inline unsigned long
modredcul_get_ul (const residueredcul_t s, 
	          const modulusredcul_t m MAYBE_UNUSED)
{
  unsigned long t;
  ASSERT_EXPENSIVE (s[0] < m[0].m);
  modredcul_frommontgomery (&t, s, m);
  return t;
}


MAYBE_UNUSED
static inline void
modredcul_get_uls (modintredcul_t r, const residueredcul_t s, 
		   const modulusredcul_t m MAYBE_UNUSED)
{
  ASSERT_EXPENSIVE (s[0] < m[0].m);
  modredcul_frommontgomery (r, s, m);
}


MAYBE_UNUSED
static inline int
modredcul_equal (const residueredcul_t a, const residueredcul_t b, 
             const modulusredcul_t m MAYBE_UNUSED)
{
  ASSERT_EXPENSIVE (a[0] < m[0].m && b[0] < m[0].m);
  return (a[0] == b[0]);
}


MAYBE_UNUSED
static inline int
modredcul_is0 (const residueredcul_t a, const modulusredcul_t m MAYBE_UNUSED)
{
  ASSERT_EXPENSIVE (a[0] < m[0].m);
  return (a[0] == 0UL);
}


MAYBE_UNUSED
static inline int
modredcul_is1 (const residueredcul_t a, const modulusredcul_t m MAYBE_UNUSED)
{
  ASSERT_EXPENSIVE (a[0] < m[0].m);
  return (a[0] == m[0].one[0]);
}


MAYBE_UNUSED
static inline void
modredcul_add (residueredcul_t r, const residueredcul_t a, 
               const residueredcul_t b, const modulusredcul_t m)
{
  ASSERT_EXPENSIVE (a[0] < m[0].m && b[0] < m[0].m);
#ifdef MODTRACE
  printf ("modul_add: a = %lu, b = %lu", a[0], b[0]);
#endif

#if (defined(__i386__) || defined(__x86_64__)) && defined(__GNUC__)
  {
    unsigned long t = a[0] - m[0].m, tr = a[0] + b[0];
    
    __asm__ (
      "add %2, %1\n\t"   /* t += b */
      "cmovc %1, %0\n\t"  /* if (cy) tr = t */
      : "+r" (tr), "+&r" (t)
      : "g" (b[0])
      : "cc"
    );
    ASSERT_EXPENSIVE (tr == ((a[0] >= m[0].m - b[0]) ? (a[0] - (m[0].m - b[0])) : (a[0] + b[0])));
    r[0] = tr;
  }
#else
  r[0] = (a[0] >= m[0].m - b[0]) ? (a[0] - (m[0].m - b[0])) : (a[0] + b[0]);
#endif

#ifdef MODTRACE
  printf (", r = %lu\n", r[0]);
#endif
  ASSERT_EXPENSIVE (r[0] < m[0].m);
}


MAYBE_UNUSED
static inline void
modredcul_add_ul (residueredcul_t r, const residueredcul_t a,
                  const unsigned long b, const modulusredcul_t m)
{
  residueredcul_t t;

  /* TODO: speed up */
  ASSERT_EXPENSIVE (a[0] < m[0].m);
  modredcul_init_noset0 (t, m);
  modredcul_set_ul (t, b, m);
  modredcul_add (r, a, t, m);
  modredcul_clear (t, m);
}


MAYBE_UNUSED
static inline void
modredcul_sub (residueredcul_t r, const residueredcul_t a, 
               const residueredcul_t b, const modulusredcul_t m)
{
  ASSERT_EXPENSIVE (a[0] < m[0].m && b[0] < m[0].m);
#ifdef MODTRACE
  printf ("submod_ul: a = %lu, b = %lu", a[0], b[0]);
#endif

#if (defined(__i386__) || defined(__x86_64__)) && defined(__GNUC__)
  {
    unsigned long t = 0, tr = a[0];
    __asm__ (
      "sub %2, %0\n\t"  /* tr -= b */
      "cmovc %3, %1\n\t" /* if (a < b) t = m */
      "add %1, %0\n"    /* tr += t. Moving this out of the asm block results
                            in slowdown!?! */
      : "+&r" (tr), "+&r" (t)
      : "g" (b[0]), "rm" (m[0].m)
      : "cc"
    );
    r[0] = tr;
  }
#elif 1
  /* Seems to be faster than the one below */
  {
    unsigned long t;
    t = a[0] - b[0];
    if (a[0] < b[0])
      t += m[0].m;
    r[0] = t;
  }
#else
  r[0] = (a[0] < b[0]) ? (a[0] - b[0] + m[0]) : (a[0] - b[0]);
#endif

#ifdef MODTRACE
  printf (", r = %lu\n", r[0]);
#endif
  ASSERT_EXPENSIVE (r[0] < m[0].m);
}


MAYBE_UNUSED
static inline void
modredcul_sub_ul (residueredcul_t r, const residueredcul_t a,
                  const unsigned long b, const modulusredcul_t m)
{
  residueredcul_t t;

  /* TODO: speed up */
  ASSERT_EXPENSIVE (a[0] < m[0].m);
  modredcul_init_noset0 (t, m);
  modredcul_set_ul (t, b, m);
  modredcul_sub (r, a, t, m);
  modredcul_clear (t, m);
}


MAYBE_UNUSED
static inline void
modredcul_neg (residueredcul_t r, const residueredcul_t a, 
               const modulusredcul_t m)
{
  ASSERT_EXPENSIVE (a[0] < m[0].m);
  if (a[0] == 0UL)
    r[0] = a[0];
  else
    r[0] = m[0].m - a[0];
}


/* Computes (a / 2^wordsize) % m, but result can be r = m. 
   Input a must not be equal 0 */
MAYBE_UNUSED
static inline void
modredcul_redcsemi_ul_not0 (residueredcul_t r, const unsigned long a, 
                            const modulusredcul_t m)
{
  unsigned long tlow, thigh;

  ASSERT (a != 0);

  tlow = a * m[0].invm; /* tlow <= 2^w-1 */
  ularith_mul_ul_ul_2ul (&tlow, &thigh, tlow, m[0].m);
  /* thigh:tlow <= (2^w-1) * m */
  r[0] = thigh + 1UL; 
  /* (thigh+1):tlow <= 2^w + (2^w-1) * m  <= 2^w + 2^w*m - m 
                    <= 2^w * (m + 1) - m */
  /* r <= floor ((2^w * (m + 1) - m) / 2^w) <= floor((m + 1) - m/2^w)
       <= m */
}


/* Computes ((a + b) / 2^wordsize) % m. a <= m is permissible */
MAYBE_UNUSED
static inline void
modredcul_addredc_ul (residueredcul_t r, const residueredcul_t a, 
                      const unsigned long b, const modulusredcul_t m)
{
  unsigned long slow, shigh, tlow, thigh;
  
  ASSERT_EXPENSIVE (a[0] <= m[0].m);
  slow = b;
  shigh = 0UL;
  ularith_add_ul_2ul (&slow, &shigh, a[0]);
  
  tlow = slow * m[0].invm;
  ularith_mul_ul_ul_2ul (&tlow, &thigh, tlow, m[0].m);
  ASSERT_EXPENSIVE (slow + tlow == 0UL);
  r[0] = thigh + shigh + ((slow != 0UL) ? 1UL : 0UL);
  
  /* r = ((a+b) + (((a+b)%2^w * invm) % 2^w) * m) / 2^w  Use a<=m-1, b<=2^w-1
     r <= (m + 2^w - 1 + (2^w - 1) * m) / 2^w
        = (m - 1 + 2^w + m*2^w - m) / 2^w
        = (- 1 + 2^w + m2^w) / 2^w
        = m + 1 - 1/2^w
     r <= m, since r is an integer
  */
  if (r[0] == m[0].m)
    r[0] = 0UL;
}


/* Computes ((a + b) / 2^wordsize) % m, but result can be == m.
   a <= m is permissible */
MAYBE_UNUSED
static inline void
modredcul_addredcsemi_ul (residueredcul_t r, const residueredcul_t a, 
                          const unsigned long b, const modulusredcul_t m)
{
  unsigned long slow, shigh, tlow;
  unsigned char sb;
  
  ASSERT_EXPENSIVE(a[0] <= m[0].m);
  slow = b;
#if defined(__x86_64__) && defined(__GNUC__)
   __asm__ ( "addq %2, %0\n\t" /* cy * 2^w + slow = a + b */
            "setne %1\n\t"     /* if (slow != 0) sb = 1 */
            "adcb $0, %1\n"    /* sb += cy */
            : "+&r" (slow), "=qm" (sb)
            : "rm" (a[0])
            : "cc");
  shigh = sb;
#elif defined(__i386__) && defined(__GNUC__)
   __asm__ ( "addl %2, %0\n\t"
            "setne %1\n\t"
            "adcb $0, %1\n"
            : "+&r" (slow), "=qm" (sb)
            : "rm" (a[0])
            : "cc");
  shigh = sb;
#else
  shigh = 0UL;
  ularith_add_ul_2ul (&slow, &shigh, a[0]);
  shigh += (slow != 0UL) ? 1UL : 0UL;
#endif

  tlow = slow * m[0].invm;
  ularith_mul_ul_ul_2ul (&tlow, r, tlow, m[0].m);
  ASSERT_EXPENSIVE (slow + tlow == 0UL);
  r[0] += shigh;

  /* r = ((a+b) + (((a+b)%2^w * invm) % 2^w) * m) / 2^w
     r <= ((a+b) + (2^w - 1) * m) / 2^w
     r <= (m + 2^w-1 + m*2^w - m) / 2^w
     r <= (2^w -1 + p2^w) / 2^w
     r <= p + 1 - 1/2^w
     r <= p
  */
}


MAYBE_UNUSED
static inline void
modredcul_mul (residueredcul_t r, const residueredcul_t a, 
               const residueredcul_t b, const modulusredcul_t m)
{
  unsigned long plow, phigh, tlow, thigh;

  ASSERT_EXPENSIVE (m[0].m % 2 != 0);
  ASSERT_EXPENSIVE (a[0] < m[0].m && b[0] < m[0].m);
#if defined(MODTRACE)
  printf ("(%lu * %lu / 2^%ld) %% %lu", 
          a[0], b[0], 8 * sizeof(unsigned long), m[0].m);
#endif
  
  ularith_mul_ul_ul_2ul (&plow, &phigh, a[0], b[0]);
  tlow = plow * m[0].invm;
  ularith_mul_ul_ul_2ul (&tlow, &thigh, tlow, m[0].m);
  /* Let w = 2^wordsize. We know (phigh * w + plow) + (thigh * w + tlow) 
     == 0 (mod w) so either plow == tlow == 0, or plow !=0 and tlow != 0. 
     In the former case we want phigh + thigh + 1, in the latter 
     phigh + thigh */
  /* Since a <= p-1 and b <= p-1, and p <= w-1, a*b <= w^2 - 4*w + 4, so
     adding 1 to phigh is safe */
#if 0
  /* Slower? */
  ularith_add_ul_2ul (&plow, &phigh, tlow);
#else
  phigh += (plow != 0UL) ? 1UL : 0UL;
#endif

  modredcul_add (r, &phigh, &thigh, m);

#if defined(MODTRACE)
  printf (" == %lu /* PARI */ \n", r[0]);
#endif
}
                         

/* Computes (a * b + c)/ 2^wordsize % m. Requires that 
   a * b + c < 2^wordsize * m */

MAYBE_UNUSED
static inline void
modredcul_muladdredc (residueredcul_t r, const residueredcul_t a, 
		      const residueredcul_t b, const residueredcul_t c, 
		      const modulusredcul_t m)
{
  unsigned long plow, phigh, tlow, thigh;
  ASSERT_EXPENSIVE (m[0].m % 2 != 0);
#if defined(MODTRACE)
  printf ("(%lu * %lu / 2^%ld) %% %lu", 
          a[0], b, 8 * sizeof(unsigned long), m[0]);
#endif
  
  ularith_mul_ul_ul_2ul (&plow, &phigh, a[0], b[0]);
  ularith_add_ul_2ul (&plow, &phigh, c[0]);
  tlow = plow * m[0].invm;
  ularith_mul_ul_ul_2ul (&tlow, &thigh, tlow, m[0].m);
  phigh += (plow != 0UL ? 1UL : 0UL);
  modredcul_add (r, &phigh, &thigh, m);
  
#if defined(MODTRACE)
  printf (" == %lu /* PARI */ \n", r[0]);
#endif
}
                         

MAYBE_UNUSED
static inline void
modredcul_div2 (residueredcul_t r, const residueredcul_t a, 
                const modulusredcul_t m)
{
  ASSERT_EXPENSIVE (m[0].m % 2UL != 0UL);
  if (a[0] % 2UL == 0UL)
    r[0] = a[0] / 2UL;
  else
    r[0] = a[0] / 2UL + m[0].m / 2UL + 1UL;
}


MAYBE_UNUSED
static inline int
modredcul_next (residueredcul_t r, const modulusredcul_t m)
{
    return (++r[0] == m[0].m);
}


MAYBE_UNUSED
static inline int
modredcul_finished (const residueredcul_t r, const modulusredcul_t m)
{
    return (r[0] == m[0].m);
}


/* prototypes of non-inline functions */
void modredcul_div3 (residueredcul_t, const residueredcul_t, 
                     const modulusredcul_t);
void modredcul_div7 (residueredcul_t, const residueredcul_t, 
                     const modulusredcul_t);
void modredcul_gcd (unsigned long *, const residueredcul_t, 
                    const modulusredcul_t);
void modredcul_pow_ul (residueredcul_t, const residueredcul_t, 
                   const unsigned long, const modulusredcul_t);
void modredcul_pow_mp (residueredcul_t, const residueredcul_t, 
                   const unsigned long *, const int, const modulusredcul_t);
void modredcul_2pow_mp (residueredcul_t, const residueredcul_t, 
                    const unsigned long *, const int, const unsigned long, 
                    const modulusredcul_t);
void modredcul_V_ul (residueredcul_t, const residueredcul_t, 
		     const residueredcul_t, const unsigned long, 
		     const modulusredcul_t);
void modredcul_V_mp (residueredcul_t, const residueredcul_t, 
		     const unsigned long *, const int, const modulusredcul_t);
int modredcul_sprp (const residueredcul_t, const modulusredcul_t);
int modredcul_inv (residueredcul_t, const residueredcul_t, 
		   const modulusredcul_t);
int modredcul_jacobi (const residueredcul_t, const modulusredcul_t);

#endif
