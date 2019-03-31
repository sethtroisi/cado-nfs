/* Some commonly used assembly helper macros for arithmetic on 
   unsigned 64-bit integers. 
   Defining U64ARITH_VERBOSE_ASM puts in the asm output a line for each
   function call that shows what registers/memory locations the operands 
   are in.
   Defining U64ARITH_NO_ASM avoids asm macros and uses the C fallback code
   where available.
*/

#ifndef U64_ARITH_H__
#define U64_ARITH_H__

#include <assert.h>
#include <stdint.h>
#include <gmp.h>
#include "macros.h"

#ifdef WANT_ASSERT_EXPENSIVE
#define ASSERT_EXPENSIVE(x) ASSERT(x)
#else
#define ASSERT_EXPENSIVE(x)
#endif


/* Let a = a1 + 2^64 * a2, b = b1 + 2^64 * b2. Return 1 if a > b,
   and 0 if a <= b. */
static inline int
u64arith_gt_2_2(uint64_t, uint64_t, uint64_t, uint64_t) ATTRIBUTE((const));
static inline int
u64arith_gt_2_2(const uint64_t a1, const uint64_t a2,
                const uint64_t b1, const uint64_t b2)
{
  return a2 > b2 || (a2 == b2 && a1 > b1);
}


/* Add two uint64_t to two uint64_t with carry propagation from 
   low word (r1) to high word (r2). Any carry out from high word is lost. */

static inline void
u64arith_add_2_2 (uint64_t *r1, uint64_t *r2, 
		  const uint64_t a1, const uint64_t a2)
{
#ifdef U64ARITH_VERBOSE_ASM
  __asm__ ("# u64arith_add_2_2 (%0, %1, %2, %3)\n" : : 
           "X" (*r1), "X" (*r2), "X" (a1), "X" (a2));
#endif
#if !defined (U64ARITH_NO_ASM) && defined(HAVE_GCC_STYLE_AMD64_INLINE_ASM)
  __asm__ __VOLATILE (
    "addq %2, %0\n\t"
    "adcq %3, %1\n"
    : "+&r" (*r1), "+r" (*r2)
    : "rme" (a1), "rme" (a2)
    : "cc");
#else
  *r1 += a1;
  *r2 += a2 + (*r1 < a1);
#endif
}

/* Add a uint64_t to two uint64_t with carry propagation from low word (r1)
   to high word (r2). Any carry out from high word is lost. */

static inline void
u64arith_add_1_2 (uint64_t *r1, uint64_t *r2, const uint64_t a)
{
  /* We have assembly only for x86_64 and on that architecture, this function
     would generate the same code as u64arith_add_2_2() with an immediate
     zero value for the high word, so we fall back to that. */
  u64arith_add_2_2(r1, r2, a, 0);
}


/* Adds two uint64_t from two uint64_t with carry propagation 
   from low word (r1) to high word (r2). Returns 1 if there was a carry out 
   from high word, otherwise returns 0. */

static inline char
u64arith_add_2_2_cy (uint64_t *r1, uint64_t *r2, 
		     const uint64_t a1, const uint64_t a2)
{
  char cy;
#ifdef U64ARITH_VERBOSE_ASM
  __asm__ ("# u64arith_add_2_2_cy (%0, %1, %2, %3)\n" : : 
           "X" (*r1), "X" (*r2), "X" (a1), "X" (a2));
#endif
#if !defined (U64ARITH_NO_ASM) && defined(HAVE_GCC_STYLE_AMD64_INLINE_ASM)
  __asm__ __VOLATILE (
    "addq %3, %0\n\t"
    "adcq %4, %1\n\t"
    "setc %2\n"
    : "+&r" (*r1), "+r" (*r2), "=r" (cy)
    : "rme" (a1), "rme" (a2)
    : "cc");
#else
  uint64_t u1 = *r1 + a1,
           u2 = *r2 + a2 + (u1 < *r1);
  /* Overflow occurred iff the sum is smaller than one of the summands */
  cy = u64arith_gt_2_2(a1, a2, u1, u2);
  *r1 = u1;
  *r2 = u2;
#endif
  return cy;
}


/* Requires a < m and b <= m, then r == a+b (mod m) and r < m */
static inline void
u64arith_addmod_1_1 (uint64_t *r, const uint64_t a,
                     const uint64_t b, const uint64_t m)
{
  ASSERT_EXPENSIVE (a < m && b <= m);

#if (defined(__i386__) && defined(__GNUC__)) || defined(HAVE_GCC_STYLE_AMD64_INLINE_ASM)
  {
    uint64_t t = a + b, tr = a - m;

    __asm__ __VOLATILE (
      "add %2, %0\n\t"   /* tr += b */
      "cmovnc %1, %0\n\t"  /* if (!cy) tr = t */
      : "+&r" (tr)
      : "rm" (t), "rme" (b)
      : "cc"
    );
    ASSERT_EXPENSIVE (tr == ((a >= m - b) ? (a - (m - b)) : (a + b)));
    r[0] = tr;
  }
#else
  r[0] = (b >= m - a) ? (b - (m - a)) : (a + b);
#endif

  ASSERT_EXPENSIVE (r[0] < m);
}


/* Subtract an uint64_t from two uint64_t with borrow propagation 
   from low word (r1) to high word (r2). Any borrow out from high word is 
   lost. */

static inline void
u64arith_sub_1_2 (uint64_t *r1, uint64_t *r2, 
                  const uint64_t a)
{
#ifdef U64ARITH_VERBOSE_ASM
  __asm__ ("# u64arith_sub_1_2  (%0, %1, %2)\n" : : 
           "X" (*r1), "X" (*r2), "X" (a));
#endif
#if !defined (U64ARITH_NO_ASM) && defined(HAVE_GCC_STYLE_AMD64_INLINE_ASM)
  __asm__ __VOLATILE (
    "subq %2, %0\n\t"
    "sbbq $0, %1\n"
    : "+&r" (*r1), "+r" (*r2)
    : "rme" (a)
    : "cc");
#else
  uint64_t u = *r1;
  *r1 -= a;
  if (*r1 > u)
    (*r2)--;
#endif
}


/* Subtract two uint64_t from two uint64_t with borrow propagation 
   from low word (r1) to high word (r2). Any borrow out from high word is 
   lost. */

static inline void
u64arith_sub_2_2 (uint64_t *r1, uint64_t *r2, 
		  const uint64_t a1, const uint64_t a2)
{
#ifdef U64ARITH_VERBOSE_ASM
  __asm__ ("# u64arith_sub_2_2 (%0, %1, %2, %3)\n" : : 
           "X" (*r1), "X" (*r2), "X" (a1), "X" (a2));
#endif
#if !defined (U64ARITH_NO_ASM) && defined(HAVE_GCC_STYLE_AMD64_INLINE_ASM)
  __asm__ __VOLATILE (
    "subq %2, %0\n\t"
    "sbbq %3, %1\n"
    : "+&r" (*r1), "+r" (*r2)
    : "rme" (a1), "rme" (a2)
    : "cc");
#else
  uint64_t u = *r1;
  *r1 -= a1;
  *r2 -= a2;
  if (*r1 > u)
    (*r2)--;
#endif
}

/* Subtract two uint64_t from two uint64_t with borrow propagation 
   from low word (r1) to high word (r2). Returns 1 if there was a borrow out 
   from high word, otherwise returns 0. */

static inline char
u64arith_sub_2_2_cy (uint64_t *r1, uint64_t *r2, 
                     const uint64_t a1, const uint64_t a2)
{
  char cy;
#ifdef U64ARITH_VERBOSE_ASM
  __asm__ ("# u64arith_sub_2_2_cy (%0, %1, %2, %3)\n" : : 
           "X" (*r1), "X" (*r2), "X" (a1), "X" (a2));
#endif
#if !defined (U64ARITH_NO_ASM) && defined(HAVE_GCC_STYLE_AMD64_INLINE_ASM)
  __asm__ __VOLATILE (
    "subq %3, %0\n\t"
    "sbbq %4, %1\n\t"
    "setc %2\n"
    : "+&r" (*r1), "+r" (*r2), "=r" (cy)
    : "rme" (a1), "rme" (a2)
    : "cc");
#else
  uint64_t u1 = *r1 - a1, u2 = *r2 - a2;
  if (a1 > *r1)
    u2--;
  cy = u64arith_gt_2_2(a1, a2, *r1, *r2);
  *r1 = u1;
  *r2 = u2;
#endif
  return cy;
}


/* Subtract only if result is non-negative */

static inline void
u64arith_sub_1_1_ge (uint64_t *r, const uint64_t a)
{
  uint64_t t = *r;
#ifdef U64ARITH_VERBOSE_ASM
  __asm__ ("# u64arith_sub_1_1_ge (%0, %1, %2, %3)\n" : : 
           "X" (*r1), "X" (*r2), "X" (a1), "X" (a2));
#endif
#if !defined (U64ARITH_NO_ASM) && defined(HAVE_GCC_STYLE_AMD64_INLINE_ASM)
  __asm__ __VOLATILE (
    "subq %2, %0\n\t" /* r -= a */
    "cmovc %1, %0\n\t" /* If there's a borrow, restore r from t */
    : "+&r" (*r)
    : "r" (t), "rme" (a)
    : "cc");
#else
  t -= a;
  if (*r >= a)
    *r = t;
#endif
}


static inline void
u64arith_sub_2_2_ge (uint64_t *r1, uint64_t *r2, 
                     const uint64_t a1, const uint64_t a2)
{
  uint64_t t1 = *r1, t2 = *r2;
#ifdef U64ARITH_VERBOSE_ASM
  __asm__ ("# u64arith_sub_2_2_ge (%0, %1, %2, %3)\n" : : 
           "X" (*r1), "X" (*r2), "X" (a1), "X" (a2));
#endif
#if !defined (U64ARITH_NO_ASM) && defined(HAVE_GCC_STYLE_AMD64_INLINE_ASM)
  __asm__ __VOLATILE (
    "subq %4, %0\n\t" /* r1 -= a1 */
    "sbbq %5, %1\n\t" /* r2 -= a2 + cy */
    "cmovc %2, %0\n\t" /* If there's a borrow, restore r1 from t1 */
    "cmovc %3, %1\n\t" /* and r2 from t2 */
    : "+&r" (*r1), "+&r" (*r2)
    : "r" (t1), "r" (t2), "rme" (a1), "rme" (a2)
    : "cc");
#else
  if (!u64arith_gt_2_2(a1, a2, *r1, *r2))
    {
      *r1 = t1 - a1;
      *r2 = t2 - a2 - (a1 > t1);
    }
#endif
}


static inline void
u64arith_submod_1_1 (uint64_t *r, const uint64_t a,
                     const uint64_t b, const uint64_t m)
{
  ASSERT_EXPENSIVE (a < m && b < m);
#if !defined (U64ARITH_NO_ASM) && defined(HAVE_GCC_STYLE_AMD64_INLINE_ASM)
  {
    uint64_t tr, t = a;
    __asm__ __VOLATILE (
      "sub %2, %1\n\t"  /* t -= b ( = a - b) */
      "lea (%1,%3,1), %0\n\t" /* tr = t + m ( = a - b + m) */
      "cmovnc %1, %0\n\t" /* if (a >= b) tr = t */
      : "=&r" (tr), "+&r" (t)
      : "rme" (b), "r" (m)
      : "cc"
    );
    r[0] = tr;
  }
#elif 1
  /* Seems to be faster than the one below */
  {
    uint64_t t = 0UL, tr;
    if ((tr = a - b) > a)
      t = m;
    r[0] = tr + t;
  }
#else
  r[0] = (a < b) ? (a - b + m) : (a - b);
#endif

  ASSERT_EXPENSIVE (r[0] < m);
}


/* Multiply two uint64_t "a" and "b" and put the result as 
   r2:r1 (r2 being the high word) */

static inline void
u64arith_mul_1_1_2 (uint64_t *r1, uint64_t *r2, 
		    const uint64_t a, const uint64_t b)
{
#ifdef U64ARITH_VERBOSE_ASM
  __asm__ ("# u64arith_mul_1_1_2 (%0, %1, %2, %3)\n" : : 
           "X" (*r1), "X" (*r2), "X" (a), "X" (b));
#endif
#if !defined (U64ARITH_NO_ASM) && defined(HAVE_GCC_STYLE_AMD64_INLINE_ASM)
  __asm__ __VOLATILE (
    "mulq %3"
    : "=a" (*r1), "=d" (*r2)
    : "%0" (a), "rm" (b)
    : "cc");
#elif 0 && !defined (U64ARITH_NO_ASM) && defined(HAVE_GCC_STYLE_ARM_INLINE_ASM)
/*
  This is the 32-bit code. Do not use.
  TODO: use correct instruction for ARM64. Need an ARM64 system to test.
  Raspberry Pi in principle has a 64-bit cpu, but Raspbian runs it in
  32-bit mode. There are experimental 64-bit OS for the Pi, though. */
*/
  __asm__ __VOLATILE(
   "umull   %[r1], %[r2], %[a], %[b]\n\t"
  : [r1] "=&r" (*r1), [r2] "=&r" (*r2)
  : [a] "r" (a), [b] "r" (b)
  );
#else
  const uint64_t mask = (1UL << 32) - 1UL;
  const uint64_t ah = a >> 32, al = a & mask, bh = b >> 32, bl = b & mask;
  uint64_t t1, t2, p1, p2;

  t1 = al * bl;
  t2 = 0UL;
  p1 = ah * bl;
  p2 = p1 >> 32;
  p1 = (p1 & mask) << 32;
  u64arith_add_2_2 (&t1, &t2, p1, p2);
  p1 = al * bh;
  p2 = p1 >> 32;
  p1 = (p1 & mask) << 32;
  u64arith_add_2_2 (&t1, &t2, p1, p2);
  t2 += ah * bh;
  *r1 = t1; 
  *r2 = t2;
#endif
}


static inline void
u64arith_sqr_1_2 (uint64_t *r1, uint64_t *r2, 
		  const uint64_t a)
{
#ifdef U64ARITH_VERBOSE_ASM
  __asm__ ("# u64arith_sqr_1_2 (%0, %1, %2)\n" : : 
           "X" (*r1), "X" (*r2), "X" (a));
#endif
#if !defined (U64ARITH_NO_ASM) && defined(HAVE_GCC_STYLE_AMD64_INLINE_ASM)
  __asm__ __VOLATILE (
    "mulq %%rax"
    : "=a" (*r1), "=d" (*r2)
    : "0" (a)
    : "cc");
#elif 0 && !defined (U64ARITH_NO_ASM) && defined(HAVE_GCC_STYLE_ARM_INLINE_ASM)
/*
  This is the 32-bit code. Do not use.
  TODO: use correct instruction for ARM64. Need an ARM64 system to test.
  Raspberry Pi in principle has a 64-bit cpu, but Raspbian runs it in
  32-bit mode. There are experimental 64-bit OS for the Pi, though. */
*/
  __asm__ __VOLATILE(
   "umull   %[r1], %[r2], %[a], %[a]\n\t"
  : [r1] "=&r" (*r1), [r2] "=&r" (*r2)
  : [a] "r" (a)
  );
#else
  const uint64_t mask = (1UL << 32) - 1UL;
  const uint64_t ah = a >> 32, al = a & mask;
  uint64_t t1, t2, p1, p2;

  t1 = al * al;
  t2 = 0UL;
  p1 = ah * al;
  p2 = p1 >> 32;
  p1 = (p1 & mask) << 32;
  u64arith_add_2_2 (&t1, &t2, p1, p2);
  u64arith_add_2_2 (&t1, &t2, p1, p2);
  t2 += ah * ah;
  *r1 = t1; 
  *r2 = t2;
#endif
}


/* Integer division of a two uint64_t value a2:a1 by a uint64_t divisor. Returns
   quotient and remainder. */

static inline void
u64arith_divqr_2_1_1 (uint64_t *q, uint64_t *r,
		      const uint64_t a1, const uint64_t a2,
		      const uint64_t b)
{
  ASSERT(a2 < b); /* Or there will be quotient overflow */
#ifdef U64ARITH_VERBOSE_ASM
  __asm__ ("# u64arith_divqr_2_1_1 (%0, %1, %2, %3, %4)\n" : : 
           "X" (*q), "X" (*r), "X" (a1), "X" (a2), "X" (b));
#endif
#if !defined (U64ARITH_NO_ASM) && defined(HAVE_GCC_STYLE_AMD64_INLINE_ASM)
  __asm__ __VOLATILE (
    "divq %4"
    : "=a" (*q), "=d" (*r)
    : "0" (a1), "1" (a2), "rm" (b)
    : "cc");
#else
  mp_limb_t A[2] = {a1, a2};
  ASSERT(sizeof(uint64_t) == sizeof(mp_limb_t));
  r[0] = mpn_divmod_1 (A, A, 2, b);
  q[0] = A[0];
#endif
}


/* Integer division of a two uint64_t value by a uint64_t divisor. Returns
   only remainder. */

static inline void
u64arith_divr_2_1_1 (uint64_t *r, uint64_t a1, const uint64_t a2,
                     const uint64_t b)
{
  ASSERT(a2 < b); /* Or there will be quotient overflow */
#ifdef U64ARITH_VERBOSE_ASM
  __asm__ ("# u64arith_divr_2_1_1 (%0, %1, %2, %3)\n" : : 
           "X" (*r), "X" (a1), "X" (a2), "X" (b));
#endif
#if !defined (U64ARITH_NO_ASM) && defined(HAVE_GCC_STYLE_AMD64_INLINE_ASM)
  __asm__ __VOLATILE (
    "divq %3"
    : "+a" (a1), "=d" (*r)
    : "1" (a2), "rm" (b)
    : "cc");
#else
  mp_limb_t A[2] = {a1, a2};
  ASSERT(sizeof(uint64_t) == sizeof(mp_limb_t));
  r[0] = mpn_divmod_1 (A, A, 2, b);
#endif
}


/* Set *r to lo shifted right by i bits, filling in the low bits from hi into the high
   bits of *r. I.e., *r = (hi * 2^64 + lo) / 2^i. Assumes 0 <= i < 64. */
static inline void
u64arith_shrd (uint64_t *r, const uint64_t hi, const uint64_t lo,
              const uint64_t i)
{
  ASSERT_EXPENSIVE (0 <= i && i < 64);
#ifdef U64ARITH_VERBOSE_ASM
/* Disable the "uninitialized" warning here, as *r is only written to and
   does not need to be initialized, but we need to write (*r) here so the
   "X" constraint can be resolved even when r does not have an address, e.g.,
   when it is passed around in a register. It seems that "X" is assumed by
   gcc as possibly referring to an input, and since "X" matches anything,
   that's probably a neccessary assumtion to make. */
#if GNUC_VERSION_ATLEAST(4,4,0)
#if GNUC_VERSION_ATLEAST(4,6,0)
#pragma GCC diagnostic push
#endif
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#endif
  __asm__ ("# u64arith_shrd (*r=%0, hi=%1, lo=%2, i=%3)\n" : : 
           "X" (*r), "X" (hi), "X" (lo), "X" (i));
#if GNUC_VERSION_ATLEAST(4,6,0)
#pragma GCC diagnostic pop
#endif
#endif

#if !defined (U64ARITH_NO_ASM) && defined(HAVE_GCC_STYLE_AMD64_INLINE_ASM)
  __asm__ __VOLATILE (
    "shrdq %3, %1, %0\n"
    : "=rm" (*r)
    : "r" (hi), "0" (lo), "cJ" (i) /* i can be in %cl or a literal constant < 64 */
    : "cc");
#else
  if (i > 0) /* shl by 64 is no-op on x86! */
    *r = (lo >> i) | (hi << (64 - i));
  else
    *r = lo;
#endif
}

/* Set *r to hi shifted left by i bits, filling in the high bits from lo into the low
   bits of *r. I.e., *r = (hi + lo*2^-64) * 2^i. Assumes 0 <= i < 64. */
static inline void
u64arith_shld (uint64_t *r, const uint64_t lo, const uint64_t hi,
              const unsigned char i)
{
  ASSERT_EXPENSIVE (0 <= i && i < 64);
#ifdef U64ARITH_VERBOSE_ASM
#if GNUC_VERSION_ATLEAST(4,4,0)
#if GNUC_VERSION_ATLEAST(4,6,0)
#pragma GCC diagnostic push
#endif
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#endif
  __asm__ ("# u64arith_shld (*r=%0, lo=%1, hi=%2, i=%3)\n" : : 
           "X" (*r), "X" (lo), "X" (hi), "X" (i));
#if GNUC_VERSION_ATLEAST(4,6,0)
#pragma GCC diagnostic pop
#endif
#endif

#if !defined (U64ARITH_NO_ASM) && defined(HAVE_GCC_STYLE_AMD64_INLINE_ASM)
  __asm__ __VOLATILE (
    "shldq %3, %1, %0\n"
    : "=rm" (*r)
    : "r" (lo), "0" (hi), "cJ" (i) /* i can be in %cl or a literal constant < 64 */
    : "cc");
#else
  if (i > 0) /* shr by 64 is no-op on x86! */
    *r = (hi << i) | (lo >> (64 - i));
  else
    *r = hi;
#endif
}

/* Returns number of trailing zeros in a. a must not be zero */
static inline int
u64arith_ctz (const uint64_t a)
{
#if !defined (U64ARITH_NO_ASM) && defined(__GNUC__) && \
    (__GNUC__ >= 4 || __GNUC__ >= 3 && __GNUC_MINOR__ >= 4)
  ASSERT_EXPENSIVE (a != 0UL);
  return __builtin_ctzll(a);
#else
  static const unsigned char trailing_zeros[256] =
    {8,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
     5,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
     6,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
     5,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
     7,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
     5,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
     6,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
     5,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0};
  char lsh, t = 0;
  uint64_t y = a;
  ASSERT_EXPENSIVE (a != 0UL);
  do {
    lsh = trailing_zeros [(unsigned char) y];
    y >>= lsh;
    t += lsh;
  } while (lsh == 8);
  return (int) t;
#endif
}

/* Returns number of leading zeros in a. a must not be zero */
static inline int
u64arith_clz (const uint64_t a)
{
#if !defined (U64ARITH_NO_ASM) && defined(__GNUC__) && \
    (__GNUC__ >= 4 || __GNUC__ >= 3 && __GNUC_MINOR__ >= 4)
  ASSERT_EXPENSIVE (a != 0UL);
  if (sizeof(uint64_t) == sizeof(unsigned long))
    return __builtin_clzl(a);
  else if (sizeof(uint64_t) == sizeof(unsigned long long))
    return __builtin_clzll(a);
  else
    assert(sizeof(uint64_t) == sizeof(unsigned long) || sizeof(uint64_t) == sizeof(unsigned long long));
#else
  uint64_t t = 1UL << (64 - 1);
  int i;
  ASSERT_EXPENSIVE (a != 0UL);
  for (i = 0; (a & t) == 0UL; i++)
    t >>= 1;
  return i;
#endif
}


/* Compute 1/n (mod 2^wordsize) */
static inline uint64_t
u64arith_invmod (const uint64_t n)
{
  uint64_t r;

  ASSERT (n % 2UL != 0UL);
  
  /* Suggestion from PLM: initing the inverse to (3*n) XOR 2 gives the
     correct inverse modulo 32, then 3 (for 32 bit) or 4 (for 64 bit) 
     Newton iterations are enough. */
  r = (3UL * n) ^ 2UL;
  /* Newton iteration */
  r = 2UL * r - (uint32_t) r * (uint32_t) r * (uint32_t) n;
  r = 2UL * r - (uint32_t) r * (uint32_t) r * (uint32_t) n;
#if 0
  r = 2UL * r - r * r * n;
#else
  /*
    r*n = k*2^32 + 1
        
    r' = 2 * r - r * r * n
    r' = 2 * r - r * (k*2^32 + 1)
    r' = 2 * r - r * k*2^32 - r
    r' = r - r * k*2^32
    r' = r - ((r * k) % 2^32) * 2^32
  */
  uint32_t k = (uint32_t)(r * n >> 32);
  k *= (uint32_t) r;
  r = r - ((uint64_t)k << 32);
#endif

  return r;
}

/* Compute n/2 (mod m), where m must be odd. */
static inline uint64_t
u64arith_div2mod (const uint64_t n, const uint64_t m)
{
#if !defined (U64ARITH_NO_ASM) && defined(HAVE_GCC_STYLE_AMD64_INLINE_ASM)
  uint64_t s = n, t = m;
  ASSERT_EXPENSIVE (m % 2UL != 0UL);

  __asm__ __VOLATILE(
	  "add %1, %0\n\t"
	  "rcr $1, %0\n\t" /* FIXME: rcr is SLOW! */
	  "shr $1, %1\n\t"
	  "cmovnc %1, %0\n"
	  : "+&r" (t), "+&r" (s)
	  : : "cc"
	  );
  return t;
#else
  ASSERT_EXPENSIVE (m % 2UL != 0UL);
  if (n % 2UL == 0UL)
    return n / 2UL;
  else
    return n / 2UL + m / 2UL + 1UL;
#endif
}


/* Integer (truncated) square root of n */
static inline uint64_t
u64arith_sqrt (const uint64_t n)
{
  unsigned int i;
  uint64_t xs, c, d, s2;
  const unsigned int l = 63 - (unsigned int)__builtin_clzl(n);

  d = n; /* d = n - x^2 */
  xs = 0UL;
  s2 = 1UL << (l - l % 2);

  for (i = l / 2; i != 0; i--)
    {
      /* Here, s2 = 1 << (2*i) */
      /* xs = x << (i + 1), the value of x shifted left i+1 bits */

      c = xs + s2; /* c = (x + 2^i) ^ 2 - x^2 = 2^(i+1) * x + 2^(2*i) */
      xs >>= 1; /* Now xs is shifted only i positions */
      if (d >= c)
        {
          d -= c;
          xs |= s2; /* x |= 1UL << i <=> xs |= 1UL << (2*i) */
        }
      s2 >>= 2;
    }

  c = xs + s2;
  xs >>= 1;   
  if (d >= c) 
    xs |= s2;
  
  return xs;
}

/* Given r = -rem/p (mod den), we want num/(den*2^k) (mod p) ==
   (ratio + rem/den)/2^k (mod p).
   Using (a variant of) Bezout's identity, we have, for some non-negative
   integer t,
   r * p - t * den = -rem, or
   r * p + rem = t * den,
   thus den | (r * p + rem), and thus
   t = (r * p + rem) / den is an integer and satisfies
   t = rem/den (mod p).

   We have 0 <= r <= den-1 and rem <= den-1, and thus
   0 <= t = p * r/den + rem/den <=
   p (1 - 1/den) + 1 - 1/den =
   p + 1 - (p + 1)/den < p + 1.
   Thus t is almost a properly reduced residue for rem/den (mod p).
   As p fits in uint64_t, so does t, and we can compute t modulo
   2^64; since den is odd, we can multiply by den^{-1} mod 2^64
   to effect division by den.

   Finally we compute (t + ratio)/2^k mod p = num/(den*2^k) mod p.  */

static inline uint64_t
u64arith_post_process_inverse(const uint64_t r, const uint64_t p,
  const uint64_t rem, const uint64_t den_inv,
  const uint64_t ratio, const uint64_t k)
{
  uint64_t t = (r * p + rem) * den_inv;
  const uint64_t ratio_p = (ratio >= p) ? ratio % p : ratio;
  ASSERT_ALWAYS(t <= p); /* Cheap and fairly strong test */
  /* u64arith_addmod_1_1() accepts third operand == p and still produces
     a properly reduced sum mod p. */
  u64arith_addmod_1_1 (&t, ratio_p, t, p);

  ASSERT_EXPENSIVE(t < p);
  ASSERT_EXPENSIVE(k == 0 || p % 2 == 1);
  for (uint64_t j = 0; j < k; j++) {
    t = u64arith_div2mod(t, p);
  }
  return t;
}


/* Computes r = ((phigh * 2^64 + plow) / 2^64) % m
   Requires phigh < m and invm = -1/m (mod 2^64). */

static inline void
u64arith_redc(uint64_t *r, const uint64_t plow,
             const uint64_t phigh, const uint64_t m,
             const uint64_t invm)
{
  uint64_t t = phigh;

  ASSERT_EXPENSIVE (phigh < m);

#if !defined (U64ARITH_NO_ASM) && defined(HAVE_GCC_STYLE_AMD64_INLINE_ASM)

  /* TODO: are the register constraints watertight?
     %rax gets modified but putting tlow as an output constraint with "+"
     will keep r from getting allocated in %rax, which is a shame
     since we often want the result in %rax for the next multiply. */

  __asm__ __VOLATILE (
    "imulq %[invm], %%rax\n\t"
    "cmpq $1, %%rax \n\t"                /* if plow != 0, increase t */
    "sbbq $-1, %[t]\n\t"
    "mulq %[m]\n\t"
    "lea (%[t],%%rdx,1), %[r]\n\t"  /* compute (rdx + thigh) mod m */
    "subq %[m], %[t]\n\t"
    "addq %%rdx, %[t]\n\t"
    "cmovcq %[t], %[r]\n\t"
    : [t] "+&r" (t), [r] "=&r" (r[0])
    : [invm] "rm" (invm), [m] "rm" (m), "a" (plow)
    : "%rdx", "cc"
  );
#else
  uint64_t tlow, thigh;
  tlow = plow * invm;
  u64arith_mul_1_1_2 (&tlow, &thigh, tlow, m);
  /* Let w = 2^wordsize. We know (phigh * w + plow) + (thigh * w + tlow)
     == 0 (mod w) so either plow == tlow == 0, or plow !=0 and tlow != 0.
     In the former case we want phigh + thigh + 1, in the latter
     phigh + thigh. Since t = phigh < m, and modredcul_add can handle the
     case where the second operand is equal to m, adding 1 is safe */

  t += (plow != 0UL) ? 1UL : 0UL; /* Does not depend on the mul */

  u64arith_addmod_1_1(r, t, thigh, m);
#endif
  ASSERT_EXPENSIVE (r[0] < m);
}


#endif /* ifndef U64_ARITH_H__ */
