/* Some commonly used assembly helper macros for arithmetic on 
   unsigned long. 
   Defining ULARITH_VERBOSE_ASM puts in the asm output a line for each
   function call that shows what registers/memory locations the operands 
   are in.
   Defining ULARITH_NO_ASM avoids asm macros and uses the C fallback code
   where available.
*/

#ifndef UL_ARITH_H__

#define UL_ARITH_H__

#include <assert.h>

/* <limits.h> defines LONG_BIT only with _XOPEN_SOURCE defined, but if 
   another header (such as <stdio.h>) already included <features.h> before 
   _XOPEN_SOURCE was set to 1, future includes of <features.h> are
   short-circuited and _XOPEN_SOURCE is ignored. */

#ifndef CHAR_BIT
#define CHAR_BIT 8
#endif
#ifndef LONG_BIT
#define LONG_BIT	((int) sizeof(long) * CHAR_BIT)
#endif

#ifndef ASSERT
#define ASSERT(x)	assert(x)
#endif

#ifndef	MAYBE_UNUSED
#if defined(__GNUC__)
#define MAYBE_UNUSED __attribute__ ((unused))
#else
#define MAYBE_UNUSED
#endif
#endif

/* Add an unsigned long to two unsigned longs with carry propagation from 
   low word (r1) to high word (r2). Any carry out from high word is lost. */

static inline void
ularith_add_ul_2ul (unsigned long *r1, unsigned long *r2, 
			const unsigned long a)
{
#ifdef ULARITH_VERBOSE_ASM
  __asm__ ("# ularith_add_ul_2ul (%0, %1, %2)\n" : : 
           "X" (*r1), "X" (*r2), "X" (a));
#endif
#if !defined (ULARITH_NO_ASM) && defined(__x86_64__) && defined(__GNUC__)
  __asm__ ( "addq %2, %0\n\t"
            "adcq $0, %1\n"
            : "+&r" (*r1), "+r" (*r2) 
            : "g" (a)
            : "cc"); /* TODO: add commutativity and alternative for add to 
                        memory */
#elif !defined (ULARITH_NO_ASM) && defined(__i386__) && defined(__GNUC__)
  __asm__ ( "addl %2, %0\n\t"
            "adcl $0, %1\n"
            : "+&r" (*r1), "+r" (*r2)
            : "g" (a)
            : "cc");
#else
  *r1 += a;
  if (*r1 < a)
    (*r2)++;
#endif
}


/* Add two unsigned longs to two unsigned longs with carry propagation from 
   low word (r1) to high word (r2). Any carry out from high word is lost. */

static inline void
ularith_add_2ul_2ul (unsigned long *r1, unsigned long *r2, 
			 const unsigned long a1, const unsigned long a2)
{
#ifdef ULARITH_VERBOSE_ASM
  __asm__ ("# ularith_add_2ul_2ul (%0, %1, %2, %3)\n" : : 
           "X" (*r1), "X" (*r2), "X" (a1), "X" (a2));
#endif
#if !defined (ULARITH_NO_ASM) && defined(__x86_64__) && defined(__GNUC__)
  __asm__ ( "addq %2, %0\n\t"
            "adcq %3, %1\n"
            : "+&r" (*r1), "+r" (*r2)
            : "g" (a1), "g" (a2)
            : "cc");
#elif !defined (ULARITH_NO_ASM) && defined(__i386__) && defined(__GNUC__)
  __asm__ ( "addl %2, %0\n\t"
            "adcl %3, %1\n"
            : "+&r" (*r1), "+r" (*r2)
            : "g" (a1), "g" (a2)
            : "cc");
#else
  *r1 += a1;
  *r2 += a2 + (*r1 < a1);
#endif
}


/* Subtract an unsigned long from two unsigned longs with borrow propagation 
   from low word (r1) to high word (r2). Any borrow out from high word is 
   lost. */

static inline void
ularith_sub_ul_2ul (unsigned long *r1, unsigned long *r2, 
			const unsigned long a)
{
#ifdef ULARITH_VERBOSE_ASM
  __asm__ ("# ularith_sub_ul_2ul  (%0, %1, %2)\n" : : 
           "X" (*r1), "X" (*r2), "X" (a));
#endif
#if !defined (ULARITH_NO_ASM) && defined(__x86_64__) && defined(__GNUC__)
  __asm__ ( "subq %2, %0\n\t"
            "sbbq $0, %1\n"
            : "+&r" (*r1), "+r" (*r2)
            : "g" (a)
            : "cc");
#elif !defined (ULARITH_NO_ASM) && defined(__i386__) && defined(__GNUC__)
  __asm__ ( "subl %2, %0\n\t"
            "sbbl $0, %1\n"
            : "+&r" (*r1), "+r" (*r2)
            : "g" (a)
            : "cc");
#else
  unsigned long u = *r1;
  *r1 -= a;
  if (*r1 > u)
    (*r2)--;
#endif
}


/* Subtract two unsigned longs from two unsigned longs with borrow propagation 
   from low word (r1) to high word (r2). Any borrow out from high word is 
   lost. */

static inline void
ularith_sub_2ul_2ul (unsigned long *r1, unsigned long *r2, 
			 const unsigned long a1, const unsigned long a2)
{
#ifdef ULARITH_VERBOSE_ASM
  __asm__ ("# ularith_sub_2ul_2ul (%0, %1, %2, %3)\n" : : 
           "X" (*r1), "X" (*r2), "X" (a1), "X" (a2));
#endif
#if !defined (ULARITH_NO_ASM) && defined(__x86_64__) && defined(__GNUC__)
  __asm__ ( "subq %2, %0\n\t"
            "sbbq %3, %1\n"
            : "+&r" (*r1), "+r" (*r2)
            : "g" (a1), "g" (a2)
            : "cc");
#elif !defined (ULARITH_NO_ASM) && defined(__i386__) && defined(__GNUC__)
  __asm__ ( "subl %2, %0\n\t"
            "sbbl %3, %1\n"
            : "+&r" (*r1), "+r" (*r2)
            : "g" (a1), "g" (a2)
            : "cc");
#else
  unsigned long u = *r1;
  *r1 -= a1;
  *r2 -= a2;
  if (*r1 > u)
    (*r2)--;
#endif
}


/* Subtract only if result is non-negative */

static inline void
ularith_sub_2ul_2ul_ge (unsigned long *r1, unsigned long *r2, 
			const unsigned long a1, const unsigned long a2)
{
  unsigned long t1 = *r1, t2 = *r2;
#ifdef ULARITH_VERBOSE_ASM
  __asm__ ("# ularith_sub_2ul_2ul_ge (%0, %1, %2, %3)\n" : : 
           "X" (*r1), "X" (*r2), "X" (a1), "X" (a2));
#endif
#if !defined (ULARITH_NO_ASM) && defined(__x86_64__) && defined(__GNUC__)
  __asm__ ( "subq %4, %0\n\t" /* r1 -= a1 */
            "sbbq %5, %1\n\t" /* r2 -= a2 + cy */
	    "cmovc %2, %0\n\t" /* If there's a borrow, restore r1 from t1 */
	    "cmovc %3, %1\n\t" /* and r2 from t2 */
            : "+&r" (*r1), "+&r" (*r2)
            : "r" (t1), "r" (t2), "g" (a1), "g" (a2)
            : "cc");
#elif !defined (ULARITH_NO_ASM) && defined(__i386__) && defined(__GNUC__)
  __asm__ ( "subl %4, %0\n\t"
            "sbbl %5, %1\n\t"
	    "cmovc %2, %0\n\t"
	    "cmovc %3, %1\n\t"
            : "+&r" (*r1), "+&r" (*r2)
            : "r" (t1), "r" (t2), "g" (a1), "g" (a2)
            : "cc");
#else
  t1 -= a1;
  t2 -= a2 + (t1 > *r1);
  if (t2 <= *r2)
    {
      *r1 = t1;
      *r2 = t2;
    }
#endif
}


/* Multiply two unsigned long "a" and "b" and put the result as 
   r2:r1 (r2 being the high word) */

static inline void
ularith_mul_ul_ul_2ul (unsigned long *r1, unsigned long *r2, 
			   const unsigned long a, const unsigned long b)
{
#ifdef ULARITH_VERBOSE_ASM
  __asm__ ("# ularith_mul_ul_ul_2ul (%0, %1, %2, %3)\n" : : 
           "X" (*r1), "X" (*r2), "X" (a), "X" (b));
#endif
#if defined(__x86_64__) && defined(__GNUC__)
  __asm__ ( "mulq %3"
	    : "=a" (*r1), "=d" (*r2)
	    : "%0" (a), "rm" (b)
	    : "cc");
#elif defined(__i386__) && defined(__GNUC__)
  __asm__ ( "mull %3"
	    : "=a" (*r1), "=d" (*r2)
	    : "%0" (a), "rm" (b)
	    : "cc");
#else
  abort ();
#endif
}


/* Integer division of a two ulong value a2:a1 by a ulong divisor. Returns
   quotient and remainder. */

static inline void
ularith_div_2ul_ul_ul (unsigned long *q, unsigned long *r, 
			   const unsigned long a1, const unsigned long a2, 
			   const unsigned long b)
{
  ASSERT(a2 < b); /* Or there will be quotient overflow */
#ifdef ULARITH_VERBOSE_ASM
  __asm__ ("# ularith_div_2ul_ul_ul (%0, %1, %2, %3, %4)\n" : : 
           "X" (*q), "X" (*r), "X" (a1), "X" (a2), "X" (b));
#endif
#if defined(__x86_64__) && defined(__GNUC__)
  __asm__ ( "divq %4"
            : "=a" (*q), "=d" (*r)
            : "0" (a1), "1" (a2), "rm" (b)
            : "cc");
#elif defined(__i386__) && defined(__GNUC__)
  __asm__ ( "divl %4"
            : "=a" (*q), "=d" (*r)
            : "0" (a1), "1" (a2), "rm" (b)
            : "cc");
#else
  abort ();
#endif
}


/* Integer division of a two longint value by a longint divisor. Returns
   only remainder. */

static inline void
ularith_div_2ul_ul_ul_r (unsigned long *r, unsigned long a1,
                 const unsigned long a2, const unsigned long b)
{
  ASSERT(a2 < b); /* Or there will be quotient overflow */
#ifdef ULARITH_VERBOSE_ASM
  __asm__ ("# ularith_div_2ul_ul_ul_r (%0, %1, %2, %3)\n" : : 
           "X" (*r), "X" (a1), "X" (a2), "X" (b));
#endif
#if defined(__x86_64__) && defined(__GNUC__)
  __asm__ ( "divq %3"
            : "+a" (a1), "=d" (*r)
            : "1" (a2), "rm" (b)
            : "cc");
#elif defined(__i386__) && defined(__GNUC__)
  __asm__ ( "divl %3"
            : "+a" (a1), "=d" (*r)
            : "1" (a2), "rm" (b)
            : "cc");
#else
  abort ();
#endif
}


/* Shift *r right by i bits, filling in the low bits from a into the high
   bits of *r */
MAYBE_UNUSED
static inline void
ularith_shrd (unsigned long *r, const unsigned long a, const int i)
{
#ifdef ULARITH_VERBOSE_ASM
  __asm__ ("# ularith_shrd (%0, %1, %2)\n" : : 
           "X" (*r), "X" (a), "X" (i));
#endif
#if !defined (ULARITH_NO_ASM) && defined(__x86_64__) && defined(__GNUC__)
  __asm__ ("shrdq %b2, %1, %0\n": 
  /* the b modifier makes gcc print the byte part of the register (%cl) */
           "+rm" (*r) :
           "r" (a), "cJ" (i) : /* i can be in %cx or a constant < 64 */
           "cc"
          );
#elif !defined (ULARITH_NO_ASM) && defined(__i386__) && defined(__GNUC__)
  __asm__ ("shrdl %b2, %1, %0\n": 
           "+rm" (*r) :
           "r" (a), "cI" (i) : /* i can be in %cx or a constant < 32 */
           "cc"
          );
#else
  *r = (*r >> i) | (a << (LONG_BIT - i));
#endif
}

/* Shift *r left by i bits, filling in the high bits from a into the low
   bits of *r */
MAYBE_UNUSED
static inline void
ularith_shld (unsigned long *r, const unsigned long a, const int i)
{
#ifdef ULARITH_VERBOSE_ASM
  __asm__ ("# ularith_shld (%0, %1, %2)\n" : : 
           "X" (*r), "X" (a), "X" (i));
#endif
#if !defined (ULARITH_NO_ASM) && defined(__x86_64__) && defined(__GNUC__)
  __asm__ ("shldq %b2, %1, %0\n": 
  /* the b modifier makes gcc print the byte part of the register (%cl) */
           "+rm" (*r) :
           "r" (a), "cJ" (i) : /* i can be in %cx or a constant < 64 */
           "cc"
          );
#elif !defined (ULARITH_NO_ASM) && defined(__i386__) && defined(__GNUC__)
  __asm__ ("shldl %b2, %1, %0\n": 
           "+rm" (*r) :
           "r" (a), "cI" (i) : /* i can be in %cx or a constant < 32 */
           "cc"
          );
#else
  *r = (*r << i) | (a >> (LONG_BIT - i));
#endif
}

/* Returns number of trailing zeros in a. a must not be zero */
MAYBE_UNUSED
static inline int
ularith_ctz (const unsigned long a)
{
#if defined(__GNUC__) && (__GNUC__ >= 4 || __GNUC__ >= 3 && __GNUC_MINOR__ >= 4)
  ASSERT (a != 0UL);
  return __builtin_ctzl(a);
#else
  unsigned long t = a;
  int i;
  ASSERT (a != 0UL);
  for (i = 0; (t & 1UL) == 0UL; i++)
    t >>= 1;
  return i;
#endif
}

/* Returns number of leading zeros in a. a must not be zero */
MAYBE_UNUSED
static inline int
ularith_clz (const unsigned long a)
{
#if defined(__GNUC__) && (__GNUC__ >= 4 || __GNUC__ >= 3 && __GNUC_MINOR__ >= 4)
  ASSERT (a != 0UL);
  return __builtin_clzl(a);
#else
  unsigned long t = 1UL << (LONG_BIT - 1);
  int i;
  ASSERT (a != 0UL);
  for (i = 0; (a & t) == 0UL; i++)
    t >>= 1;
  return i;
#endif
}

#endif
