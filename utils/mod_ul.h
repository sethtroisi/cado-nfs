/* Some functions for modular arithmetic with residues and modulus in
   unsigned long variables. Due to inlining, this file must be included
   in the caller's source code with #include */

/* Naming convention: all function start with modul, for 
   MODulus Unsigned Long, followed by underscore, functionality of function
  (add, mul, etc), and possibly underscore and specification of what argument
  types the function takes (_ul, etc).
  There are typedef's that rename all functions to mod_* instead of 
  modul_*, which one day might become an automatic renaming scheme so that
  different modulus sizes can be used simply by #including different mod_*.c
  files, but without changing anything else in the source code. */

#ifndef MOD_UL_H__

#define MOD_UL_H__

/**********************************************************************/
#include <assert.h>

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

#ifdef  mod_init
#warning "mod_ul.h included after modredc_ul.h ; bindings overwritten"
#endif

/* Define the default names and types for arithmetic with these functions */
#undef mod_init
#undef mod_init_noset0
#undef mod_clear
#undef mod_set
#undef mod_set_ul
#undef mod_set_ul_reduced
#undef mod_swap
#undef mod_initmod_ul
#undef mod_getmod_ul
#undef mod_clearmod
#undef mod_get_ul
#undef mod_equal
#undef mod_is0
#undef mod_add
#undef mod_add_ul
#undef mod_sub
#undef mod_sub_ul
#undef mod_neg
#undef mod_mul
#undef mod_addredc_ul
#undef mod_mulredc_ul
#undef mod_muladdredc_ul
#undef mod_div2
#undef mod_div3
#undef mod_pow_ul
#undef mod_pow_mp
#undef mod_2pow_mp
#undef mod_V_ul
#undef mod_V_mp
#undef mod_sprp
#undef mod_gcd
#undef mod_inv
#undef mod_jacobi
#undef mod_set0
#undef mod_set1
#undef mod_next
#undef mod_finished
#undef residue_t
#undef modulus_t

#define mod_init             modul_init
#define mod_init_noset0      modul_init_noset0
#define mod_clear            modul_clear
#define mod_set              modul_set
#define mod_set_ul           modul_set_ul
#define mod_set_ul_reduced   modul_set_ul_reduced
#define mod_swap             modul_swap
#define mod_initmod_ul       modul_initmod_ul
#define mod_getmod_ul        modul_getmod_ul
#define mod_clearmod         modul_clearmod
#define mod_get_ul           modul_get_ul
#define mod_equal            modul_equal
#define mod_is0              modul_is0
#define mod_add              modul_add
#define mod_add_ul           modul_add_ul
#define mod_sub              modul_sub
#define mod_sub_ul           modul_sub_ul
#define mod_neg              modul_neg
#define mod_mul              modul_mul
#define mod_tomontgomery     modul_tomontgomery
#define mod_frommontgomery   modul_frommontgomery
#define mod_addredc_ul       modul_addredc_ul
#define mod_mulredc          modul_mulredc
#define mod_mulredc_ul       modul_mulredc_ul
#define mod_muladdredc_ul    modul_muladdredc_ul
#define mod_div2             modul_div2
#define mod_div3             modul_div3
#define mod_invmodlong       modul_invmodlong
#define mod_powredc_ul       modul_powredc_ul
#define mod_powredc_mp       modul_powredc_mp
#define mod_2powredc_mp      modul_2powredc_mp
#define mod_Vredc_ul         modul_Vredc_ul
#define mod_Vredc_mp         modul_Vredc_mp
#define mod_sprp             modul_sprp
#define mod_gcd              modul_gcd
#define mod_inv              modul_inv
#define mod_jacobi           modul_jacobi
#define mod_set0          modul_set0
#define mod_set1          modul_set1
#define mod_next          modul_next
#define mod_finished          modul_finished

typedef unsigned long residueul_t[1];
typedef unsigned long modulusul_t[1];
typedef residueul_t residue_t;
typedef modulusul_t modulus_t;

/* Initialises a residue_t type and sets it to zero */
MAYBE_UNUSED
static inline void
modul_init (residueul_t r, const modulusul_t m MAYBE_UNUSED)
{
  r[0] = 0UL;
  return;
}


/* Initialises a residue_t type, but does not set it to zero. For fixed length
   residue_t types, that leaves nothing to do at all. */
MAYBE_UNUSED
static inline void
modul_init_noset0 (residueul_t r MAYBE_UNUSED, 
                   const modulusul_t m MAYBE_UNUSED)
{
  return;
}


MAYBE_UNUSED
static inline void
modul_clear (residueul_t r MAYBE_UNUSED, 
             const modulusul_t m MAYBE_UNUSED)
{
  return;
}


MAYBE_UNUSED
static inline void
modul_set (residueul_t r, const residueul_t s, const 
	   modulusul_t m MAYBE_UNUSED)
{
  ASSERT_EXPENSIVE (s[0] < m[0]);
  r[0] = s[0];
}


MAYBE_UNUSED
static inline void
modul_set_ul (residueul_t r, const unsigned long s, const modulusul_t m)
{
  r[0] = s % m[0];
}

/* Sets the residue_t to the class represented by the integer s. Assumes that
   s is reduced (mod m), i.e. 0 <= s < m */

MAYBE_UNUSED
static inline void
modul_set_ul_reduced (residueul_t r, const unsigned long s, 
                      const modulusul_t m MAYBE_UNUSED)
{
  ASSERT (s < m[0]);
  r[0] = s;
}

/* These two are so trivial that we don't really require m in the
 * interface. For 1 we might, as the internal representation might 
 * not use "1" for 1 (e.g. when using Montgomery's REDC.)
 * For interface homogeneity we make even modul_set0 take the m parameter.
 */
MAYBE_UNUSED 
static inline void 
modul_set0 (residueul_t r, MAYBE_UNUSED const modulusul_t m) 
{ 
  r[0] = 0UL; 
}


MAYBE_UNUSED 
static inline void 
modul_set1 (residueul_t r, MAYBE_UNUSED const modulusul_t m) 
{ 
  r[0] = 1UL; 
}


/* Exchanges the values of the two arguments */

MAYBE_UNUSED
static inline void
modul_swap (residueul_t a, residueul_t b, 
            const modulusul_t m MAYBE_UNUSED)
{
  unsigned long t;
  ASSERT_EXPENSIVE (a[0] < m[0] && b[0] < m[0]);
  t = a[0];
  a[0] = b[0];
  b[0] = t;
}
                          

MAYBE_UNUSED
static inline void
modul_initmod_ul (modulusul_t r, const unsigned long s)
{
  r[0] = s;
}

MAYBE_UNUSED
static inline unsigned long
modul_getmod_ul (const modulusul_t m)
{
  return m[0];
}

MAYBE_UNUSED
static inline void
modul_clearmod (modulusul_t m MAYBE_UNUSED)
{
  return;
}

MAYBE_UNUSED
static inline unsigned long
modul_get_ul (const residueul_t s, 
	      const modulusul_t m MAYBE_UNUSED)
{
  ASSERT_EXPENSIVE (s[0] < m[0]);
  return s[0];
}

MAYBE_UNUSED
static inline int
modul_equal (const residueul_t a, const residueul_t b, 
             const modulusul_t m MAYBE_UNUSED)
{
  ASSERT_EXPENSIVE (a[0] < m[0] && b[0] < m[0]);
  return (a[0] == b[0]);
}

MAYBE_UNUSED
static inline int
modul_is0 (const residueul_t a, const modulusul_t m MAYBE_UNUSED)
{
  ASSERT_EXPENSIVE (a[0] < m[0]);
  return (a[0] == 0UL);
}

MAYBE_UNUSED
static inline int
modul_is1 (const residueul_t a, const modulusul_t m MAYBE_UNUSED)
{
  ASSERT_EXPENSIVE (a[0] < m[0]);
  return (a[0] == 1UL);
}

MAYBE_UNUSED
static inline void
modul_add (residueul_t r, const residueul_t a, const residueul_t b, 
	   const modulusul_t m)
{
  ASSERT_EXPENSIVE (a[0] < m[0] && b[0] < m[0]);
#ifdef MODTRACE
  printf ("modul_add: a = %lu, b = %lu", a[0], b[0]);
#endif

#if (defined(__i386__) || defined(__x86_64__)) && defined(__GNUC__)
  {
    unsigned long t = a[0] - m[0], tr = a[0] + b[0];
    
    __asm__ (
      "add %2, %1\n\t"   /* t += b */
      "cmovc %1, %0\n\t"  /* if (cy) tr = t */
      : "+r" (tr), "+&r" (t)
      : "g" (b[0])
      : "cc"
    );
    ASSERT_EXPENSIVE (tr == (a[0] >= m[0] - b[0]) ? (a[0] - (m[0] - b[0])) : (a[0] + b[0]));
    r[0] = tr;
  }
#else
  r[0] = (a[0] >= m[0] - b[0]) ? (a[0] - (m[0] - b[0])) : (a[0] + b[0]);
#endif

#ifdef MODTRACE
  printf (", r = %lu\n", r[0]);
#endif
  ASSERT_EXPENSIVE (r[0] < m[0]);
}

/* FIXME: This function is really modul_add_ul_reduced */
MAYBE_UNUSED
static inline void
modul_add_ul (residueul_t r, const residueul_t a, const unsigned long b, 
	      const modulusul_t m)
{
  ASSERT_EXPENSIVE (a[0] < m[0] && b < m[0]);
#ifdef MODTRACE
  printf ("modul_add_ul: a = %lu, b = %lu", a[0], b);
#endif
  r[0] = (a[0] >= m[0] - b) ? (a[0] - (m[0] - b)) : (a[0] + b);
#ifdef MODTRACE
  printf (", r = %lu\n", r[0]);
#endif
}

MAYBE_UNUSED
static inline void
modul_sub (residueul_t r, const residueul_t a, const residueul_t b, 
	   const modulusul_t m)
{
  ASSERT_EXPENSIVE (a[0] < m[0] && b[0] < m[0]);
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
      : "g" (b[0]), "rm" (m[0])
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
      t += m[0];
    r[0] = t;
  }
#else
  r[0] = (a[0] < b[0]) ? (a[0] - b[0] + m[0]) : (a[0] - b[0]);
#endif

#ifdef MODTRACE
  printf (", r = %lu\n", r[0]);
#endif
  ASSERT_EXPENSIVE (r[0] < m[0]);
}


/* FIXME: This function is really modul_sub_ul_reduced */
MAYBE_UNUSED
static inline void
modul_sub_ul (residueul_t r, const residueul_t a, const unsigned long b, 
	      const modulusul_t m)
{
  ASSERT_EXPENSIVE (a[0] < m[0]);
  ASSERT(b < m[0]);
#ifdef MODTRACE
  printf ("modul_sub_ul: a = %lu, b = %lu", a[0], b);
#endif
  r[0] = (a[0] >= b) ? (a[0] - b) : (a[0] + (m[0] - b));
#ifdef MODTRACE
  printf (", r = %lu\n", r[0]);
#endif
}

MAYBE_UNUSED
static inline void
modul_neg (residueul_t r, const residueul_t a, const modulusul_t m)
{
  ASSERT_EXPENSIVE (a[0] < m[0]);
  if (a[0] == 0UL)
    r[0] = a[0];
  else
    r[0] = m[0] - a[0];
}

/* Add an unsigned long to two unsigned longs with carry propagation from 
   low word to high word. Any carry out from high word is lost. */

static inline void
modul_add_ul_2ul (unsigned long *r1, unsigned long *r2, const unsigned long a)
{
#if defined(__x86_64__) && defined(__GNUC__)
  __asm__ ( "addq %2, %0\n\t"
            "adcq $0, %1\n"
            : "+&r" (*r1), "+r" (*r2)
            : "rm" (a)
            : "cc");
#elif defined(__i386__) && defined(__GNUC__)
  __asm__ ( "addl %2, %0\n\t"
            "adcl $0, %1\n"
            : "+&r" (*r1), "+r" (*r2)
            : "rm" (a)
            : "cc");
#else
  abort ();
#endif
}

/* Add two unsigned longs to two unsigned longs with carry propagation from 
   low word to high word. Any carry out from high word is lost. */

static inline void
modul_add_2ul_2ul (unsigned long *r1, unsigned long *r2, 
                   const unsigned long a1, const unsigned long a2)
{
#if defined(__x86_64__) && defined(__GNUC__)
  __asm__ ( "addq %2, %0\n\t"
            "adcq %3, %1\n"
            : "+&r" (*r1), "+r" (*r2)
            : "rm" (a1), "rm" (a2)
            : "cc");
#elif defined(__i386__) && defined(__GNUC__)
  __asm__ ( "addl %2, %0\n\t"
            "adcl %3, %1\n"
            : "+&r" (*r1), "+r" (*r2)
            : "rm" (a1), "rm" (a2)
            : "cc");
#else
  abort ();
#endif
}


/* Multiply two unsigned long "a" and "b" and put the result as 
   r2:r1 (r2 being the high word) */

static inline void
mul_ul_ul_2ul (unsigned long *r1, unsigned long *r2, const unsigned long a,
               const unsigned long b)
{
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
div_2ul_ul_ul (unsigned long *q, unsigned long *r, const unsigned long a1,
               const unsigned long a2, const unsigned long b)
{
  ASSERT(a2 < b); /* Or there will be quotient overflow */
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
div_2ul_ul_ul_r (unsigned long *r, unsigned long a1,
                 const unsigned long a2, const unsigned long b)
{
  ASSERT(a2 < b); /* Or there will be quotient overflow */
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

MAYBE_UNUSED
static inline void
modul_mul (residueul_t r, const residueul_t a, const residueul_t b, 
           const modulusul_t m)
{
  unsigned long t1, t2, dummy;
#if defined(MODTRACE)
  printf ("mulmod_ul: a = %lu, b = %lu", a[0], b[0]);
#endif
  
  ASSERT_EXPENSIVE (a[0] < m[0] && b[0] < m[0]);
  mul_ul_ul_2ul (&t1, &t2, a[0], b[0]);
  div_2ul_ul_ul (&dummy, r, t1, t2, m[0]);
  
#if defined(MODTRACE)
  printf (", r = %lu\n", _r);
#endif
}

/* Computes (a * 2^wordsize) % m */
MAYBE_UNUSED
static inline void
modul_tomontgomery (residueul_t r, const residueul_t a, const modulusul_t m)
{
  ASSERT_EXPENSIVE (a[0] < m[0]);
  div_2ul_ul_ul_r (r, 0UL, a[0], m[0]);
}


/* Computes (a / 2^wordsize) % m */
MAYBE_UNUSED
static inline void
modul_frommontgomery (residueul_t r, const residueul_t a, 
                      const unsigned long invm, const modulusul_t m)
{
  unsigned long tlow, thigh;
  mul_ul_ul_2ul (&tlow, &thigh, a[0], invm);
  mul_ul_ul_2ul (&tlow, &thigh, tlow, m[0]);
  r[0] = thigh + (a[0] != 0UL ? 1UL : 0UL);
}


/* Computes (a / 2^wordsize) % m, but result can be r = m. 
   Input a must not be equal 0 */
MAYBE_UNUSED
static inline void
modul_redcsemi_ul_not0 (residueul_t r, const unsigned long a, 
                        const unsigned long invm, const modulusul_t m)
{
  unsigned long tlow, thigh;

  ASSERT (a != 0);

  tlow = a * invm; /* tlow <= 2^w-1 */
  mul_ul_ul_2ul (&tlow, &thigh, tlow, m[0]);
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
modul_addredc_ul (residueul_t r, const residueul_t a, const unsigned long b, 
                  const unsigned long invm, const modulusul_t m)
{
  unsigned long slow, shigh, tlow, thigh;
  
  ASSERT_EXPENSIVE (a[0] <= m[0]);
  slow = b;
  shigh = 0UL;
  modul_add_ul_2ul (&slow, &shigh, a[0]);
  
  tlow = slow * invm;
  mul_ul_ul_2ul (&tlow, &thigh, tlow, m[0]);
  ASSERT_EXPENSIVE (slow + tlow == 0UL);
  r[0] = thigh + shigh + (slow != 0UL ? 1UL : 0UL);
  
  /* r = ((a+b) + (((a+b)%2^w * invm) % 2^w) * m) / 2^w  Use a<=m-1, b<=2^w-1
     r <= (m + 2^w - 1 + (2^w - 1) * m) / 2^w
        = (m - 1 + 2^w + m*2^w - m) / 2^w
        = (- 1 + 2^w + m2^w) / 2^w
        = m + 1 - 1/2^w
     r <= m, since r is an integer
  */
  if (r[0] == m[0])
    r[0] = 0UL;
}


/* Computes ((a + b) / 2^wordsize) % m, but result can be == m.
   a <= m is permissible */
MAYBE_UNUSED
static inline void
modul_addredcsemi_ul (residueul_t r, const residueul_t a, 
                      const unsigned long b, const unsigned long invm, 
                      const modulusul_t m)
{
  unsigned long slow, shigh, tlow;
  unsigned char sb;
  
  ASSERT_EXPENSIVE(a[0] <= m[0]);
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
  modul_add_ul_2ul (&slow, &shigh, a[0]);
  shigh += (slow != 0UL ? 1UL : 0UL);
#endif

  tlow = slow * invm;
  mul_ul_ul_2ul (&tlow, r, tlow, m[0]);
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
modul_mulredc (residueul_t r, const residueul_t a, const residueul_t b,
               const unsigned long invm, const modulusul_t m)
{
  unsigned long plow, phigh, tlow, thigh;

  ASSERT_EXPENSIVE (m[0] % 2 != 0);
  ASSERT_EXPENSIVE (a[0] < m[0] && b[0] < m[0]);
#if defined(MODTRACE)
  printf ("(%lu * %lu / 2^%ld) %% %lu", 
          a[0], b[0], 8 * sizeof(unsigned long), m[0]);
#endif
  
  mul_ul_ul_2ul (&plow, &phigh, a[0], b[0]);
  tlow = plow * invm;
  mul_ul_ul_2ul (&tlow, &thigh, tlow, m[0]);
  /* Let w = 2^wordsize. We know (phigh * w + plow) + (thigh * w + tlow) 
     == 0 (mod w) so either plow == tlow == 0, or plow !=0 and tlow != 0. 
     In the former case we want phigh + thigh + 1, in the latter 
     phigh + thigh */
  /* Since a <= p-1 and b <= p-1, and p <= w-1, a*b <= w^2 - 4*w + 4, so
     adding 1 to phigh is safe */
#if 0
  /* Slower? */
  modul_add_ul_2ul (&plow, &phigh, tlow);
#else
  phigh += (plow != 0UL ? 1UL : 0UL);
#endif

  modul_add (r, &phigh, &thigh, m);

#if defined(MODTRACE)
  printf (" == %lu /* PARI */ \n", r[0]);
#endif
}
                         
/* FIXME: check for overflow if b > m */
MAYBE_UNUSED
static inline void
modul_mulredc_ul (residueul_t r, const residueul_t a, const unsigned long b,
                  const unsigned long invm, const modulusul_t m)
{
  unsigned long plow, phigh, tlow, thigh;
  ASSERT_EXPENSIVE (m[0] % 2 != 0);
#if defined(MODTRACE)
  printf ("(%lu * %lu / 2^%ld) %% %lu", 
          a[0], b, 8 * sizeof(unsigned long), m[0]);
#endif
  
  mul_ul_ul_2ul (&plow, &phigh, a[0], b);
  tlow = plow * invm;
  mul_ul_ul_2ul (&tlow, &thigh, tlow, m[0]);
  phigh += (plow != 0UL ? 1UL : 0UL);
  r[0] = (phigh >= m[0] - thigh) ? (phigh - (m[0] - thigh)) : (phigh + thigh);
  
#if defined(MODTRACE)
  printf (" == %lu /* PARI */ \n", r[0]);
#endif
}
                         

/* Computes (a * b + c)/ 2^wordsize % m. Requires that 
   a * b + c < 2^wordsize * m */

MAYBE_UNUSED
static inline void
modul_muladdredc_ul (residueul_t r, const residueul_t a, const unsigned long b,
                     const unsigned long c, const unsigned long invm, 
                     const modulusul_t m)
{
  unsigned long plow, phigh, tlow, thigh;
  ASSERT_EXPENSIVE (m[0] % 2 != 0);
#if defined(MODTRACE)
  printf ("(%lu * %lu / 2^%ld) %% %lu", 
          a[0], b, 8 * sizeof(unsigned long), m[0]);
#endif
  
  mul_ul_ul_2ul (&plow, &phigh, a[0], b);
  modul_add_ul_2ul (&plow, &phigh, c);
  tlow = plow * invm;
  mul_ul_ul_2ul (&tlow, &thigh, tlow, m[0]);
  phigh += (plow != 0UL ? 1UL : 0UL);
  r[0] = (phigh >= m[0] - thigh) ? (phigh - (m[0] - thigh)) : (phigh + thigh);
  
#if defined(MODTRACE)
  printf (" == %lu /* PARI */ \n", r[0]);
#endif
}
                         

MAYBE_UNUSED
static inline void
modul_div2 (residueul_t r, const residueul_t a, const modulusul_t m)
{
  if (a[0] % 2UL == 0UL)
    r[0] = a[0] / 2UL;
  else
    {
      ASSERT(m[0] % 2UL != 0UL);
      r[0] = a[0] / 2UL + m[0] / 2UL + 1UL;
    }
}

MAYBE_UNUSED
static inline int
modul_next(residueul_t r, modulus_t p)
{
    return (++r[0] == p[0]);
}

MAYBE_UNUSED
static inline int
modul_finished(residueul_t r, modulus_t p)
{
    return (r[0] == p[0]);
}

/* prototypes of non-inline functions */
void modul_div3 (residueul_t, const residueul_t, const modulusul_t);
void modul_gcd (residueul_t, const residueul_t, const modulusul_t);
unsigned long modul_invmodlong (modulusul_t);
void modul_powredc_ul (residueul_t, const residueul_t, const unsigned long, 
		       const unsigned long, const modulusul_t);
void modul_powredc_mp (residueul_t, const residueul_t, const unsigned long *,
                       const int, const unsigned long, const modulusul_t);
void modul_2powredc_mp (residueul_t, const residueul_t, const unsigned long *, 
                        const int, const unsigned long, const unsigned long, 
                        const modulusul_t);
void modul_Vredc_ul (residueul_t, const residueul_t, const unsigned long, 
                     const unsigned long, const modulusul_t);
void modul_Vredc_mp (residueul_t, const residueul_t, const unsigned long *,
                     const int, const unsigned long, const modulusul_t);
int modul_sprp (const residueul_t, const unsigned long, const modulusul_t);
int modul_inv (residueul_t, const residueul_t, const modulusul_t);
int modul_jacobi (const residueul_t, const modulusul_t);

#endif  /* MOD_UL_H__ */
