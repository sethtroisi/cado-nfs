/* Some functions for modular arithmetic with residues and modulus in
   unsigned long variables. Residues are stored in Montgomery form,
   rediction after multiplication is done with REDC. Due to inlining, 
   this file must be included in the caller's source code with #include */

/* Naming convention: all function start with modredcul, for 
   MODulus REDC Unsigned Long, followed by underscore, functionality of 
   function (add, mul, etc), and possibly underscore and specification of 
   what argument types the function takes (_ul, etc).
   There are typedef's that rename all functions to mod_* instead of 
   modredcul_*, which one day might become an automatic renaming scheme so 
   that different modulus sizes can be used simply by #including different 
   mod_*.h files, but without changing anything else in the source code. */

#ifndef MODREDC_UL_H__

#define MODREDC_UL_H__

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

/* Define the default names and types for arithmetic with these functions */
#define mod_init             modredcul_init
#define mod_init_noset0      modredcul_init_noset0
#define mod_clear            modredcul_clear
#define mod_set              modredcul_set
#define mod_set_ul           modredcul_set_ul
#define mod_set_ul_reduced   modredcul_set_ul_reduced
#define mod_swap             modredcul_swap
#define mod_initmod_ul       modredcul_initmod_ul
#define mod_getmod_ul        modredcul_getmod_ul
#define mod_clearmod         modredcul_clearmod
#define mod_get_ul           modredcul_get_ul
#define mod_equal            modredcul_equal
#define mod_is0              modredcul_is0
#define mod_add              modredcul_add
#define mod_add_ul           modredcul_add_ul
#define mod_sub              modredcul_sub
#define mod_sub_ul           modredcul_sub_ul
#define mod_neg              modredcul_neg
#define mod_mul              modredcul_mul
#define mod_addredc_ul       modredcul_addredc_ul
#define mod_mulredc_ul       modredcul_mulredc_ul
#define mod_muladdredc_ul    modredcul_muladdredc_ul
#define mod_div2             modredcul_div2
#define mod_div3             modredcul_div3
#define mod_pow_ul           modredcul_pow_ul
#define mod_pow_mp           modredcul_pow_mp
#define mod_2pow_mp          modredcul_2pow_mp
#define mod_V_ul             modredcul_V_ul
#define mod_V_mp             modredcul_V_mp
#define mod_sprp             modredcul_sprp
#define mod_gcd              modredcul_gcd
#define mod_inv              modredcul_inv
#define mod_jacobi           modredcul_jacobi
#define mod_set0             modredcul_set0
#define mod_set1             modredcul_set1
#define mod_next             modredcul_next
#define mod_finished         modredcul_finished

typedef unsigned long residueredcul_t[1];
typedef struct { 
  unsigned long m;
  unsigned long invm;
} __modulusredcul_t;
typedef __modulusredcul_t modulusredcul_t[1];
typedef residueredcul_t residue_t;
typedef modulusredcul_t modulus_t;


static inline void 
modredcul_tomontgomery (residueredcul_t, const residueredcul_t, 
                        const modulusredcul_t m);
static inline void
modredcul_frommontgomery (residueredcul_t, const residueredcul_t, 
                          const modulusredcul_t);


/* Initialises a residue_t type and sets it to zero */
MAYBE_UNUSED
static inline void
modredcul_init (residueredcul_t r, const modulusredcul_t m MAYBE_UNUSED)
{
  r[0] = 0UL;
  return;
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
  r[0] = 1UL;
  modredcul_tomontgomery (r, r, m);
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
                          

/* Compute 1/t (mod 2^wordsize) */
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

  ASSERT (r * n == 1UL);

  return r;
}


MAYBE_UNUSED
static inline void
modredcul_initmod_ul (modulusredcul_t m, const unsigned long s)
{
  m[0].m = s;
  m[0].invm = -modredcul_invmodul (s);
}


MAYBE_UNUSED
static inline unsigned long
modredcul_getmod_ul (const modulusredcul_t m)
{
  return m[0].m;
}


MAYBE_UNUSED
static inline void
modredcul_clearmod (modulusredcul_t m MAYBE_UNUSED)
{
  return;
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
  return (modredcul_get_ul (a, m) == 1UL);
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
modredcul_neg (residueredcul_t r, const residueredcul_t a, 
               const modulusredcul_t m)
{
  ASSERT_EXPENSIVE (a[0] < m[0].m);
  if (a[0] == 0UL)
    r[0] = a[0];
  else
    r[0] = m[0].m - a[0];
}


/* Add an unsigned long to two unsigned longs with carry propagation from 
   low word to high word. Any carry out from high word is lost. */

static inline void
add_ul_2ul (unsigned long *r1, unsigned long *r2, const unsigned long a)
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
add_2ul_2ul (unsigned long *r1, unsigned long *r2, 
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


/* Computes (a * 2^wordsize) % m */
MAYBE_UNUSED
static inline void
modredcul_tomontgomery (residueredcul_t r, const residueredcul_t a, 
                        const modulusredcul_t m)
{
  ASSERT_EXPENSIVE (a[0] < m[0].m);
  div_2ul_ul_ul_r (r, 0UL, a[0], m[0].m);
}


/* Computes (a / 2^wordsize) % m */
MAYBE_UNUSED
static inline void
modredcul_frommontgomery (residueredcul_t r, const residueredcul_t a, 
                          const modulusredcul_t m)
{
  unsigned long tlow, thigh;
  tlow = a[0] * m[0].invm;
  mul_ul_ul_2ul (&tlow, &thigh, tlow, m[0].m);
  r[0] = thigh + ((a[0] != 0UL) ? 1UL : 0UL);
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
  mul_ul_ul_2ul (&tlow, &thigh, tlow, m[0].m);
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
  add_ul_2ul (&slow, &shigh, a[0]);
  
  tlow = slow * m[0].invm;
  mul_ul_ul_2ul (&tlow, &thigh, tlow, m[0].m);
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
  add_ul_2ul (&slow, &shigh, a[0]);
  shigh += (slow != 0UL) ? 1UL : 0UL;
#endif

  tlow = slow * m[0].invm;
  mul_ul_ul_2ul (&tlow, r, tlow, m[0].m);
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
  
  mul_ul_ul_2ul (&plow, &phigh, a[0], b[0]);
  tlow = plow * m[0].invm;
  mul_ul_ul_2ul (&tlow, &thigh, tlow, m[0].m);
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
  phigh += (plow != 0UL) ? 1UL : 0UL;
#endif

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
