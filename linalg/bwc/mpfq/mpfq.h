#ifndef MPFQ_H_
#define MPFQ_H_

/* This header contains common declarations used by mpfq modules */

#include <gmp.h>

#ifdef __cplusplus
extern "C" {
#endif



/*** Constants for field_specify ***/

#define MPFQ_DONE 0             /* At the end of the variadic option functions */
#define MPFQ_PRIME_MPN 1        /* mp_limb_t *, size depending on implementation. Prefer MPFQ_PRIME_MPZ */
#define MPFQ_POLYNOMIAL 2       /* this expects an mpfq polynomial */
#define MPFQ_DEGREE 3           /* int */
#define MPFQ_IO_TYPE 4          /* for setopt */
#define MPFQ_GROUPSIZE 5        /* int (SIMD group size) */
#define MPFQ_PRIME_MPZ 6        /* mpz_t */

/***  Some useful macros ***/

#define MALLOC_FAILED()                                                 \
        do {                                                            \
                fprintf(stderr, "malloc failed in %s\n", __func__);     \
                abort();                                                \
        } while (0)

#define BUILD_BITMASK(x) ((x) == GMP_LIMB_BITS ? ((mp_limb_t) - 1) : (~ - ((mp_limb_t) 1 << (x))))

#define LEXGE2(X,Y,A,B) (X>A || (X == A && Y >= B))
#define LEXGE3(X,Y,Z,A,B,C) (X>A || (X == A && LEXGE2(Y,Z,B,C)))
#define LEXLE2(X,Y,A,B) LEXGE2(A,B,X,Y)
#define LEXLE3(X,Y,Z,A,B,C) LEXGE3(A,B,C,X,Y,Z)

#ifndef GNUC_VERSION
#define GNUC_VERSION(X,Y,Z)     \
    defined(__GNUC__) &&        \
(__GNUC__ == X && __GNUC_MINOR__ == Y && __GNUC_PATCHLEVEL__ == Z)
#endif

#ifndef GNUC_VERSION_ATLEAST
#define GNUC_VERSION_ATLEAST(X,Y,Z)     \
    defined(__GNUC__) &&        \
LEXGE3(__GNUC__,__GNUC_MINOR__,__GNUC_PATCHLEVEL__,X,Y,Z)
#endif

#ifndef GNUC_VERSION_ATMOST
#define GNUC_VERSION_ATMOST(X,Y,Z)     \
    defined(__GNUC__) &&        \
LEXLE3(__GNUC__,__GNUC_MINOR__,__GNUC_PATCHLEVEL__,X,Y,Z)
#endif


/* typedef unsigned long ulong; */

#ifndef MAX
#define MAX(h,i) ((h) > (i) ? (h) : (i))
#endif

#ifndef MIN
#define MIN(l,o) ((l) < (o) ? (l) : (o))
#endif

#ifndef __cplusplus
#ifndef max
#define max(a,b) ((a) > (b) ? (a) : (b))
#endif

#ifndef min
#define min(a,b) ((a) < (b) ? (a) : (b))
#endif
#endif	/* __cplusplus */


#ifndef ABS
#define ABS(x) ((x) >= 0 ? (x) : -(x))
#endif

#if GNUC_VERSION_ATLEAST(3,4,0)
/* according to
 * http://gcc.gnu.org/onlinedocs/gcc-3.1.1/gcc/Variable-Attributes.html#Variable%20Attributes
 * the 'unused' attribute already existed in 3.1.1 ; however the rules
 * for its usage remained quirky until 3.4.0, so we prefer to stick to
 * the more modern way of using the unused attribute, and recommend
 * setting the -Wno-unused flag for pre-3.4 versions of gcc
 */
#ifndef	MAYBE_UNUSED
#define	MAYBE_UNUSED	__attribute__((unused))
#endif
#endif

#if GNUC_VERSION_ATLEAST(3,1,0) /* apparently */
#ifndef MPFQ_EXPECT
#define MPFQ_EXPECT(x,val)   __builtin_expect(x,val)
#endif
#endif

#if GNUC_VERSION_ATLEAST(3,4,0)
#define clzl(x)         __builtin_clzl(x)
#define HAVE_clzl
#endif

#if GNUC_VERSION_ATLEAST(3,4,0)
#define ctzl(x)         __builtin_ctzl(x)
#define HAVE_ctzl
#endif

#if GNUC_VERSION_ATLEAST(3,4,0)
#define parityl(x)      __builtin_parityl(x)
#define HAVE_parityl
#endif

#ifndef	MAYBE_UNUSED
#define	MAYBE_UNUSED	/**/
#endif

#ifndef MPFQ_EXPECT
#define MPFQ_EXPECT(x,val)   (x)
#endif

#ifndef	MPFQ_UNLIKELY
#define	MPFQ_UNLIKELY(x)	MPFQ_EXPECT(x, 0)
#endif
#ifndef	MPFQ_LIKELY
#define	MPFQ_LIKELY(x)	MPFQ_EXPECT(x, 1)
#endif


#ifndef HAVE_clzl
/* provide slow fallbacks */
static inline int clzl(unsigned long x)
{
        static const int t[4] = { 2, 1, 0, 0 };
        int a = 0;
        int res;
#if (GMP_LIMB_BITS == 64)
        if (x >> 32) { a += 32; x >>= 32; }
#endif  
        if (x >> 16) { a += 16; x >>= 16; }
        if (x >>  8) { a +=  8; x >>=  8; }
        if (x >>  4) { a +=  4; x >>=  4; }
        if (x >>  2) { a +=  2; x >>=  2; }
        res = GMP_LIMB_BITS - 2 - a + t[x];
        return res;
}
#define HAVE_clzl
#define HAVE_clzl_fallback
#endif

#ifndef HAVE_ctzl
static inline int ctzl(unsigned long x)
{
	return GMP_LIMB_BITS - clzl(x & - x);
}
#define HAVE_ctzl
#define HAVE_ctzl_fallback
#endif

#ifndef HAVE_parityl
static inline int parityl(unsigned long x)
{
	static const int t[4] = { 0, 1, 1, 0, };
#if (GMP_LIMB_BITS == 64)
	x ^= (x >> 32);
#endif
	x ^= (x >> 16);
	x ^= (x >>  8);
	x ^= (x >>  4);
	x ^= (x >>  2);
	return t[x & 3UL];
}
#define HAVE_parityl
#define HAVE_parityl_fallback
#endif

static inline int clzlx(unsigned long * x, int n)
{
	int r = 0;
	for( ; n > 0 && MPFQ_UNLIKELY(!x[n-1]) ; --n) r+=GMP_LIMB_BITS;
	if (n == 0) return r;
	r += clzl(x[n-1]);
	return r;
}

static inline int ctzlx(unsigned long * x, int n)
{
	int r = 0;
	for( ; n > 0 && MPFQ_UNLIKELY(!*x) ; --n,++x) r+=GMP_LIMB_BITS;
	if (n == 0) return r;
	r += ctzl(*x);
	return r;
}

#ifdef __cplusplus
}
#endif

#endif	/* MPFQ_H_ */
