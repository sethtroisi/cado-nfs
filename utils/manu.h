#ifndef MANU_H_
#define MANU_H_

#ifdef	__cplusplus
#include <cstdlib>
#include <cassert>
#else
#include <stdlib.h>
#include <assert.h>
#endif

#ifdef	__cplusplus
extern "C" {
#endif

/* This stuff requires assert.h and stdio.h to be included first 
 * (or cassert and cstdio) */

/* Note: If we use the following, the full iostream stuff must be
 * included, which is a pain. We prefer to get around with cstdio alone
 */
#if defined(FORCE_CROAK_ON_CERR) && defined(__cplusplus)
#define manu_croaks(x,y)					\
	std::cout << flush;					\
	std::cerr	<< (x) << " in " << __func__ 		\
		<< " at " << __FILE__ << ':' << __LINE__	\
		<< " -- " << (y) << endl;

#else
#define manu_croaks(x,y)					\
	fprintf(stderr,"%s in %s at %s:%d -- %s\n",		\
			(x),__func__,__FILE__,__LINE__,(y))
#endif

#ifndef BEGIN_BLOCK
#define BEGIN_BLOCK	do {
#define END_BLOCK	} while (0)
#endif

/**********************************************************************/
/* A couple of helper macros that I like to define */
#ifdef ASSERT
#undef ASSERT
#endif
#ifdef ASSERT_ALWAYS
#undef ASSERT_ALWAYS
#endif
#ifdef BUG
#undef BUG
#endif
#ifdef BUG_ON
#undef BUG_ON
#endif
#ifdef BUG_IF
#undef BUG_IF
#endif
#ifdef MISSING
#undef MISSING
#endif
#ifdef WARNING
#undef WARNING
#endif
#ifdef BASH_USER
#undef BASH_USER
#endif
#ifdef DIE
#undef DIE
#endif

#define ASSERT(x)	assert(x)
#define BUG()	BEGIN_BLOCK				\
		manu_croaks("code BUG()", "Abort");	\
		abort();				\
	END_BLOCK
#define BUG_ON_MSG(x,msg)	BEGIN_BLOCK		\
	if (x) {					\
		manu_croaks("code BUG() : " msg,	\
			"Abort");			\
		abort();				\
	}						\
	END_BLOCK
#define	BUG_ON(x)	BUG_ON_MSG(x,"condition " #x " reached")
#define	ASSERT_ALWAYS(x)	BEGIN_BLOCK		\
	if (!(x)) {					\
		manu_croaks("code BUG() : "		\
			"condition " #x " failed",	\
			"Abort");			\
		abort();				\
	}						\
	END_BLOCK
#define MISSING()					\
	BEGIN_BLOCK					\
		manu_croaks("missing code", "Abort");	\
		BUG();					\
	END_BLOCK
#define WARNING(x)	manu_croaks("Warning",(x))
#define BAD_CODE_WARNING				\
	BEGIN_BLOCK					\
		static int i;				\
		if (!i++) {				\
			WARNING("This is bad code");	\
		}					\
	END_BLOCK
#define	BASH_USER(x)					\
	BEGIN_BLOCK					\
		manu_croaks("impolite user", (x));	\
		BUG();					\
	END_BLOCK
#define	DIE(x)						\
	BEGIN_BLOCK					\
		manu_croaks("die", (x));		\
		BUG();					\
	END_BLOCK

/**********************************************************************/
#ifndef TRUE
#define TRUE (0==0)
#endif

#ifndef FALSE
#define FALSE (0==1)
#endif

#ifndef ABS
#define ABS(x) ((x) >= 0 ? (x) : -(x))
#endif
	
#ifndef MIN
#define MIN(l,o) ((l) < (o) ? (l) : (o))
#endif

#ifndef MAX
#define MAX(h,i) ((h) > (i) ? (h) : (i))
#endif

/**********************************************************************/
#if defined(__GNUC__)
	/* && !defined(__cplusplus) */
#define try_inline __inline__
#ifndef	UNUSED_VARIABLE
#define UNUSED_VARIABLE __attribute__ ((unused))
#endif
#ifndef EXPECT
#define EXPECT(x,val)	__builtin_expect(x,val)
#endif
#else
#define try_inline
#ifndef	UNUSED_VARIABLE
#define UNUSED_VARIABLE
#endif
#ifndef EXPECT
#define	EXPECT(x,val)	(x)
#endif
#endif

#ifndef	EXPECT_TRUE
#define EXPECT_TRUE(x)	EXPECT(x,1)
#endif
#ifndef	EXPECT_FALSE
#define EXPECT_FALSE(x)	EXPECT(x,0)
#endif

/**********************************************************************/
#ifndef DISABLE_ALLOCA
/* This turns out to work fairly often */
#include <stdlib.h>
#define FAST_ALLOC(x)   alloca(x)
#define FAST_FREE(x)    
#else
#define FAST_ALLOC(x)   my_malloc(x)
#define FAST_FREE(x)    free(x)
#endif

/**********************************************************************/
/* That's dirty, but I really don't want to link libm in */
#define iceildiv(x,y)	(((x)+(y)-1)/(y))

#ifdef  __GNUC__
#define clzl(x)         __builtin_clzl(x)
#define ctzl(x)         __builtin_ctzl(x)
#define parityl(x)      __builtin_parityl(x)
#else
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
static inline int ctzl(unsigned long x)
{
	return GMP_LIMB_BITS - clzl(x & - x);
}
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
#endif

#ifdef	__cplusplus
}
#endif

#endif	/* MANU_H_ */
