#ifndef MANU_H_
#define MANU_H_

#include <limits.h>
#include "cado.h"

#ifdef	__cplusplus
extern "C" {
#endif

#ifndef BEGIN_BLOCK
#define BEGIN_BLOCK	do {
#define END_BLOCK	} while (0)
#endif

/**********************************************************************/
/* A couple of helper macros that I like to define */
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

#define BUG()	BEGIN_BLOCK				\
		croak__("code BUG()", "Abort");	\
		abort();				\
	END_BLOCK
#define BUG_ON_MSG(x,msg)	BEGIN_BLOCK		\
	if (x) {					\
		croak__("code BUG() : " msg,	\
			"Abort");			\
		abort();				\
	}						\
	END_BLOCK
#define	BUG_ON(x)	BUG_ON_MSG(x,"condition " #x " reached")
#define MISSING()					\
	BEGIN_BLOCK					\
		croak__("missing code", "Abort");	\
		BUG();					\
	END_BLOCK
#define WARNING(x)	croak__("Warning",(x))
#define BAD_CODE_WARNING				\
	BEGIN_BLOCK					\
		static int i;				\
		if (!i++) {				\
			WARNING("This is bad code");	\
		}					\
	END_BLOCK
#define	BASH_USER(x)					\
	BEGIN_BLOCK					\
		croak__("impolite user", (x));	\
		BUG();					\
	END_BLOCK
#define	DIE(x)						\
	BEGIN_BLOCK					\
		croak__("die", (x));		\
		BUG();					\
	END_BLOCK

/**********************************************************************/
#ifndef TRUE
#define TRUE (0==0)
#endif

#ifndef FALSE
#define FALSE (0==1)
#endif

/**********************************************************************/
#ifndef DISABLE_ALLOCA
/* This turns out to work fairly often for defining alloca */
#include <stdlib.h>
#define FAST_ALLOC(x)   alloca(x)
#define FAST_FREE(x)    
#else
#define FAST_ALLOC(x)   my_malloc(x)
#define FAST_FREE(x)    free(x)
#endif

/**********************************************************************/

/* Handy, and does not require libm */
#ifndef iceildiv
#define iceildiv(x,y)	(((x)+(y)-1)/(y))
#endif

/* Number of words holding B bits ; better naming sought. */
#ifndef  BITS_TO_WORDS
#define	BITS_TO_WORDS(B,W)	iceildiv((B),(W))
#endif


/* those builtins seem to have appeared in 3.4 (april 2004) */
#if defined(__GNUC__) && (__GNUC__ >= 4 || __GNUC__ >= 3 && __GNUC_MINOR__ >= 4)
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

#if defined(__cplusplus) && GNUC_VERSION_ATLEAST(4,3,0)
/* Starting with gcc 4.3, -Wempty-body moans for loops like
 * for(;(x=x->next)!=NULL;y++);
 * It must shut up.
 */
#pragma GCC diagnostic ignored "-Wempty-body"
#endif

#ifdef	__cplusplus
}
#endif

#endif	/* MANU_H_ */
