#ifndef CADO_MACROS_H_
#define CADO_MACROS_H_

/**********************************************************************/
/* Common asserting/debugging defines */
/* See README.macro_usage */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#define ASSERT(x)	assert(x)

#define croak__(x,y) do {						\
        fprintf(stderr,"%s in %s at %s:%d -- %s\n",			\
                (x),__func__,__FILE__,__LINE__,(y));			\
    } while (0)

#define ASSERT_ALWAYS(x)						\
    do {								\
        if (!(x)) {							\
            croak__("code BUG() : condition " #x " failed",		\
                    "Abort");						\
            abort();							\
        }								\
    } while (0)
#define FATAL_ERROR_CHECK(cond, msg)					\
    do {								\
      if (UNLIKELY((cond))) {                                           \
            croak__("Fatal error", msg);				\
          abort();                                                      \
        }								\
    } while (0)

#define DIE_ERRNO_DIAG(tst, func, arg) do {				\
    if ((tst)) {					        	\
        fprintf(stderr, func "(%s): %s\n", arg, strerror(errno));       \
        exit(1);					        	\
    }							        	\
} while (0)


/*********************************************************************/
/* Helper macros */
/* See README.macro_usage */

#ifndef ABS
#define ABS(x) ((x) >= 0 ? (x) : -(x))
#endif

#ifndef MIN
#define MIN(l,o) ((l) < (o) ? (l) : (o))
#endif

#ifndef MAX
#define MAX(h,i) ((h) > (i) ? (h) : (i))
#endif

#if defined(__GNUC__)
#ifndef	MAYBE_UNUSED
#define MAYBE_UNUSED __attribute__ ((unused))
#endif
#ifndef EXPECT
#define EXPECT(x,val)	__builtin_expect(x,val)
#endif
#else
#ifndef	MAYBE_UNUSED
#define MAYBE_UNUSED
#endif
#ifndef EXPECT
#define	EXPECT(x,val)	(x)
#endif
#endif

#ifndef	LIKELY
#define LIKELY(x)	EXPECT(x,1)
#endif
#ifndef	UNLIKELY
#define UNLIKELY(x)	EXPECT(x,0)
#endif

#define LEXGE2(X,Y,A,B) (X>A || (X == A && Y >= B))
#define LEXGE3(X,Y,Z,A,B,C) (X>A || (X == A && LEXGE2(Y,Z,B,C)))
#define LEXLE2(X,Y,A,B) LEXGE2(A,B,X,Y)
#define LEXLE3(X,Y,Z,A,B,C) LEXGE3(A,B,C,X,Y,Z)

#define GNUC_VERSION(X,Y,Z)     \
    defined(__GNUC__) &&        \
(__GNUC__ == X && __GNUC_MINOR__ == Y && __GNUC_PATCHLEVEL__ == Z)
#define GNUC_VERSION_ATLEAST(X,Y,Z)     \
    defined(__GNUC__) &&        \
LEXGE3(__GNUC__,__GNUC_MINOR__,__GNUC_PATCHLEVEL__,X,Y,Z)
#define GNUC_VERSION_ATMOST(X,Y,Z)     \
    defined(__GNUC__) &&        \
LEXLE3(__GNUC__,__GNUC_MINOR__,__GNUC_PATCHLEVEL__,X,Y,Z)


#endif	/* CADO_MACROS_H_ */
