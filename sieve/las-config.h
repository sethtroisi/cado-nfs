#ifndef LAS_CONFIG_H_
#define LAS_CONFIG_H_

#include "cado_config.h"

#ifdef HAVE_SSE2
#define SSE_NORM_INIT
#endif

#ifdef SSE_NORM_INIT
#include <emmintrin.h>
#endif

/* As its name says, this is a ugly hack that initializes all lognorms to the
   maximal value (255) on the rational side. But it seems to work well, and to
   miss only about 7% to 8% relations wrt a more accurate estimation. */
//#define UGLY_HACK

/* number of bits used to estimate the norms */
#define NORM_BITS 10

/* define PROFILE to keep certain function from being inlined, in order to
   make them show up on profiler output */
//#define PROFILE

/* Trick to discard lognorms that will probably lead to non L-smooth
   cofactors. Disabled for now since it requires accurate sieving
   (in particular by prime powers and bad primes). */
// #define COFACTOR_TRICK

/* This triggers code which fills bucket in several passes, one for each
 * congruence class mod 2 (three such, the trivial one leading to
 * spurious reports). It's currently only part of the story, and at the
 * moment it is almost neutral in terms of efficiency (small slowdown
 * because of access pattern).  But it's the way to go if one wants to
 * support I=16. There are many other places where changes must be made.
 * This particular flag affects only the treatment of the ``bucket
 * sieved'' primes, not the pattern-sieved, or small-sieved.
 */
#define MOD2_CLASSES_BS 0       /* define to 0 or 1 */

/* un-sieving of locations where gcd(i,j)>1 instead of testing gcd for
 * each survivor. Appears slower than default. This code has always been
 * #ifdef'd out, but maybe can be improved enough to make it worthwhile
 */
#define xxxUNSIEVE_NOT_COPRIME  /* see las-unsieve.c */


/* default sieve region side is 2^DEFAULT_I */
#define DEFAULT_I 12

/* default bucket region: 2^16 = 64K == close to L1 size, but this is the
   (current) largest possible value, otherwise bucket.h must be changed,
   since it stores positions on 16 bits */
#ifndef LOG_BUCKET_REGION
#define LOG_BUCKET_REGION 16
#endif

#if LOG_BUCKET_REGION > 16
#error "Too large LOG_BUCKET_REGION, please adapt bucket.h first"
#endif

/* This flag is necessary to support I=16. Otherwise it's a useless
 * burden
 */
#define SUPPORT_I16

/* Define SKIP_GCD3 to skip updates where 3 divides gcd(i,j) in the
   bucket sieving phase. Slightly slower than not skipping them
   in single-thread mode, but might be useful for multi-threading,
   or when memory is tight */
// #define SKIP_GCD3

/* These parameters control the size of the buckets. 
 * The number of updates that a bucket can accumulate is estimated as
 *   (loglog(factor base bound) - loglog(bucket sieving threshold)) 
 *     * BUCKET_LIMIT_FACTOR * I * J + BUCKET_LIMIT_ADD 
 * We don't store updates where 2 divides gcd(i,j) which reduces the number 
 * of updates to about (1-1/4)=3/4, so 0.8 should be safe.
 * If we don't store updates where 3 divides gcd(i,j) either, their number
 * is reduced to about (1-1/4)*(1-1/9)=2/3, then 0.7 should be ok
 *
 * A significant part of the inaccuracy in predicting the bucket sizes
 * stems from the use of the Mertens estimate in lieu of proper sums. Now
 * that we _do_ compute this sum, we gain some precision.
 */

#ifndef BUCKET_LIMIT_FACTOR
#ifdef SKIP_GCD3
#define BUCKET_LIMIT_FACTOR (2.0/3.0 * 1.05)
#else
#define BUCKET_LIMIT_FACTOR (3.0/4.0 * 1.05)
#endif
#endif

#ifndef BUCKET_LIMIT_ADD
#define BUCKET_LIMIT_ADD 0
#endif

/* Guard for the logarithms of norms, so that the value does not wrap around
   zero due to roundoff errors. */
#define GUARD 4

/* GUARD+LOG_MAX should be as near as possible from 256, to get more accuracy
   in the norm computations, but not too much, otherwise a norm might be
   rounded to zero. */
#define LOG_MAX (255.9 - (double) GUARD)

/* See PROFILE flag above */
/* Some functions should not be inlined when we profile or it's hard or
   impossible to tell them apart from the rest in the profiler output */
#ifdef PROFILE
#define NOPROFILE_INLINE
#define NOPROFILE_STATIC
#else
#define NOPROFILE_INLINE static inline
#define NOPROFILE_STATIC static
#endif

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

void las_display_config_flags(FILE * stream);

#ifdef __cplusplus
}
#endif

#endif	/* LAS_CONFIG_H_ */
