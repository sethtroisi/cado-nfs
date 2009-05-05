#ifndef CADO_UTILS_MISC_H_
#define CADO_UTILS_MISC_H_

#include <stddef.h>
#include "macros.h"

#ifdef __cplusplus
extern "C" {
#endif

extern char * cado_strndup(const char * a, size_t n);

extern void * malloc_check(const size_t x);

extern void * malloc_aligned(size_t size, size_t alignment);
extern void free_aligned(void * ptr, size_t size, size_t alignment);

extern void * malloc_pagealigned(size_t sz); 
extern void free_pagealigned(void * ptr, size_t sz); 

/* k must be a power of 2. Returns the smallest multiple of k greater
 * than or equal to n
 */
static inline unsigned long
next_multiple_of_powerof2(unsigned long n, unsigned long k)
{
    ASSERT((k & (k-1)) == 0);
    return ((n-1)|(k-1)) + 1;
}
/* those builtins seem to have appeared in 3.4 (april 2004) */
#if defined(__GNUC__) && (__GNUC__ >= 4 || __GNUC__ >= 3 && __GNUC_MINOR__ >= 4)
#define clzl(x)         __builtin_clzl(x)
#define ctzl(x)         __builtin_ctzl(x)
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
#endif

#ifdef __cplusplus
}
#endif

#endif	/* CADO_UTILS_MISC_H_ */
