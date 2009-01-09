#ifndef ABASE_U128_H_
#define ABASE_U128_H_

/* This is used in order to enable GF(2)-only code */
#define ABASE_BINARY

#include <stddef.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

/* I know it's ugly */
#include <alloca.h>

#include "macros.h"
#include "random_generation.h"
// #include "electric_alloc.h"


/* provides an interface for the arithmetic base of computations.
 * Typically, we have plain old datatypes here, e.g. uint128_t's
 *
 * The rather contrived interface must be regarded while having in mind
 * how all this could carry over to data types stored within mpns, even
 * possibly when the type width is going to be determined at runtime.
 */
typedef void * abase_u128_obj_t;
typedef void * abase_u128_obj_ptr;
typedef const void * abase_u128_obj_srcptr;

#define abase_u128_obj_init(x) (x=NULL)
#define abase_u128_obj_clear(x)

#include <emmintrin.h>

typedef __v2di abase_u128_base_type;
#define abase_u128_nbits(x)  128

// repeating is just the fact of putting several base type elements next
// to each other.
#define abase_u128_repeat(x)  1
#define abase_u128_needs_scrap(x)    false

#define abase_u128_max_accumulate(x) UINT_MAX
#define abase_u128_max_accumulate_wide(x) UINT_MAX

static inline
abase_u128_base_type * abase_u128_init(abase_u128_obj_srcptr x MAYBE_UNUSED,
        unsigned int n)
{
    return (abase_u128_base_type *) malloc(n * abase_u128_repeat(x) * sizeof(abase_u128_base_type));
    // return (abase_u128_base_type *) electric_alloc(n * abase_u128_repeat(x) * sizeof(abase_u128_base_type));
}

static inline
abase_u128_base_type * abase_u128_clear(abase_u128_obj_srcptr x MAYBE_UNUSED,
        abase_u128_base_type * p,
        unsigned int n MAYBE_UNUSED)
{
    free(p);
    // electric_free(p, n * abase_u128_repeat(x) * sizeof(abase_u128_base_type));
    return NULL;
}

static inline
void abase_u128_zero(abase_u128_obj_srcptr x MAYBE_UNUSED,
        abase_u128_base_type * p,
        unsigned int n MAYBE_UNUSED)
{
    memset(p, 0, n * abase_u128_repeat(x) * sizeof(abase_u128_base_type));
}

static inline
void abase_u128_random(abase_u128_obj_srcptr x MAYBE_UNUSED,
        abase_u128_base_type * p,
        unsigned int n MAYBE_UNUSED)
{
    myrand_area(p, n * abase_u128_repeat(x) * sizeof(abase_u128_base_type));
}

static inline
void abase_u128_copy(abase_u128_obj_srcptr x MAYBE_UNUSED,
        abase_u128_base_type * q,
        const abase_u128_base_type * p,
        unsigned int n MAYBE_UNUSED)
{
    memcpy(q, p, n * abase_u128_repeat(x) * sizeof(abase_u128_base_type));
}

/*
static inline
abase_u128_base_type *
abase_u128_next(abase_u128_obj_srcptr * x MAYBE_UNUSED,
        abase_u128_base_type * p,
        unsigned int k)
{
    return p + k * abase_u128_repeat(x);
}

... and a const counterpart. In fact, I'd rather do the addition from the caller side, it's not much different.
*/

static inline
ptrdiff_t
abase_u128_offset(abase_u128_obj_srcptr x MAYBE_UNUSED, unsigned int k)
{
    return k * abase_u128_repeat(x);
}

static inline
size_t
abase_u128_bytes(abase_u128_obj_srcptr x MAYBE_UNUSED, unsigned int k)
{
    return abase_u128_offset(x, k) * sizeof(abase_u128_base_type);
}

static inline
int abase_u128_is_zero(abase_u128_obj_srcptr x MAYBE_UNUSED,
        const abase_u128_base_type * p,
        unsigned int n MAYBE_UNUSED)
{
    for(unsigned int i = 0 ; i < n * abase_u128_repeat(x) * sizeof(abase_u128_base_type) ; i++) {
        if (((unsigned char *)p)[i])
            return 0;
    }
    return 1;
}

static inline
void abase_u128_add(abase_u128_obj_srcptr x MAYBE_UNUSED,
        abase_u128_base_type * dst, const abase_u128_base_type * src)
{
    for(unsigned int i = 0 ; i < abase_u128_repeat(x) ; i++) {
        dst[i] ^= src[i];
    }
}

#if 0
    // we need an out-of-place reduction routine. The binary one is
    // stupid of course. But note that the output traits are in
    // accordance with the ``wide version'' of the input traits.
    inline void reduce(t * dst, const t * src, unsigned int n = 1) const {
        memcpy(dst, src, n * repeat() * sizeof(t));
    }
#endif

/* The initf versions do alloca -- it's meant for temporaries
 */
static inline
abase_u128_base_type * abase_u128_initf(abase_u128_obj_srcptr x MAYBE_UNUSED)
{
    return (abase_u128_base_type *) alloca(abase_u128_repeat(x) * sizeof(abase_u128_base_type));
}

static inline
abase_u128_base_type * abase_u128_clearf(abase_u128_obj_srcptr x MAYBE_UNUSED,
        abase_u128_base_type * p MAYBE_UNUSED)
{
    return NULL;
}

static inline
unsigned int
abase_u128_read(abase_u128_obj_srcptr x MAYBE_UNUSED,
        FILE * f, abase_u128_base_type * p, unsigned int n)
{
    return fread(p, abase_u128_repeat(x) * sizeof(abase_u128_base_type), n, f);
}

static inline
unsigned int
abase_u128_write(abase_u128_obj_srcptr x MAYBE_UNUSED,
        FILE * f, abase_u128_base_type * p, unsigned int n)
{
    return fread(p, abase_u128_repeat(x) * sizeof(abase_u128_base_type), n, f);
}

/* Define our convenience macros */
#define abobj_t   abase_u128_obj_t
#define abobj_ptr   abase_u128_obj_ptr
#define abobj_srcptr   abase_u128_obj_srcptr
#define abobj_init(x)   abase_u128_obj_init(x)
#define abobj_clear(x)  abase_u128_obj_clear(x)

#define abt     abase_u128_base_type
#define abnbits(x)      abase_u128_nbits(x)
#define abrepeat(x)     abase_u128_repeat(x)
#define abneeds_scrap(x)        abase_u128_needs_scrap(x)
#define abmax_accumulate(x)     abase_u128_max_accumulate(x)
#define abmax_accumulate_wide(x)        abase_u128_max_accumulate_wide(x)
#define abinit(x,n)     abase_u128_init(x,n)
#define abinitf(x,n)    abase_u128_initf(x,n)
#define abclear(x,p,n)  abase_u128_clear(x,p,n)
#define abclearf(x,p,n) abase_u128_clearf(x,p,n)
#define abzero(x,p,n)   abase_u128_zero(x,p,n)
#define abis_zero(x,p,n)        abase_u128_is_zero(x,p,n)
#define abrandom(x,p,n) abase_u128_random(x,p,n)
#define aboffset(x,k)   abase_u128_offset(x,k)
#define abbytes(x,k)   abase_u128_bytes(x,k)
#define abcopy(x,q,p,n) abase_u128_copy(x,q,p,n)
#define abadd(x,q,p) abase_u128_add(x,q,p)
#define abread(x,f,p,n) abase_u128_read(x,f,p,n)
#define abwrite(x,f,p,n) abase_u128_write(x,f,p,n)

#endif	/* ABASE_U128_H_ */
