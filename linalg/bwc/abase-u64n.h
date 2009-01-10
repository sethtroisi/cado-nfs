#ifndef ABASE_U64N_H_
#define ABASE_U64N_H_

/* This provides a version of abase_u64, with a runtime-determined repeat
 * value. The runtime value is the only thing to be seen in the abobj
 * data. Note that this structure is bound to be pretty slow, as the
 * runtime-specified repeat count is checked in all tight loops !!
 *
 * For another version which does the same for a constant repeat count,
 * hence not suffering from this speed penalty, see abase-u64k.h
 */

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
 * Typically, we have plain old datatypes here, e.g. uint64_t's
 *
 * The rather contrived interface must be regarded while having in mind
 * how all this could carry over to data types stored within mpns, even
 * possibly when the type width is going to be determined at runtime.
 */
typedef unsigned int abase_u64n_obj_t[1];
typedef unsigned int * abase_u64n_obj_ptr;
typedef const unsigned int * abase_u64n_obj_srcptr;

#define abase_u64n_obj_init(x) (x[0] = 2)
#define abase_u64n_obj_clear(x)
#define abase_u64n_obj_init_set(y,x) (y[0]=x[0])

typedef uint64_t abase_u64n_base_type;

#define abase_u64n_repeat(x)  (x[0])

// nbits takes into account the repeat count. So it's repeat(x) times the
// number of bits in abase_u64n_base_type.
#define abase_u64n_nbits(x)  (64*x[0])

#define abase_u64n_max_accumulate(x) UINT_MAX
#define abase_u64n_max_accumulate_wide(x) UINT_MAX

static inline
abase_u64n_base_type * abase_u64n_init(abase_u64n_obj_srcptr x MAYBE_UNUSED,
        unsigned int n)
{
    return (abase_u64n_base_type *) malloc(n * abase_u64n_repeat(x) * sizeof(abase_u64n_base_type));
    // return (abase_u64n_base_type *) electric_alloc(n * abase_u64n_repeat(x) * sizeof(abase_u64n_base_type));
}

static inline
abase_u64n_base_type * abase_u64n_clear(abase_u64n_obj_srcptr x MAYBE_UNUSED,
        abase_u64n_base_type * p,
        unsigned int n MAYBE_UNUSED)
{
    free(p);
    // electric_free(p, n * abase_u64n_repeat(x) * sizeof(abase_u64n_base_type));
    return NULL;
}

static inline
void abase_u64n_zero(abase_u64n_obj_srcptr x MAYBE_UNUSED,
        abase_u64n_base_type * p,
        unsigned int n MAYBE_UNUSED)
{
    memset(p, 0, n * abase_u64n_repeat(x) * sizeof(abase_u64n_base_type));
}

static inline
void abase_u64n_random(abase_u64n_obj_srcptr x MAYBE_UNUSED,
        abase_u64n_base_type * p,
        unsigned int n MAYBE_UNUSED)
{
    myrand_area(p, n * abase_u64n_repeat(x) * sizeof(abase_u64n_base_type));
}

static inline
void abase_u64n_copy(abase_u64n_obj_srcptr x MAYBE_UNUSED,
        abase_u64n_base_type * q,
        const abase_u64n_base_type * p,
        unsigned int n MAYBE_UNUSED)
{
    memcpy(q, p, n * abase_u64n_repeat(x) * sizeof(abase_u64n_base_type));
}

/*
static inline
abase_u64n_base_type *
abase_u64n_next(abase_u64n_obj_srcptr * x MAYBE_UNUSED,
        abase_u64n_base_type * p,
        unsigned int k)
{
    return p + k * abase_u64n_repeat(x);
}

... and a const counterpart. In fact, I'd rather do the addition from the caller side, it's not much different.
*/

static inline
ptrdiff_t
abase_u64n_offset(abase_u64n_obj_srcptr x MAYBE_UNUSED, unsigned int k)
{
    return k * abase_u64n_repeat(x);
}

static inline
size_t
abase_u64n_bytes(abase_u64n_obj_srcptr x MAYBE_UNUSED, unsigned int k)
{
    return abase_u64n_offset(x, k) * sizeof(abase_u64n_base_type);
}

static inline
int abase_u64n_is_zero(abase_u64n_obj_srcptr x MAYBE_UNUSED,
        const abase_u64n_base_type * p,
        unsigned int n MAYBE_UNUSED)
{
    for(unsigned int i = 0 ; i < n * abase_u64n_repeat(x) ; i++) {
        if (p[i])
            return 0;
    }
    return 1;
}

static inline
void abase_u64n_add(abase_u64n_obj_srcptr x MAYBE_UNUSED,
        abase_u64n_base_type * dst, const abase_u64n_base_type * src)
{
    for(unsigned int i = 0 ; i < abase_u64n_repeat(x) ; i++) {
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
abase_u64n_base_type * abase_u64n_initf(abase_u64n_obj_srcptr x MAYBE_UNUSED)
{
    return (abase_u64n_base_type *) alloca(abase_u64n_repeat(x) * sizeof(abase_u64n_base_type));
}

static inline
abase_u64n_base_type * abase_u64n_clearf(abase_u64n_obj_srcptr x MAYBE_UNUSED,
        abase_u64n_base_type * p MAYBE_UNUSED)
{
    return NULL;
}

static inline
unsigned int
abase_u64n_read(abase_u64n_obj_srcptr x MAYBE_UNUSED,
        FILE * f, abase_u64n_base_type * p, unsigned int n)
{
    return fread(p, abase_u64n_repeat(x) * sizeof(abase_u64n_base_type), n, f);
}

static inline
unsigned int
abase_u64n_write(abase_u64n_obj_srcptr x MAYBE_UNUSED,
        FILE * f, abase_u64n_base_type * p, unsigned int n)
{
    return fread(p, abase_u64n_repeat(x) * sizeof(abase_u64n_base_type), n, f);
}


/* This is the only non trivial function from this file. Given the number
 * K of entries represented by the abase type (which is
 * abase_u64n_nbits(x), here n * 64), write
 * a K times K matrix as K elements in the area pointed to by w.  */
#include <emmintrin.h>
static inline void
abase_u64n_dotprod(abase_u64n_obj_srcptr x MAYBE_UNUSED,
        abase_u64n_base_type * w,
        const abase_u64n_base_type * u,
        const abase_u64n_base_type * v,
        unsigned int n)
{
    memset(w, 0, abase_u64n_nbits(x) * abase_u64n_repeat(x) * sizeof(abase_u64n_base_type));
    for(unsigned int i = 0 ; i < n ; i++) {
        __v2di * w0 = (__v2di*) w;
        for(unsigned int l = 0 ; l < abase_u64n_repeat(x) ; l++) {
            // TODO: It's possible to expand more, and use a __v2di
            // mb[4][2], or even [4]. This wouldn't change the code much
            // (see the u128 version), and is likely to speed things up a
            // wee bit maybe.
            __v2di mb[4] = {
                (__v2di) {0, 0},
                (__v2di) {*v, 0},
                (__v2di) {0, *v},
                (__v2di) {*v, *v},
            };
            v++;
            __v2di *sw = w0;
            const abase_u64n_base_type * ut = u;
            for(unsigned int k = 0 ; k < abase_u64n_repeat(x) ; k++) {
                uint64_t a = *ut++;
                for (unsigned int j = 0; j < 64; j += 2) {
                    *sw ^= mb[a & 3];
                    a >>= 2;
                    sw += abase_u64n_repeat(x);
                }
            }
            w0++;
        }
        u += abase_u64n_repeat(x);
    }
}


/* Define our convenience macros */
#define abobj_t   abase_u64n_obj_t
#define abobj_ptr   abase_u64n_obj_ptr
#define abobj_srcptr   abase_u64n_obj_srcptr
#define abobj_init(x)   abase_u64n_obj_init(x)
#define abobj_init_set(y,x)   abase_u64n_obj_init_set(y,x)
#define abobj_clear(x)  abase_u64n_obj_clear(x)

#define abt     abase_u64n_base_type
#define abnbits(x)      abase_u64n_nbits(x)
#define abrepeat(x)     abase_u64n_repeat(x)
#define abmax_accumulate(x)     abase_u64n_max_accumulate(x)
#define abmax_accumulate_wide(x)        abase_u64n_max_accumulate_wide(x)
#define abinit(x,n)     abase_u64n_init(x,n)
#define abinitf(x,n)    abase_u64n_initf(x,n)
#define abclear(x,p,n)  abase_u64n_clear(x,p,n)
#define abclearf(x,p,n) abase_u64n_clearf(x,p,n)
#define abzero(x,p,n)   abase_u64n_zero(x,p,n)
#define abis_zero(x,p,n)        abase_u64n_is_zero(x,p,n)
#define abrandom(x,p,n) abase_u64n_random(x,p,n)
#define aboffset(x,k)   abase_u64n_offset(x,k)
#define abbytes(x,k)   abase_u64n_bytes(x,k)
#define abcopy(x,q,p,n) abase_u64n_copy(x,q,p,n)
#define abadd(x,q,p) abase_u64n_add(x,q,p)
#define abread(x,f,p,n) abase_u64n_read(x,f,p,n)
#define abwrite(x,f,p,n) abase_u64n_write(x,f,p,n)
#define abdotprod(x,w,u,v,n) abase_u64n_dotprod(x,w,u,v,n)

#endif	/* ABASE_U64N_H_ */
