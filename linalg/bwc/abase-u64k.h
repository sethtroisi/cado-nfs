#ifndef ABASE_U64K_H_
#define ABASE_U64K_H_

/* This provides a version of abase_u64, with a constant repeat
 * value. It's similar to ABASE_U64 in just about every respect.
 */

#define ABASE_U64K_REPEAT_COUNT 2

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
typedef void * abase_u64k_obj_t;
typedef void * abase_u64k_obj_ptr;
typedef const void * abase_u64k_obj_srcptr;

#define abase_u64k_obj_init(x) (x=NULL)
#define abase_u64k_obj_clear(x)
#define abase_u64k_obj_init_set(y,x) (y=x)

typedef uint64_t abase_u64k_base_type;

#define abase_u64k_repeat(x)  ABASE_U64K_REPEAT_COUNT

// nbits takes into account the repeat count. So it's repeat(x) times the
// number of bits in abase_u64k_base_type.
#define abase_u64k_nbits(x)  (64*ABASE_U64K_REPEAT_COUNT)

#define abase_u64k_max_accumulate(x) UINT_MAX
#define abase_u64k_max_accumulate_wide(x) UINT_MAX

#include "pad.h"
#define P(X)    PAD(abase_u64k,X)
#define ABASE_F(t,n,a) static inline t P(n) a
#include "abase-binary-generic.h"
                
/* This one is in the C file. */
extern void
P(dotprod)(P(obj_srcptr) x,
        P(base_type) * w,
        const P(base_type) * u,
        const P(base_type) * v,
        unsigned int n);

#ifndef ABASE_DONTBIND_u64k
/* Bind our interface as the default one */
#define ABASE_BIND(X)   PAD(abase_u64k,X)
#include "abase-api.h"
#endif

#undef  P
#undef  ABASE_F

#endif	/* ABASE_U64K_H_ */
