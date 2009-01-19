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

#include "pad.h"
#define P(X)    PAD(abase_u64n,X)
#define ABASE_F(t,n,a) static inline t P(n) a

#define AVOID_GENERIC_obj_set_nbys
/* This generic version assumes compile-time-constant stride */
ABASE_F(void,obj_set_nbys,(P(obj_ptr) x MAYBE_UNUSED, unsigned int nbys))
{
    unsigned int q = nbys / 64;
    unsigned int r = nbys % 64;
    ASSERT_ALWAYS(r == 0);
    x[0] = q;
}

#include "abase-binary-generic.h"

/* This one is in the C file. */
extern void
P(dotprod)(P(obj_srcptr) x,
        P(base_type) * w,
        const P(base_type) * u,
        const P(base_type) * v,
        unsigned int n);

#ifndef ABASE_DONTBIND_u64n
/* Bind our interface as the default one */
#define ABASE_BIND(X)    PAD(abase_u64n,X)
#include "abase-api.h"
#endif

#undef  P
#undef  ABASE_F

#endif	/* ABASE_U64N_H_ */
