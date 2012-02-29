#ifndef __TYPES_H__
#define __TYPES_H__

#include "fppol.h"
#include "macros.h"
#include "cppmeta.h"



/* Field-specific size of polynomials according to their role in the sieve.
 *
 * - sq:      the type that contains a special q (and the corresponding root).
 *
 * - ai:      the type that contains coordinates of the (skewed) reduced basis
 *            of the q-lattice.
 *
 * - fbprime: the type for the polynomials generators of the prime ideals in
 *            the factor base; the size of such an element will probably never
 *            be larger than 2^32.
 *
 * - ij:      the type for a polynomial i (or j).
 *
 * /!\ Sizes must be either 16, 32, or 64.
 *****************************************************************************/

// Over GF(2).
#if   defined(USE_F2)
# define      __sq_SIZE 64
# define      __ai_SIZE 32
# define __fbprime_SIZE 32
# define      __ij_SIZE 32

// Over GF(3).
#elif defined(USE_F3)
# define      __sq_SIZE 64
# define      __ai_SIZE 32
# define __fbprime_SIZE 32
# define      __ij_SIZE 32
#endif

// Multiprecision has no fixed size.
#define __fppol_SIZE



/* Type aliases for sieve-related types.
 *****************************************************************************/

#define __DECL_ALIAS_TYPE(type)                                    \
  typedef CAT(CAT(fppol, __##type##_SIZE), _t)      type##_t;      \
  typedef CAT(CAT(fppol, __##type##_SIZE), _ptr)    type##_ptr;    \
  typedef CAT(CAT(fppol, __##type##_SIZE), _srcptr) type##_srcptr;

__DECL_ALIAS_TYPE(sq)
__DECL_ALIAS_TYPE(ai)
__DECL_ALIAS_TYPE(fbprime)
__DECL_ALIAS_TYPE(ij)

#undef __DECL_ALIAS_TYPE



/* Function aliases for sieve-related types.
 *****************************************************************************/

// Helper constants for tests.
#define __void_NO_RETURN ,
#define __char_ptr_TYPE  , char *

// Helper macros for lists of arguments.
#define __ARGS_PROTO(type, n) type arg##n
#define __ARGS_CALL( type, n) arg##n


// Generic alias for fppol<sz>_<fun> functions.
#define __DECL_ALIAS_FUN(ret, type, fun, ...)           \
  static inline                                         \
  SWITCH(ret, TYPE, ret)                                \
  type##_##fun(FOR_ALL_ARGS(__ARGS_PROTO, __VA_ARGS__)) \
  { IF(ret, NO_RETURN, , return)                        \
    CAT(CAT(fppol, __##type##_SIZE), _##fun)            \
    (FOR_ALL_ARGS(__ARGS_CALL, __VA_ARGS__)); }


// Generic alias for fppol<sz>_<fun>_<sp> functions.
#define __DECL_ALIAS_FUN1(ret, type, fun, tp, ...)                  \
  static inline                                                     \
  SWITCH(ret, TYPE, ret)                                            \
  type##_##fun##_##tp(FOR_ALL_ARGS(__ARGS_PROTO, __VA_ARGS__))      \
  { IF(ret, NO_RETURN, , return)                                    \
    CAT(CAT(fppol, __##type##_SIZE), CAT(_##fun##_, __##tp##_SIZE)) \
    (FOR_ALL_ARGS(__ARGS_CALL, __VA_ARGS__)); }


// Generic alias for fppol<sz>_<fun>_<sp>x<sq> functions.
#define __DECL_ALIAS_FUN2(ret, type, fun, tp, tq, ...)                \
  static inline                                                       \
  SWITCH(ret, TYPE, ret)                                              \
  type##_##fun##_##tp##x##tq(FOR_ALL_ARGS(__ARGS_PROTO, __VA_ARGS__)) \
  { IF(ret, NO_RETURN, , return)                                      \
    CAT(CAT(fppol, __##type##_SIZE),                                  \
        CAT(CAT(_##fun##_, __##tp##_SIZE), CAT(x, __##tq##_SIZE)))    \
    (FOR_ALL_ARGS(__ARGS_CALL, __VA_ARGS__)); }


#define __DECL_ALIAS_MUL(type)                                          \
    __DECL_ALIAS_FUN (void, type, mul,    type##_ptr,                   \
                                          type##_srcptr, type##_srcptr) \
    __DECL_ALIAS_FUN (void, type, addmul, type##_ptr,                   \
                                          type##_srcptr, type##_srcptr) \
    __DECL_ALIAS_FUN (void, type, submul, type##_ptr,                   \
                                          type##_srcptr, type##_srcptr)

#define __DECL_ALIAS_MUL_yy(type, tp)                                     \
    __DECL_ALIAS_FUN1(void, type, mul,    tp, type##_ptr,                 \
                                              type##_srcptr, tp##_srcptr) \
    __DECL_ALIAS_FUN1(void, type, addmul, tp, type##_ptr,                 \
                                              type##_srcptr, tp##_srcptr) \
    __DECL_ALIAS_FUN1(void, type, submul, tp, type##_ptr,                 \
                                              type##_srcptr, tp##_srcptr)

#define __DECL_ALIAS_MUL_yyxzz(type, tp, tq)                                \
    __DECL_ALIAS_FUN2(void, type, mul,    tp, tq, type##_ptr,               \
                                                  tp##_srcptr, tq##_srcptr) \
    __DECL_ALIAS_FUN2(void, type, addmul, tp, tq, type##_ptr,               \
                                                  tp##_srcptr, tq##_srcptr) \
    __DECL_ALIAS_FUN2(void, type, submul, tp, tq, type##_ptr,               \
                                                  tp##_srcptr, tq##_srcptr)


// All function aliases bundled up into a single macro.
#define __DECL_ALIAS_ALL(type)                                                \
  __DECL_ALIAS_FUN (void, type, set_zero,       type##_ptr)                   \
  __DECL_ALIAS_FUN (void, type, set_one,        type##_ptr)                   \
  __DECL_ALIAS_FUN (void, type, set_ti,         type##_ptr,    unsigned)      \
  __DECL_ALIAS_FUN (void, type, set,            type##_ptr,    type##_srcptr) \
  __DECL_ALIAS_FUN1(int,  type, set,    sq,     type##_ptr,    sq_srcptr)     \
  __DECL_ALIAS_FUN1(int,  type, set,    ai,     type##_ptr,    ai_srcptr)     \
  __DECL_ALIAS_FUN1(int,  type, set,    fbprime,type##_ptr,    fbprime_srcptr)\
  __DECL_ALIAS_FUN1(int,  type, set,    ij,     type##_ptr,    ij_srcptr)     \
  __DECL_ALIAS_FUN (int,  type, set_mp,         type##_ptr,    fppol_srcptr)  \
  __DECL_ALIAS_FUN1(void, fppol,set,    type,   fppol_ptr,     type##_srcptr) \
  __DECL_ALIAS_FUN (void, type, get_coeff,      fp_ptr,                       \
                                                type##_srcptr, unsigned)      \
  __DECL_ALIAS_FUN (void, type, set_coeff,      type##_ptr,                   \
                                                fp_srcptr,     unsigned)      \
  __DECL_ALIAS_FUN (int,  type, is_zero,        type##_srcptr)                \
  __DECL_ALIAS_FUN (int,  type, eq,             type##_srcptr, type##_srcptr) \
  __DECL_ALIAS_FUN (int,  type, is_monic,       type##_srcptr)                \
  __DECL_ALIAS_FUN (int,  type, is_valid,       type##_srcptr)                \
  __DECL_ALIAS_FUN (void, type, opp,            type##_ptr,    type##_srcptr) \
  __DECL_ALIAS_FUN (void, type, add,            type##_ptr,                   \
                                                type##_srcptr, type##_srcptr) \
  __DECL_ALIAS_FUN (void, type, sub,            type##_ptr,                   \
                                                type##_srcptr, type##_srcptr) \
  __DECL_ALIAS_FUN (void, type, shl,            type##_ptr,                   \
                                                type##_srcptr, unsigned)      \
  __DECL_ALIAS_FUN (void, type, shr,            type##_ptr,                   \
                                                type##_srcptr, unsigned)      \
  __DECL_ALIAS_FUN (int,      type, deg,        type##_srcptr)                \
  __DECL_ALIAS_FUN (char_ptr, type, get_str,    char *,        type##_srcptr) \
  __DECL_ALIAS_FUN (int,      type, set_str,    type##_ptr,    const char *)  \
  __DECL_ALIAS_FUN (void,     type, out,        FILE *,        type##_srcptr) \
  __DECL_ALIAS_FUN (int,      type, inp,        type##_ptr,    FILE *)        \
  __DECL_ALIAS_MUL      (type)                                                \
  __DECL_ALIAS_MUL_yy   (fppol, type)                                         \
  __DECL_ALIAS_MUL_yyxzz(fppol, type, sq)                                     \
  __DECL_ALIAS_MUL_yyxzz(fppol, type, ai)                                     \
  __DECL_ALIAS_MUL_yyxzz(fppol, type, fbprime)                                \
  __DECL_ALIAS_MUL_yyxzz(fppol, type, ij)                                     \
  __DECL_ALIAS_FUN (int,  type, divrem,         type##_ptr,    type##_ptr,    \
                                                type##_srcptr, type##_srcptr) \
  __DECL_ALIAS_FUN (int,  type, div,            type##_ptr,                   \
                                                type##_srcptr, type##_srcptr) \
  __DECL_ALIAS_FUN (int,  type, rem,            type##_ptr,                   \
                                                type##_srcptr, type##_srcptr) \
  __DECL_ALIAS_FUN (void, type, shl1mod,        type##_ptr,                   \
                                                type##_srcptr, type##_srcptr) \
  __DECL_ALIAS_FUN (void, type, mulmod,         type##_ptr,    type##_srcptr, \
                                                type##_srcptr, type##_srcptr) \
  __DECL_ALIAS_FUN (int,  type, invmod,         type##_ptr,                   \
                                                type##_srcptr, type##_srcptr) \
  __DECL_ALIAS_FUN (void, type, gcd,            type##_ptr,                   \
                                                type##_srcptr, type##_srcptr)

__DECL_ALIAS_ALL(sq)
__DECL_ALIAS_ALL(ai)
__DECL_ALIAS_ALL(fbprime)
__DECL_ALIAS_ALL(ij)

#undef __void_NO_RETURN
#undef __char_ptr_TYPE
#undef __ARGS_PROTO
#undef __ARGS_CALL
#undef __DECL_ALIAS_FUN
#undef __DECL_ALIAS_FUN1
#undef __DECL_ALIAS_FUN2
#undef __DECL_ALIAS_MUL_yy
#undef __DECL_ALIAS_MUL_yyxzz
#undef __DECL_ALIAS_ALL




///////////////////////////////////////////////////////////////////
// Other types that are built from the others

// Information about the q-lattice.
// This includes the reduced basis.
typedef struct {
    sq_t q;
    sq_t rho;
    ai_t a0;
    ai_t a1;
    ai_t b0;
    ai_t b1;
    int side;
} qlat_struct_t;

typedef qlat_struct_t qlat_t[1];
typedef qlat_struct_t* qlat_ptr;
typedef const qlat_struct_t* qlat_srcptr;

// TODO: with alignements, we loose a lot of space
// See the types in cado-nfs and try to imitate ?
typedef struct {
    fbprime_t p;
    fbprime_t r;
    unsigned char degp;
} fbideal_t;

// a vector of (i,j) polynomials in the sieve.
typedef struct {
    ij_t i;
    ij_t j;
} ijvec_struct;

typedef ijvec_struct ijvec_t[1];

static inline void ijvec_add(ijvec_t W, ijvec_t V, ijvec_t U) {
    ij_add(W->i, V->i, U->i);
    ij_add(W->j, V->j, U->j);
}

// corresponding position in S. (might return 2 integers at some point)
typedef unsigned int ijpos_t;

// conversion: to be adapted for each case.
static inline ijpos_t ijvec2pos(ijvec_t V,
        MAYBE_UNUSED int I, MAYBE_UNUSED int J) {
    return (ijpos_t)(V->j[0])*(1U<<I) + (ijpos_t)V->i[0];
}
static inline void ijpos2vec(ijvec_t V, ijpos_t P, 
        MAYBE_UNUSED int I, MAYBE_UNUSED int J) {
    V->i[0] = (P) % (1U<<I);    // TODO: bit-fiddling, here
    V->j[0] = (P) / (1U<<I);
}



#endif   /* __TYPES_H__ */
