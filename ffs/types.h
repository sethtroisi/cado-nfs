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



/* Type aliases for sieve-related types.
 *****************************************************************************/

#define __ALIAS_TYPE(type)                                         \
  typedef CAT(CAT(fppol, __##type##_SIZE), _t)      type##_t;      \
  typedef CAT(CAT(fppol, __##type##_SIZE), _ptr)    type##_ptr;    \
  typedef CAT(CAT(fppol, __##type##_SIZE), _srcptr) type##_srcptr;

__ALIAS_TYPE(sq)
__ALIAS_TYPE(ai)
__ALIAS_TYPE(fbprime)
__ALIAS_TYPE(ij)

#undef __ALIAS_TYPE



/* Function aliases for sieve-related types.
 *****************************************************************************/

// Size of know polynomial types.
#define __fppol16_SIZE 16
#define __fppol32_SIZE 32
#define __fppol64_SIZE 64
#define   __fppol_SIZE

// Get type and size for function arguments.
#undef  IS_NUM
#undef  IS_MP
#define  __8_IS_NUM ,
#define __16_IS_NUM ,
#define __32_IS_NUM ,
#define __64_IS_NUM ,
#define __mp_IS_MP  ,
#define __ARG_TYPE(t) IF(t, IS_MP, fppol, IF(t, IS_NUM, fppol##t, t))
#define __ARG_SIZE(t) IF(t, IS_MP, mp,    IF(t, IS_NUM, t, __##t##_SIZE))

// No return statement when return type is void.
#undef  NO_RETURN
#define __void_NO_RETURN ,

// Make type name out of some common abbreviations.
#undef  TYPE
#define    __FILE_ptr_TYPE(...)   , FILE *
#define    __char_ptr_TYPE(...)   , char *
#define __char_srcptr_TYPE(...)   , const char *
#define    ___ptr_TYPE(t, ...)    ,  t##_ptr
#define ___srcptr_TYPE(t, ...)    ,  t##_srcptr
#define     ___tp_TYPE(t, tp, tq) , tp##_srcptr
#define     ___tq_TYPE(t, tp, tq) , tq##_srcptr
#define __GET_TYPE(type, ...) SWITCH(type, TYPE(__VA_ARGS__), type)

// Helper macros for lists of arguments.
#define __ARGS_PROTO(ctx, type, n) __GET_TYPE(type, EXPAND(ctx)) arg##n
#define __ARGS_CALL( ctx, type, n) arg##n


// Generic alias for fppol<sz>_<fun> functions.
#define __ALIAS_FUN(ret, type, fun, ...)                       \
  static inline                                                \
  __GET_TYPE(ret, type, )                                      \
  type##_##fun(FOR_ALL(__ARGS_PROTO, (type, , ), __VA_ARGS__)) \
  { IF(ret, NO_RETURN, , return)                               \
    CAT(CAT(fppol, __##type##_SIZE), _##fun)                   \
      (FOR_ALL(__ARGS_CALL, , __VA_ARGS__)); }

// Generic alias for fppol<sz>_<fun>_<sp> functions.
#define __ALIAS_FUN1(ret, type, fun, tp, ...)                        \
  static inline                                                      \
  __GET_TYPE(ret, type, __ARG_TYPE(tp), )                            \
  type##_##fun##_##tp                                                \
    (FOR_ALL(__ARGS_PROTO, (type, __ARG_TYPE(tp), ), __VA_ARGS__))   \
  { IF(ret, NO_RETURN, , return)                                     \
    CAT(CAT(fppol, __##type##_SIZE), CAT(_##fun##_, __ARG_SIZE(tp))) \
      (FOR_ALL(__ARGS_CALL, , __VA_ARGS__)); }

// Generic alias for fppol<sz>_<fun>_<sp>x<sq> functions.
#define __ALIAS_FUN2(ret, type, fun, tp, tq, ...)                    \
  static inline                                                      \
  __GET_TYPE(ret, type, __ARG_TYPE(tp), __ARG_TYPE(tq))              \
  type##_##fun##_##tp##x##tq                                         \
    (FOR_ALL(__ARGS_PROTO, (type, __ARG_TYPE(tp), __ARG_TYPE(tq)),   \
             __VA_ARGS__))                                           \
  { IF(ret, NO_RETURN, , return)                                     \
    CAT(CAT(fppol, __##type##_SIZE),                                 \
        CAT(CAT(_##fun##_, __ARG_SIZE(tp)), CAT(x, __ARG_SIZE(tq)))) \
      (FOR_ALL(__ARGS_CALL, , __VA_ARGS__)); }

// Aliases for multiplication.
#define __ALIAS_MUL(type)                                    \
    __ALIAS_FUN (void, type, mul,    _ptr, _srcptr, _srcptr) \
    __ALIAS_FUN (void, type, addmul, _ptr, _srcptr, _srcptr) \
    __ALIAS_FUN (void, type, submul, _ptr, _srcptr, _srcptr)

#define __ALIAS_MUL1(type, tp)                               \
    __ALIAS_FUN1(void, type, mul,    tp, _ptr, _srcptr, _tp) \
    __ALIAS_FUN1(void, type, addmul, tp, _ptr, _srcptr, _tp) \
    __ALIAS_FUN1(void, type, submul, tp, _ptr, _srcptr, _tp)

#define __ALIAS_MUL2(type, tp, tq)                           \
    __ALIAS_FUN2(void, type, mul,    tp, tq, _ptr, _tp, _tq) \
    __ALIAS_FUN2(void, type, addmul, tp, tq, _ptr, _tp, _tq) \
    __ALIAS_FUN2(void, type, submul, tp, tq, _ptr, _tp, _tq)


// All function aliases bundled up into a single macro.
#define __ALIAS_FUN_ALL(type)                                                 \
  __ALIAS_FUN (void,     type, set_zero,       _ptr)                          \
  __ALIAS_FUN (void,     type, set_one,        _ptr)                          \
  __ALIAS_FUN (void,     type, set_ti,         _ptr, unsigned)                \
  __ALIAS_FUN (void,     type, set,            _ptr, _srcptr)                 \
  __ALIAS_FUN1(int,      type, set, sq,        _ptr, _tp)                     \
  __ALIAS_FUN1(int,      type, set, ai,        _ptr, _tp)                     \
  __ALIAS_FUN1(int,      type, set, fbprime,   _ptr, _tp)                     \
  __ALIAS_FUN1(int,      type, set, ij,        _ptr, _tp)                     \
  __ALIAS_FUN1(int,      type, set,  8,        _ptr, _tp)                     \
  __ALIAS_FUN1(int,      type, set, 16,        _ptr, _tp)                     \
  __ALIAS_FUN1(int,      type, set, 32,        _ptr, _tp)                     \
  __ALIAS_FUN1(int,      type, set, 64,        _ptr, _tp)                     \
  __ALIAS_FUN1(int,      type, set, mp,        _ptr, _tp)                     \
  __ALIAS_FUN1(int,   fppol16, set, type,      _ptr, _tp)                     \
  __ALIAS_FUN1(int,   fppol32, set, type,      _ptr, _tp)                     \
  __ALIAS_FUN1(int,   fppol64, set, type,      _ptr, _tp)                     \
  __ALIAS_FUN1(void,    fppol, set, type,      _ptr, _tp)                     \
  __ALIAS_FUN (void,     type, swap,           _ptr, _ptr)                    \
  __ALIAS_FUN (void,     type, get_coeff,      fp_ptr, _srcptr, unsigned)     \
  __ALIAS_FUN (void,     type, set_coeff,      _ptr, fp_srcptr, unsigned)     \
  __ALIAS_FUN (int,      type, set_next,       _ptr, _srcptr, unsigned)       \
  __ALIAS_FUN (int,      type, monic_set_next, _ptr, _srcptr, unsigned)       \
  __ALIAS_FUN (uint64_t, type, get_ui,         _srcptr,  unsigned, unsigned)  \
  __ALIAS_FUN (int,      type, set_ui,   _ptr, uint64_t, unsigned, unsigned)  \
  __ALIAS_FUN (int,      type, deg,            _srcptr)                       \
  __ALIAS_FUN (int,      type, is_zero,        _srcptr)                       \
  __ALIAS_FUN (int,      type, in_fp,          _srcptr)                       \
  __ALIAS_FUN (int,      type, eq,             _srcptr, _srcptr)              \
  __ALIAS_FUN (int,      type, cmp,            _srcptr, _srcptr)              \
  __ALIAS_FUN (int,      type, is_monic,       _srcptr)                       \
  __ALIAS_FUN (int,      type, is_valid,       _srcptr)                       \
  __ALIAS_FUN (void,     type, opp,            _ptr, _srcptr)                 \
  __ALIAS_FUN (void,     type, add,            _ptr, _srcptr, _srcptr)        \
  __ALIAS_FUN (void,     type, sub,            _ptr, _srcptr, _srcptr)        \
  __ALIAS_FUN (void,     type, smul,           _ptr, _srcptr, fp_srcptr)      \
  __ALIAS_FUN (void,     type, sdiv,           _ptr, _srcptr, fp_srcptr)      \
  __ALIAS_FUN (void,     type, mul_ti,         _ptr, _srcptr, unsigned)       \
  __ALIAS_FUN (void,     type, div_ti,         _ptr, _srcptr, unsigned)       \
  __ALIAS_FUN (void,     type, mod_ti,         _ptr, _srcptr, unsigned)       \
  __ALIAS_FUN (void,     type, add_disjoint,   _ptr, _srcptr, _srcptr)        \
  __ALIAS_FUN (void,     type, diff,           _ptr, _srcptr, _srcptr)        \
  __ALIAS_FUN (int,      type, divrem,   _ptr, _ptr, _srcptr, _srcptr)        \
  __ALIAS_FUN (int,      type, div,            _ptr, _srcptr, _srcptr)        \
  __ALIAS_FUN (int,      type, rem,            _ptr, _srcptr, _srcptr)        \
  __ALIAS_FUN (void,     type, multmod,        _ptr, _srcptr, _srcptr)        \
  __ALIAS_FUN (void,     type, mulmod,         _ptr, _srcptr, _srcptr,_srcptr)\
  __ALIAS_FUN (int,      type, invmod,         _ptr, _srcptr, _srcptr)        \
  __ALIAS_FUN (void,     type, gcd,            _ptr, _srcptr, _srcptr)        \
  __ALIAS_FUN (char_ptr, type, get_str,        char_ptr, _srcptr)             \
  __ALIAS_FUN (int,      type, set_str,        _ptr,     char_srcptr)         \
  __ALIAS_FUN (void,     type, out,            FILE_ptr, _srcptr)             \
  __ALIAS_FUN (int,      type, inp,            _ptr,     FILE_ptr)            \
  __ALIAS_MUL (          type)                                                \
  __ALIAS_MUL1(         fppol, type)                                          \
  __ALIAS_MUL2(         fppol, type, sq)                                      \
  __ALIAS_MUL2(         fppol, type, ai)                                      \
  __ALIAS_MUL2(         fppol, type, fbprime)                                 \
  __ALIAS_MUL2(         fppol, type, ij)

__ALIAS_FUN_ALL(sq)
__ALIAS_FUN_ALL(ai)
__ALIAS_FUN_ALL(fbprime)
__ALIAS_FUN_ALL(ij)

#undef  __8_IS_NUM
#undef __16_IS_NUM
#undef __32_IS_NUM
#undef __64_IS_NUM
#undef __mp_IS_MP
#undef __void_NO_RETURN
#undef    __FILE_ptr_TYPE
#undef    __char_ptr_TYPE
#undef __char_srcptr_TYPE
#undef    ___ptr_TYPE
#undef ___srcptr_TYPE
#undef     ___tp_TYPE
#undef     ___tq_TYPE
#undef __ARG_TYPE
#undef __ARG_SIZE
#undef __GET_TYPE
#undef __ARGS_PROTO
#undef __ARGS_CALL
#undef __ALIAS_FUN
#undef __ALIAS_FUN1
#undef __ALIAS_FUN2
#undef __ALIAS_MUL
#undef __ALIAS_MUL1
#undef __ALIAS_MUL2
#undef __ALIAS_FUN_ALL



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


/* Vector in the reduced q-lattice as a pair of polynomials (i,j).
 * Note that, by convention, j should be kept monic.
 *****************************************************************************/

// Structure definition.
typedef struct {
  ij_t i;
  ij_t j;
} __ijvec_struct;

// Type and pointer shorthands.
typedef       __ijvec_struct  ijvec_t[1];
typedef       __ijvec_struct *ijvec_ptr;
typedef const __ijvec_struct *ijvec_srcptr;


/* Basis of the p-lattice seen as a GF(p)-vector space of (i,j)-vectors.
 *****************************************************************************/

// Structure definition.
typedef struct {
  unsigned I, J;
  unsigned dim;
  ijvec_t *v;
} __ijbasis_struct;

// Type and pointer shorthands.
typedef       __ijbasis_struct  ijbasis_t[1];
typedef       __ijbasis_struct *ijbasis_ptr;
typedef const __ijbasis_struct *ijbasis_srcptr;


// Element of the factor base as an ideal gothp = (p, r).
// Given the two vectors (a0, b0) and (a1, b1) of the reduced q-lattice,
// lambda is precomputed as lambda = - (a1 - r*b1) / (a0 - r*b0) mod p.
// In the case of a projective root (i.e., a0 - r*b0 = 0 mod p),
// proj is set to 1.
// TODO: with alignements, we loose a lot of space
// See the types in cado-nfs and try to imitate ?
typedef struct {
  fbprime_t p;
  fbprime_t r;
  fbprime_t lambda;
  ij_t      tildep;    // 1/p mod sublat_info.modulus
  ij_t      i0;
  //ij_t      j0, j;   // FIXME: projective-root sieving.
  uint8_t   degp;
  _Bool     proj;
  _Bool     power;
} __fbideal_struct;

typedef       __fbideal_struct  fbideal_t[1];
typedef       __fbideal_struct *fbideal_ptr;
typedef const __fbideal_struct *fbideal_srcptr;

typedef struct {
  unsigned   n;  // nb of entries in the factor base
  fbideal_t *elts;
} __factor_base_struct;

typedef       __factor_base_struct  factor_base_t[1];
typedef       __factor_base_struct *factor_base_ptr;
typedef const __factor_base_struct *factor_base_srcptr;

/* ffspol is a type to store the defining polynomial of the function
 * field.
 ****************************************************************************/

// Polynomials over GF(p)[x,t], represented as polynomials in x whose
// coefficients lie in GF(p)[t] (using the multiprecision type fppol_t).
// By convention, deg(0) = -1.
typedef struct {
  int       deg;
  unsigned  alloc;
  fppol_t  *coeffs;
} __ffspol_struct;

// Type and pointer shorthands.
// - ffspol_t:      type of polynomials.
// - ffspol_ptr:    read/write pointer (internal type)
// - ffspol_srcptr: read-only  pointer (internal type)
typedef       __ffspol_struct  ffspol_t[1];
typedef       __ffspol_struct *ffspol_ptr;
typedef const __ffspol_struct *ffspol_srcptr;



#endif   /* __TYPES_H__ */
