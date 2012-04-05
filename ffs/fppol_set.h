// Always include "fppol.h" first so that everything gets loaded in the correct
// order.
#include "fppol.h"

#ifndef __FPPOL_SET_H__
#define __FPPOL_SET_H__

#include "cppmeta.h"



/* Fixed-size polynomials.
 *****************************************************************************/

// Bitwise OR of all bit vectors: coefficient-wise non-zero test.
// Generic prototype:
//   uint<sz>_t fppol<sz>_fold_or(fppol<sz>_srcptr p);
#define __DECL_FPPOLxx_FOLD_OR(sz)                        \
  static inline                                           \
  uint##sz##_t fppol##sz##_fold_or(fppol##sz##_srcptr p);


// All declarations bundled up into a single macro.
#define __DECL_FPPOLxx_ARITH_ALL(sz) \
        __DECL_FPPOLxx_FOLD_OR  (sz) \

__DECL_FPPOLxx_ARITH_ALL(16)
__DECL_FPPOLxx_ARITH_ALL(32)
__DECL_FPPOLxx_ARITH_ALL(64)

#undef __DECL_FPPOLxx_FOLD_OR
#undef __DECL_FPPOLxx_DEG
#undef __DECL_FPPOLxx_ARITH_ALL


// Set to zero.
// Generic prototype:
//   void fppol<sz>_set_zero(fppol<sz>_ptr r);
#define __DEF_FPPOLxx_SET_ZERO(sz)                       \
  static inline                                          \
  void fppol##sz##_set_zero(fppol##sz##_ptr r)           \
  { for (unsigned k = 0; k < __FP_BITS; ++k) r[k] = 0; }


// Set to one.
// Generic prototype:
//   void fppol<sz>_set_one(fppol<sz>_ptr r);
#define __DEF_FPPOLxx_SET_ONE(sz)             \
  static inline                               \
  void fppol##sz##_set_one(fppol##sz##_ptr r) \
  { __FP_ONE(sz, r); }


// Set to t^i.
// Generic prototype:
//   void fppol<sz>_set_ti(fppol<sz>_ptr r, unsigned i);
#define __DEF_FPPOLxx_SET_TI(sz)                           \
  static inline                                            \
  void fppol##sz##_set_ti(fppol##sz##_ptr r, unsigned i)   \
  { __FP_ONE(sz, r);                                       \
    for (unsigned k = 0; k < __FP_BITS; ++k) r[k] <<= i; }


// Set to another polynomial.
// Generic prototype:
//   void fppol<sz>_set(fppol<sz>_ptr r, fppol<sz>_srcptr p);
#define __DEF_FPPOLxx_SET(sz)                                   \
  static inline                                                 \
  void fppol##sz##_set(fppol##sz##_ptr r, fppol##sz##_srcptr p) \
  { for (unsigned k = 0; k < __FP_BITS; ++k) r[k] = p[k]; }


// Set to an <sp>-term polynomial.
// Return 1 if successful.
// Generic prototype:
//   int fppol<sz>_set_<sp>(fppol<sz>_ptr r, fppol<sp>_srcptr p);
#define __DEF_FPPOLxx_SET_yy(sz, sp)                                \
  static inline                                                     \
  int fppol##sz##_set_##sp(fppol##sz##_ptr r, fppol##sp##_srcptr p) \
  { for (unsigned k = 0; k < __FP_BITS; ++k)                        \
      r[k] = (uint##sz##_t)p[k];                                    \
    return (sz >= sp) || !(fppol##sp##_fold_or(p) >> sz); }


// Set to a multiprecision polynomial.
// Return 1 if successful.
// Generic prototype:
//   int fppol<sz>_set_mp(fppol<sz>_ptr r, fppol_srcptr p);
#define __DEF_FPPOLxx_SET_MP(sz)                            \
  static inline                                             \
  int fppol##sz##_set_mp(fppol##sz##_ptr r, fppol_srcptr p) \
  { fppol##sz##_set_64(r, p->limbs[0]);                     \
    return p->deg < sz; }


// Swap two polynomials.
// Generic prototype:
//   void fppol<sz>_swap(fppol<sz>_ptr p, fppol<sz>_ptr q);
#define __DEF_FPPOLxx_SWAP(sz)                                \
  static inline                                               \
  void fppol##sz##_swap(fppol##sz##_ptr p, fppol##sz##_ptr q) \
  { fppol##sz##_t t;       fppol##sz##_set(t, p);             \
    fppol##sz##_set(p, q); fppol##sz##_set(q, t); }


// Get degree-i coefficient.
// Generic prototype:
//   void fppol<sz>_get_coeff(fp_ptr r, fppol<sz>_srcptr p, unsigned i);
#define __DEF_FPPOLxx_GET_COEFF(sz)                                      \
  static inline                                                          \
  void fppol##sz##_get_coeff(fp_ptr r, fppol##sz##_srcptr p, unsigned i) \
  { for (unsigned k = 0; k < __FP_BITS; ++k)                             \
      r[k] = (p[k] >> i) & 1; }


// Set degree-i coefficient.
// Generic prototype:
//   void fppol<sz>_set_coeff(fppol<sz>_ptr r, fp_srcptr x, unsigned i);
#define __DEF_FPPOLxx_SET_COEFF(sz)                                      \
  static inline                                                          \
  void fppol##sz##_set_coeff(fppol##sz##_ptr r, fp_srcptr x, unsigned i) \
  { uint##sz##_t m = ~(1ul<<i);                                          \
    for (unsigned k = 0; k < __FP_BITS; ++k)                             \
      r[k] = (r[k] & m) | ((uint##sz##_t)x[k] << i); }


// Set to next n-term (monic) polynomial, in lexicographical order.
// Return zero in case of overflow.
// Generic prototype:
//   int fppol<sz>_<monic>_set_next(fppol<sz>_ptr r, fppol<sz>_srcptr p,
//                                  unsigned n);
#define __DEF_FPPOLxx_SET_NEXT(sz, monic)                                 \
  static inline                                                           \
  int CAT(CAT(fppol##sz, SWITCH(monic, EMPTY, _monic)), _set_next)        \
    (fppol##sz##_ptr r, fppol##sz##_srcptr p, unsigned n)                 \
  { CAT(CAT(__FP, SWITCH(monic, EMPTY, _MONIC)), _SET_NEXT)(sz, r, p, n); \
    return fppol##sz##_fold_or(r) & (((uint##sz##_t)1 << n) - 1); }


// Conversion of an (n+m)-term polynomial, whose m most significant
// coefficients form a monic polynomial, to an unsigned int.
// /!\ Assume that deg(p) < n+m.
// Generic prototype:
//   unsigned fppol<sz>_get_ui(fppol<sz>_srcptr p, unsigned n, unsigned m);
#define __DEF_FPPOLxx_GET_UI(sz)                                             \
  static inline                                                              \
  unsigned fppol##sz##_get_ui(fppol##sz##_srcptr p, MAYBE_UNUSED unsigned n, \
                                                    MAYBE_UNUSED unsigned m) \
  { unsigned r; __FP_GET_UI(sz, r, p, n, m); return r; }


// Conversion of an (n+m)-term polynomial, whose m most significant
// coefficients form a monic polynomial, from an unsigned int.
// Return 1 if successful.
// /!\ Assume that deg(p) < n+m.
// Generic prototype:
//   int fppol<sz>_set_ui(fppol<sz>_ptr r, unsigned x, unsigned n, unsigned m);
#define __DEF_FPPOLxx_SET_UI(sz)                                           \
  static inline                                                            \
  int fppol##sz##_set_ui(fppol##sz##_ptr r, unsigned x,                    \
                         MAYBE_UNUSED unsigned n, MAYBE_UNUSED unsigned m) \
  { __FP_SET_UI(sz, r, x, n, m); return 1; }


// All definitions bundled up into a single macro.
#define __DEF_FPPOLxx_SET_ALL(sz)       \
        __DEF_FPPOLxx_SET_ZERO (sz)     \
        __DEF_FPPOLxx_SET_ONE  (sz)     \
        __DEF_FPPOLxx_SET_TI   (sz)     \
        __DEF_FPPOLxx_SET      (sz)     \
        __DEF_FPPOLxx_SET_yy   (sz, 16) \
        __DEF_FPPOLxx_SET_yy   (sz, 32) \
        __DEF_FPPOLxx_SET_yy   (sz, 64) \
        __DEF_FPPOLxx_SET_MP   (sz)     \
        __DEF_FPPOLxx_SWAP     (sz)     \
        __DEF_FPPOLxx_GET_COEFF(sz)     \
        __DEF_FPPOLxx_SET_COEFF(sz)     \
        __DEF_FPPOLxx_SET_NEXT (sz,  )  \
        __DEF_FPPOLxx_SET_NEXT (sz, 1)  \
        __DEF_FPPOLxx_GET_UI   (sz)     \
        __DEF_FPPOLxx_SET_UI   (sz)

__DEF_FPPOLxx_SET_ALL(16)
__DEF_FPPOLxx_SET_ALL(32)
__DEF_FPPOLxx_SET_ALL(64)

#undef __DEF_FPPOLxx_SET_ZERO
#undef __DEF_FPPOLxx_SET_ONE
#undef __DEF_FPPOLxx_SET_TI
#undef __DEF_FPPOLxx_SET
#undef __DEF_FPPOLxx_SET_yy
#undef __DEF_FPPOLxx_SET_yy
#undef __DEF_FPPOLxx_SET_yy
#undef __DEF_FPPOLxx_SET_MP
#undef __DEF_FPPOLxx_SWAP
#undef __DEF_FPPOLxx_GET_COEFF
#undef __DEF_FPPOLxx_SET_COEFF
#undef __DEF_FPPOLxx_SET_NEXT
#undef __DEF_FPPOLxx_GET_UI
#undef __DEF_FPPOLxx_SET_UI
#undef __DEF_FPPOLxx_SET_ALL



/* Multiprecision polynomials.
 *****************************************************************************/

// Set to zero.
void fppol_set_zero(fppol_ptr r);

// Set to one.
void fppol_set_one(fppol_ptr r);

// Set to t^i.
void fppol_set_ti(fppol_ptr r, unsigned i);

// Set to another polynomial.
void fppol_set(fppol_ptr r, fppol_srcptr p);

// Set to an <sz>-term polynomial.
// Generic prototype:
//   int fppol_set_<sz>(fppol_ptr r, fppol<sz>_srcptr p);
#define __DECL_FPPOL_SET_xx(sz)                           \
  void fppol_set_##sz(fppol_ptr r, fppol##sz##_srcptr p);

__DECL_FPPOL_SET_xx(16)
__DECL_FPPOL_SET_xx(32)
__DECL_FPPOL_SET_xx(64)

#undef __DECL_FPPOL_SET_xx

// Swap two polynomials.
// /!\ Shallow copy only.
static inline
void fppol_swap(fppol_ptr p, fppol_ptr q)
{ __fppol_struct t = *p; *p = *q; *q = t; }

// Get degree-i coefficient.
static inline
void fppol_get_coeff(fp_ptr r, fppol_srcptr p, unsigned i)
{ if ((int)i > p->deg) fp_set_zero(r);
  else                 fppol64_get_coeff(r, p->limbs[i>>6], i&0x3f); }

// Set degree-i coefficient.
void fppol_set_coeff(fppol_ptr r, fp_srcptr x, unsigned i);

#endif  /* __FPPOL_SET_H__ */
