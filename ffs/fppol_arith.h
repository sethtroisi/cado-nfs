// Always include "fppol.h" first so that everything gets loaded in the correct
// order.
#include "fppol.h"

#ifndef __FPPOL_ARITH_H__
#define __FPPOL_ARITH_H__

#include "macros.h"
#include "cppmeta.h"



/* Base field elements.
 *****************************************************************************/

// Set to zero.
static inline
void fp_set_zero(fp_ptr r)
{ for (unsigned k = 0; k < __FP_BITS; ++k) r[k] = 0; }


// Set to one.
static inline
void fp_set_one(fp_ptr r)
{ __FP_ONE(8, r); }


// Set to largest element.
static inline
void fp_set_max(fp_ptr r)
{ __FP_MAX(8, r); }


// Set to another element.
static inline
void fp_set(fp_ptr r, fp_srcptr p)
{ for (unsigned k = 0; k < __FP_BITS; ++k) r[k] = p[k]; }


// Test if zero.
static inline
int fp_is_zero(fp_srcptr p)
{ uint8_t t = p[0];
  for (unsigned k = 1; k < __FP_BITS; ++k) t |= p[k];
  return !t; }


// Test if one.
static inline
int fp_is_one(fp_srcptr p)
{ uint8_t t = 0;
  for (unsigned k = 1; k < __FP_BITS; ++k) t |= p[k];
  return p[0] == 1 && !t; }


// Inverse.
static inline
void fp_inv(fp_ptr r, fp_srcptr p)
{ __FP_SINV(8, r, p); }


// Multiplication.
static inline
void fp_mul(fp_ptr r, fp_srcptr p, fp_srcptr q)
{ __FP_SMUL(8, r, p, q); }


// Division.
static inline
void fp_div(fp_ptr r, fp_srcptr p, fp_srcptr q)
{ __FP_SDIV(8, r, p, q); }



/* Fixed-size polynomials.
 *****************************************************************************/

// Bitwise OR of all bit vectors: coefficient-wise non-zero test.
// Generic prototype:
//   uint<sz>_t fppol<sz>_fold_or(fppol<sz>_srcptr p);
#define __DECL_FPPOLxx_FOLD_OR(sz)                        \
  static inline                                           \
  uint##sz##_t fppol##sz##_fold_or(fppol##sz##_srcptr p);


// Degree.
// By convention, deg(0) = -1.
// Generic prototype:
//   int fppol<sz>_deg(fppol<sz>_srcptr p);
#define __DECL_FPPOLxx_DEG(sz)               \
  int fppol##sz##_deg(fppol##sz##_srcptr p);


// All declarations bundled up into a single macro.
#define __DECL_FPPOLxx_ARITH_ALL(sz) \
        __DECL_FPPOLxx_FOLD_OR  (sz) \
        __DECL_FPPOLxx_DEG      (sz)

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


// Set to largest n-term polynomial.
// Generic prototype:
//   void fppol<sz>_set_max(fppol<sz>_ptr r, unsigned n);
#define __DEF_FPPOLxx_SET_MAX(sz)                                  \
  static inline                                                    \
  void fppol##sz##_set_max(fppol##sz##_ptr r, unsigned n)          \
  { uint##sz##_t m = ((uint##sz##_t)1 << n) - 1;                   \
    __FP_MAX(sz, r);                                               \
    for (unsigned k = 0; k < __FP_BITS; ++k) r[k] = (-r[k]) & m; }


// Set to largest n-term monic polynomial.
// Generic prototype:
//   void fppol<sz>_monic_set_max(fppol<sz>_ptr r, unsigned n);
#define __DEF_FPPOLxx_MONIC_SET_MAX(sz)                         \
  static inline                                                 \
  void fppol##sz##_monic_set_max(fppol##sz##_ptr r, unsigned n) \
  { fppol##sz##_set_max(r, n-1);                                \
    fp_t t; fp_set_one(t);                                      \
    for (unsigned k = 0; k < __FP_BITS; ++k)                    \
      r[k] |= (uint##sz##_t)t[k] << (n-1); }


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


// Bitwise OR of all bit vectors: coefficient-wise non-zero test.
// Generic prototype:
//   uint<sz>_t fppol<sz>_fold_or(fppol<sz>_srcptr p);
#define __DEF_FPPOLxx_FOLD_OR(sz)                        \
  static inline                                          \
  uint##sz##_t fppol##sz##_fold_or(fppol##sz##_srcptr p) \
  { uint##sz##_t t = p[0];                               \
    for (unsigned k = 1; k < __FP_BITS; ++k) t |= p[k];  \
    return t; }


// Test if zero.
// Generic prototype:
//   int fppol<sz>_is_zero(fppol<sz>_srcptr p);
#define __DEF_FPPOLxx_IS_ZERO(sz)               \
  static inline                                 \
  int fppol##sz##_is_zero(fppol##sz##_srcptr p) \
  { return !fppol##sz##_fold_or(p); }


// Test if in GF(p) (i.e., deg <= 0).
// Generic prototype:
//   int fppol<sz>_in_fp(fppol<sz>_srcptr p);
#define __DEF_FPPOLxx_IN_FP(sz)               \
  static inline                               \
  int fppol##sz##_in_fp(fppol##sz##_srcptr p) \
  { return !(fppol##sz##_fold_or(p) >> 1); }


// Test if equal.
// Generic prototype:
//   int fppol<sz>_eq(fppol<sz>_srcptr p, fppol<sz>_srcptr q);
#define __DEF_FPPOLxx_EQ(sz)                                     \
  static inline                                                  \
  int fppol##sz##_eq(fppol##sz##_srcptr p, fppol##sz##_srcptr q) \
  { int rc = p[0] == q[0];                                       \
    for (unsigned k = 1; k < __FP_BITS; ++k)                     \
      rc = rc && p[k] == q[k];                                   \
    return rc; }


// Test if monic.
// Generic prototype:
//   int fppol<sz>_is_monic(fppol<sz>_srcptr p);
#define __DEF_FPPOLxx_IS_MONIC(sz)               \
  static inline                                  \
  int fppol##sz##_is_monic(fppol##sz##_srcptr p) \
  { int d = fppol##sz##_deg(p);                  \
    if (UNLIKELY(d == -1)) return 0;             \
    fp_t lc;                                     \
    fppol##sz##_get_coeff(lc, p, d);             \
    return fp_is_one(lc); }


// Test if valid representation.
// Generic prototype:
//   int fppol<sz>_is_valid(fppol<sz>_srcptr p);
#define __DEF_FPPOLxx_IS_VALID(sz)                            \
  static inline                                               \
  int fppol##sz##_is_valid(MAYBE_UNUSED fppol##sz##_srcptr p) \
  { return __FP_IS_VALID(sz, p); }


// Opposite.
// Generic prototype:
//   void fppol<sz>_opp(fppol<sz>_ptr r, fppol<sz>_srcptr p);
#define __DEF_FPPOLxx_OPP(sz)                                   \
  static inline                                                 \
  void fppol##sz##_opp(fppol##sz##_ptr r, fppol##sz##_srcptr p) \
  { __FP_OPP(sz, r, p); }


// Addition.
// Generic prototype:
//   void fppol<sz>_add(fppol<sz>_ptr    r,
//                      fppol<sz>_srcptr p, fppol<sz>_srcptr q);
#define __DEF_FPPOLxx_ADD(sz)                                      \
  static inline                                                    \
  void fppol##sz##_add(fppol##sz##_ptr    r,                       \
                       fppol##sz##_srcptr p, fppol##sz##_srcptr q) \
  { __FP_ADD(sz, r, p, q); }


// Subtraction.
// Generic prototype:
//   void fppol<sz>_sub(fppol<sz>_ptr    r,
//                      fppol<sz>_srcptr p, fppol<sz>_srcptr q);
#define __DEF_FPPOLxx_SUB(sz)                                      \
  static inline                                                    \
  void fppol##sz##_sub(fppol##sz##_ptr    r,                       \
                       fppol##sz##_srcptr p, fppol##sz##_srcptr q) \
  { __FP_SUB(sz, r, p, q); }


// Scalar multiplication.
// Generic prototype:
//   void fppol<sz>_smul(fppol<sz>_ptr    r,
//                       fppol<sz>_srcptr p, fp_srcptr x);
#define __DEF_FPPOLxx_SMUL(sz)                             \
  static inline                                            \
  void fppol##sz##_smul(fppol##sz##_ptr r,                 \
                        fppol##sz##_srcptr p, fp_srcptr x) \
  { fppol##sz##_t t;                                       \
    for (unsigned k = 0; k < __FP_BITS; ++k)               \
      t[k] = -(uint##sz##_t)x[k];                          \
    __FP_SMUL(sz, r, p, t); }


// Scalar division.
// Generic prototype:
//   void fppol<sz>_sdiv(fppol<sz>_ptr r,
//                       fppol<sz>_srcptr p, fp_srcptr x);
#define __DEF_FPPOLxx_SDIV(sz)                             \
  static inline                                            \
  void fppol##sz##_sdiv(fppol##sz##_ptr    r,              \
                        fppol##sz##_srcptr p, fp_srcptr x) \
  { fppol##sz##_t t;                                       \
    for (unsigned k = 0; k < __FP_BITS; ++k)               \
      t[k] = -(uint##sz##_t)x[k];                          \
    __FP_SDIV(sz, r, p, t); }


// Coefficient-wise inversion.
// Generic prototype:
//   void fppol<sz>_sinv(fppol<sz>_ptr r, fppol<sz>_srcptr p);
#define __DEF_FPPOLxx_SINV(sz)                                   \
  static inline                                                  \
  void fppol##sz##_sinv(fppol##sz##_ptr r, fppol##sz##_srcptr p) \
  { __FP_SINV(sz, r, p); }


// Addition without overlapping coefficients (i.e., coefficient-wise OR).
// Generic prototype:
//   void fppol<sz>_add_disjoint(fppol<sz>_ptr    r,
//                               fppol<sz>_srcptr p, fppol<sz>_srcptr q);
#define __DEF_FPPOLxx_ADD_DISJOINT(sz)                                      \
  static inline                                                             \
  void fppol##sz##_add_disjoint(fppol##sz##_ptr r,                          \
                                fppol##sz##_srcptr p, fppol##sz##_srcptr q) \
  { for (unsigned k = 0; k < __FP_BITS; ++k) r[k] = p[k] | q[k]; }


// Multiplication by t^i (i.e., left shift).
// Generic prototype:
//   void fppol<sz>_mul_ti(fppol<sz>_ptr r, fppol<sz>_srcptr p, unsigned i);
#define __DEF_FPPOLxx_MUL_TI(sz)                                              \
  static inline                                                               \
  void fppol##sz##_mul_ti(fppol##sz##_ptr r, fppol##sz##_srcptr p, unsigned i)\
  { for (unsigned k = 0; k < __FP_BITS; ++k) r[k] = i < sz ? p[k] << i : 0; }


// Quotient of division by t^i (i.e., right shift).
// Generic prototype:
//   void fppol<sz>_shr(fppol<sz>_ptr r, fppol<sz>_srcptr p, unsigned i);
#define __DEF_FPPOLxx_DIV_TI(sz)                                              \
  static inline                                                               \
  void fppol##sz##_div_ti(fppol##sz##_ptr r, fppol##sz##_srcptr p, unsigned i)\
  { for (unsigned k = 0; k < __FP_BITS; ++k) r[k] = i < sz ? p[k] >> i : 0; }


// Remainder modulo t^i (i.e., i first coefficients).
// Generic prototype:
//   void fppol<sz>_mod_ti(fppol<sz>_ptr r, fppol<sz>_srcptr p, unsigned i);
#define __DEF_FPPOLxx_MOD_TI(sz)                                              \
  static inline                                                               \
  void fppol##sz##_mod_ti(fppol##sz##_ptr r, fppol##sz##_srcptr p, unsigned i)\
  { uint##sz##_t m = ((uint##sz##_t)1 << i) - 1;                              \
    for (unsigned k = 0; k < __FP_BITS; ++k) r[k] = p[k] & m; }


// Conversion of a (monic) polynomial to an unsigned int, after a preliminary
// multiplication by t^i.
// Generic prototype:
//   unsigned fppol<sz>_<monic>_get_ui_mul_ti(fppol<sz>_srcptr p, unsigned i);
#define __DEF_FPPOLxx_GET_UI_MUL_TI(sz, monic)                               \
  static inline                                                              \
  unsigned CAT(CAT(fppol##sz, SWITCH(monic, EMPTY, _monic)), _get_ui_mul_ti) \
    (fppol##sz##_srcptr p, unsigned i)                                       \
  { unsigned r;                                                              \
    CAT(CAT(__FP, SWITCH(monic, EMPTY, _MONIC)), _GET_UI)(sz, r, p, i);      \
    return r; }


// Conversion of a (monic) polynomial to an unsigned int.
// Generic prototype:
//   unsigned fppol<sz>_<monic>_get_ui(fppol<sz>_srcptr p);
#define __DEF_FPPOLxx_GET_UI(sz, monic)                                       \
  static inline                                                               \
  unsigned CAT(CAT(fppol##sz, SWITCH(monic, EMPTY, _monic)), _get_ui)         \
    (fppol##sz##_srcptr p)                                                    \
  { return                                                                    \
    CAT(CAT(fppol##sz, SWITCH(monic, EMPTY, _monic)), _get_ui_mul_ti)(p, 0); }


// Return the largest unsigned int corresponding to an n-term (monic)
// polynomial.
// Generic prototype:
//   unsigned fppol<sz>_<monic>_get_ui_max(unsigned n);
#define __DEF_FPPOLxx_GET_UI_MAX(sz, monic)                                   \
  static inline                                                               \
  unsigned CAT(CAT(fppol##sz, SWITCH(monic, EMPTY, _monic)), _get_ui_max)     \
    (unsigned n)                                                              \
  { fppol##sz##_t t;                                                          \
    CAT(CAT(       fppol##sz, SWITCH(monic, EMPTY, _monic)), _set_max)(t, n); \
    return CAT(CAT(fppol##sz, SWITCH(monic, EMPTY, _monic)), _get_ui) (t); }


// Conversion of a (monic) polynomial from an unsigned int.
// Return 1 if successful.
// Generic prototype:
//   int fppol<sz>_<monic>_set_ui(fppol<sz>_ptr r, unsigned x);
#define __DEF_FPPOLxx_SET_UI(sz, monic)                                \
  static inline                                                        \
  int CAT(CAT(fppol##sz, SWITCH(monic, EMPTY, _monic)), _set_ui)       \
    (fppol##sz##_ptr r, unsigned x)                                    \
  { CAT(CAT(__FP, SWITCH(monic, EMPTY, _MONIC)), _SET_UI)(sz, , r, x); \
    return 1; }


// Conversion of a (monic) polynomial from the next valid unsigned int
// representation directly following x.
// Return the unsigned int representation used for the conversion.
// The behavior is undefined if the next valid representation is out of bounds.
// Generic prototype:
//   unsigned fppol<sz>_<monic>_set_next_ui(fppol<sz>_ptr r, unsigned x);
#define __DEF_FPPOLxx_SET_NEXT_UI(sz, monic)                               \
  static inline                                                            \
  unsigned CAT(CAT(fppol##sz, SWITCH(monic, EMPTY, _monic)), _set_next_ui) \
    (fppol##sz##_ptr r, unsigned x)                                        \
  { CAT(CAT(__FP, SWITCH(monic, EMPTY, _MONIC)), _SET_UI)(sz, 1, r, x);    \
    return x; }


// All definitions bundled up into a single macro.
#define __DEF_FPPOLxx_ARITH_ALL(sz)         \
        __DEF_FPPOLxx_SET_ZERO     (sz)     \
        __DEF_FPPOLxx_SET_ONE      (sz)     \
        __DEF_FPPOLxx_SET_TI       (sz)     \
        __DEF_FPPOLxx_SET_MAX      (sz)     \
        __DEF_FPPOLxx_MONIC_SET_MAX(sz)     \
        __DEF_FPPOLxx_SET          (sz)     \
        __DEF_FPPOLxx_SET_yy       (sz, 16) \
        __DEF_FPPOLxx_SET_yy       (sz, 32) \
        __DEF_FPPOLxx_SET_yy       (sz, 64) \
        __DEF_FPPOLxx_SET_MP       (sz)     \
        __DEF_FPPOLxx_GET_COEFF    (sz)     \
        __DEF_FPPOLxx_SET_COEFF    (sz)     \
        __DEF_FPPOLxx_FOLD_OR      (sz)     \
        __DEF_FPPOLxx_IS_ZERO      (sz)     \
        __DEF_FPPOLxx_IN_FP        (sz)     \
        __DEF_FPPOLxx_EQ           (sz)     \
        __DEF_FPPOLxx_IS_MONIC     (sz)     \
        __DEF_FPPOLxx_IS_VALID     (sz)     \
        __DEF_FPPOLxx_OPP          (sz)     \
        __DEF_FPPOLxx_ADD          (sz)     \
        __DEF_FPPOLxx_SUB          (sz)     \
        __DEF_FPPOLxx_SMUL         (sz)     \
        __DEF_FPPOLxx_SDIV         (sz)     \
        __DEF_FPPOLxx_SINV         (sz)     \
        __DEF_FPPOLxx_ADD_DISJOINT (sz)     \
        __DEF_FPPOLxx_MUL_TI       (sz)     \
        __DEF_FPPOLxx_DIV_TI       (sz)     \
        __DEF_FPPOLxx_MOD_TI       (sz)     \
        __DEF_FPPOLxx_GET_UI_MUL_TI(sz,  )  \
        __DEF_FPPOLxx_GET_UI_MUL_TI(sz, 1)  \
        __DEF_FPPOLxx_GET_UI       (sz,  )  \
        __DEF_FPPOLxx_GET_UI       (sz, 1)  \
        __DEF_FPPOLxx_GET_UI_MAX   (sz,  )  \
        __DEF_FPPOLxx_GET_UI_MAX   (sz, 1)  \
        __DEF_FPPOLxx_SET_UI       (sz,  )  \
        __DEF_FPPOLxx_SET_UI       (sz, 1)  \
        __DEF_FPPOLxx_SET_NEXT_UI  (sz,  )  \
        __DEF_FPPOLxx_SET_NEXT_UI  (sz, 1)

__DEF_FPPOLxx_ARITH_ALL(16)
__DEF_FPPOLxx_ARITH_ALL(32)
__DEF_FPPOLxx_ARITH_ALL(64)

#undef __DEF_FPPOLxx_SET_ZERO
#undef __DEF_FPPOLxx_SET_ONE
#undef __DEF_FPPOLxx_SET_TI
#undef __DEF_FPPOLxx_SET_MAX
#undef __DEF_FPPOLxx_MONIC_SET_MAX
#undef __DEF_FPPOLxx_SET
#undef __DEF_FPPOLxx_SET_yy
#undef __DEF_FPPOLxx_SET_yy
#undef __DEF_FPPOLxx_SET_yy
#undef __DEF_FPPOLxx_SET_MP
#undef __DEF_FPPOLxx_GET_COEFF
#undef __DEF_FPPOLxx_SET_COEFF
#undef __DEF_FPPOLxx_FOLD_OR
#undef __DEF_FPPOLxx_IS_ZERO
#undef __DEF_FPPOLxx_IN_FP
#undef __DEF_FPPOLxx_EQ
#undef __DEF_FPPOLxx_IS_MONIC
#undef __DEF_FPPOLxx_IS_VALID
#undef __DEF_FPPOLxx_OPP
#undef __DEF_FPPOLxx_ADD
#undef __DEF_FPPOLxx_SUB
#undef __DEF_FPPOLxx_SMUL
#undef __DEF_FPPOLxx_SDIV
#undef __DEF_FPPOLxx_SINV
#undef __DEF_FPPOLxx_ADD_DISJOINT
#undef __DEF_FPPOLxx_MUL_TI
#undef __DEF_FPPOLxx_DIV_TI
#undef __DEF_FPPOLxx_MOD_TI
#undef __DEF_FPPOLxx_GET_UI
#undef __DEF_FPPOLxx_GET_UI_MAX
#undef __DEF_FPPOLxx_SET_UI
#undef __DEF_FPPOLxx_SET_NEXT_UI
#undef __DEF_FPPOLxx_ARITH_ALL



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

// Get degree-i coefficient.
static inline
void fppol_get_coeff(fp_ptr r, fppol_srcptr p, unsigned i)
{ if ((int)i > p->deg) fp_set_zero(r);
  else                 fppol64_get_coeff(r, p->limbs[i>>6], i&0x3f); }

// Set degree-i coefficient.
void fppol_set_coeff(fppol_ptr r, fp_srcptr x, unsigned i);

// Test if zero.
static inline
int fppol_is_zero(fppol_srcptr p)
{ return p->deg == -1; }

// Test if in GF(p) (i.e., deg <= 0).
static inline
int fppol_in_fp(fppol_srcptr p)
{ return p->deg <= 0; }

// Test if equal.
int fppol_eq(fppol_srcptr p, fppol_srcptr q);

// Test if monic.
int fppol_is_monic(fppol_srcptr p);

// Test if valid representation.
int fppol_is_valid(fppol_srcptr p);

// Opposite.
void fppol_opp(fppol_ptr r, fppol_srcptr p);

// Addition.
void fppol_add(fppol_ptr r, fppol_srcptr p, fppol_srcptr q);

// Subtraction.
void fppol_sub(fppol_ptr r, fppol_srcptr p, fppol_srcptr q);

// Scalar multiplication.
void fppol_smul(fppol_ptr r, fppol_srcptr p, fp_srcptr x);

// Scalar division.
void fppol_sdiv(fppol_ptr r, fppol_srcptr p, fp_srcptr x);

// Coefficient-wise inversion.
void fppol_sinv(fppol_ptr r, fppol_srcptr p);

// Addition without overlapping coefficients (i.e., coefficient-wise OR).
void fppol_add_disjoint(fppol_ptr r, fppol_srcptr p, fppol_srcptr q);

// Multiplication by t^i (i.e., left shift).
void fppol_mul_ti(fppol_ptr r, fppol_srcptr p, unsigned i);

// Quotient of division by t^i (i.e., right shift).
void fppol_div_ti(fppol_ptr r, fppol_srcptr p, unsigned i);

// Remainder modulo t^i (i.e., i first coefficients).
void fppol_mod_ti(fppol_ptr r, fppol_srcptr p, unsigned i);

// Degree.
static inline
int fppol_deg(fppol_srcptr p)
{ return p->deg; }

#endif  /* __FPPOL_ARITH_H__ */
