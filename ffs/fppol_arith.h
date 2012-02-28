// Always include "fppol.h" first so that everything gets loaded in the correct
// order.
#include "fppol.h"

#ifndef __FPPOL_ARITH_H__
#define __FPPOL_ARITH_H__

#include "macros.h"



/* Base field elements.
 *****************************************************************************/

// Set to zero.
static inline
void fp_set_zero(fp_ptr r)
{ for (unsigned k = 0; k < __FP_BITS; ++k) r[k] = 0; }


// Set to one.
static inline
void fp_set_one(fp_ptr r)
{ r[0] = 1;
  for (unsigned k = 1; k < __FP_BITS; ++k) r[k] = 0; }


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
#define __DEF_FPPOLxx_SET_ONE(sz)                        \
  static inline                                          \
  void fppol##sz##_set_one(fppol##sz##_ptr r)            \
  { r[0] = 1;                                            \
    for (unsigned k = 1; k < __FP_BITS; ++k) r[k] = 0; }


// Set to t^i.
// Generic prototype:
//   void fppol<sz>_set_ti(fppol<sz>_ptr r, unsigned i);
#define __DEF_FPPOLxx_SET_TI(sz)                         \
  static inline                                          \
  void fppol##sz##_set_ti(fppol##sz##_ptr r, unsigned i) \
  { r[0] = (uint##sz##_t)1 << i;                         \
    for (unsigned k = 1; k < __FP_BITS; ++k) r[k] = 0; }


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
  { fppol##sz##_set_64(r, p->limb[0]);                      \
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


// Left shift.
// Generic prototype:
//   void fppol<sz>_shl(fppol<sz>_ptr r, fppol<sz>_srcptr p, unsigned i);
#define __DEF_FPPOLxx_SHL(sz)                                                 \
  static inline                                                               \
  void fppol##sz##_shl(fppol##sz##_ptr r, fppol##sz##_srcptr p, unsigned i)   \
  { for (unsigned k = 0; k < __FP_BITS; ++k) r[k] = i < sz ? p[k] << i : 0; }


// Right shift.
// Generic prototype:
//   void fppol<sz>_shr(fppol<sz>_ptr r, fppol<sz>_srcptr p, unsigned i);
#define __DEF_FPPOLxx_SHR(sz)                                                 \
  static inline                                                               \
  void fppol##sz##_shr(fppol##sz##_ptr r, fppol##sz##_srcptr p, unsigned i)   \
  { for (unsigned k = 0; k < __FP_BITS; ++k) r[k] = i < sz ? p[k] >> i : 0; }


// All definitions bundled up into a single macro.
#define __DEF_FPPOLxx_ARITH_ALL(sz)     \
        __DEF_FPPOLxx_SET_ZERO (sz)     \
        __DEF_FPPOLxx_SET_ONE  (sz)     \
        __DEF_FPPOLxx_SET_TI   (sz)     \
        __DEF_FPPOLxx_SET      (sz)     \
        __DEF_FPPOLxx_SET_yy   (sz, 16) \
        __DEF_FPPOLxx_SET_yy   (sz, 32) \
        __DEF_FPPOLxx_SET_yy   (sz, 64) \
        __DEF_FPPOLxx_SET_MP   (sz)     \
        __DEF_FPPOLxx_GET_COEFF(sz)     \
        __DEF_FPPOLxx_SET_COEFF(sz)     \
        __DEF_FPPOLxx_FOLD_OR  (sz)     \
        __DEF_FPPOLxx_IS_ZERO  (sz)     \
        __DEF_FPPOLxx_EQ       (sz)     \
        __DEF_FPPOLxx_IS_MONIC (sz)     \
        __DEF_FPPOLxx_IS_VALID (sz)     \
        __DEF_FPPOLxx_OPP      (sz)     \
        __DEF_FPPOLxx_ADD      (sz)     \
        __DEF_FPPOLxx_SUB      (sz)     \
        __DEF_FPPOLxx_SMUL     (sz)     \
        __DEF_FPPOLxx_SDIV     (sz)     \
        __DEF_FPPOLxx_SHL      (sz)     \
        __DEF_FPPOLxx_SHR      (sz)

__DEF_FPPOLxx_ARITH_ALL(16)
__DEF_FPPOLxx_ARITH_ALL(32)
__DEF_FPPOLxx_ARITH_ALL(64)

#undef __DEF_FPPOLxx_SET_ZERO
#undef __DEF_FPPOLxx_SET_ONE
#undef __DEF_FPPOLxx_SET_TI
#undef __DEF_FPPOLxx_SET
#undef __DEF_FPPOLxx_SET_yy
#undef __DEF_FPPOLxx_SET_yy
#undef __DEF_FPPOLxx_SET_yy
#undef __DEF_FPPOLxx_SET_MP
#undef __DEF_FPPOLxx_GET_COEFF
#undef __DEF_FPPOLxx_SET_COEFF
#undef __DEF_FPPOLxx_FOLD_OR
#undef __DEF_FPPOLxx_IS_ZERO
#undef __DEF_FPPOLxx_EQ
#undef __DEF_FPPOLxx_IS_MONIC
#undef __DEF_FPPOLxx_IS_VALID
#undef __DEF_FPPOLxx_OPP
#undef __DEF_FPPOLxx_ADD
#undef __DEF_FPPOLxx_SUB
#undef __DEF_FPPOLxx_SMUL
#undef __DEF_FPPOLxx_SDIV
#undef __DEF_FPPOLxx_SHL
#undef __DEF_FPPOLxx_SHR
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

// Set to something coded in an uint64_t.
void fppol_set_ui(fppol_ptr r, uint64_t x);

// Get degree-i coefficient.
static inline
void fppol_get_coeff(fp_ptr r, fppol_srcptr p, unsigned i)
{ if ((int)i > p->deg) fp_set_zero(r);
  else                 fppol64_get_coeff(r, p->limb[i>>6], i&0x3f); }

// Set degree-i coefficient.
void fppol_set_coeff(fppol_ptr r, fp_srcptr x, unsigned i);

// Test if zero.
static inline
int fppol_is_zero(fppol_srcptr p)
{ return p->deg == -1; }

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

// Left shift.
void fppol_shl(fppol_ptr r, fppol_srcptr p, unsigned i);

// Right shift.
void fppol_shr(fppol_ptr r, fppol_srcptr p, unsigned i);

// Degree.
static inline
int fppol_deg(fppol_srcptr p)
{ return p->deg; }

#endif  /* __FPPOL_ARITH_H__ */
