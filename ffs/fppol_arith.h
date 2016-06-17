// Always include "fppol.h" first so that everything gets loaded in the correct
// order.
#include "fppol.h"

#ifndef __FPPOL_ARITH_H__
#define __FPPOL_ARITH_H__



/* Fixed-size polynomials.
 *****************************************************************************/

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


// Difference between two polynomials (i.e., coefficient-wise XOR).
// Generic prototype:
//   void fppol<sz>_diff(fppol<sz>_ptr    r,
//                       fppol<sz>_srcptr p, fppol<sz>_srcptr q);
#define __DEF_FPPOLxx_DIFF(sz)                                      \
  static inline                                                     \
  void fppol##sz##_diff(fppol##sz##_ptr r,                          \
                        fppol##sz##_srcptr p, fppol##sz##_srcptr q) \
  { for (unsigned k = 0; k < __FP_BITS; ++k) r[k] = p[k] ^ q[k]; }


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


// All definitions bundled up into a single macro.
#define __DEF_FPPOLxx_ARITH_ALL(sz)    \
        __DEF_FPPOLxx_OPP         (sz) \
        __DEF_FPPOLxx_ADD         (sz) \
        __DEF_FPPOLxx_SUB         (sz) \
        __DEF_FPPOLxx_SMUL        (sz) \
        __DEF_FPPOLxx_SDIV        (sz) \
        __DEF_FPPOLxx_SINV        (sz) \
        __DEF_FPPOLxx_ADD_DISJOINT(sz) \
        __DEF_FPPOLxx_DIFF        (sz) \
        __DEF_FPPOLxx_MUL_TI      (sz) \
        __DEF_FPPOLxx_DIV_TI      (sz) \
        __DEF_FPPOLxx_MOD_TI      (sz)

__DEF_FPPOLxx_ARITH_ALL( 8)
__DEF_FPPOLxx_ARITH_ALL(16)
__DEF_FPPOLxx_ARITH_ALL(32)
__DEF_FPPOLxx_ARITH_ALL(64)

#undef __DEF_FPPOLxx_OPP
#undef __DEF_FPPOLxx_ADD
#undef __DEF_FPPOLxx_SUB
#undef __DEF_FPPOLxx_SMUL
#undef __DEF_FPPOLxx_SDIV
#undef __DEF_FPPOLxx_SINV
#undef __DEF_FPPOLxx_ADD_DISJOINT
#undef __DEF_FPPOLxx_DIFF
#undef __DEF_FPPOLxx_MUL_TI
#undef __DEF_FPPOLxx_DIV_TI
#undef __DEF_FPPOLxx_MOD_TI
#undef __DEF_FPPOLxx_ARITH_ALL



/* Multiprecision polynomials.
 *****************************************************************************/

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

#endif  /* __FPPOL_ARITH_H__ */
