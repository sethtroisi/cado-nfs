#ifndef __F2POL_H__
#define __F2POL_H__

#include <stdint.h>

#include "macros.h"
#include "cppmeta.h"



/* Basic bitsliced arithmetic.
 *****************************************************************************/

// Number of elements in the field.
#define __FP_SIZE 2

// Field characteristic.
#define __FP_CHAR 2

// Number of bits per element.
#define __FP_BITS 1

// Test if valid representation.
#define __FP_IS_VALID(sz, p) 1


// One.
#define __FP_ONE_0 1

// Largest element in the base field.
#define __FP_MAX_0 1

// Opposite.
#define __FP_OPP_0(p0) (p0)

// Addition.
#define __FP_ADD_0(p0, q0) ((p0) ^ (q0))

// Subtraction.
#define __FP_SUB_0(p0, q0) __FP_ADD_0(p0, __FP_OPP_0(q0))

// Coefficient-wise inverse.
#define __FP_SINV_0(p0) (p0)

// Coefficient-wise multiplication.
#define __FP_SMUL_0(p0, q0) ((p0) & (q0))

// Coefficient-wise division.
#define __FP_SDIV_0(p0, q0) __FP_SMUL_0(p0, __FP_SINV_0(q0))


// Generic constant definition.
#define __FP_CST(cst, sz, r) \
  do { r[0] = __FP_##cst##_0; } while (0)

// Generic unary operation.
#define __FP_UOP(op, sz, r, p) \
  do { r[0] = __FP_##op##_0(p[0]); } while (0)

// Generic binary operation.
#define __FP_BOP(op, sz, r, p, q) \
  do { r[0] = __FP_##op##_0(p[0], q[0]); } while (0)

// Definition of all coefficient-wise operations.
#define __FP_ONE( sz, r)       __FP_CST(ONE,  sz, r)
#define __FP_MAX( sz, r)       __FP_CST(MAX,  sz, r)
#define __FP_OPP( sz, r, p)    __FP_UOP(OPP,  sz, r, p)
#define __FP_ADD( sz, r, p, q) __FP_BOP(ADD,  sz, r, p, q)
#define __FP_SUB( sz, r, p, q) __FP_BOP(SUB,  sz, r, p, q)
#define __FP_SINV(sz, r, p)    __FP_UOP(SINV, sz, r, p)
#define __FP_SMUL(sz, r, p, q) __FP_BOP(SMUL, sz, r, p, q)
#define __FP_SDIV(sz, r, p, q) __FP_BOP(SDIV, sz, r, p, q)



/* Integer conversions.
 *****************************************************************************/

// Conversion of an n-term polynomial to an unsigned int.
#define __FP_GET_UI(sz, r, p, n) \
  do { r = (unsigned)p[0]; } while (0)

// Conversion of an n-term polynomial from an unsigned int.
#define __FP_SET_UI(sz, next, r, x, n) \
  do { SWITCH(next, EMPTY, ++x;)       \
       r[0] = (uint##sz##_t)x; } while (0)

// Conversions to/from an unsigned int in the case of monic polynomials.
#define __FP_MONIC_GET_UI __FP_GET_UI
#define __FP_MONIC_SET_UI __FP_SET_UI



/* Multiplications.
 *****************************************************************************/

// Naive <sr> <- <sp> x <sq> multiplication.
#define __FP_MUL_xx_yyxzz_NAIVE(sr, sp, sq, r, p, q, op)    \
  do {                                                      \
    fppol##sp##_t __m;                                      \
    fppol##sr##_t __t = {0};                                \
    for (unsigned __i = sq; __i--; ) {                      \
      __m[0] = -(uint##sp##_t)((q[0] >> __i) & 1);          \
      __t[0] = __FP_ADD_0(__t[0] << 1, __m[0] & p[0]);      \
    }                                                       \
    CAT(fppol##sr##_, SWITCH(op, OP(r, __t), set(r, __t))); \
  } while (0)

#define __FP_MUL_16_16x16_NAIVE(            r, p, q, op) \
        __FP_MUL_xx_yyxzz_NAIVE(16, 16, 16, r, p, q, op)
#define __FP_MUL_32_16x16_NAIVE(            r, p, q, op) \
        __FP_MUL_xx_yyxzz_NAIVE(32, 16, 16, r, p, q, op)
#define __FP_MUL_32_32x16_NAIVE(            r, p, q, op) \
        __FP_MUL_xx_yyxzz_NAIVE(32, 32, 16, r, p, q, op)
#define __FP_MUL_32_32x32_NAIVE(            r, p, q, op) \
        __FP_MUL_xx_yyxzz_NAIVE(32, 32, 32, r, p, q, op)
#define __FP_MUL_64_32x16_NAIVE(            r, p, q, op) \
        __FP_MUL_xx_yyxzz_NAIVE(64, 32, 16, r, p, q, op)
#define __FP_MUL_64_32x32_NAIVE(            r, p, q, op) \
        __FP_MUL_xx_yyxzz_NAIVE(64, 32, 32, r, p, q, op)
#define __FP_MUL_64_64x16_NAIVE(            r, p, q, op) \
        __FP_MUL_xx_yyxzz_NAIVE(64, 64, 16, r, p, q, op)
#define __FP_MUL_64_64x32_NAIVE(            r, p, q, op) \
        __FP_MUL_xx_yyxzz_NAIVE(64, 64, 32, r, p, q, op)
#define __FP_MUL_64_64x64_NAIVE(            r, p, q, op) \
        __FP_MUL_xx_yyxzz_NAIVE(64, 64, 64, r, p, q, op)


// Naive 128 <- 64 x <sq> multiplication.
#define __FP_MUL_128_64xxx_NAIVE(sq, rh, rl, p, q, op)    \
  do {                                                    \
    fppol64_t __m;                                        \
    fppol64_t __h = {0}, __l = {0};                       \
    for (unsigned __i = sq; __i--; ) {                    \
      __m[0] = -(uint64_t)((q[0] >> __i) & 1);            \
      __h[0] = (__h[0] << 1) | (__l[0] >> 63);            \
      __l[0] = __FP_ADD_0(__l[0] << 1, __m[0] & p[0]);    \
    }                                                     \
    CAT(fppol64_, SWITCH(op, OP(rh, __h), set(rh, __h))); \
    CAT(fppol64_, SWITCH(op, OP(rl, __l), set(rl, __l))); \
  } while (0)

#define __FP_MUL_128_64x16_NAIVE(    rh, rl, p, q, op) \
        __FP_MUL_128_64xxx_NAIVE(16, rh, rl, p, q, op)
#define __FP_MUL_128_64x32_NAIVE(    rh, rl, p, q, op) \
        __FP_MUL_128_64xxx_NAIVE(32, rh, rl, p, q, op)
#define __FP_MUL_128_64x64_NAIVE(    rh, rl, p, q, op) \
        __FP_MUL_128_64xxx_NAIVE(64, rh, rl, p, q, op)


// Select multiplication algorithms.
#define __FP_MUL_16_16x16  __FP_MUL_16_16x16_NAIVE
#define __FP_MUL_32_16x16  __FP_MUL_32_16x16_NAIVE
#define __FP_MUL_32_32x16  __FP_MUL_32_32x16_NAIVE
#define __FP_MUL_32_32x32  __FP_MUL_32_32x32_NAIVE
#define __FP_MUL_64_32x16  __FP_MUL_64_32x16_NAIVE
#define __FP_MUL_64_32x32  __FP_MUL_64_32x32_NAIVE
#define __FP_MUL_64_64x16  __FP_MUL_64_64x16_NAIVE
#define __FP_MUL_64_64x32  __FP_MUL_64_64x32_NAIVE
#define __FP_MUL_64_64x64  __FP_MUL_64_64x64_NAIVE
#define __FP_MUL_128_64x16 __FP_MUL_128_64x16_NAIVE
#define __FP_MUL_128_64x32 __FP_MUL_128_64x32_NAIVE
#define __FP_MUL_128_64x64 __FP_MUL_128_64x64_NAIVE

#endif  /* __F2POL_H__ */
