#ifndef __F2POL_H__
#define __F2POL_H__

#include <stdint.h>

#include "macros.h"
#include "cppmeta.h"

#ifdef HAVE_GF2X
#include <gf2x/gf2x_mul1.h>
#endif



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

// Conversion from integer.
#define __FP_SET_Z_0(x) ((x) & 0x1)

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


// Generic operation.
#define __FP_OP(op, sz, r, args) \
  do { r[0] = __FP_##op##_0 args; } while (0)

// Definition of all coefficient-wise operations.
#define __FP_ONE(  sz, r)       __FP_OP(ONE,   sz, r, )
#define __FP_SET_Z(sz, r, x)    __FP_OP(SET_Z, sz, r, (x))
#define __FP_OPP(  sz, r, p)    __FP_OP(OPP,   sz, r, (p[0]))
#define __FP_ADD(  sz, r, p, q) __FP_OP(ADD,   sz, r, (p[0], q[0]))
#define __FP_SUB(  sz, r, p, q) __FP_OP(SUB,   sz, r, (p[0], q[0]))
#define __FP_SINV( sz, r, p)    __FP_OP(SINV,  sz, r, (p[0]))
#define __FP_SMUL( sz, r, p, q) __FP_OP(SMUL,  sz, r, (p[0], q[0]))
#define __FP_SDIV( sz, r, p, q) __FP_OP(SDIV,  sz, r, (p[0], q[0]))



/* Iteration and integer conversions.
 *****************************************************************************/

// Next polynomial, in lexicographical order.
#define __FP_SET_NEXT(sz, r, p, n) \
  do { r[0] = p[0] + 1; } while (0)

// Next monic polynomial, in lexicographical order.
#define __FP_MONIC_SET_NEXT __FP_SET_NEXT


// Conversion of an (n+m)-term polynomial, whose m most significant
// coefficients form a monic polynomial, to an unsigned int.
#define __FP_GET_UI(sz, r, p, n, m) \
  do { r = (uint64_t)p[0]; } while (0)

// Conversion of an (n+m)-term polynomial, whose m most significant
// coefficients form a monic polynomial, from an unsigned int.
#define __FP_SET_UI(sz, r, x, n, m) \
  do { r[0] = (uint##sz##_t)x; } while (0)

// Conversion of an polynomial to an uint64_t (evaluation in 2^__FP_BITS)
#define __FP_GET_UI_SPARSE(sz, r, p) \
  do { r = (uint64_t)p[0]; } while (0)

// Conversion of an uint64_t (evaluation in 2^__FP_BITS) in polynomial
#define __FP_SET_UI_SPARSE(sz, r, x) \
  do { ASSERT(sz==64 || x>>sz == 0); r[0] = (uint##sz##_t)x; } while (0)



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

#define __FP_MUL_8_8x8_NAIVE(               r, p, q, op) \
        __FP_MUL_xx_yyxzz_NAIVE( 8,  8,  8, r, p, q, op)
#define __FP_MUL_16_8x8_NAIVE(              r, p, q, op) \
        __FP_MUL_xx_yyxzz_NAIVE(16,  8,  8, r, p, q, op)
#define __FP_MUL_16_16x8_NAIVE(             r, p, q, op) \
        __FP_MUL_xx_yyxzz_NAIVE(16, 16,  8, r, p, q, op)
#define __FP_MUL_16_16x16_NAIVE(            r, p, q, op) \
        __FP_MUL_xx_yyxzz_NAIVE(16, 16, 16, r, p, q, op)
#define __FP_MUL_32_16x8_NAIVE(             r, p, q, op) \
        __FP_MUL_xx_yyxzz_NAIVE(32, 16,  8, r, p, q, op)
#define __FP_MUL_32_16x16_NAIVE(            r, p, q, op) \
        __FP_MUL_xx_yyxzz_NAIVE(32, 16, 16, r, p, q, op)
#define __FP_MUL_32_32x8_NAIVE(             r, p, q, op) \
        __FP_MUL_xx_yyxzz_NAIVE(32, 32,  8, r, p, q, op)
#define __FP_MUL_32_32x16_NAIVE(            r, p, q, op) \
        __FP_MUL_xx_yyxzz_NAIVE(32, 32, 16, r, p, q, op)
#define __FP_MUL_32_32x32_NAIVE(            r, p, q, op) \
        __FP_MUL_xx_yyxzz_NAIVE(32, 32, 32, r, p, q, op)
#define __FP_MUL_64_32x8_NAIVE(             r, p, q, op) \
        __FP_MUL_xx_yyxzz_NAIVE(64, 32,  8, r, p, q, op)
#define __FP_MUL_64_32x16_NAIVE(            r, p, q, op) \
        __FP_MUL_xx_yyxzz_NAIVE(64, 32, 16, r, p, q, op)
#define __FP_MUL_64_32x32_NAIVE(            r, p, q, op) \
        __FP_MUL_xx_yyxzz_NAIVE(64, 32, 32, r, p, q, op)
#define __FP_MUL_64_64x8_NAIVE(             r, p, q, op) \
        __FP_MUL_xx_yyxzz_NAIVE(64, 64,  8, r, p, q, op)
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

#define __FP_MUL_128_64x8_NAIVE(     rh, rl, p, q, op) \
        __FP_MUL_128_64xxx_NAIVE( 8, rh, rl, p, q, op)
#define __FP_MUL_128_64x16_NAIVE(    rh, rl, p, q, op) \
        __FP_MUL_128_64xxx_NAIVE(16, rh, rl, p, q, op)
#define __FP_MUL_128_64x32_NAIVE(    rh, rl, p, q, op) \
        __FP_MUL_128_64xxx_NAIVE(32, rh, rl, p, q, op)
#define __FP_MUL_128_64x64_NAIVE(    rh, rl, p, q, op) \
        __FP_MUL_128_64xxx_NAIVE(64, rh, rl, p, q, op)

#ifdef HAVE_GF2X
#define __FP_MUL_128_64x64_GF2X(rh,rl,p,q,op)  \
do {    \
    fppol64_t c[2];\
    gf2x_mul1((unsigned long *)&c[0], p[0], q[0]);\
    CAT(fppol64_, SWITCH(op, OP(rh, c[1]), set(rh, c[1])));\
    CAT(fppol64_, SWITCH(op, OP(rl, c[0]), set(rl, c[0])));\
} while (0)

#define __FP_MUL_64_32x32_GF2X(r,p,q,op)  \
do {    \
    fppol64_t c[2];\
    unsigned long a = p[0], b = q[0]; \
    gf2x_mul1((unsigned long *)&c[0], a, b);\
    CAT(fppol64_, SWITCH(op, OP(r, c[0]), set(r, c[0])));\
} while (0)

#define __FP_MUL_32_32x32_GF2X(r,p,q,op)  \
do {    \
    fppol64_t c[2];\
    unsigned long a = p[0], b = q[0]; \
    gf2x_mul1((unsigned long *)&c[0], a, b);\
    fppol32_t cc; \
    fppol32_set_64(cc, c[0]); \
    CAT(fppol32_, SWITCH(op, OP(r, cc), set(r, cc)));\
} while (0)
#endif

// Select multiplication algorithms.
#define __FP_MUL_8_8x8     __FP_MUL_8_8x8_NAIVE
#define __FP_MUL_16_8x8    __FP_MUL_16_8x8_NAIVE
#define __FP_MUL_16_16x8   __FP_MUL_16_16x8_NAIVE
#define __FP_MUL_16_16x16  __FP_MUL_16_16x16_NAIVE
#define __FP_MUL_32_16x8   __FP_MUL_32_16x8_NAIVE
#define __FP_MUL_32_16x16  __FP_MUL_32_16x16_NAIVE
#define __FP_MUL_32_32x8   __FP_MUL_32_32x8_NAIVE
#define __FP_MUL_32_32x16  __FP_MUL_32_32x16_NAIVE
#define __FP_MUL_64_32x8   __FP_MUL_64_32x8_NAIVE
#define __FP_MUL_64_32x16  __FP_MUL_64_32x16_NAIVE
#define __FP_MUL_64_64x8   __FP_MUL_64_64x8_NAIVE
#define __FP_MUL_64_64x16  __FP_MUL_64_64x16_NAIVE
#define __FP_MUL_64_64x32  __FP_MUL_64_64x32_NAIVE
#define __FP_MUL_64_64x64  __FP_MUL_64_64x64_NAIVE
#define __FP_MUL_128_64x8  __FP_MUL_128_64x8_NAIVE
#define __FP_MUL_128_64x16 __FP_MUL_128_64x16_NAIVE
#define __FP_MUL_128_64x32 __FP_MUL_128_64x32_NAIVE
#ifdef HAVE_GF2X
#define __FP_MUL_128_64x64 __FP_MUL_128_64x64_GF2X
#define __FP_MUL_64_32x32  __FP_MUL_64_32x32_GF2X
#define __FP_MUL_32_32x32  __FP_MUL_32_32x32_GF2X
#else
#define __FP_MUL_128_64x64 __FP_MUL_128_64x64_NAIVE
#define __FP_MUL_64_32x32  __FP_MUL_64_32x32_NAIVE
#define __FP_MUL_32_32x32  __FP_MUL_32_32x32_NAIVE
#endif

#endif  /* __F2POL_H__ */
