#ifndef __F3POL_H__
#define __F3POL_H__

#include <stdint.h>

#include "cppmeta.h"



/* Basic bitsliced arithmetic.
 *****************************************************************************/

// Number of elements in the field.
#define __FP_SIZE 3

// Field characteristic.
#define __FP_CHAR 3

// Number of bits per element.
#define __FP_BITS 2

// Test if valid representation.
#define __FP_IS_VALID(sz, p) (!(p[0] & p[1]))


// Opposite.
#define __FP_OPP_0(p0, p1) (p1)
#define __FP_OPP_1(p0, p1) (p0)

#define __FP_OPP(sz, r, p)                              \
  do { uint##sz##_t __t[] = { __FP_OPP_0(p[0], p[1]),   \
                              __FP_OPP_1(p[0], p[1]) }; \
       r[0] = __t[0]; r[1] = __t[1]; } while (0)


// Addition.
#define __FP_ADD_0(p0, p1, q0, q1) (((p1)&(q1)) | (~((p1)|(q1)) & ((p0)^(q0))))
#define __FP_ADD_1(p0, p1, q0, q1) (((p0)&(q0)) | (~((p0)|(q0)) & ((p1)^(q1))))

#define __FP_ADD(sz, r, p, q)                                       \
  do { uint##sz##_t __t[] = { __FP_ADD_0(p[0], p[1], q[0], q[1]),   \
                              __FP_ADD_1(p[0], p[1], q[0], q[1]) }; \
       r[0] = __t[0]; r[1] = __t[1]; } while (0)


// Subtraction.
#define __FP_SUB_0(p0, p1, q0, q1) __FP_ADD_0(p0, p1, __FP_OPP_0(q0, q1), \
                                                      __FP_OPP_1(q0, q1))
#define __FP_SUB_1(p0, p1, q0, q1) __FP_ADD_1(p0, p1, __FP_OPP_0(q0, q1), \
                                                      __FP_OPP_1(q0, q1))

#define __FP_SUB(sz, r, p, q)                                       \
  do { uint##sz##_t __t[] = { __FP_SUB_0(p[0], p[1], q[0], q[1]),   \
                              __FP_SUB_1(p[0], p[1], q[0], q[1]) }; \
       r[0] = __t[0]; r[1] = __t[1]; } while (0)


// Coefficient-wise inverse.
#define __FP_SINV_0(p0, p1) (p0)
#define __FP_SINV_1(p0, p1) (p1)

#define __FP_SINV(sz, r, p)                              \
  do { uint##sz##_t __t[] = { __FP_SINV_0(p[0], p[1]),   \
                              __FP_SINV_1(p[0], p[1]) }; \
       r[0] = __t[0]; r[1] = __t[1]; } while (0)


// Coefficient-wise multiplication.
#define __FP_SMUL_0(p0, p1, q0, q1) (((p0)&(q0)) | ((p1)&(q1)))
#define __FP_SMUL_1(p0, p1, q0, q1) (((p0)&(q1)) | ((p1)&(q0)))

#define __FP_SMUL(sz, r, p, q)                                       \
  do { uint##sz##_t __t[] = { __FP_SMUL_0(p[0], p[1], q[0], q[1]),   \
                              __FP_SMUL_1(p[0], p[1], q[0], q[1]) }; \
       r[0] = __t[0]; r[1] = __t[1]; } while (0)


// Coefficient-wise division.
#define __FP_SDIV_0(p0, p1, q0, q1) __FP_SMUL_0(p0, p1, __FP_SINV_0(q0, q1), \
                                                        __FP_SINV_1(q0, q1))
#define __FP_SDIV_1(p0, p1, q0, q1) __FP_SMUL_1(p0, p1, __FP_SINV_0(q0, q1), \
                                                        __FP_SINV_1(q0, q1))

#define __FP_SDIV(sz, r, p, q)                                       \
  do { uint##sz##_t __t[] = { __FP_SDIV_0(p[0], p[1], q[0], q[1]),   \
                              __FP_SDIV_1(p[0], p[1], q[0], q[1]) }; \
       r[0] = __t[0]; r[1] = __t[1]; } while (0)



/* Multiplications.
 *****************************************************************************/

// Naive <sr> <- <sp> x <sq> multiplication.
#define __FP_MUL_xx_yyxzz_NAIVE(sr, sp, sq, r, p, q, op)                 \
  do {                                                                   \
    fppol##sp##_t __m;                                                   \
    fppol##sr##_t __t = {0,0}, __s;                                      \
    for (unsigned __i = sq; __i--; ) {                                   \
      __m[0] = -(uint##sp##_t)((q[0] >> __i) & 1); __t[0] <<= 1;         \
      __m[1] = -(uint##sp##_t)((q[1] >> __i) & 1); __t[1] <<= 1;         \
      __s[0] = __FP_ADD_0(__t[0], __t[1], __m[0] & p[0], __m[0] & p[1]); \
      __s[1] = __FP_ADD_1(__t[0], __t[1], __m[0] & p[0], __m[0] & p[1]); \
      __t[0] = __FP_SUB_0(__s[0], __s[1], __m[1] & p[0], __m[1] & p[1]); \
      __t[1] = __FP_SUB_1(__s[0], __s[1], __m[1] & p[0], __m[1] & p[1]); \
    }                                                                    \
    CAT(fppol##sr##_, SWITCH(op, OP(r, __t), set(r, __t)));              \
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
#define __FP_MUL_128_64xxx_NAIVE(sq, rh, rl, p, q, op)                   \
  do {                                                                   \
    fppol64_t __m;                                                       \
    fppol64_t __h = {0,0}, __l = {0,0}, __s;                             \
    for (unsigned __i = sq; __i--; ) {                                   \
      __m[0] = -(uint64_t)((q[0] >> __i) & 1);                           \
      __m[1] = -(uint64_t)((q[1] >> __i) & 1);                           \
      __h[0] = (__h[0] << 1) | (__l[0] >> 63); __l[0] <<= 1;             \
      __h[1] = (__h[1] << 1) | (__l[1] >> 63); __l[1] <<= 1;             \
      __s[0] = __FP_ADD_0(__l[0], __l[1], __m[0] & p[0], __m[0] & p[1]); \
      __s[1] = __FP_ADD_1(__l[0], __l[1], __m[0] & p[0], __m[0] & p[1]); \
      __l[0] = __FP_SUB_0(__s[0], __s[1], __m[1] & p[0], __m[1] & p[1]); \
      __l[1] = __FP_SUB_1(__s[0], __s[1], __m[1] & p[0], __m[1] & p[1]); \
    }                                                                    \
    CAT(fppol64_, SWITCH(op, OP(rh, __h), set(rh, __h)));                \
    CAT(fppol64_, SWITCH(op, OP(rl, __l), set(rl, __l)));                \
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

#endif  /* __F3POL_H__ */
