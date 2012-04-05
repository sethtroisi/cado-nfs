#ifndef __F3POL_H__
#define __F3POL_H__

#include <stdint.h>

#include "macros.h"
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


// One.
#define __FP_ONE_0 1
#define __FP_ONE_1 0

// Conversion from integer.
#define __FP_SET_Z_0(x) ( (x)     & 0x1)
#define __FP_SET_Z_1(x) (((x)>>1) & 0x1)

// Opposite.
#define __FP_OPP_0(p0, p1) (p1)
#define __FP_OPP_1(p0, p1) (p0)

// Addition.
#define __FP_ADD_0(p0, p1, q0, q1) (((p0)|(q0)) ^ (((p0)|(p1)) & ((q0)|(q1))))
#define __FP_ADD_1(p0, p1, q0, q1) (((p1)|(q1)) ^ (((p0)|(p1)) & ((q0)|(q1))))

// Subtraction.
#define __FP_SUB_0(p0, p1, q0, q1) \
        __FP_ADD_0(p0, p1, __FP_OPP_0(q0, q1), __FP_OPP_1(q0, q1))
#define __FP_SUB_1(p0, p1, q0, q1) \
        __FP_ADD_1(p0, p1, __FP_OPP_0(q0, q1), __FP_OPP_1(q0, q1))

// Coefficient-wise inverse.
#define __FP_SINV_0(p0, p1) (p0)
#define __FP_SINV_1(p0, p1) (p1)

// Coefficient-wise multiplication.
#define __FP_SMUL_0(p0, p1, q0, q1) (((p0)&(q0)) | ((p1)&(q1)))
#define __FP_SMUL_1(p0, p1, q0, q1) (((p0)&(q1)) | ((p1)&(q0)))

// Coefficient-wise division.
#define __FP_SDIV_0(p0, p1, q0, q1) \
        __FP_SMUL_0(p0, p1, __FP_SINV_0(q0, q1), __FP_SINV_1(q0, q1))
#define __FP_SDIV_1(p0, p1, q0, q1) \
        __FP_SMUL_1(p0, p1, __FP_SINV_0(q0, q1), __FP_SINV_1(q0, q1))


// Generic operation.
#define __FP_OP(op, sz, r, args)                                    \
  do { uint##sz##_t __t[] = { (uint##sz##_t)(__FP_##op##_0 args),   \
                              (uint##sz##_t)(__FP_##op##_1 args) }; \
       r[0] = __t[0]; r[1] = __t[1]; } while (0)

// Definition of all coefficient-wise operations.
#define __FP_ONE(  sz, r)       __FP_OP(ONE,   sz, r, )
#define __FP_SET_Z(sz, r, x)    __FP_OP(SET_Z, sz, r, (x))
#define __FP_OPP(  sz, r, p)    __FP_OP(OPP,   sz, r, (p[0], p[1]))
#define __FP_ADD(  sz, r, p, q) __FP_OP(ADD,   sz, r, (p[0], p[1], q[0], q[1]))
#define __FP_SUB(  sz, r, p, q) __FP_OP(SUB,   sz, r, (p[0], p[1], q[0], q[1]))
#define __FP_SINV( sz, r, p)    __FP_OP(SINV,  sz, r, (p[0], p[1]))
#define __FP_SMUL( sz, r, p, q) __FP_OP(SMUL,  sz, r, (p[0], p[1], q[0], q[1]))
#define __FP_SDIV( sz, r, p, q) __FP_OP(SDIV,  sz, r, (p[0], p[1], q[0], q[1]))



/* Iteration and integer conversions.
 *****************************************************************************/

// Next polynomial, in lexicographical order.
#define __FP_SET_NEXT(sz, r, p, n)   \
  do {                               \
    uint##sz##_t __t0 = p[1] + 1;    \
    uint##sz##_t __t1 = p[1] ^ __t0; \
    __t1 = (__t1 >> 1) ^ __t1;       \
    r[0] =  p[0] ^ __t1;             \
    r[1] = (r[0] & __t1) ^ __t0;     \
  } while (0)

// Next monic polynomial, in lexicographical order.
#define __FP_MONIC_SET_NEXT(sz, r, p, n)     \
  do {                                       \
    uint##sz##_t __t0 = p[1] + 1;            \
    uint##sz##_t __t1 = p[1] ^ __t0;         \
    if (~__t1 & p[0]) { /* Not last coeff */ \
      __t1 = (__t1 >> 1) ^ __t1;             \
      r[0] =  p[0] ^ __t1;                   \
      r[1] = (r[0] & __t1) ^ __t0;           \
    }                                        \
    else {                  /* Last coeff */ \
      r[0] = p[0] ? __t1 + 1 : 1;            \
      r[1] = 0;                              \
    }                                        \
  } while (0)


// Internal look-up tables for conversion.
extern const uint8_t        __f3_get_ui_conv[];
extern const uint8_t  __f3_monic_get_ui_conv[];
extern const uint16_t       __f3_set_ui_conv[];
extern const uint16_t __f3_monic_set_ui_conv[];

// Conversion of an (n+m)-term polynomial, whose m most significant
// coefficients form a monic polynomial, to an unsigned int.
#define __FP_GET_UI(sz, r, p, n, m)                      \
  do {                                                   \
    r = 0;                                               \
    uint##sz##_t __t[2] = { p[0], p[1] };                \
    unsigned __i = 0, __j = 0, __b = MAX(n, n+m-4), __w; \
    while (__j < __b && __t[0] | __t[1]) {               \
      __w = (__t[0] & 0x1f) | ((__t[1] & 0x1f) << 5);    \
      r |= (uint64_t)__f3_get_ui_conv[__w] << __i;       \
      __i += 8; __j += 5; __t[0] >>= 5; __t[1] >>= 5;    \
    }                                                    \
    if (__t[0] | __t[1]) {                               \
      __w = (__t[0] & 0x1f) | ((__t[1] & 0x1f) << 5);    \
      r |= (uint64_t)__f3_monic_get_ui_conv[__w] << __i; \
    }                                                    \
  } while (0)

// Conversion of an (n+m)-term polynomial, whose m most significant
// coefficients form a monic polynomial, from an unsigned int.
#define __FP_SET_UI(sz, r, x, n, m)                       \
  do {                                                    \
    r[0] = r[1] = 0;                                      \
    unsigned __i = 0, __j = 0, __b = MAX(n, n+m-4), __w;  \
    while (__j < __b && x) {                              \
      __w = (x >> 8) & 0xff;                              \
      if (UNLIKELY(__w >= 243))           /* 243 = 3^5 */ \
        return 0;                                         \
      __w = __f3_set_ui_conv[__w];                        \
      r[0] |= ((uint##sz##_t)(__w & 0x1f) << __j);        \
      r[1] |= ((uint##sz##_t)(__w >> 5  ) << __j);        \
      __i += 8; __j += 5; x >>= 8;                        \
    }                                                     \
    if (x) {                                              \
      __w = (x >> 8) & 0xff;                              \
      if (UNLIKELY(__w >= 122)) /* 122 = (3^5-1)/2 + 1 */ \
        return 0;                                         \
      __w = __f3_monic_set_ui_conv[__w];                  \
      r[0] |= ((uint##sz##_t)(__w & 0x1f) << __j);        \
      r[1] |= ((uint##sz##_t)(__w >> 5  ) << __j);        \
      __i += 8; __j += 5; x >>= 8;                        \
    }                                                     \
  } while (0)



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
