/* This file is part of the gf2x library.

   Copyright 2010, 2012, 2013, 2015
   Richard Brent, Pierrick Gaudry, Emmanuel Thome', Paul Zimmermann
   Nicolas Estibals (for this file)

   This program is free software; you can redistribute it and/or modify it
   under the terms of either:
    - If the archive contains a file named toom-gpl.c (not a trivial
    placeholder), the GNU General Public License as published by the Free
    Software Foundation; either version 3 of the License, or (at your
    option) any later version.
    - If the archive contains a file named toom-gpl.c which is a trivial
    placeholder, the GNU Lesser General Public License as published by
    the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.
   
   This program is distributed in the hope that it will be useful, but WITHOUT
   ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
   FITNESS FOR A PARTICULAR PURPOSE.  See the license text for more details.
   
   You should have received a copy of the GNU General Public License as
   well as the GNU Lesser General Public License along with this program;
   see the files COPYING and COPYING.LIB.  If not, write to the Free
   Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
   02110-1301, USA.
*/
#ifndef GF2X_MUL7_H_
#define GF2X_MUL7_H_

#include "gf2x.h"
/* All gf2x source files for lowlevel functions must include gf2x-small.h
 * This is mandatory for the tuning mechanism. */
#include "gf2x/gf2x-small.h"

#if GF2X_WORDSIZE != 64
#error "This code is for 64-bit only"
#endif

#ifndef GF2X_HAVE_PCLMUL_SUPPORT
#error "This code needs pclmul support"
#endif

/* TODO: if somebody comes up with a neat way to improve the interface so
 * as to remove the false dependency on pclmul, that would be nice.
 */
static inline __m128i
GF2X_FUNC(mul7cl_mul1) (unsigned long a, unsigned long b)
{
    return _mm_clmulepi64_si128(_gf2x_mm_setr_epi64(a,0), _gf2x_mm_setr_epi64(b,0), 0);
}


/* variant with 22 multiplications */
GF2X_STORAGE_CLASS_mul7
void gf2x_mul7 (unsigned long *c, const unsigned long *a, const unsigned long *b)
{
#define PXOR(lop, rop) _mm_xor_si128((lop), (rop))
#define PXOR3(op1, op2, op3) PXOR(op1, PXOR(op2, op3))
#define PXOR4(op1, op2, op3, op4) PXOR(op1, PXOR3(op2, op3, op4))
#define PXOR5(op1, op2, op3, op4, op5) PXOR(op1, PXOR4(op2, op3, op4, op5))
#define PXOR6(op1, op2, op3, op4, op5, op6) PXOR(op1, PXOR5(op2, op3, op4, op5, op6))
#define PZERO    _mm_setzero_si128()
    /* Montgomery formulae with 22 multiplications, see
     Five, Six, and Seven-Term {K}aratsuba-Like Formulae,
     IEEE Transactions on Computers, volume 54, number 3, p. 362-369, 2005 */
    unsigned long ta[5], tb[5], pa[22], pb[22];
    ta[0] = a[0] ^ a[4];
    ta[1] = a[3] ^ a[5];
    ta[2] = a[2] ^ a[6];
    ta[3] = a[1] ^ ta[0];
    ta[4] = ta[1] ^ ta[2];
    pa[0] = a[6];
    pa[1] = a[5];
    pa[2] = a[5] ^ a[6];
    pa[3] = a[4];
    pa[4] = a[4] ^ a[6];
    pa[5] = a[3];
    pa[6] = ta[1];
    pa[7] = a[2];
    pa[8] = ta[2];
    pa[9] = a[1];
    pa[10]= a[1] ^ a[3];
    pa[11]= a[1] ^ a[2] ^ a[4] ^ a[5];
    pa[12]= a[1] ^ ta[4];
    pa[13]= a[0];
    pa[14]= ta[0];
    pa[15]= a[0] ^ a[2];
    pa[16]= a[0] ^ ta[4];
    pa[17]= a[3] ^ ta[0] ^ ta[2];
    pa[18]= a[0] ^ a[1];
    pa[19]= a[3] ^ a[6] ^ ta[3];
    pa[20]= ta[1] ^ ta[3];
    pa[21]= ta[3] ^ ta[4];

    tb[0] = b[0] ^ b[4];
    tb[1] = b[3] ^ b[5];
    tb[2] = b[2] ^ b[6];
    tb[3] = b[1] ^ tb[0];
    tb[4] = tb[1] ^ tb[2];
    pb[0] = b[6];
    pb[1] = b[5];
    pb[2] = b[5] ^ b[6];
    pb[3] = b[4];
    pb[4] = b[4] ^ b[6];
    pb[5] = b[3];
    pb[6] = tb[1];
    pb[7] = b[2];
    pb[8] = tb[2];
    pb[9] = b[1];
    pb[10]= b[1] ^ b[3];
    pb[11]= b[1] ^ b[2] ^ b[4] ^ b[5];
    pb[12]= b[1] ^ tb[4];
    pb[13]= b[0];
    pb[14]= tb[0];
    pb[15]= b[0] ^ b[2];
    pb[16]= b[0] ^ tb[4];
    pb[17]= b[3] ^ tb[0] ^ tb[2];
    pb[18]= b[0] ^ b[1];
    pb[19]= b[3] ^ b[6] ^ tb[3];
    pb[20]= tb[1] ^ tb[3];
    pb[21]= tb[3] ^ tb[4];

    __m128i p[22];

    p[0] = GF2X_FUNC(mul7cl_mul1)(pa[0], pb[0]);
    p[1] = GF2X_FUNC(mul7cl_mul1)(pa[1], pb[1]);
    p[2] = GF2X_FUNC(mul7cl_mul1)(pa[2], pb[2]);
    p[3] = GF2X_FUNC(mul7cl_mul1)(pa[3], pb[3]);
    p[4] = GF2X_FUNC(mul7cl_mul1)(pa[4], pb[4]);
    p[5] = GF2X_FUNC(mul7cl_mul1)(pa[5], pb[5]);
    p[6] = GF2X_FUNC(mul7cl_mul1)(pa[6], pb[6]);
    p[7] = GF2X_FUNC(mul7cl_mul1)(pa[7], pb[7]);
    p[8] = GF2X_FUNC(mul7cl_mul1)(pa[8], pb[8]);
    p[9] = GF2X_FUNC(mul7cl_mul1)(pa[9], pb[9]);
    p[10]= GF2X_FUNC(mul7cl_mul1)(pa[10], pb[10]);
    p[11]= GF2X_FUNC(mul7cl_mul1)(pa[11], pb[11]);
    p[12]= GF2X_FUNC(mul7cl_mul1)(pa[12], pb[12]);
    p[13]= GF2X_FUNC(mul7cl_mul1)(pa[13], pb[13]);
    p[14]= GF2X_FUNC(mul7cl_mul1)(pa[14], pb[14]);
    p[15]= GF2X_FUNC(mul7cl_mul1)(pa[15], pb[15]);
    p[16]= GF2X_FUNC(mul7cl_mul1)(pa[16], pb[16]);
    p[17]= GF2X_FUNC(mul7cl_mul1)(pa[17], pb[17]);
    p[18]= GF2X_FUNC(mul7cl_mul1)(pa[18], pb[18]);
    p[19]= GF2X_FUNC(mul7cl_mul1)(pa[19], pb[19]);
    p[20]= GF2X_FUNC(mul7cl_mul1)(pa[20], pb[20]);
    p[21]= GF2X_FUNC(mul7cl_mul1)(pa[21], pb[21]);

    __m128i t[13];

    t[0]  = PXOR(p[0], p[1]);
    t[1]  = PXOR(p[9], p[13]);
    t[2]  = PXOR(p[3], p[6]);
    t[3]  = PXOR(p[7], p[10]);
    t[4]  = PXOR(p[11], p[18]);
    t[5]  = PXOR(p[4], t[3]);
    t[6]  = PXOR(p[15], t[2]);
    t[7]  = PXOR(p[20], t[5]);
    t[8]  = PXOR(p[5], p[14]);
    t[9]  = PXOR(p[2], p[17]);
    t[10] = PXOR(p[5], p[8]);
    t[11] = PXOR(p[21], t[6]);
    t[12] = PXOR(p[16], t[4]);

    __m128i cc[13];

    cc[0]  = p[13];
    cc[2]  = PXOR3(p[7], p[15], t[1]);
    cc[4]  = PXOR4(p[3], t[1], t[3], t[8]);
    cc[6]  = PXOR5(p[0], p[12], p[13], t[7], t[11]);
    cc[8]  = PXOR4(p[7], t[0], t[2], t[10]);
    cc[10] = PXOR3(p[3], p[4], t[0]);
    cc[12] = p[0];

    cc[1]  = PXOR(p[18], t[1]);
    cc[3]  = PXOR3(t[7], t[9], t[12]);
    cc[5]  = PXOR6(p[2], p[11], p[19], t[1], t[10], t[11]);
    cc[7]  = PXOR5(p[21], t[0], t[5], t[8], t[12]);
    cc[9]  = PXOR5(p[12], p[19], t[4], t[6], t[9]);
    cc[11] = PXOR(p[2], t[0]);

    _mm_storeu_si128((__m128i*)(c),    PXOR(cc[0] ,                           _mm_slli_si128(cc[1], 8)));
    _mm_storeu_si128((__m128i*)(c+2),  PXOR(cc[2] , PXOR(_mm_srli_si128(cc[1], 8), _mm_slli_si128(cc[3], 8))));
    _mm_storeu_si128((__m128i*)(c+4),  PXOR(cc[4] , PXOR(_mm_srli_si128(cc[3], 8), _mm_slli_si128(cc[5], 8))));
    _mm_storeu_si128((__m128i*)(c+6),  PXOR(cc[6] , PXOR(_mm_srli_si128(cc[5], 8), _mm_slli_si128(cc[7], 8))));
    _mm_storeu_si128((__m128i*)(c+8),  PXOR(cc[8] , PXOR(_mm_srli_si128(cc[7], 8), _mm_slli_si128(cc[9], 8))));
    _mm_storeu_si128((__m128i*)(c+10), PXOR(cc[10], PXOR(_mm_srli_si128(cc[9], 8), _mm_slli_si128(cc[11], 8))));
    _mm_storeu_si128((__m128i*)(c+12), PXOR(cc[12],  _mm_srli_si128(cc[11], 8)));
#undef PXOR
#undef PXOR3
#undef PXOR4
#undef PXOR5
#undef PXOR6
#undef PZERO
}

#endif  /* GF2X_MUL7_H_ */
