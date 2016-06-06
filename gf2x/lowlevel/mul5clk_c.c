/* This file is part of the gf2x library.

   Copyright 2010, 2013, 2015
   Richard Brent, Pierrick Gaudry, Emmanuel Thome', Paul Zimmermann,
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
#ifndef GF2X_MUL5_H_
#define GF2X_MUL5_H_

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
GF2X_FUNC(mul5clk_c_mul1) (unsigned long a, unsigned long b)
{
    return _mm_clmulepi64_si128(_gf2x_mm_setr_epi64(a,0), _gf2x_mm_setr_epi64(b,0), 0);
}


GF2X_STORAGE_CLASS_mul5
void gf2x_mul5 (unsigned long *c, const unsigned long *a,
        const unsigned long *b)
{
#define PXOR(lop, rop) _mm_xor_si128((lop), (rop))
  /* Montgomery formulae with 13 multiplications, see
     Five, Six, and Seven-Term {K}aratsuba-Like Formulae,
     IEEE Transactions on Computers, volume 54, number 3, p. 362-369, 2005 */
  unsigned long ta[3], tb[3];

  ta[0] = a[0]  ^ a[4]         ; tb[0] = b[0]  ^ b[4];
  ta[1] = a[1]  ^ a[2]         ; tb[1] = b[1]  ^ b[2];
  ta[2] = a[3]  ^ ta[0]        ; tb[2] = b[3]  ^ tb[0];

  unsigned long pa[8], pb[8];
  pa[0] = ta[1] ^ ta[2]        ; pb[0] = tb[1] ^ tb[2];
  pa[1] = a[2]  ^ ta[2]        ; pb[1] = b[2]  ^ tb[2];
  pa[2] = ta[0] ^ ta[1]        ; pb[2] = tb[0] ^ tb[1];
  pa[3] = a[1]  ^ ta[2]        ; pb[3] = b[1]  ^ tb[2];
  pa[4] = a[0]  ^ a[2]  ^ a[3] ; pb[4] = b[0]  ^ b[2]  ^ b[3];
  pa[5] = a[4]  ^ ta[1]        ; pb[5] = b[4]  ^ tb[1];
  pa[6] = a[3]  ^ a[4]         ; pb[6] = b[3]  ^ b[4];
  pa[7] = a[0]  ^ a[1]         ; pb[7] = b[0]  ^ b[1];

  __m128i p0  = GF2X_FUNC(mul5clk_c_mul1)(pa[0], pb[0]);
  __m128i p2  = GF2X_FUNC(mul5clk_c_mul1)(pa[1], pb[1]);
  __m128i p4  = GF2X_FUNC(mul5clk_c_mul1)(pa[2], pb[2]);
  __m128i p6  = GF2X_FUNC(mul5clk_c_mul1)(pa[3], pb[3]);
  __m128i p8  = GF2X_FUNC(mul5clk_c_mul1)(pa[4], pb[4]);
  __m128i p10 = GF2X_FUNC(mul5clk_c_mul1)(pa[5], pb[5]);
  __m128i p12 = GF2X_FUNC(mul5clk_c_mul1)(pa[6], pb[6]);
  __m128i p14 = GF2X_FUNC(mul5clk_c_mul1)(pa[7], pb[7]);
  __m128i p16 = GF2X_FUNC(mul5clk_c_mul1)(ta[0], tb[0]);
  __m128i p18 = GF2X_FUNC(mul5clk_c_mul1)(a[4],  b[4]);
  __m128i p20 = GF2X_FUNC(mul5clk_c_mul1)(a[3],  b[3]);
  __m128i p22 = GF2X_FUNC(mul5clk_c_mul1)(a[1],  b[1]);
  __m128i p24 = GF2X_FUNC(mul5clk_c_mul1)(a[0],  b[0]);
  __m128i t0  = PXOR(p14, p24);
  __m128i t2  = PXOR(p12, p18);
  __m128i t4  = PXOR(p2, p16);
  __m128i t6  = PXOR(p0, p6);
  __m128i t8  = PXOR(p4, p16);
  __m128i t10 = PXOR(p10, t0);
  __m128i t12 = PXOR(p8, t2);

  __m128i ce0 = p24;
  __m128i ce2 = PXOR(p18, PXOR(t8, t10));
  __m128i ce4 = PXOR(p0, PXOR(p20, PXOR(p22, PXOR(t10, t12))));
  __m128i ce6 = PXOR(p24, PXOR(t4, t12));
  __m128i ce8 = p18;

  __m128i co1 = PXOR(p22, t0);
  __m128i co3 = PXOR(t2, PXOR(t4, t6));
  __m128i co5 = PXOR(t0, PXOR(t6, t8));
  __m128i co7 = PXOR(p20, t2);

  _mm_storeu_si128((__m128i*)(c),   PXOR(ce0, _mm_slli_si128(co1, 8)));
  _mm_storeu_si128((__m128i*)(c+2), PXOR(ce2, PXOR(_mm_srli_si128(co1, 8), _mm_slli_si128(co3, 8))));
  _mm_storeu_si128((__m128i*)(c+4), PXOR(ce4, PXOR(_mm_srli_si128(co3, 8), _mm_slli_si128(co5, 8))));
  _mm_storeu_si128((__m128i*)(c+6), PXOR(ce6, PXOR(_mm_srli_si128(co5, 8), _mm_slli_si128(co7, 8))));
  _mm_storeu_si128((__m128i*)(c+8), PXOR(ce8, _mm_srli_si128(co7, 8)));
#undef PXOR
}


#endif  /* GF2X_MUL5_H_ */
