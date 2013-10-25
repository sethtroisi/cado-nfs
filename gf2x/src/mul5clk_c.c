/* This file is part of the gf2x library.

   Copyright 2007, 2008, 2009, 2010
   Richard Brent, Pierrick Gaudry, Emmanuel Thome', Paul Zimmermann,
   Nicolas Estibals (for this file)

   This program is free software; you can redistribute it and/or modify it
   under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation; either version 2.1 of the License, or (at
   your option) any later version.
   
   This program is distributed in the hope that it will be useful, but WITHOUT
   ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
   FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
   License for more details.
   
   You should have received a copy of the GNU Lesser General Public
   License along with CADO-NFS; see the file COPYING.  If not, write to
   the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
   Boston, MA 02110-1301, USA.
*/
#ifndef GF2X_MUL5_H_
#define GF2X_MUL5_H_

#include "gf2x.h"
/* All gf2x source files for lowlevel functions must include gf2x-small.h
 * This is mandatory for the tuning mechanism. */
#include "gf2x/gf2x-small.h"

#include <emmintrin.h>
#include <wmmintrin.h>

#if GF2X_WORDSIZE != 64
#error "This code is for 64-bit only"
#endif

#include "gf2x/gf2x-config.h"

#ifndef HAVE_PCLMUL_SUPPORT
#error "This code needs pclmul support"
#endif

/* TODO: if somebody comes up with a neat way to improve the interface so
 * as to remove the false dependency on pclmul, that would be nice.
 */
static inline __v2di
GF2X_FUNC(mul5clk_c_mul1) (unsigned long a, unsigned long b)
{   
    __v2di aa = (__v2di) { (long long)a, 0 };
    __v2di bb = (__v2di) { (long long)b, 0 };
    return _mm_clmulepi64_si128(aa, bb, 0);
}   
GF2X_STORAGE_CLASS_mul5
void gf2x_mul5 (unsigned long *c, const unsigned long *a,
        const unsigned long *b)
{
  /* Montgomery formulae with 13 multiplications, see
     Five, Six, and Seven-Term {K}aratsuba-Like Formulae,
     IEEE Transactions on Computers, volume 54, number 3, p. 362-369, 2005 */
  unsigned long ta[3], tb[3], pa[8], pb[8];
  __v2di p0, p2, p4, p6, p8, p10, p12, p14, p16, p18, p20, p22, p24;
  __v2di t0, t2, t4, t6, t8, t10, t12;
  ta[0] = a[0]  ^ a[4]         ; tb[0] = b[0]  ^ b[4];
  ta[1] = a[1]  ^ a[2]         ; tb[1] = b[1]  ^ b[2];
  ta[2] = a[3]  ^ ta[0]        ; tb[2] = b[3]  ^ tb[0];
  pa[0] = ta[1] ^ ta[2]        ; pb[0] = tb[1] ^ tb[2];
  pa[1] = a[2]  ^ ta[2]        ; pb[1] = b[2]  ^ tb[2];
  pa[2] = ta[0] ^ ta[1]        ; pb[2] = tb[0] ^ tb[1];
  pa[3] = a[1]  ^ ta[2]        ; pb[3] = b[1]  ^ tb[2];
  pa[4] = a[0]  ^ a[2]  ^ a[3] ; pb[4] = b[0]  ^ b[2]  ^ b[3];
  pa[5] = a[4]  ^ ta[1]        ; pb[5] = b[4]  ^ tb[1];
  pa[6] = a[3]  ^ a[4]         ; pb[6] = b[3]  ^ b[4];
  pa[7] = a[0]  ^ a[1]         ; pb[7] = b[0]  ^ b[1];
  p0  = GF2X_FUNC(mul5clk_c_mul1)(pa[0], pb[0]);
  p2  = GF2X_FUNC(mul5clk_c_mul1)(pa[1], pb[1]);
  p4  = GF2X_FUNC(mul5clk_c_mul1)(pa[2], pb[2]);
  p6  = GF2X_FUNC(mul5clk_c_mul1)(pa[3], pb[3]);
  p8  = GF2X_FUNC(mul5clk_c_mul1)(pa[4], pb[4]);
  p10 = GF2X_FUNC(mul5clk_c_mul1)(pa[5], pb[5]);
  p12 = GF2X_FUNC(mul5clk_c_mul1)(pa[6], pb[6]);
  p14 = GF2X_FUNC(mul5clk_c_mul1)(pa[7], pb[7]);
  p16 = GF2X_FUNC(mul5clk_c_mul1)(ta[0], tb[0]);
  p18 = GF2X_FUNC(mul5clk_c_mul1)(a[4],  b[4]);
  p20 = GF2X_FUNC(mul5clk_c_mul1)(a[3],  b[3]);
  p22 = GF2X_FUNC(mul5clk_c_mul1)(a[1],  b[1]);
  p24 = GF2X_FUNC(mul5clk_c_mul1)(a[0],  b[0]);
  t0  = p14 ^ p24;
  t2  = p12 ^ p18;
  t4  = p2  ^ p16;
  t6  = p0  ^ p6;
  t8  = p4  ^ p16;
  t10 = p10 ^ t0;
  t12 = p8  ^ t2;

  __v2di ce0 = p24;
  __v2di ce2 = p18 ^ t8  ^ t10;
  __v2di ce4 = p0  ^ p20 ^ p22 ^ t10 ^ t12;
  __v2di ce6 = p24 ^ t4  ^ t12;
  __v2di ce8 = p18;

  __v2di co1 = p22 ^ t0;
  __v2di co3 = t2  ^ t4  ^ t6;
  __v2di co5 = t0  ^ t6  ^ t8;
  __v2di co7 = p20 ^ t2;

  _mm_storeu_si128((__v2di*)(c),   ce0 ^ _mm_slli_si128(co1, 8));
  _mm_storeu_si128((__v2di*)(c+2), ce2 ^ _mm_srli_si128(co1, 8) ^ _mm_slli_si128(co3, 8));
  _mm_storeu_si128((__v2di*)(c+4), ce4 ^ _mm_srli_si128(co3, 8) ^ _mm_slli_si128(co5, 8));
  _mm_storeu_si128((__v2di*)(c+6), ce6 ^ _mm_srli_si128(co5, 8) ^ _mm_slli_si128(co7, 8));
  _mm_storeu_si128((__v2di*)(c+8), ce8 ^ _mm_srli_si128(co7, 8));
}


#endif  /* GF2X_MUL5_H_ */
