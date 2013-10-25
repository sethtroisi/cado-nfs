/* This file is part of the gf2x library.

   Copyright 2007, 2008, 2009, 2010
   Richard Brent, Pierrick Gaudry, Emmanuel Thome', Paul Zimmermann

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
#ifndef GF2X_MUL9_H_
#define GF2X_MUL9_H_

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
GF2X_FUNC(mul9k3_mul1) (unsigned long a, unsigned long b)
{   
    __v2di aa = (__v2di) { (long long)a, 0 };
    __v2di bb = (__v2di) { (long long)b, 0 };
    return _mm_clmulepi64_si128(aa, bb, 0);
}   

/* same as mul5, but stores a[4]*b[4] into {d,2} */
GF2X_STORAGE_CLASS_mul5
void mul9k3_gf2x_mul5c (unsigned long *c, const unsigned long *a,
                        const unsigned long *b, unsigned long *d)
{
  /* Montgomery formulae with 13 multiplications */
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
  p0  = GF2X_FUNC(mul9k3_mul1)(pa[0], pb[0]);
  p2  = GF2X_FUNC(mul9k3_mul1)(pa[1], pb[1]);
  p4  = GF2X_FUNC(mul9k3_mul1)(pa[2], pb[2]);
  p6  = GF2X_FUNC(mul9k3_mul1)(pa[3], pb[3]);
  p8  = GF2X_FUNC(mul9k3_mul1)(pa[4], pb[4]);
  p10 = GF2X_FUNC(mul9k3_mul1)(pa[5], pb[5]);
  p12 = GF2X_FUNC(mul9k3_mul1)(pa[6], pb[6]);
  p14 = GF2X_FUNC(mul9k3_mul1)(pa[7], pb[7]);
  p16 = GF2X_FUNC(mul9k3_mul1)(ta[0], tb[0]);
  p18 = GF2X_FUNC(mul9k3_mul1)(a[4],  b[4]);
  p20 = GF2X_FUNC(mul9k3_mul1)(a[3],  b[3]);
  p22 = GF2X_FUNC(mul9k3_mul1)(a[1],  b[1]);
  p24 = GF2X_FUNC(mul9k3_mul1)(a[0],  b[0]);
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

  _mm_storeu_si128((__v2di*)(d),   p18); /* a[4] * b[4] */
  _mm_storeu_si128((__v2di*)(c),   ce0 ^ _mm_slli_si128(co1, 8));
  _mm_storeu_si128((__v2di*)(c+2), ce2 ^ _mm_srli_si128(co1, 8) ^ _mm_slli_si128(co3, 8));
  _mm_storeu_si128((__v2di*)(c+4), ce4 ^ _mm_srli_si128(co3, 8) ^ _mm_slli_si128(co5, 8));
  _mm_storeu_si128((__v2di*)(c+6), ce6 ^ _mm_srli_si128(co5, 8) ^ _mm_slli_si128(co7, 8));
  _mm_storeu_si128((__v2di*)(c+8), ce8 ^ _mm_srli_si128(co7, 8));
}

/* same as mul5, but assumes {d,2} contains a[4]*b[4] */
GF2X_STORAGE_CLASS_mul5
void mul9k3_gf2x_mul5b (unsigned long *c, const unsigned long *a,
                        const unsigned long *b, const unsigned long *d)
{
  /* Montgomery formulae with 13 multiplications */
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
  p0  = GF2X_FUNC(mul9k3_mul1)(pa[0], pb[0]);
  p2  = GF2X_FUNC(mul9k3_mul1)(pa[1], pb[1]);
  p4  = GF2X_FUNC(mul9k3_mul1)(pa[2], pb[2]);
  p6  = GF2X_FUNC(mul9k3_mul1)(pa[3], pb[3]);
  p8  = GF2X_FUNC(mul9k3_mul1)(pa[4], pb[4]);
  p10 = GF2X_FUNC(mul9k3_mul1)(pa[5], pb[5]);
  p12 = GF2X_FUNC(mul9k3_mul1)(pa[6], pb[6]);
  p14 = GF2X_FUNC(mul9k3_mul1)(pa[7], pb[7]);
  p16 = GF2X_FUNC(mul9k3_mul1)(ta[0], tb[0]);
  /* p18 = GF2X_FUNC(mul9k3_mul1)(a[4],  b[4]); */
  p18 = _mm_loadu_si128((__v2di *) d);
  p20 = GF2X_FUNC(mul9k3_mul1)(a[3],  b[3]);
  p22 = GF2X_FUNC(mul9k3_mul1)(a[1],  b[1]);
  p24 = GF2X_FUNC(mul9k3_mul1)(a[0],  b[0]);
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

/* (based on mul9k.c)
   1 call to mul4 and 2 calls to mul5 minus one multiplication, i.e.,
   34 multiplications with mul5clk_c */
GF2X_STORAGE_CLASS_mul9
void gf2x_mul9 (unsigned long *c, const unsigned long *a, const unsigned long *b)
{
  unsigned long aa[5], bb[5], ab[10], ab5, ab6, ab7, ab8, ab9, d[2];

    gf2x_mul4 (c+10, a+5, b+5);
    mul9k3_gf2x_mul5c (c, a, b, d); /* a[4]*b[4] is cached in d */
    aa[0] = a[0] ^ a[5];
    aa[1] = a[1] ^ a[6];
    aa[2] = a[2] ^ a[7];
    aa[3] = a[3] ^ a[8];
    aa[4] = a[4];
    bb[0] = b[0] ^ b[5];
    bb[1] = b[1] ^ b[6];
    bb[2] = b[2] ^ b[7];
    bb[3] = b[3] ^ b[8];
    bb[4] = b[4];
    mul9k3_gf2x_mul5b (ab, aa, bb, d);
    ab5 = ab[5] ^ c[5];
    ab6 = ab[6] ^ c[6];
    ab7 = ab[7] ^ c[7];
    ab8 = ab[8] ^ c[8];
    ab9 = ab[9] ^ c[9];
    c[5] ^= ab[0] ^ c[0] ^ c[10];
    c[6] ^= ab[1] ^ c[1] ^ c[11];
    c[7] ^= ab[2] ^ c[2] ^ c[12];
    c[8] ^= ab[3] ^ c[3] ^ c[13];
    c[9] ^= ab[4] ^ c[4] ^ c[14];
    c[10] ^= ab5 ^ c[15];
    c[11] ^= ab6 ^ c[16];
    c[12] ^= ab7 ^ c[17];
    c[13] ^= ab8;
    c[14] ^= ab9;
}

#endif  /* GF2X_MUL9_H_ */
