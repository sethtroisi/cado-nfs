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
#ifndef GF2X_MUL5_H_
#define GF2X_MUL5_H_

#include "gf2x.h"
/* All gf2x source files for lowlevel functions must include gf2x-small.h
 * This is mandatory for the tuning mechanism. */
#include "gf2x/gf2x-small.h"

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
GF2X_FUNC(mul3clk_c_gf2x_mul1) (unsigned long a, unsigned long b)
{
    __v2di aa = (__v2di) { (long long)a, 0 };
    __v2di bb = (__v2di) { (long long)b, 0 };
    return _mm_clmulepi64_si128(aa, bb, 0);
}

/* uses the variant of Karatsuba with 6 multiplications
   {d, 2} <- {a+2,1} * {b+2,1} */
GF2X_STORAGE_CLASS_mul3 void
gf2x_mul3c (unsigned long *c, const unsigned long *a, const unsigned long *b,
            unsigned long *d)
{
  unsigned long aa[3], bb[3];
  __v2di p0, p1, p2;
  __v2di pp0, pp1, pp2;

  aa[0] = a[1]^a[2];
  aa[1] = a[0]^a[2];
  aa[2] = a[0]^a[1];
  bb[0] = b[1]^b[2];
  bb[1] = b[0]^b[2];
  bb[2] = b[0]^b[1];
  p0  = GF2X_FUNC(mul3clk_c_gf2x_mul1)(a[0], b[0]);
  p1  = GF2X_FUNC(mul3clk_c_gf2x_mul1)(a[1], b[1]);
  p2  = GF2X_FUNC(mul3clk_c_gf2x_mul1)(a[2], b[2]);
  pp0 = GF2X_FUNC(mul3clk_c_gf2x_mul1)(aa[0], bb[0]);
  pp1 = GF2X_FUNC(mul3clk_c_gf2x_mul1)(aa[1], bb[1]);
  pp2 = GF2X_FUNC(mul3clk_c_gf2x_mul1)(aa[2], bb[2]);

  __v2di ce0, ce2, ce4;
  __v2di co1, co3;

  ce0 = p0;
  ce2 = p0^p1^p2^pp1;
  ce4 = p2;

  co1 = p0^p1^pp2;
  co3 = pp0^p1^p2;

  _mm_storeu_si128((__v2di*)(d),   p2);
  _mm_storeu_si128((__v2di*)(c),   ce0 ^ _mm_slli_si128(co1, 8));
  _mm_storeu_si128((__v2di*)(c+2), ce2 ^ _mm_srli_si128(co1, 8) ^ _mm_slli_si128(co3, 8));
  _mm_storeu_si128((__v2di*)(c+4), ce4 ^ _mm_srli_si128(co3, 8));
}

/* uses the variant of Karatsuba with 6 multiplications,
   assumes {d, 2} = {a+2,1} * {b+2,1} */
GF2X_STORAGE_CLASS_mul3 void
gf2x_mul3b (unsigned long *c, const unsigned long *a, const unsigned long *b,
            const unsigned long *d)
{
  unsigned long aa[3], bb[3];
  __v2di p0, p1, p2;
  __v2di pp0, pp1, pp2;

  aa[0] = a[1]^a[2];
  aa[1] = a[0]^a[2];
  aa[2] = a[0]^a[1];
  bb[0] = b[1]^b[2];
  bb[1] = b[0]^b[2];
  bb[2] = b[0]^b[1];
  p0  = GF2X_FUNC(mul3clk_c_gf2x_mul1)(a[0], b[0]);
  p1  = GF2X_FUNC(mul3clk_c_gf2x_mul1)(a[1], b[1]);
  p2  = _mm_loadu_si128((__v2di *) d);
  pp0 = GF2X_FUNC(mul3clk_c_gf2x_mul1)(aa[0], bb[0]);
  pp1 = GF2X_FUNC(mul3clk_c_gf2x_mul1)(aa[1], bb[1]);
  pp2 = GF2X_FUNC(mul3clk_c_gf2x_mul1)(aa[2], bb[2]);

  __v2di ce0, ce2, ce4;
  __v2di co1, co3;

  ce0 = p0;
  ce2 = p0^p1^p2^pp1;
  ce4 = p2;

  co1 = p0^p1^pp2;
  co3 = pp0^p1^p2;

  _mm_storeu_si128((__v2di*)(c),   ce0 ^ _mm_slli_si128(co1, 8));
  _mm_storeu_si128((__v2di*)(c+2), ce2 ^ _mm_srli_si128(co1, 8) ^ _mm_slli_si128(co3, 8));
  _mm_storeu_si128((__v2di*)(c+4), ce4 ^ _mm_srli_si128(co3, 8));
}

/* based on mul5k_b.c, version with M(5)=M(2)+2M(3)-1=14 multiplications */
GF2X_STORAGE_CLASS_mul5
void gf2x_mul5 (unsigned long *c, const unsigned long *a, const unsigned long *b)
{
  unsigned long aa[3], bb[3], ab[6], ab3, ab4, ab5, d[2];

    gf2x_mul2 (c+6, a+3, b+3);
    gf2x_mul3c (c, a, b, d);
    aa[0] = a[0] ^ a[3];
    aa[1] = a[1] ^ a[4];
    aa[2] = a[2];
    bb[0] = b[0] ^ b[3];
    bb[1] = b[1] ^ b[4];
    bb[2] = b[2];
    gf2x_mul3b (ab, aa, bb, d);
    ab3 = ab[3] ^ c[3];
    ab4 = ab[4] ^ c[4];
    ab5 = ab[5] ^ c[5];
    c[3] ^= ab[0] ^ c[0] ^ c[6];
    c[4] ^= ab[1] ^ c[1] ^ c[7];
    c[5] ^= ab[2] ^ c[2] ^ c[8];
    c[6] ^= ab3 ^ c[9];
    c[7] ^= ab4;
    c[8] ^= ab5;
}

#endif  /* GF2X_MUL5_H_ */
