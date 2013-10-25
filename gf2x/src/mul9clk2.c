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
GF2X_FUNC(mul9clk_mul1) (unsigned long a, unsigned long b)
{
    __v2di aa = (__v2di) { (long long)a, 0 };
    __v2di bb = (__v2di) { (long long)b, 0 };
    return _mm_clmulepi64_si128(aa, bb, 0);
}

/* uses the variant of Karatsuba with 6 multiplications */
static void GF2X_FUNC(mul9clk_mul3) (__v2di *ce, __v2di *co, const unsigned long *a, const unsigned long *b)
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
  p0  = GF2X_FUNC(mul9clk_mul1)(a[0], b[0]);
  p1  = GF2X_FUNC(mul9clk_mul1)(a[1], b[1]);
  p2  = GF2X_FUNC(mul9clk_mul1)(a[2], b[2]);
  pp0 = GF2X_FUNC(mul9clk_mul1)(aa[0], bb[0]);
  pp1 = GF2X_FUNC(mul9clk_mul1)(aa[1], bb[1]);
  pp2 = GF2X_FUNC(mul9clk_mul1)(aa[2], bb[2]);

  ce[0] = p0;
  ce[1] = p0^p1^p2^pp1;
  ce[2] = p2;

  co[0] = p0^p1^pp2;
  co[1] = pp0^p1^p2;
}

/* recursive application of the Karatsuba-3 algorithm with 6 multiplies,
   i.e., with a total of 6*6=36 multiplies */
GF2X_STORAGE_CLASS_mul9
void gf2x_mul9 (unsigned long *c, const unsigned long *a, const unsigned long *b)
{
  unsigned long aa[9], bb[9];
  aa[0] = a[3]^a[6];
  aa[1] = a[4]^a[7];

  aa[2] = a[5]^a[8];

  aa[3] = a[0]^a[6];
  aa[4] = a[1]^a[7];

  aa[5] = a[2]^a[8];

  aa[6] = a[0]^a[3];
  aa[7] = a[1]^a[4];

  aa[8] = a[2]^a[5];

  bb[0] = b[3]^b[6];
  bb[1] = b[4]^b[7];
  bb[2] = b[5]^b[8];

  bb[3] = b[0]^b[6];
  bb[4] = b[1]^b[7];

  bb[5] = b[2]^b[8];
  bb[6] = b[0]^b[3];
  bb[7] = b[1]^b[4];
  bb[8] = b[2]^b[5];

  __v2di p0e[3], p0o[2];
  __v2di p1e[3], p1o[2];
  __v2di p2e[3], p2o[2];
  __v2di q0e[3], q0o[2];
  __v2di q1e[3], q1o[2];
  __v2di q2e[3], q2o[2];
  GF2X_FUNC(mul9clk_mul3) (p0e, p0o, a+0, b+0);
  GF2X_FUNC(mul9clk_mul3) (p1e, p1o, a+3, b+3);
  GF2X_FUNC(mul9clk_mul3) (p2e, p2o, a+6, b+6);
  GF2X_FUNC(mul9clk_mul3) (q0e, q0o, aa+0, bb+0);
  GF2X_FUNC(mul9clk_mul3) (q1e, q1o, aa+3, bb+3);
  GF2X_FUNC(mul9clk_mul3) (q2e, q2o, aa+6, bb+6);
 
  __v2di e,h,l;
  e = p0e[0];
  l = p0o[0];
  _mm_storeu_si128((__v2di*)(c), e ^ _mm_slli_si128(l, 8));

  e = p0e[1];                      h=l; l=p0o[1]^p1e[0]^q2e[0]^p0e[0];
  _mm_storeu_si128((__v2di*)(c+2), e ^ _mm_slli_si128(l, 8) ^ _mm_srli_si128(h, 8));

  e = p0e[2]^p0o[0]^p1o[0]^q2o[0]; h=l; l=p1e[1]^q2e[1]^p0e[1];
  _mm_storeu_si128((__v2di*)(c+4), e ^ _mm_slli_si128(l, 8) ^ _mm_srli_si128(h, 8));

  e = p0e[0]^q1e[0]^p2e[0]^p1e[0]^p0o[1]^p1o[1]^q2o[1]; h=l; l=p0o[0]^q1o[0]^p2o[0]^p1o[0]^p1e[2]^q2e[2]^p0e[2];
  _mm_storeu_si128((__v2di*)(c+6), e ^ _mm_slli_si128(l, 8) ^ _mm_srli_si128(h, 8));

  e = p0e[1]^q1e[1]^p2e[1]^p1e[1]; h=l; l=p0o[1]^q1o[1]^p2o[1]^p1o[1]^q0e[0]^p1e[0]^p2e[0];
  _mm_storeu_si128((__v2di*)(c+8), e ^ _mm_slli_si128(l, 8) ^ _mm_srli_si128(h, 8));

  e = p0e[2]^q1e[2]^p2e[2]^p1e[2]^q0o[0]^p1o[0]^p2o[0]; h=l; l=q0e[1]^p1e[1]^p2e[1];
  _mm_storeu_si128((__v2di*)(c+10), e ^ _mm_slli_si128(l, 8) ^ _mm_srli_si128(h, 8));

  e = p2e[0]^q0o[1]^p1o[1]^p2o[1]; h=l; l=p2o[0]^q0e[2]^p1e[2]^p2e[2];
  _mm_storeu_si128((__v2di*)(c+12), e ^ _mm_slli_si128(l, 8) ^ _mm_srli_si128(h, 8));

  e = p2e[1]; h=l; l=p2o[1];
  _mm_storeu_si128((__v2di*)(c+14), e ^ _mm_slli_si128(l, 8) ^ _mm_srli_si128(h, 8));

  e = p2e[2];
  _mm_storeu_si128((__v2di*)(c+16), e ^ _mm_srli_si128(h, 8));

  /*
  c[0]  =                   p0e[0]                   ^                                           ;
  c[1]  =                   p0e[1]                   ^                   p0o[0]                  ;
  c[2]  =                   p0e[2]                   ^                   p0o[1]                  ;
  c[3]  =                   p0e[3]                   ^                   p0o[2]p1e[0]q2e[0]p0e[0];
  c[4]  =                   p0e[4]p0o[0]p1o[0]q2o[0] ^                   p0o[3]p1e[1]q2e[1]p0e[1];
  c[5]  =                   p0e[5]p0o[1]p1o[1]q2o[1] ^                         p1e[2]q2e[2]p0e[2];
  c[6]  = p0e[0]q1e[0]p2e[0]p1e[0]p0o[2]p1o[2]q2o[2] ^                         p1e[3]q2e[3]p0e[3];
  c[7]  = p0e[1]q1e[1]p2e[1]p1e[1]p0o[3]p1o[3]q2o[3] ^ p0o[0]q1o[0]p2o[0]p1o[0]p1e[4]q2e[4]p0e[4];
  c[8]  = p0e[2]q1e[2]p2e[2]p1e[2]                   ^ p0o[1]q1o[1]p2o[1]p1o[1]p1e[5]q2e[5]p0e[5];
  c[9]  = p0e[3]q1e[3]p2e[3]p1e[3]                   ^ p0o[2]q1o[2]p2o[2]p1o[2]q0e[0]p1e[0]p2e[0];
  c[10] = p0e[4]q1e[4]p2e[4]p1e[4]q0o[0]p1o[0]p2o[0] ^ p0o[3]q1o[3]p2o[3]p1o[3]q0e[1]p1e[1]p2e[1];
  c[11] = p0e[5]q1e[5]p2e[5]p1e[5]q0o[1]p1o[1]p2o[1] ^                         q0e[2]p1e[2]p2e[2];
  c[12] = p2e[0]                  q0o[2]p1o[2]p2o[2] ^                         q0e[3]p1e[3]p2e[3];
  c[13] = p2e[1]                  q0o[3]p1o[3]p2o[3] ^                   p2o[0]q0e[4]p1e[4]p2e[4];
  c[14] = p2e[2]                                     ^                   p2o[1]q0e[5]p1e[5]p2e[5];
  c[15] = p2e[3]                                     ^                   p2o[2]                  ;
  c[16] = p2e[4]                                     ^                   p2o[3]                  ;
  c[17] = p2e[5]                                     ^                                           ;
  */
}                          

#endif  /* GF2X_MUL9_H_ */
