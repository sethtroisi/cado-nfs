/* This file is part of the gf2x library.

   Copyright 2007, 2008, 2009, 2010, 2012, 2013, 2015
   Richard Brent, Pierrick Gaudry, Emmanuel Thome', Paul Zimmermann

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
GF2X_FUNC(mul5k3_mul1) (unsigned long a, unsigned long b)
{
    return _mm_clmulepi64_si128(_gf2x_mm_setr_epi64(a,0), _gf2x_mm_setr_epi64(b,0), 0);
}

#define PXOR(lop, rop) _mm_xor_si128((lop), (rop))

/* uses the variant of Karatsuba with 6 multiplications
   {d, 2} <- {a+2,1} * {b+2,1} */
GF2X_STORAGE_CLASS_mul3 void
GF2X_FUNC(mul5k3_mul3c) (unsigned long *c, const unsigned long *a, const unsigned long *b,
            unsigned long *d)
{
  unsigned long aa[3], bb[3];
  __m128i p0, p1, p2;
  __m128i pp0, pp1, pp2;

  aa[0] = a[1]^a[2];
  aa[1] = a[0]^a[2];
  aa[2] = a[0]^a[1];
  bb[0] = b[1]^b[2];
  bb[1] = b[0]^b[2];
  bb[2] = b[0]^b[1];
  p0  = GF2X_FUNC(mul5k3_mul1)(a[0], b[0]);
  p1  = GF2X_FUNC(mul5k3_mul1)(a[1], b[1]);
  p2  = GF2X_FUNC(mul5k3_mul1)(a[2], b[2]);
  pp0 = GF2X_FUNC(mul5k3_mul1)(aa[0], bb[0]);
  pp1 = GF2X_FUNC(mul5k3_mul1)(aa[1], bb[1]);
  pp2 = GF2X_FUNC(mul5k3_mul1)(aa[2], bb[2]);

  __m128i ce0 = p0;
  __m128i ce2 = PXOR(p0, PXOR(p1, PXOR(p2, pp1)));
  __m128i ce4 = p2;

  __m128i co1 = PXOR(p0, PXOR(p1, pp2));
  __m128i co3 = PXOR(pp0, PXOR(p1, p2));

  _mm_storeu_si128((__m128i*)(d),   p2);
  _mm_storeu_si128((__m128i*)(c),   PXOR(ce0, _mm_slli_si128(co1, 8)));
  _mm_storeu_si128((__m128i*)(c+2), PXOR(ce2, PXOR(_mm_srli_si128(co1, 8), _mm_slli_si128(co3, 8))));
  _mm_storeu_si128((__m128i*)(c+4), PXOR(ce4, _mm_srli_si128(co3, 8)));
}

/* uses the variant of Karatsuba with 6 multiplications,
   assumes {d, 2} = {a+2,1} * {b+2,1} */
GF2X_STORAGE_CLASS_mul3 void
GF2X_FUNC(mul5k3_mul3b) (unsigned long *c, const unsigned long *a, const unsigned long *b,
            const unsigned long *d)
{
  unsigned long aa[3], bb[3];
  __m128i p0, p1, p2;
  __m128i pp0, pp1, pp2;

  aa[0] = a[1]^a[2];
  aa[1] = a[0]^a[2];
  aa[2] = a[0]^a[1];
  bb[0] = b[1]^b[2];
  bb[1] = b[0]^b[2];
  bb[2] = b[0]^b[1];
  p0  = GF2X_FUNC(mul5k3_mul1)(a[0], b[0]);
  p1  = GF2X_FUNC(mul5k3_mul1)(a[1], b[1]);
  p2  = _mm_loadu_si128((__m128i *) d);
  pp0 = GF2X_FUNC(mul5k3_mul1)(aa[0], bb[0]);
  pp1 = GF2X_FUNC(mul5k3_mul1)(aa[1], bb[1]);
  pp2 = GF2X_FUNC(mul5k3_mul1)(aa[2], bb[2]);

  __m128i ce0 = p0;
  __m128i ce2 = PXOR(p0, PXOR(p1, PXOR(p2, pp1)));
  __m128i ce4 = p2;

  __m128i co1 = PXOR(p0, PXOR(p1, pp2));
  __m128i co3 = PXOR(pp0, PXOR(p1, p2));

  _mm_storeu_si128((__m128i*)(c),   PXOR(ce0, _mm_slli_si128(co1, 8)));
  _mm_storeu_si128((__m128i*)(c+2), PXOR(ce2, PXOR(_mm_srli_si128(co1, 8), _mm_slli_si128(co3, 8))));
  _mm_storeu_si128((__m128i*)(c+4), PXOR(ce4, _mm_srli_si128(co3, 8)));
}
#undef PXOR

/* based on mul5k_b.c, version with M(5)=M(2)+2M(3)-1=14 multiplications */
GF2X_STORAGE_CLASS_mul5
void gf2x_mul5 (unsigned long *c, const unsigned long *a, const unsigned long *b)
{
  unsigned long aa[3], bb[3], ab[6], ab3, ab4, ab5, d[2];

    gf2x_mul2 (c+6, a+3, b+3);
    GF2X_FUNC(mul5k3_mul3c) (c, a, b, d);
    aa[0] = a[0] ^ a[3];
    aa[1] = a[1] ^ a[4];
    aa[2] = a[2];
    bb[0] = b[0] ^ b[3];
    bb[1] = b[1] ^ b[4];
    bb[2] = b[2];
    GF2X_FUNC(mul5k3_mul3b) (ab, aa, bb, d);
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
