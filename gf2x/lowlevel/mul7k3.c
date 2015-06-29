/* This file is part of the gf2x library.

   Copyright 2010, 2012, 2013, 2015
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

#define PXOR(lop, rop) _mm_xor_si128((lop), (rop))
#define PXOR3(op1, op2, op3) PXOR(op1, PXOR(op2, op3))
#define PZERO    _mm_setzero_si128()

/* TODO: if somebody comes up with a neat way to improve the interface so
 * as to remove the false dependency on pclmul, that would be nice.
 */
static inline void GF2X_FUNC(mul7k3_mul2)(__m128i *t, __m128i ss1,
	       __m128i ss2)
{
    __m128i t00 = _mm_clmulepi64_si128(ss1, ss2, 0);
    __m128i t11 = _mm_clmulepi64_si128(ss1, ss2, 0x11);
    ss1 = PXOR(ss1, _mm_shuffle_epi32(ss1, _MM_SHUFFLE(1, 0, 3, 2)));
    ss2 = PXOR(ss2, _mm_shuffle_epi32(ss2, _MM_SHUFFLE(1, 0, 3, 2)));
    __m128i tk = PXOR(PXOR(t00, t11), _mm_clmulepi64_si128(ss1, ss2, 0));
    t[0] = PXOR(t00, _mm_unpacklo_epi64(PZERO, tk));
    t[1] = PXOR(t11, _mm_unpackhi_epi64(tk, PZERO));
}
static inline void
GF2X_FUNC(mul7k3_mul2b)(__m128i *t, __m128i ss1, __m128i ss2, const unsigned long * sc)
{
    __m128i t00 = _mm_clmulepi64_si128(ss1, ss2, 0);
    ss1 = PXOR(ss1, _mm_shuffle_epi32(ss1, _MM_SHUFFLE(1, 0, 3, 2)));
    ss2 = PXOR(ss2, _mm_shuffle_epi32(ss2, _MM_SHUFFLE(1, 0, 3, 2)));
    __m128i c = _mm_loadu_si128((__m128i*)sc);
    __m128i tk = PXOR(PXOR(t00, c), _mm_clmulepi64_si128(ss1, ss2, 0));
    t[0] = PXOR(t00, _mm_unpacklo_epi64(PZERO, tk));
    t[1] = PXOR(c, _mm_unpackhi_epi64(tk, PZERO));
}
static inline void
GF2X_FUNC(mul7k3_mul2c)(__m128i *t, __m128i ss1, __m128i ss2, unsigned long *sc)
{
    __m128i t00 = _mm_clmulepi64_si128(ss1, ss2, 0);
    __m128i t11 = _mm_clmulepi64_si128(ss1, ss2, 0x11);
    ss1 = PXOR(ss1, _mm_shuffle_epi32(ss1, _MM_SHUFFLE(1, 0, 3, 2)));
    ss2 = PXOR(ss2, _mm_shuffle_epi32(ss2, _MM_SHUFFLE(1, 0, 3, 2)));
    __m128i tk = PXOR(PXOR(t00, t11), _mm_clmulepi64_si128(ss1, ss2, 0));
    _mm_storeu_si128((__m128i*)sc, t11);
    t[0] = PXOR(t00, _mm_unpacklo_epi64(PZERO, tk));
    t[1] = PXOR(t11, _mm_unpackhi_epi64(tk, PZERO));
}

/* specialized Karatsuba with 3 calls to mul2, i.e., 9 multiplications
   {d,2} <- {a+3,1} * {b+3,1} */
GF2X_STORAGE_CLASS_mul4
void GF2X_FUNC(mul7k3_mul4c) (unsigned long *c, const unsigned long *a, const unsigned long *b, unsigned long *d)
{
  __m128i ab[2];
  __m128i lo[2], hi[2];
  __m128i a0 = _mm_loadu_si128((__m128i*)a);
  __m128i a2 = _mm_loadu_si128((__m128i*)(a+2));
  __m128i b0 = _mm_loadu_si128((__m128i*)b);
  __m128i b2 = _mm_loadu_si128((__m128i*)(b+2));
  GF2X_FUNC(mul7k3_mul2)(lo, a0, b0);
  GF2X_FUNC(mul7k3_mul2c)(hi, a2, b2, d);
  __m128i middle = PXOR(lo[1], hi[0]);
  GF2X_FUNC(mul7k3_mul2)(ab, PXOR(a0, a2), PXOR(b0, b2));
  _mm_storeu_si128((__m128i*)(c + 0), lo[0]);
  _mm_storeu_si128((__m128i*)(c + 2), PXOR3(ab[0], lo[0], middle));
  _mm_storeu_si128((__m128i*)(c + 4), PXOR3(ab[1], hi[1], middle));
  _mm_storeu_si128((__m128i*)(c + 6), hi[1]);
}

/* specialized Karatsuba with 3 calls to mul2, i.e., 9 multiplications,
   assume {d,2} = {a+3,1} * {b+3,1} */
GF2X_STORAGE_CLASS_mul4
void GF2X_FUNC(mul7k3_mul4b) (unsigned long *c, const unsigned long *a, const unsigned long *b, unsigned long *d)
{
  __m128i ab[2];
  __m128i lo[2], hi[2];
  __m128i a0 = _mm_loadu_si128((__m128i*)a);
  __m128i a2 = _mm_loadu_si128((__m128i*)(a+2));
  __m128i b0 = _mm_loadu_si128((__m128i*)b);
  __m128i b2 = _mm_loadu_si128((__m128i*)(b+2));
  GF2X_FUNC(mul7k3_mul2)(lo, a0, b0);
  GF2X_FUNC(mul7k3_mul2b)(hi, a2, b2, d);
  __m128i middle = PXOR(lo[1], hi[0]);
  GF2X_FUNC(mul7k3_mul2)(ab, PXOR(a0, a2), PXOR(b0, b2));
  _mm_storeu_si128((__m128i*)(c + 0), lo[0]);
  _mm_storeu_si128((__m128i*)(c + 2), PXOR3(ab[0], lo[0], middle));
  _mm_storeu_si128((__m128i*)(c + 4), PXOR3(ab[1], hi[1], middle));
  _mm_storeu_si128((__m128i*)(c + 6), hi[1]);
}

#undef PXOR
#undef PXOR3
#undef PZERO

/* based on mul7k.c, version with M(3)+2M(4)-1=23 multiplications */
GF2X_STORAGE_CLASS_mul7
void gf2x_mul7 (unsigned long *c, const unsigned long *a, const unsigned long *b)
{
    unsigned long aa[4], bb[4], ab[8], ab4, ab5, ab6, ab7, d[2];

    gf2x_mul3 (c+8, a+4, b+4);
    GF2X_FUNC(mul7k3_mul4c) (c, a, b, d);
    aa[0] = a[0] ^ a[4];
    aa[1] = a[1] ^ a[5];
    aa[2] = a[2] ^ a[6];
    aa[3] = a[3];
    bb[0] = b[0] ^ b[4];
    bb[1] = b[1] ^ b[5];
    bb[2] = b[2] ^ b[6];
    bb[3] = b[3];
    GF2X_FUNC(mul7k3_mul4b) (ab, aa, bb, d);
    ab4 = ab[4] ^ c[4];
    ab5 = ab[5] ^ c[5];
    ab6 = ab[6] ^ c[6];
    ab7 = ab[7] ^ c[7];
    c[4] ^= ab[0] ^ c[0] ^ c[8];
    c[5] ^= ab[1] ^ c[1] ^ c[9];
    c[6] ^= ab[2] ^ c[2] ^ c[10];
    c[7] ^= ab[3] ^ c[3] ^ c[11];
    c[8] ^= ab4 ^ c[12];
    c[9] ^= ab5 ^ c[13];
    c[10] ^= ab6;
    c[11] ^= ab7;
}

#endif  /* GF2X_MUL7_H_ */
