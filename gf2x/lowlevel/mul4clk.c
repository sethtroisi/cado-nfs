/* This file is part of the gf2x library.

   Copyright 2010, 2013, 2015
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

#ifndef GF2X_MUL4_H_
#define GF2X_MUL4_H_

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
#define PXOR(lop, rop) _mm_xor_si128((lop), (rop))
#define PZERO    _mm_setzero_si128()
static inline void
GF2X_FUNC(mul4clk_mul2)(__m128i * t, __m128i ss1, __m128i ss2)
{
    __m128i t00 = _mm_clmulepi64_si128(ss1, ss2, 0);
    __m128i t11 = _mm_clmulepi64_si128(ss1, ss2, 0x11);
    ss1 = PXOR(ss1, _mm_shuffle_epi32(ss1, _MM_SHUFFLE(1,0,3,2)));
    ss2 = PXOR(ss2, _mm_shuffle_epi32(ss2, _MM_SHUFFLE(1,0,3,2)));
    __m128i tk = PXOR(PXOR(t00, t11), _mm_clmulepi64_si128(ss1, ss2, 0));
    t[0] = PXOR(t00, _mm_unpacklo_epi64(PZERO, tk));
    t[1] = PXOR(t11, _mm_unpackhi_epi64(tk, PZERO));
}
/* specialized Karatsuba with 3 calls to mul2, i.e., 9 multiplications */
GF2X_STORAGE_CLASS_mul4
void gf2x_mul4 (unsigned long *c, const unsigned long *a, const unsigned long *b)
{
  __m128i ab[2];
  __m128i lo[2], hi[2];
  __m128i a0 = _mm_loadu_si128((__m128i*)a);
  __m128i a2 = _mm_loadu_si128((__m128i*)(a+2));
  __m128i b0 = _mm_loadu_si128((__m128i*)b);
  __m128i b2 = _mm_loadu_si128((__m128i*)(b+2));
  GF2X_FUNC(mul4clk_mul2)(lo, a0, b0);
  GF2X_FUNC(mul4clk_mul2)(hi, a2, b2);
  __m128i middle = PXOR(lo[1], hi[0]);
  GF2X_FUNC(mul4clk_mul2)(ab, PXOR(a0, a2), PXOR(b0, b2));
  _mm_storeu_si128((__m128i*)(c + 0), lo[0]);
  _mm_storeu_si128((__m128i*)(c + 2), PXOR(ab[0], PXOR(lo[0], middle)));
  _mm_storeu_si128((__m128i*)(c + 4), PXOR(ab[1], PXOR(hi[1], middle)));
  _mm_storeu_si128((__m128i*)(c + 6), hi[1]);
}
#undef PXOR
#undef PZERO

#endif  /* GF2X_MUL4_H_ */
