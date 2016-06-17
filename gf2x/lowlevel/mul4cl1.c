/* This file is part of the gf2x library.

   Copyright 2007, 2008, 2009, 2010, 2011, 2012, 2013, 2015
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
static inline __m128i
GF2X_FUNC(mul4cl1_mul1) (unsigned long a, unsigned long b)
{
    __m128i aa = _gf2x_mm_setr_epi64(a, 0);
    __m128i bb = _gf2x_mm_setr_epi64(b, 0);
    return _mm_clmulepi64_si128(aa, bb, 0);
}

#define PXOR(lop, rop) _mm_xor_si128((lop), (rop))
#define PXOR3(op1, op2, op3) PXOR(op1, PXOR(op2, op3))
#define PXOR4(op1, op2, op3, op4) PXOR(op1, PXOR3(op2, op3, op4))
#define PXOR5(op1, op2, op3, op4, op5) PXOR(op1, PXOR4(op2, op3, op4, op5))
#define PXOR6(op1, op2, op3, op4, op5, op6) PXOR(op1, PXOR5(op2, op3, op4, op5, op6))
#define PXOR7(op1, op2, op3, op4, op5, op6, op7) PXOR(op1, PXOR6(op2, op3, op4, op5, op6, op7))
#define PZERO    _mm_setzero_si128()


GF2X_STORAGE_CLASS_mul4
void gf2x_mul4 (unsigned long *c, const unsigned long *a, const unsigned long *b)
{
  // <http://eprint.iacr.org/2011/540>
  __m128i m0, m1, m2, m3, m4, m5, m6, m7, m8, t0, t1;
  unsigned long aa4, aa5, aa6, aa7, aa8;
  unsigned long bb4, bb5, bb6, bb7, bb8;

  aa4 = a[0] ^ a[1]; bb4 = b[0] ^ b[1];
  aa5 = a[2] ^ a[3]; bb5 = b[2] ^ b[3];
  aa6 = a[0] ^ a[2]; bb6 = b[0] ^ b[2];
  aa7 = a[1] ^ a[3]; bb7 = b[1] ^ b[3];
  aa8 = aa4 ^ aa5;   bb8 = bb4 ^ bb5;

  m0 = GF2X_FUNC(mul4cl1_mul1)( a[0],  b[0]);
  m1 = GF2X_FUNC(mul4cl1_mul1)( a[1],  b[1]);
  m2 = GF2X_FUNC(mul4cl1_mul1)( a[2],  b[2]);
  m3 = GF2X_FUNC(mul4cl1_mul1)( a[3],  b[3]);
  m4 = GF2X_FUNC(mul4cl1_mul1)(aa4, bb4);
  m5 = GF2X_FUNC(mul4cl1_mul1)(aa5, bb5);
  m6 = GF2X_FUNC(mul4cl1_mul1)(aa6, bb6);
  m7 = GF2X_FUNC(mul4cl1_mul1)(aa7, bb7);
  m8 = GF2X_FUNC(mul4cl1_mul1)(aa8, bb8);

  t0 = PXOR(m0, m1);
  t1 = PXOR(m2, m3);

  __m128i ce0 = m0;
  __m128i ce2 = PXOR3(t0, m2, m6);
  __m128i ce4 = PXOR3(t1, m1, m7);
  __m128i ce6 = m3;

  __m128i co1 = PXOR(t0, m4);
  __m128i co5 = PXOR(t1, m5);
  __m128i co3 = PXOR5(co1, co5, m6, m7, m8);

  _mm_storeu_si128((__m128i*)(c),   PXOR(ce0, _mm_slli_si128(co1, 8)));
  _mm_storeu_si128((__m128i*)(c+2), PXOR3(ce2, _mm_srli_si128(co1, 8), _mm_slli_si128(co3, 8)));
  _mm_storeu_si128((__m128i*)(c+4), PXOR3(ce4, _mm_srli_si128(co3, 8), _mm_slli_si128(co5, 8)));
  _mm_storeu_si128((__m128i*)(c+6), PXOR(ce6, _mm_srli_si128(co5, 8)));
}

#undef PXOR
#undef PXOR3
#undef PXOR4
#undef PXOR5
#undef PXOR6
#undef PXOR7
#undef PZERO
#endif  /* GF2X_MUL4_H_ */
