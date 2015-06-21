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

/* Implements 128x128 -> 256 bit product using pclmulqdq instruction. */

#ifndef GF2X_MUL2_H_
#define GF2X_MUL2_H_

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

/* Karatsuba with 3 multiplications */
GF2X_STORAGE_CLASS_mul2
void gf2x_mul2(unsigned long * t, unsigned long const * s1,
        unsigned long const * s2)
{
#define PXOR(lop, rop) _mm_xor_si128((lop), (rop))
    __m128i ss1 = _mm_loadu_si128((__m128i *)s1);
    __m128i ss2 = _mm_loadu_si128((__m128i *)s2);


    __m128i t00 = _mm_clmulepi64_si128(ss1, ss2, 0);
    __m128i t11 = _mm_clmulepi64_si128(ss1, ss2, 0x11);

    ss1 = PXOR(ss1, _mm_shuffle_epi32(ss1, _MM_SHUFFLE(1,0,3,2)));
    ss2 = PXOR(ss2, _mm_shuffle_epi32(ss2, _MM_SHUFFLE(1,0,3,2)));

    __m128i tk = PXOR(t00, PXOR(t11, _mm_clmulepi64_si128(ss1, ss2, 0)));

    /* mul2cl.c is essentially identical, just replaces srli and srli by
     * unpacklo and unpackhi */
    _mm_storeu_si128((__m128i *)t,     PXOR(t00, _mm_slli_si128(tk, 8)));
    _mm_storeu_si128((__m128i *)(t+2), PXOR(t11, _mm_srli_si128(tk, 8)));
#undef PXOR
}
#endif  /* GF2X_MUL2_H_ */
