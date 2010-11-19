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
/* Implements 128x128 -> 256 bit product using pclmulqdq instruction. */

#ifndef GF2X_MUL2_H_
#define GF2X_MUL2_H_

#include "gf2x.h"
#include "gf2x/gf2x-impl.h"
/* All gf2x source files for lowlevel functions must include gf2x-small.h
 * This is mandatory for the tuning mechanism. */
#include "gf2x/gf2x-small.h"

#include <stdint.h>
#include <wmmintrin.h>

#if GF2X_WORDSIZE != 64
#error "This code is for 64-bit only"
#endif

#include "gf2x/gf2x-config.h"

#ifndef HAVE_PCLMUL_SUPPORT
#error "This code needs pclmul support"
#endif

/* Karatsuba with 3 multiplications */
GF2X_STORAGE_CLASS_mul2
void gf2x_mul2(unsigned long * t, unsigned long const * s1,
        unsigned long const * s2)
{
    __v2di ss1, ss2, s1s, s2s;
    __v2di t00, t11, tk;
    ss1 = _mm_loadu_si128((__v2di *)s1);
    ss2 = _mm_loadu_si128((__v2di *)s2);


    t00 = _mm_clmulepi64_si128(ss1, ss2, 0);
    t11 = _mm_clmulepi64_si128(ss1, ss2, 17);
    
    s1s = _mm_shuffle_epi32(ss1, 78);
    ss1 ^= s1s;
    s2s = _mm_shuffle_epi32(ss2, 78);
    ss2 ^= s2s;
    
    tk = t00 ^ t11 ^ _mm_clmulepi64_si128(ss1, ss2, 0);

    _mm_storeu_si128((__v2di *)t, t00 ^ _mm_slli_si128(tk, 8));
    _mm_storeu_si128((__v2di *)(t+2), t11 ^ _mm_srli_si128(tk, 8));
}
#endif  /* GF2X_MUL2_H_ */
