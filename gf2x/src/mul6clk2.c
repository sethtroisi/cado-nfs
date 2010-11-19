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
#ifndef GF2X_MUL6_H_
#define GF2X_MUL6_H_

#include "gf2x.h"
/* All gf2x source files for lowlevel functions must include gf2x-small.h
 * This is mandatory for the tuning mechanism. */
#include "gf2x/gf2x-small.h"

#include <stdint.h>
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
/* This specialized version avoids loads, and relies on the destination
 * being aligned, so that aligned stores are possible */
static inline void
GF2X_FUNC(mul6clk2_mul2)(__v2di * t, __v2di ss1, __v2di ss2)
{
    typedef union {
        __v2di s;
        unsigned long x[2];
    } __v2di_proxy;

    __v2di_proxy t00, t11, tk;
    t00.s = _mm_clmulepi64_si128(ss1, ss2, 0);
    t11.s = _mm_clmulepi64_si128(ss1, ss2, 17);
    ss1 ^= _mm_shuffle_epi32(ss1, 78);  // 78 == 0b01001110 (swap)
    ss2 ^= _mm_shuffle_epi32(ss2, 78);  // 78 == 0b01001110 (swap)
    tk.s = _mm_clmulepi64_si128(ss1, ss2, 0);
    tk.s ^= t00.s ^ t11.s;
    t00.x[1] ^= tk.x[0];
    t11.x[0] ^= tk.x[1];
    t[0] = t00.s;
    t[1] = t11.s;
}

/* variant with 6 calls to mul2, i.e., 18 multiplications */
GF2X_STORAGE_CLASS_mul6
void gf2x_mul6 (unsigned long *c, const unsigned long *a, const unsigned long *b)
{
    __v2di aa[3], bb[3];
    __v2di p0[2], p1[2], p2[2];
    __v2di pp0[2], pp1[2], pp2[2];
    __v2di a0 = _mm_loadu_si128((__v2di*)(a));
    __v2di a1 = _mm_loadu_si128((__v2di*)(a+2));
    __v2di a2 = _mm_loadu_si128((__v2di*)(a+4));
    __v2di b0 = _mm_loadu_si128((__v2di*)(b));
    __v2di b1 = _mm_loadu_si128((__v2di*)(b+2));
    __v2di b2 = _mm_loadu_si128((__v2di*)(b+4));
    aa[0] = a1^a2;
    aa[1] = a0^a2;
    aa[2] = a0^a1;
    bb[0] = b1^b2;
    bb[1] = b0^b2;
    bb[2] = b0^b1;
    GF2X_FUNC(mul6clk2_mul2)(p0, a0, b0);
    GF2X_FUNC(mul6clk2_mul2)(p1, a1, b1);
    GF2X_FUNC(mul6clk2_mul2)(p2, a2, b2);
    GF2X_FUNC(mul6clk2_mul2)(pp0, aa[0], bb[0]);
    GF2X_FUNC(mul6clk2_mul2)(pp1, aa[1], bb[1]);
    GF2X_FUNC(mul6clk2_mul2)(pp2, aa[2], bb[2]);
    _mm_storeu_si128((__v2di*)(c + 2*0), p0[0]);
    _mm_storeu_si128((__v2di*)(c + 2*1), p0[0]^p1[0]^pp2[0]       ^ p0[1]);
    _mm_storeu_si128((__v2di*)(c + 2*2), p0[0]^p1[0]^p2[0]^pp1[0] ^ p0[1]^p1[1]^pp2[1]);
    _mm_storeu_si128((__v2di*)(c + 2*3), pp0[0]^p1[0]^p2[0]       ^ p0[1]^p1[1]^p2[1]^pp1[1]);
    _mm_storeu_si128((__v2di*)(c + 2*4), p2[0]                    ^ pp0[1]^p1[1]^p2[1]);
    _mm_storeu_si128((__v2di*)(c + 2*5),                            p2[1]);
}

#endif  /* GF2X_MUL6_H_ */
