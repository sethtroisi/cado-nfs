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



#ifndef GF2X_MUL1_H_
#define GF2X_MUL1_H_

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

GF2X_STORAGE_CLASS_mul1 void
gf2x_mul1 (unsigned long *c, unsigned long a, unsigned long b)
{
    __m128i aa = _gf2x_mm_setr_epi64(a, 0);
    __m128i bb = _gf2x_mm_setr_epi64(b, 0);
    _mm_storeu_si128((__m128i*)c, _mm_clmulepi64_si128(aa, bb, 0));
}

GF2X_STORAGE_CLASS_mul_1_n unsigned long
gf2x_mul_1_n (unsigned long *cp, const unsigned long *bp, long sb, unsigned long a)
{
    long i;
    typedef union {
        __m128i s;
        unsigned long x[2];
    } __m128i_proxy;

    __m128i y = _gf2x_mm_set1_epi64(a);
    __m128i x;
    __m128i_proxy cc;


    // do two at a time
    for (i = 0; i + 2 < sb; i += 2) {
        x = _gf2x_mm_setr_epi64(bp[i], bp[i+1]);
        cc.s = _mm_clmulepi64_si128(x, y, 0);
        if (i == 0)
            cp[i] = cc.x[0];
        else
            cp[i] ^= cc.x[0];
        cp[i+1] = cc.x[1];
        cc.s = _mm_clmulepi64_si128(x, y, 1);
        cp[i+1] ^= cc.x[0];
        cp[i+2] = cc.x[1];
    }
    // last is different, to handle carry out
    unsigned long cy;
    if (i == sb - 2) {  // case bp is even
        x = _gf2x_mm_setr_epi64(bp[i], bp[i+1]);
        cc.s = _mm_clmulepi64_si128(x, y, 0);
        if (i == 0)
            cp[i] = cc.x[0];
        else
            cp[i] ^= cc.x[0];
        cp[i+1] = cc.x[1];
        cc.s = _mm_clmulepi64_si128(x, y, 1);
        cp[i+1] ^= cc.x[0];
        cy = cc.x[1];
    } else { //case bp is odd
        x = _gf2x_mm_setr_epi64(bp[i], 0);
        cc.s = _mm_clmulepi64_si128(x, y, 0);
        if (i == 0)
            cp[i] = cc.x[0];
        else
            cp[i] ^= cc.x[0];
        cy = cc.x[1];
    }
    return cy;
}

GF2X_STORAGE_CLASS_addmul_1_n unsigned long
gf2x_addmul_1_n (unsigned long *dp, const unsigned long *cp, const unsigned long* bp, long sb, unsigned long a)
{
    long i;
    typedef union {
        __m128i s;
        unsigned long x[2];
    } __m128i_proxy;

    __m128i y = _gf2x_mm_set1_epi64(a);
    __m128i x;
    __m128i_proxy dd;

    // do two at a time
    for (i = 0; i + 2 < sb; i += 2) {
        x = _gf2x_mm_setr_epi64(bp[i], bp[i+1]);
        dd.s = _mm_clmulepi64_si128(x, y, 0);
        if (i == 0)
            dp[i] = cp[i] ^ dd.x[0];
        else
            dp[i] ^= dd.x[0];
        dp[i+1] = cp[i+1] ^ dd.x[1];
        dd.s = _mm_clmulepi64_si128(x, y, 1);
        dp[i+1] ^= dd.x[0];
        dp[i+2] = cp[i+2] ^ dd.x[1];
    }
    unsigned long cy;
    if (i == sb - 2) {  // case bp is even
        x = _gf2x_mm_setr_epi64(bp[i], bp[i+1]);
        dd.s = _mm_clmulepi64_si128(x, y, 0);
        if (i == 0)
            dp[i] = cp[i] ^ dd.x[0];
        else
            dp[i] ^= dd.x[0];
        dp[i+1] = cp[i+1] ^ dd.x[1];
        dd.s = _mm_clmulepi64_si128(x, y, 1);
        dp[i+1] ^= dd.x[0];
        cy = dd.x[1];
    } else {
        x = _gf2x_mm_setr_epi64(bp[i], 0);
        dd.s = _mm_clmulepi64_si128(x, y, 0);
        if (i == 0)
            dp[i] = cp[i] ^ dd.x[0];
        else
            dp[i] ^= dd.x[0];
        cy = dd.x[1];
    }
    return cy;
}

#endif   /* GF2X_MUL1_H_ */
