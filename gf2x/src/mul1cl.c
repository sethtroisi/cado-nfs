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


#ifndef GF2X_MUL1_H_
#define GF2X_MUL1_H_

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

GF2X_STORAGE_CLASS_mul1 void
gf2x_mul1 (unsigned long *c, unsigned long a, unsigned long b)
{
    __v2di aa = (__v2di) { (long long)a, 0 };
    __v2di bb = (__v2di) { (long long)b, 0 };
    _mm_storeu_si128((__v2di*)c, _mm_clmulepi64_si128(aa, bb, 0));
}

GF2X_STORAGE_CLASS_mul_1_n unsigned long
gf2x_mul_1_n (unsigned long *cp, const unsigned long *bp, long sb, unsigned long a)
{
    long i;
    typedef union {
        __v2di s;
        unsigned long x[2];
    } __v2di_proxy;

    __v2di y = (__v2di) { (long long)a, (long long)a };
    __v2di x;
    __v2di_proxy cc;


    // do two at a time
    for (i = 0; i + 2 < sb; i += 2) { 
        x = (__v2di) { (long long)bp[i], (long long)bp[i+1] };
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
        x = (__v2di) { (long long)bp[i], (long long)bp[i+1] }; 
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
        x = (__v2di) { (long long)bp[i], 0 }; 
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
        __v2di s;
        unsigned long x[2];
    } __v2di_proxy;

    __v2di y = (__v2di) { (long long)a, (long long)a };
    __v2di x;
    __v2di_proxy dd;

    // do two at a time
    for (i = 0; i + 2 < sb; i += 2) { 
        x = (__v2di) { (long long)bp[i], (long long)bp[i+1] };
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
        x = (__v2di) { (long long)bp[i], (long long)bp[i+1] }; 
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
        x = (__v2di) { (long long)bp[i], 0 };
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
