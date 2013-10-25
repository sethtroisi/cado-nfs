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
/* Implements 128x128 -> 256 bit product using SSE2 instructions. */

#ifndef GF2X_MUL2_H_
#define GF2X_MUL2_H_

#include "gf2x.h"
/* All gf2x source files for lowlevel functions must include gf2x-small.h
 * This is mandatory for the tuning mechanism. */
#include "gf2x/gf2x-small.h"

#include "gf2x/gf2x-impl.h"

#include <stdint.h>
#include <emmintrin.h>

#if GF2X_WORDSIZE != 64
#error "This code is for 64-bit only"
#endif

#include "gf2x/gf2x-config.h"

#ifndef HAVE_SSE2_SUPPORT
#error "This code needs sse-2 support"
#endif

#ifndef	GNUC_VERSION
#define GNUC_VERSION(X,Y,Z)     \
    (defined(__GNUC__) &&        \
    (__GNUC__ == X && __GNUC_MINOR__ == Y && __GNUC_PATCHLEVEL__ == Z))
#endif
#if (GNUC_VERSION(4,3,0) || GNUC_VERSION(4,3,1))
#warning "Your GCC version is buggy. Binary fields may fail randomly"
/* Gcc bug reports 37101 and 37340 -- the only convenient fix is to
 * upgrade to 4.3.2 */
#endif

GF2X_STORAGE_CLASS_mul2
void gf2x_mul2(unsigned long * t, unsigned long const * s1,
        unsigned long const * s2)
{
#define SHL(x, r) _mm_slli_epi64((x), (r))
#define SHR(x, r) _mm_srli_epi64((x), (r))
#define SHLD(x, r) _mm_slli_si128((x), (r) >> 3)
#define SHRD(x, r) _mm_srli_si128((x), (r) >> 3)
    __v2di u;
    __v2di t0;
    __v2di t1;
    __v2di t2;

    __v2di g[16];
    __v2di w;
    __v2di m = (__v2di) { 0xeeeeeeeeeeeeeeeeLL, 0xeeeeeeeeeeeeeeeeLL, };
    /* sequence update walk */
    g[ 0] = (__v2di) { 0, };
    __v2di b0 = _mm_loadu_si128((__v2di*) s2);
    g[ 1] = b0;
    __v2di v1 = (__v2di) { (long long)s1[0], (long long)s1[0], };
    w = -SHR(b0,63);
    __v2di v2 = (__v2di) { (long long)s1[1], (long long)s1[1], };
    v1 = SHR(v1 & m, 1); t1 = v1 & w;
    g[ 2] = SHL(b0, 1); g[ 3] = g[ 2] ^ b0;
    v2 = SHR(v2 & m, 1); t2 = v2 & w;
    g[ 4] = SHL(g[ 2], 1); g[ 5] = g[ 4] ^ b0;
    w = -SHR(g[ 2],63);
    g[ 6] = SHL(g[ 3], 1); g[ 7] = g[ 6] ^ b0;
    v1 = SHR(v1 & m, 1); t1 ^= v1 & w;
    g[ 8] = SHL(g[ 4], 1); g[ 9] = g[ 8] ^ b0;
    v2 = SHR(v2 & m, 1); t2 ^= v2 & w;
    g[10] = SHL(g[ 5], 1); g[11] = g[10] ^ b0;
    w = -SHR(g[4],63);
    g[12] = SHL(g[ 6], 1); g[13] = g[12] ^ b0;
    v1 = SHR(v1 & m, 1); t1 ^= v1 & w;
    g[14] = SHL(g[ 7], 1); g[15] = g[14] ^ b0;
    v2 = SHR(v2 & m, 1); t2 ^= v2 & w;



    /* round 0 */
    u = g[s1[0]       & 15]; t0  = u;
    u = g[s1[0] >>  4 & 15]; t0 ^= SHL(u,  4); t1 ^= SHR(u, 60);
    u = g[s1[0] >>  8 & 15]; t0 ^= SHL(u,  8); t1 ^= SHR(u, 56);
    u = g[s1[0] >> 12 & 15]; t0 ^= SHL(u, 12); t1 ^= SHR(u, 52);
    u = g[s1[0] >> 16 & 15]; t0 ^= SHL(u, 16); t1 ^= SHR(u, 48);
    u = g[s1[0] >> 20 & 15]; t0 ^= SHL(u, 20); t1 ^= SHR(u, 44);
    u = g[s1[0] >> 24 & 15]; t0 ^= SHL(u, 24); t1 ^= SHR(u, 40);
    u = g[s1[0] >> 28 & 15]; t0 ^= SHL(u, 28); t1 ^= SHR(u, 36);
    u = g[s1[0] >> 32 & 15]; t0 ^= SHL(u, 32); t1 ^= SHR(u, 32);
    u = g[s1[0] >> 36 & 15]; t0 ^= SHL(u, 36); t1 ^= SHR(u, 28);
    u = g[s1[0] >> 40 & 15]; t0 ^= SHL(u, 40); t1 ^= SHR(u, 24);
    u = g[s1[0] >> 44 & 15]; t0 ^= SHL(u, 44); t1 ^= SHR(u, 20);
    u = g[s1[0] >> 48 & 15]; t0 ^= SHL(u, 48); t1 ^= SHR(u, 16);
    u = g[s1[0] >> 52 & 15]; t0 ^= SHL(u, 52); t1 ^= SHR(u, 12);
    u = g[s1[0] >> 56 & 15]; t0 ^= SHL(u, 56); t1 ^= SHR(u,  8);
    u = g[s1[0] >> 60 & 15]; t0 ^= SHL(u, 60); t1 ^= SHR(u,  4);

    /* round 1 */
    u = g[s1[1]       & 15]; t1 ^= u;
    u = g[s1[1] >>  4 & 15]; t1 ^= SHL(u,  4); t2 ^= SHR(u, 60);
    u = g[s1[1] >>  8 & 15]; t1 ^= SHL(u,  8); t2 ^= SHR(u, 56);
    u = g[s1[1] >> 12 & 15]; t1 ^= SHL(u, 12); t2 ^= SHR(u, 52);
    u = g[s1[1] >> 16 & 15]; t1 ^= SHL(u, 16); t2 ^= SHR(u, 48);
    u = g[s1[1] >> 20 & 15]; t1 ^= SHL(u, 20); t2 ^= SHR(u, 44);
    u = g[s1[1] >> 24 & 15]; t1 ^= SHL(u, 24); t2 ^= SHR(u, 40);
    u = g[s1[1] >> 28 & 15]; t1 ^= SHL(u, 28); t2 ^= SHR(u, 36);
    u = g[s1[1] >> 32 & 15]; t1 ^= SHL(u, 32); t2 ^= SHR(u, 32);
    u = g[s1[1] >> 36 & 15]; t1 ^= SHL(u, 36); t2 ^= SHR(u, 28);
    u = g[s1[1] >> 40 & 15]; t1 ^= SHL(u, 40); t2 ^= SHR(u, 24);
    u = g[s1[1] >> 44 & 15]; t1 ^= SHL(u, 44); t2 ^= SHR(u, 20);
    u = g[s1[1] >> 48 & 15]; t1 ^= SHL(u, 48); t2 ^= SHR(u, 16);
    u = g[s1[1] >> 52 & 15]; t1 ^= SHL(u, 52); t2 ^= SHR(u, 12);
    u = g[s1[1] >> 56 & 15]; t1 ^= SHL(u, 56); t2 ^= SHR(u,  8);
    u = g[s1[1] >> 60 & 15]; t1 ^= SHL(u, 60); t2 ^= SHR(u,  4);
    /* end */

    /* store result */
    _mm_storeu_si128((__v2di*)t, t0 ^ SHLD(t1, 64));
    _mm_storeu_si128((__v2di*)(t+2), t2 ^ SHRD(t1, 64));
}

#endif  /* GF2X_MUL2_H_ */
