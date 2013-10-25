/* This file is part of the gf2x library.

   Copyright 2007, 2008, 2009
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

#ifndef GF2X_MUL4_H_
#define GF2X_MUL4_H_

#include "gf2x.h"
#include "gf2x/gf2x-impl.h"
/* All gf2x source files for lowlevel functions must include gf2x-small.h
 * This is mandatory for the tuning mechanism. */
#include "gf2x/gf2x-small.h"

#include <stdint.h>
#include <emmintrin.h>

#if GF2X_WORDSIZE != 32
#error "This code is for 32-bit only"
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

GF2X_STORAGE_CLASS_mul4
void gf2x_mul4(unsigned long *t, unsigned long const *s1,
               unsigned long const *s2)
{
    typedef union {
        __v2di s;
        unsigned long x[4];
        uint64_t x64[2];
    } __v2di_proxy;
#define SHL(x, r) _mm_slli_epi64((x), (r))
#define SHR(x, r) _mm_srli_epi64((x), (r))
#define SHLD(x, r) _mm_slli_si128((x), (r) >> 3)
#define SHRD(x, r) _mm_srli_si128((x), (r) >> 3)
    __v2di u;
    __v2di t0;
    __v2di t1;
    __v2di t2;

    __v2di g[16];
    /* sequence update walk */
    g[0] = (__v2di) { 0, };
    g[1] = (__v2di) (__v4si) { (long)s2[0], (long)s2[1], (long)s2[2], (long)s2[3], };
    g[2] = SHL(g[1], 1);
    g[3] = g[2] ^ g[1];
    g[4] = SHL(g[2], 1);
    g[5] = g[4] ^ g[1];
    g[6] = SHL(g[3], 1);
    g[7] = g[6] ^ g[1];
    g[8] = SHL(g[4], 1);
    g[9] = g[8] ^ g[1];
    g[10] = SHL(g[5], 1);
    g[11] = g[10] ^ g[1];
    g[12] = SHL(g[6], 1);
    g[13] = g[12] ^ g[1];
    g[14] = SHL(g[7], 1);
    g[15] = g[14] ^ g[1];

    /* round 0 */
    u = g[s1[0]       & 15];
    t0  = u;
    u = g[s1[0] >>  4 & 15];
    t0 ^= SHL(u,  4); t1  = SHR(u, 60);
    u = g[s1[0] >>  8 & 15];
    t0 ^= SHL(u,  8); t1 ^= SHR(u, 56);
    u = g[s1[0] >> 12 & 15];
    t0 ^= SHL(u, 12); t1 ^= SHR(u, 52);
    u = g[s1[0] >> 16 & 15];
    t0 ^= SHL(u, 16); t1 ^= SHR(u, 48);
    u = g[s1[0] >> 20 & 15];
    t0 ^= SHL(u, 20); t1 ^= SHR(u, 44);
    u = g[s1[0] >> 24 & 15];
    t0 ^= SHL(u, 24); t1 ^= SHR(u, 40);
    u = g[s1[0] >> 28 & 15];
    t0 ^= SHL(u, 28); t1 ^= SHR(u, 36);
    u = g[s1[1]       & 15];
    t0 ^= SHL(u, 32); t1 ^= SHR(u, 32);
    u = g[s1[1] >>  4 & 15];
    t0 ^= SHL(u, 36); t1 ^= SHR(u, 28);
    u = g[s1[1] >>  8 & 15];
    t0 ^= SHL(u, 40); t1 ^= SHR(u, 24);
    u = g[s1[1] >> 12 & 15];
    t0 ^= SHL(u, 44); t1 ^= SHR(u, 20);
    u = g[s1[1] >> 16 & 15];
    t0 ^= SHL(u, 48); t1 ^= SHR(u, 16);
    u = g[s1[1] >> 20 & 15];
    t0 ^= SHL(u, 52); t1 ^= SHR(u, 12);
    u = g[s1[1] >> 24 & 15];
    t0 ^= SHL(u, 56); t1 ^= SHR(u,  8);
    u = g[s1[1] >> 28 & 15];
    t0 ^= SHL(u, 60); t1 ^= SHR(u,  4);

    /* round 1 */
    u = g[s1[2]       & 15];
    t1 ^= u;
    u = g[s1[2] >>  4 & 15];
    t1 ^= SHL(u,  4); t2  = SHR(u, 60);
    u = g[s1[2] >>  8 & 15];
    t1 ^= SHL(u,  8); t2 ^= SHR(u, 56);
    u = g[s1[2] >> 12 & 15];
    t1 ^= SHL(u, 12); t2 ^= SHR(u, 52);
    u = g[s1[2] >> 16 & 15];
    t1 ^= SHL(u, 16); t2 ^= SHR(u, 48);
    u = g[s1[2] >> 20 & 15];
    t1 ^= SHL(u, 20); t2 ^= SHR(u, 44);
    u = g[s1[2] >> 24 & 15];
    t1 ^= SHL(u, 24); t2 ^= SHR(u, 40);
    u = g[s1[2] >> 28 & 15];
    t1 ^= SHL(u, 28); t2 ^= SHR(u, 36);
    u = g[s1[3]       & 15];
    t1 ^= SHL(u, 32); t2 ^= SHR(u, 32);
    u = g[s1[3] >>  4 & 15];
    t1 ^= SHL(u, 36); t2 ^= SHR(u, 28);
    u = g[s1[3] >>  8 & 15];
    t1 ^= SHL(u, 40); t2 ^= SHR(u, 24);
    u = g[s1[3] >> 12 & 15];
    t1 ^= SHL(u, 44); t2 ^= SHR(u, 20);
    u = g[s1[3] >> 16 & 15];
    t1 ^= SHL(u, 48); t2 ^= SHR(u, 16);
    u = g[s1[3] >> 20 & 15];
    t1 ^= SHL(u, 52); t2 ^= SHR(u, 12);
    u = g[s1[3] >> 24 & 15];
    t1 ^= SHL(u, 56); t2 ^= SHR(u,  8);
    u = g[s1[3] >> 28 & 15];
    t1 ^= SHL(u, 60); t2 ^= SHR(u,  4);
    /* end */

    /* repair steps */
    /* repair section 200711-200803 */
    __v2di v1 = (__v2di) (__v4si) { (long)s1[0], (long)s1[1], (long)s1[0], (long)s1[1], };
    v1 = SHR(v1, 1);
    __v2di v2 = (__v2di) (__v4si) { (long)s1[2], (long)s1[3], (long)s1[2], (long)s1[3], };
    v2 = SHR(v2, 1);
    __v2di w;
    __v2di m = (__v2di) (__v4si) { 0x77777777L, 0x77777777L, 0x77777777L, 0x77777777L, };
    w = -SHR(g[1],63);
    v1 = v1 & m;
    t1 ^= v1 & w;
    v2 = v2 & m;
    t2 ^= v2 & w;
    w = -SHR(g[2],63);
    v1 = SHR(v1, 1) & m;
    t1 ^= v1 & w;
    v2 = SHR(v2, 1) & m;
    t2 ^= v2 & w;
    w = -SHR(g[4],63);
    v1 = SHR(v1, 1) & m;
    t1 ^= v1 & w;
    v2 = SHR(v2, 1) & m;
    t2 ^= v2 & w;

    /* store result */
    {
        {
            __v2di_proxy r;
            r.s = t0 ^ SHLD(t1, 64);
            t[0] = r.x[0];
            t[1] = r.x[1];
            t[2] = r.x[2];
            t[3] = r.x[3];
        }
        {
            __v2di_proxy r;
            r.s = t2 ^ SHRD(t1, 64);
            t[4] = r.x[0];
            t[5] = r.x[1];
            t[6] = r.x[2];
            t[7] = r.x[3];
        }
    }
#undef SHL
#undef SHR
#undef SHLD
#undef SHRD
}
#endif  /* GF2X_MUL4_H_ */
