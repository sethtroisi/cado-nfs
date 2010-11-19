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
#ifndef GF2X_MUL7_H_
#define GF2X_MUL7_H_

#include "gf2x.h"
/* All gf2x source files for lowlevel functions must include gf2x-small.h
 * This is mandatory for the tuning mechanism. */
#include "gf2x/gf2x-small.h"

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
static inline void
GF2X_FUNC(mul4clk_mul2)(__v2di * t, __v2di ss1, __v2di ss2)
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

static inline void
GF2X_FUNC(mul4clk_mul2c)(__v2di * t, __v2di ss1, __v2di ss2, unsigned long *d)
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
    d[0] = t11.x[0];
    d[1] = t11.x[1];
    tk.s ^= t00.s ^ t11.s;
    t00.x[1] ^= tk.x[0];
    t11.x[0] ^= tk.x[1];
    t[0] = t00.s;
    t[1] = t11.s;
}

static inline void
GF2X_FUNC(mul4clk_mul2b)(__v2di * t, __v2di ss1, __v2di ss2, unsigned long *d)
{
    typedef union {
        __v2di s;
        unsigned long x[2];
    } __v2di_proxy;

    __v2di_proxy t00, t11, tk;
    t00.s = _mm_clmulepi64_si128(ss1, ss2, 0);
    // t11.s = _mm_clmulepi64_si128(ss1, ss2, 17);
    t11.x[0] = d[0];
    t11.x[1] = d[1];
    ss1 ^= _mm_shuffle_epi32(ss1, 78);  // 78 == 0b01001110 (swap)
    ss2 ^= _mm_shuffle_epi32(ss2, 78);  // 78 == 0b01001110 (swap)
    tk.s = _mm_clmulepi64_si128(ss1, ss2, 0);
    tk.s ^= t00.s ^ t11.s;
    t00.x[1] ^= tk.x[0];
    t11.x[0] ^= tk.x[1];
    t[0] = t00.s;
    t[1] = t11.s;
}

/* specialized Karatsuba with 3 calls to mul2, i.e., 9 multiplications
   {d,2} <- {a+3,1} * {b+3,1} */
GF2X_STORAGE_CLASS_mul4
void gf2x_mul4c (unsigned long *c, const unsigned long *a, const unsigned long *b, unsigned long *d)
{
  __v2di ab[2];
  __v2di lo[2], hi[2];
  __v2di a0 = _mm_loadu_si128((__v2di*)a);
  __v2di a2 = _mm_loadu_si128((__v2di*)(a+2));
  __v2di b0 = _mm_loadu_si128((__v2di*)b);
  __v2di b2 = _mm_loadu_si128((__v2di*)(b+2));
  GF2X_FUNC(mul4clk_mul2)(lo, a0, b0);
  GF2X_FUNC(mul4clk_mul2c)(hi, a2, b2, d);
  __v2di middle = lo[1] ^ hi[0];
  GF2X_FUNC(mul4clk_mul2)(ab, a0 ^ a2, b0 ^ b2);
  _mm_storeu_si128((__v2di*)(c + 0), lo[0]);
  _mm_storeu_si128((__v2di*)(c + 2), ab[0] ^ lo[0] ^ middle);
  _mm_storeu_si128((__v2di*)(c + 4), ab[1] ^ hi[1] ^ middle);
  _mm_storeu_si128((__v2di*)(c + 6), hi[1]);
}

/* specialized Karatsuba with 3 calls to mul2, i.e., 9 multiplications,
   assume {d,2} = {a+3,1} * {b+3,1} */
GF2X_STORAGE_CLASS_mul4
void gf2x_mul4b (unsigned long *c, const unsigned long *a, const unsigned long *b, unsigned long *d)
{
  __v2di ab[2];
  __v2di lo[2], hi[2];
  __v2di a0 = _mm_loadu_si128((__v2di*)a);
  __v2di a2 = _mm_loadu_si128((__v2di*)(a+2));
  __v2di b0 = _mm_loadu_si128((__v2di*)b);
  __v2di b2 = _mm_loadu_si128((__v2di*)(b+2));
  GF2X_FUNC(mul4clk_mul2)(lo, a0, b0);
  GF2X_FUNC(mul4clk_mul2b)(hi, a2, b2, d);
  __v2di middle = lo[1] ^ hi[0];
  GF2X_FUNC(mul4clk_mul2)(ab, a0 ^ a2, b0 ^ b2);
  _mm_storeu_si128((__v2di*)(c + 0), lo[0]);
  _mm_storeu_si128((__v2di*)(c + 2), ab[0] ^ lo[0] ^ middle);
  _mm_storeu_si128((__v2di*)(c + 4), ab[1] ^ hi[1] ^ middle);
  _mm_storeu_si128((__v2di*)(c + 6), hi[1]);
}

/* based on mul7k.c, version with M(3)+2M(4)-1=23 multiplications */
GF2X_STORAGE_CLASS_mul7
void gf2x_mul7 (unsigned long *c, const unsigned long *a, const unsigned long *b)
{
    unsigned long aa[4], bb[4], ab[8], ab4, ab5, ab6, ab7, d[2];

    gf2x_mul3 (c+8, a+4, b+4);
    gf2x_mul4c (c, a, b, d);
    aa[0] = a[0] ^ a[4];
    aa[1] = a[1] ^ a[5];
    aa[2] = a[2] ^ a[6];
    aa[3] = a[3];
    bb[0] = b[0] ^ b[4];
    bb[1] = b[1] ^ b[5];
    bb[2] = b[2] ^ b[6];
    bb[3] = b[3];
    gf2x_mul4b (ab, aa, bb, d);
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
