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
#ifndef GF2X_MUL4_H_
#define GF2X_MUL4_H_

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
/* specialized Karatsuba with 3 calls to mul2, i.e., 9 multiplications */
GF2X_STORAGE_CLASS_mul4
void gf2x_mul4 (unsigned long *c, const unsigned long *a, const unsigned long *b)
{
  __v2di ab[2];
  __v2di lo[2], hi[2];
  __v2di a0 = _mm_loadu_si128((__v2di*)a);
  __v2di a2 = _mm_loadu_si128((__v2di*)(a+2));
  __v2di b0 = _mm_loadu_si128((__v2di*)b);
  __v2di b2 = _mm_loadu_si128((__v2di*)(b+2));
  GF2X_FUNC(mul4clk_mul2)(lo, a0, b0);
  GF2X_FUNC(mul4clk_mul2)(hi, a2, b2);
  __v2di middle = lo[1] ^ hi[0];
  GF2X_FUNC(mul4clk_mul2)(ab, a0 ^ a2, b0 ^ b2);
  _mm_storeu_si128((__v2di*)(c + 0), lo[0]);
  _mm_storeu_si128((__v2di*)(c + 2), ab[0] ^ lo[0] ^ middle);
  _mm_storeu_si128((__v2di*)(c + 4), ab[1] ^ hi[1] ^ middle);
  _mm_storeu_si128((__v2di*)(c + 6), hi[1]);
}

#endif  /* GF2X_MUL4_H_ */
