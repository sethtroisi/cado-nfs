/* This file is part of the gf2x library.

   Copyright 2007, 2008, 2009
   Richard Brent, Pierrick Gaudry, Emmanuel Thome', Paul Zimmermann

   This program is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by the
   Free Software Foundation; either version 2 of the License, or (at your
   option) any later version.

   This program is distributed in the hope that it will be useful, but WITHOUT
   ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
   FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
   more details.

   You should have received a copy of the GNU General Public License along
   with this program; see the file COPYING.  If not, write to the Free
   Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
   02111-1307, USA.
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
static inline __v2di
GF2X_FUNC(mul4clk_c_mul1) (unsigned long a, unsigned long b)
{
    __v2di aa = (__v2di) { (long long)a, 0 };
    __v2di bb = (__v2di) { (long long)b, 0 };
    return _mm_clmulepi64_si128(aa, bb, 0);
}
GF2X_STORAGE_CLASS_mul4
void gf2x_mul4 (unsigned long *c, const unsigned long *a, const unsigned long *b)
{
  // <http://eprint.iacr.org/2011/540>
  __v2di m0, m1, m2, m3, m4, m5, m6, m7, m8, t0, t1;
  unsigned long aa4, aa5, aa6, aa7, aa8;
  unsigned long bb4, bb5, bb6, bb7, bb8;

  aa4 = a[0] ^ a[1]; bb4 = b[0] ^ b[1];
  aa5 = a[2] ^ a[3]; bb5 = b[2] ^ b[3];
  aa6 = a[0] ^ a[2]; bb6 = b[0] ^ b[2];
  aa7 = a[1] ^ a[3]; bb7 = b[1] ^ b[3];
  aa8 = aa4 ^ aa5;   bb8 = bb4 ^ bb5;

  m0 = GF2X_FUNC(mul4clk_c_mul1)( a[0],  b[0]);
  m1 = GF2X_FUNC(mul4clk_c_mul1)( a[1],  b[1]);
  m2 = GF2X_FUNC(mul4clk_c_mul1)( a[2],  b[2]);
  m3 = GF2X_FUNC(mul4clk_c_mul1)( a[3],  b[3]);
  m4 = GF2X_FUNC(mul4clk_c_mul1)(aa4, bb4);
  m5 = GF2X_FUNC(mul4clk_c_mul1)(aa5, bb5);
  m6 = GF2X_FUNC(mul4clk_c_mul1)(aa6, bb6);
  m7 = GF2X_FUNC(mul4clk_c_mul1)(aa7, bb7);
  m8 = GF2X_FUNC(mul4clk_c_mul1)(aa8, bb8);

  t0 = m0 ^ m1;
  t1 = m2 ^ m3;

  __v2di ce0 = m0;
  __v2di ce2 = t0 ^ m2 ^ m6;
  __v2di ce4 = t1 ^ m1 ^ m7;
  __v2di ce6 = m3;

  __v2di co1 = t0 ^ m4;
  __v2di co5 = t1 ^ m5;
  __v2di co3 = co1 ^ co5 ^ m6 ^ m7 ^ m8;

  _mm_storeu_si128((__v2di*)(c),   ce0 ^ _mm_slli_si128(co1, 8));
  _mm_storeu_si128((__v2di*)(c+2), ce2 ^ _mm_srli_si128(co1, 8) ^ _mm_slli_si128(co3, 8));
  _mm_storeu_si128((__v2di*)(c+4), ce4 ^ _mm_srli_si128(co3, 8) ^ _mm_slli_si128(co5, 8));
  _mm_storeu_si128((__v2di*)(c+6), ce6 ^ _mm_srli_si128(co5, 8));
}

#endif  /* GF2X_MUL4_H_ */
