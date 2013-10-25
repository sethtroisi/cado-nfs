/* This file is part of the gf2x library.

   Copyright 2007, 2008, 2009, 2010
   Richard Brent, Pierrick Gaudry, Emmanuel Thome', Paul Zimmermann
   Nicolas Estibals (for this file)

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
GF2X_FUNC(mul7clk_c_mul1) (unsigned long a, unsigned long b)
{   
    __v2di aa = (__v2di) { (long long)a, 0 };
    __v2di bb = (__v2di) { (long long)b, 0 };
    return _mm_clmulepi64_si128(aa, bb, 0);
}

/* variant with 22 multiplications */
GF2X_STORAGE_CLASS_mul7
void gf2x_mul7 (unsigned long *c, const unsigned long *a, const unsigned long *b)
{
    /* Montgomery formulae with 22 multiplications, see
     Five, Six, and Seven-Term {K}aratsuba-Like Formulae,
     IEEE Transactions on Computers, volume 54, number 3, p. 362-369, 2005 */
    unsigned long ta[5], tb[5], pa[22], pb[22];
    ta[0] = a[0] ^ a[4];
    ta[1] = a[3] ^ a[5];
    ta[2] = a[2] ^ a[6];
    ta[3] = a[1] ^ ta[0];
    ta[4] = ta[1] ^ ta[2];
    pa[0] = a[6];
    pa[1] = a[5];
    pa[2] = a[5] ^ a[6];
    pa[3] = a[4];
    pa[4] = a[4] ^ a[6];
    pa[5] = a[3];
    pa[6] = ta[1];
    pa[7] = a[2];
    pa[8] = ta[2];
    pa[9] = a[1];
    pa[10]= a[1] ^ a[3];
    pa[11]= a[1] ^ a[2] ^ a[4] ^ a[5];
    pa[12]= a[1] ^ ta[4];
    pa[13]= a[0];
    pa[14]= ta[0];
    pa[15]= a[0] ^ a[2];
    pa[16]= a[0] ^ ta[4];
    pa[17]= a[3] ^ ta[0] ^ ta[2];
    pa[18]= a[0] ^ a[1];
    pa[19]= a[3] ^ a[6] ^ ta[3];
    pa[20]= ta[1] ^ ta[3];
    pa[21]= ta[3] ^ ta[4];

    tb[0] = b[0] ^ b[4];
    tb[1] = b[3] ^ b[5];
    tb[2] = b[2] ^ b[6];
    tb[3] = b[1] ^ tb[0];
    tb[4] = tb[1] ^ tb[2];
    pb[0] = b[6];
    pb[1] = b[5];
    pb[2] = b[5] ^ b[6];
    pb[3] = b[4];
    pb[4] = b[4] ^ b[6];
    pb[5] = b[3];
    pb[6] = tb[1];
    pb[7] = b[2];
    pb[8] = tb[2];
    pb[9] = b[1];
    pb[10]= b[1] ^ b[3];
    pb[11]= b[1] ^ b[2] ^ b[4] ^ b[5];
    pb[12]= b[1] ^ tb[4];
    pb[13]= b[0];
    pb[14]= tb[0];
    pb[15]= b[0] ^ b[2];
    pb[16]= b[0] ^ tb[4];
    pb[17]= b[3] ^ tb[0] ^ tb[2];
    pb[18]= b[0] ^ b[1];
    pb[19]= b[3] ^ b[6] ^ tb[3];
    pb[20]= tb[1] ^ tb[3];
    pb[21]= tb[3] ^ tb[4];

    __v2di p[22];
    
    p[0] = GF2X_FUNC(mul7clk_c_mul1)(pa[0], pb[0]);
    p[1] = GF2X_FUNC(mul7clk_c_mul1)(pa[1], pb[1]);
    p[2] = GF2X_FUNC(mul7clk_c_mul1)(pa[2], pb[2]);
    p[3] = GF2X_FUNC(mul7clk_c_mul1)(pa[3], pb[3]);
    p[4] = GF2X_FUNC(mul7clk_c_mul1)(pa[4], pb[4]);
    p[5] = GF2X_FUNC(mul7clk_c_mul1)(pa[5], pb[5]);
    p[6] = GF2X_FUNC(mul7clk_c_mul1)(pa[6], pb[6]);
    p[7] = GF2X_FUNC(mul7clk_c_mul1)(pa[7], pb[7]);
    p[8] = GF2X_FUNC(mul7clk_c_mul1)(pa[8], pb[8]);
    p[9] = GF2X_FUNC(mul7clk_c_mul1)(pa[9], pb[9]);
    p[10]= GF2X_FUNC(mul7clk_c_mul1)(pa[10], pb[10]);
    p[11]= GF2X_FUNC(mul7clk_c_mul1)(pa[11], pb[11]);
    p[12]= GF2X_FUNC(mul7clk_c_mul1)(pa[12], pb[12]);
    p[13]= GF2X_FUNC(mul7clk_c_mul1)(pa[13], pb[13]);
    p[14]= GF2X_FUNC(mul7clk_c_mul1)(pa[14], pb[14]);
    p[15]= GF2X_FUNC(mul7clk_c_mul1)(pa[15], pb[15]);
    p[16]= GF2X_FUNC(mul7clk_c_mul1)(pa[16], pb[16]);
    p[17]= GF2X_FUNC(mul7clk_c_mul1)(pa[17], pb[17]);
    p[18]= GF2X_FUNC(mul7clk_c_mul1)(pa[18], pb[18]);
    p[19]= GF2X_FUNC(mul7clk_c_mul1)(pa[19], pb[19]);
    p[20]= GF2X_FUNC(mul7clk_c_mul1)(pa[20], pb[20]);
    p[21]= GF2X_FUNC(mul7clk_c_mul1)(pa[21], pb[21]);

    __v2di t[13];

    t[0]  = p[0] ^ p[1];
    t[1]  = p[9] ^ p[13];
    t[2]  = p[3] ^ p[6];
    t[3]  = p[7] ^ p[10];
    t[4]  = p[11] ^ p[18];
    t[5]  = p[4] ^ t[3];
    t[6]  = p[15] ^ t[2];
    t[7]  = p[20] ^ t[5];
    t[8]  = p[5] ^ p[14];
    t[9]  = p[2] ^ p[17];
    t[10] = p[5] ^ p[8];
    t[11] = p[21] ^ t[6];
    t[12] = p[16] ^ t[4];

    __v2di cc[13];

    cc[0]  = p[13];
    cc[2]  = p[7] ^ p[15] ^ t[1];
    cc[4]  = p[3] ^ t[1] ^ t[3] ^ t[8];
    cc[6]  = p[0] ^ p[12] ^ p[13] ^ t[7] ^ t[11];
    cc[8]  = p[7] ^ t[0] ^ t[2] ^ t[10];
    cc[10] = p[3] ^ p[4] ^ t[0];
    cc[12] = p[0];

    cc[1]  = p[18] ^ t[1];
    cc[3]  = t[7] ^ t[9] ^ t[12];
    cc[5]  = p[2] ^ p[11] ^ p[19] ^ t[1] ^ t[10] ^ t[11];
    cc[7]  = p[21] ^ t[0] ^ t[5] ^ t[8] ^ t[12];
    cc[9]  = p[12] ^ p[19] ^ t[4] ^ t[6] ^ t[9];
    cc[11] = p[2] ^ t[0];

    _mm_storeu_si128((__v2di*)(c),    cc[0]                           ^ _mm_slli_si128(cc[1], 8));
    _mm_storeu_si128((__v2di*)(c+2),  cc[2]  ^ _mm_srli_si128(cc[1], 8) ^ _mm_slli_si128(cc[3], 8));
    _mm_storeu_si128((__v2di*)(c+4),  cc[4]  ^ _mm_srli_si128(cc[3], 8) ^ _mm_slli_si128(cc[5], 8));
    _mm_storeu_si128((__v2di*)(c+6),  cc[6]  ^ _mm_srli_si128(cc[5], 8) ^ _mm_slli_si128(cc[7], 8));
    _mm_storeu_si128((__v2di*)(c+8),  cc[8]  ^ _mm_srli_si128(cc[7], 8) ^ _mm_slli_si128(cc[9], 8));
    _mm_storeu_si128((__v2di*)(c+10), cc[10] ^ _mm_srli_si128(cc[9], 8) ^ _mm_slli_si128(cc[11], 8));
    _mm_storeu_si128((__v2di*)(c+12), cc[12] ^ _mm_srli_si128(cc[11], 8));
}

#endif  /* GF2X_MUL7_H_ */
