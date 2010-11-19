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
#ifndef GF2X_MUL5_H_
#define GF2X_MUL5_H_

#include "gf2x.h"
/* All gf2x source files for lowlevel functions must include gf2x-small.h
 * This is mandatory for the tuning mechanism. */
#include "gf2x/gf2x-small.h"

GF2X_STORAGE_CLASS_mul5
void gf2x_mul5 (unsigned long *c, const unsigned long *a, const unsigned long *b)
{
    /* adapted from gf2x_mul6, with a[5]=b[5]=c[10]=c[11]=0 */
    unsigned long d01[4], d1[4], d12[4], aa[2], bb[2];
    gf2x_mul2 (c,         a,     b); /* D0 */
    gf2x_mul2 (d1,    a + 2, b + 2); /* D1 */
    gf2x_mul1 (c + 8, a[4], b[4]);   /* D2 has only two words */
    aa[0] = a[0] ^ a[2]; aa[1] = a[1] ^ a[3];
    bb[0] = b[0] ^ b[2]; bb[1] = b[1] ^ b[3];
    gf2x_mul2 (d01, aa, bb);         /* D01 */
    aa[0] = a[0] ^ a[4]; aa[1] = a[1];
    bb[0] = b[0] ^ b[4]; bb[1] = b[1];
    gf2x_mul2 (c + 4, aa, bb);       /* D02 */
    aa[0] = a[2] ^ a[4]; aa[1] = a[3];
    bb[0] = b[2] ^ b[4]; bb[1] = b[3];
    gf2x_mul2 (d12, aa, bb);         /* D12 */
    /* low(D1) + high(D0) is used three times */
    c[2] ^= d1[0]; c[3] ^= d1[1];   /* low(D1) + high(D0) */
    c[4] ^= c[2];  c[5] ^= c[3];    /* low(D02) + low(D1) + high(D0) */
    d12[0] ^= c[2]; d12[1] ^= c[3]; /* low(D12) + low(D1) + high(D0) */
    /* low(D2) + high(D1) is used three times */
    c[8] ^= d1[2]; c[9] ^= d1[3];   /* low(D2) + high(D1) */
    c[6] ^= c[8];  c[7] ^= c[9];    /* high(D02) + low(D2) + high(D1) */
    d01[2] ^= c[8]; d01[3] ^= c[9]; /* high(D01) + low(D2) + high(D1) */
    c[2] ^= d01[0] ^ c[0]; c[3] ^= d01[1] ^ c[1]; /* l(D1)+h(D0)+l(D01)+l(D0) */
    c[4] ^= c[0] ^ d01[2]; c[5] ^= c[1] ^ d01[3];
    c[6] ^= d12[0]; c[7] ^= d12[1];
    c[8] ^= d12[2]; c[9] ^= d12[3];
}

#endif  /* GF2X_MUL5_H_ */
