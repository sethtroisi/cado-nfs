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
#ifndef GF2X_MUL9_H_
#define GF2X_MUL9_H_

#include "gf2x.h"
/* All gf2x source files for lowlevel functions must include gf2x-small.h
 * This is mandatory for the tuning mechanism. */
#include "gf2x/gf2x-small.h"

/* Karatsuba-3 applied to chunks of size 3 */
GF2X_STORAGE_CLASS_mul9
void gf2x_mul9 (unsigned long *c, const unsigned long *a, const unsigned long *b)
{
  unsigned long aa[9], bb[9];
  unsigned long p0[6], p1[6], p2[6];
  unsigned long pp0[6], pp1[6], pp2[6];
  aa[0] = a[3]^a[6];
  aa[1] = a[4]^a[7];
  aa[2] = a[5]^a[8];
  aa[3] = a[0]^a[6];
  aa[4] = a[1]^a[7];
  aa[5] = a[2]^a[8];
  aa[6] = a[0]^a[3];
  aa[7] = a[1]^a[4];
  aa[8] = a[2]^a[5];
  bb[0] = b[3]^b[6];
  bb[1] = b[4]^b[7];
  bb[2] = b[5]^b[8];
  bb[3] = b[0]^b[6];
  bb[4] = b[1]^b[7];
  bb[5] = b[2]^b[8];
  bb[6] = b[0]^b[3];
  bb[7] = b[1]^b[4];
  bb[8] = b[2]^b[5];
  gf2x_mul3 (p0, a+0, b+0);
  gf2x_mul3 (p1, a+3, b+3);
  gf2x_mul3 (p2, a+6, b+6);
  gf2x_mul3 (pp0, aa+0, bb+0);
  gf2x_mul3 (pp1, aa+3, bb+3);
  gf2x_mul3 (pp2, aa+6, bb+6);
  c[0]  = p0[0];
  c[1]  = p0[1];
  c[2]  = p0[2];
  c[3]  = p0[0]^p1[0]^pp2[0]       ^ p0[3];
  c[4]  = p0[1]^p1[1]^pp2[1]       ^ p0[4];
  c[5]  = p0[2]^p1[2]^pp2[2]       ^ p0[5];
  c[6]  = p0[0]^p1[0]^p2[0]^pp1[0] ^ p0[3]^p1[3]^pp2[3];
  c[7]  = p0[1]^p1[1]^p2[1]^pp1[1] ^ p0[4]^p1[4]^pp2[4];
  c[8]  = p0[2]^p1[2]^p2[2]^pp1[2] ^ p0[5]^p1[5]^pp2[5];
  c[9]  = pp0[0]^p1[0]^p2[0]       ^ p0[3]^p1[3]^p2[3]^pp1[3];
  c[10] = pp0[1]^p1[1]^p2[1]       ^ p0[4]^p1[4]^p2[4]^pp1[4];
  c[11] = pp0[2]^p1[2]^p2[2]       ^ p0[5]^p1[5]^p2[5]^pp1[5];
  c[12] = p2[0]                    ^ pp0[3]^p1[3]^p2[3];
  c[13] = p2[1]                    ^ pp0[4]^p1[4]^p2[4];
  c[14] = p2[2]                    ^ pp0[5]^p1[5]^p2[5];
  c[15] =                            p2[3];
  c[16] =                            p2[4];
  c[17] =                            p2[5];
}

#endif  /* GF2X_MUL9_H_ */
