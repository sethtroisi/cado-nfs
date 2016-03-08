/* This file is part of the gf2x library.

   Copyright 2007, 2008, 2009, 2010, 2013, 2015
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

#ifndef GF2X_MUL5_H_
#define GF2X_MUL5_H_

#include "gf2x.h"
/* All gf2x source files for lowlevel functions must include gf2x-small.h
 * This is mandatory for the tuning mechanism. */
#include "gf2x/gf2x-small.h"

    GF2X_STORAGE_CLASS_mul5
void gf2x_mul5 (unsigned long *c, const unsigned long *a, const unsigned long *b)
{
  unsigned long aa[6], bb[6];
  unsigned long p0[4], p1[4], p2[2];
  unsigned long pp0[4], pp1[4], pp2[4];
  aa[0] = a[2]^a[4];
  aa[1] = a[3];
  aa[2] = a[0]^a[4];
  aa[3] = a[1];
  aa[4] = a[0]^a[2];
  aa[5] = a[1]^a[3];
  bb[0] = b[2]^b[4];
  bb[1] = b[3];
  bb[2] = b[0]^b[4];
  bb[3] = b[1];
  bb[4] = b[0]^b[2];
  bb[5] = b[1]^b[3];
  gf2x_mul2 (p0, a, b);
  gf2x_mul2 (p1, a+2, b+2);
  gf2x_mul1 (p2, a[4], b[4]);
  gf2x_mul2 (pp0, aa, bb);
  gf2x_mul2 (pp1, aa+2, bb+2);
  gf2x_mul2 (pp2, aa+4, bb+4);
  c[0]  = p0[0];
  c[1]  = p0[1];
  c[2]  = p0[0]^p1[0]^pp2[0]       ^ p0[2];
  c[3]  = p0[1]^p1[1]^pp2[1]       ^ p0[3];
  c[4]  = p0[0]^p1[0]^p2[0]^pp1[0] ^ p0[2]^p1[2]^pp2[2];
  c[5]  = p0[1]^p1[1]^p2[1]^pp1[1] ^ p0[3]^p1[3]^pp2[3];
  c[6]  = pp0[0]^p1[0]^p2[0]       ^ p0[2]^p1[2]^pp1[2];
  c[7]  = pp0[1]^p1[1]^p2[1]       ^ p0[3]^p1[3]^pp1[3];
  c[8]  = p2[0]                    ^ pp0[2]^p1[2];
  c[9]  = p2[1]                    ^ pp0[3]^p1[3];
}

#endif  /* GF2X_MUL5_H_ */
