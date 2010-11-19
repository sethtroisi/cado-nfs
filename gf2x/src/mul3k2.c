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

#ifndef GF2X_MUL3_H_
#define GF2X_MUL3_H_

#include "gf2x.h"
/* All gf2x source files for lowlevel functions must include gf2x-small.h
 * This is mandatory for the tuning mechanism. */
#include "gf2x/gf2x-small.h"

/* uses the variant of Karatsuba with 6 multiplications */
GF2X_STORAGE_CLASS_mul3
void gf2x_mul3 (unsigned long *c, const unsigned long *a, const unsigned long *b)
{
  unsigned long aa[3], bb[3];
  unsigned long p0[2], p1[2], p2[2];
  unsigned long pp0[2], pp1[2], pp2[2];
  aa[0] = a[1]^a[2];
  aa[1] = a[0]^a[2];
  aa[2] = a[0]^a[1];
  bb[0] = b[1]^b[2];
  bb[1] = b[0]^b[2];
  bb[2] = b[0]^b[1];
  gf2x_mul1 (p0, a[0], b[0]);
  gf2x_mul1 (p1, a[1], b[1]);
  gf2x_mul1 (p2, a[2], b[2]);
  gf2x_mul1 (pp0, aa[0], bb[0]);
  gf2x_mul1 (pp1, aa[1], bb[1]);
  gf2x_mul1 (pp2, aa[2], bb[2]);
  c[0] = p0[0];
  c[1] = p0[0]^p1[0]^pp2[0]         ^  p0[1];
  c[2] = p0[0]^p1[0]^p2[0]^pp1[0]   ^  p0[1]^p1[1]^pp2[1];
  c[3] = pp0[0]^p1[0]^p2[0]         ^  p0[1]^p1[1]^p2[1]^pp1[1];
  c[4] = p2[0]                      ^  pp0[1]^p1[1]^p2[1];
  c[5] =                               p2[1];
}

#endif  /* GF2X_MUL3_H_ */
