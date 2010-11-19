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

#include <wmmintrin.h>

#if GF2X_WORDSIZE != 64
#error "This code is for 64-bit only"
#endif

#include "gf2x/gf2x-config.h"

#ifndef HAVE_PCLMUL_SUPPORT
#error "This code needs pclmul support"
#endif

#undef GF2X_MUL2_H_
#define CARRY
#include "mul2cl.c"
#undef CARRY

#undef GF2X_MUL2_H_
#define BORROW
#include "mul2cl.c"
#undef BORROW

/* uses Montgomery's variant of Karatsuba for n=2k+1 odd,
   with M(2k+1) = M(k) + 2M(k+1) - 1, see
   Five, Six, and Seven-Term {K}aratsuba-Like Formulae,
   IEEE Transactions on Computers, volume 54, number 3, pages 362-369, 2005 */
GF2X_STORAGE_CLASS_mul3 void
gf2x_mul3 (unsigned long *c, const unsigned long *a, const unsigned long *b)
{
  unsigned long d[2], aa[2], bb[2], p[4];

  /* let A0 = {a, 2}, A1 = {a+2, 1}, B0 = {b, 2}, B1 = {b+2, 1} */

  aa[0] = a[0] ^ a[2];
  aa[1] = a[1];
  bb[0] = b[0] ^ b[2];
  bb[1] = b[1];

  gf2x_mul2c (c, a, b, d);
  /* {c, 4} = A0 * B0 and {d, 2} = {a+1, 1} * {b+1, 1} */

  gf2x_mul2b (p, aa, bb, d);
  /* {p, 4} = (A0 + A1) * (B0 + B1) */

  gf2x_mul1 (c + 4, a[2], b[2]);

  p[0] ^= c[0] ^ c[4];
  p[1] ^= c[1] ^ c[5];

  c[5] ^= p[3] ^ c[3];
  c[4] ^= p[2] ^ c[2];
  c[2] ^= p[0];
  c[3] ^= p[1];
}

#endif  /* GF2X_MUL3_H_ */
