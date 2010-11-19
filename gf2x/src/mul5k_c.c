/* This file is part of the gf2x library.

   Copyright 2007, 2008, 2009, 2010
   Richard Brent, Pierrick Gaudry, Emmanuel Thome', Paul Zimmermann,
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
#ifndef GF2X_MUL5_H_
#define GF2X_MUL5_H_

#include "gf2x.h"
/* All gf2x source files for lowlevel functions must include gf2x-small.h
 * This is mandatory for the tuning mechanism. */
#include "gf2x/gf2x-small.h"

GF2X_STORAGE_CLASS_mul5
void gf2x_mul5 (unsigned long *c, const unsigned long *a,
        const unsigned long *b)
{
  /* Montgomery formulae with 13 multiplications, see
     Five, Six, and Seven-Term {K}aratsuba-Like Formulae,
     IEEE Transactions on Computers, volume 54, number 3, p. 362-369, 2005 */
  unsigned long ta[3], tb[3], pa[8], pb[8], p[26], t[14];
  ta[0] = a[0]  ^ a[4]         ; tb[0] = b[0]  ^ b[4];
  ta[1] = a[1]  ^ a[2]         ; tb[1] = b[1]  ^ b[2];
  ta[2] = a[3]  ^ ta[0]        ; tb[2] = b[3]  ^ tb[0];
  pa[0] = ta[1] ^ ta[2]        ; pb[0] = tb[1] ^ tb[2];
  pa[1] = a[2]  ^ ta[2]        ; pb[1] = b[2]  ^ tb[2];
  pa[2] = ta[0] ^ ta[1]        ; pb[2] = tb[0] ^ tb[1];
  pa[3] = a[1]  ^ ta[2]        ; pb[3] = b[1]  ^ tb[2];
  pa[4] = a[0]  ^ a[2]  ^ a[3] ; pb[4] = b[0]  ^ b[2]  ^ b[3];
  pa[5] = a[4]  ^ ta[1]        ; pb[5] = b[4]  ^ tb[1];
  pa[6] = a[3]  ^ a[4]         ; pb[6] = b[3]  ^ b[4];
  pa[7] = a[0]  ^ a[1]         ; pb[7] = b[0]  ^ b[1];
  gf2x_mul1(p + 0,  pa[0], pb[0]);
  gf2x_mul1(p + 2,  pa[1], pb[1]);
  gf2x_mul1(p + 4,  pa[2], pb[2]);
  gf2x_mul1(p + 6,  pa[3], pb[3]);
  gf2x_mul1(p + 8,  pa[4], pb[4]);
  gf2x_mul1(p + 10, pa[5], pb[5]);
  gf2x_mul1(p + 12, pa[6], pb[6]);
  gf2x_mul1(p + 14, pa[7], pb[7]);
  gf2x_mul1(p + 16, ta[0], tb[0]);
  gf2x_mul1(p + 18, a[4],  b[4]);
  gf2x_mul1(p + 20, a[3],  b[3]);
  gf2x_mul1(p + 22, a[1],  b[1]);
  gf2x_mul1(p + 24, a[0],  b[0]);
  t[0]  = p[14] ^ p[24];   t[1]  = p[15] ^ p[25]; 
  t[2]  = p[12] ^ p[18];   t[3]  = p[13] ^ p[19];
  t[4]  = p[2]  ^ p[16];   t[5]  = p[3]  ^ p[17];
  t[6]  = p[0]  ^ p[6];    t[7]  = p[1]  ^ p[7];
  t[8]  = p[4]  ^ p[16];   t[9]  = p[5]  ^ p[17];
  t[10] = p[10] ^ t[0];    t[11] = p[11] ^ t[1];
  t[12] = p[8]  ^ t[2];    t[13] = p[9]  ^ t[3];
  c[0] = p[24];
  c[1] = p[22] ^ t[0]                          ^ p[25];
  c[2] = p[18] ^ t[8]  ^ t[10]                 ^ p[23] ^ t[1]; 
  c[3] = t[2]  ^ t[4]  ^ t[6]                  ^ p[19] ^ t[9]  ^ t[11];
  c[4] = p[0]  ^ p[20] ^ p[22] ^ t[10] ^ t[12] ^ t[3]  ^ t[5]  ^ t[7];               
  c[5] = t[0]  ^ t[6]  ^ t[8]                  ^ p[1]  ^ p[21] ^ p[23] ^ t[11] ^ t[13];
  c[6] = p[24] ^ t[4]  ^ t[12]                 ^ t[1]  ^ t[7]  ^ t[9];
  c[7] = p[20] ^ t[2]                          ^ p[25] ^ t[5]  ^ t[13];                
  c[8] = p[18]                                 ^ p[21] ^ t[3];               
  c[9] =                                        p[19];                        
}


#endif  /* GF2X_MUL5_H_ */
