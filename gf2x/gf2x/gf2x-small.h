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

/* Main file for Karatsuba and Toom-Cook multiplication over GF(2)[x]. */

#ifndef GF2X_SMALL_H_
#define GF2X_SMALL_H_

#include "gf2x.h"

#include "gf2x/gf2x-thresholds.h"

/* functions here will end up as static functions, therefore we prefer to
 * avoid the warning relative to the fact that they are unused. */

#ifndef	MAYBE_UNUSED
#if defined(__GNUC__)
#define MAYBE_UNUSED __attribute__ ((unused))
#else
#define MAYBE_UNUSED
#endif
#endif

#ifdef __cplusplus
extern "C" {
#endif

GF2X_STORAGE_CLASS_mul1 void
gf2x_mul1(unsigned long *c, unsigned long a, unsigned long b)
        MAYBE_UNUSED;
GF2X_STORAGE_CLASS_mul_1_n unsigned long
gf2x_mul_1_n(unsigned long *cp, const unsigned long *bp, long sb, unsigned long a)
        MAYBE_UNUSED;
GF2X_STORAGE_CLASS_addmul_1_n unsigned long
gf2x_addmul_1_n(unsigned long *dp,
        const unsigned long *cp, const unsigned long* bp,
        long sb, unsigned long a)
        MAYBE_UNUSED;
GF2X_STORAGE_CLASS_mul2 void
gf2x_mul2(unsigned long *c, const unsigned long *a, const unsigned long *b)
        MAYBE_UNUSED;

GF2X_STORAGE_CLASS_mul3 void
gf2x_mul3(unsigned long *c, const unsigned long *a, const unsigned long *b)
        MAYBE_UNUSED;
GF2X_STORAGE_CLASS_mul4 void
gf2x_mul4(unsigned long *c, const unsigned long *a, const unsigned long *b)
        MAYBE_UNUSED;

/* larger sizes are never tuned so far, nor do we consider it worthwhile
 * to inline them. So no GF2X_STORAGE_CLASS_* for them.
 */
static void
gf2x_mul5(unsigned long *c, const unsigned long *a, const unsigned long *b)
        MAYBE_UNUSED;
static void
gf2x_mul6(unsigned long *c, const unsigned long *a, const unsigned long *b)
        MAYBE_UNUSED;
static void
gf2x_mul7(unsigned long *c, const unsigned long *a, const unsigned long *b)
        MAYBE_UNUSED;
static void
gf2x_mul8(unsigned long *c, const unsigned long *a, const unsigned long *b)
        MAYBE_UNUSED;
static void
gf2x_mul9(unsigned long *c, const unsigned long *a, const unsigned long *b)
        MAYBE_UNUSED;
#ifdef __cplusplus
}
#endif

/* This file provides all the small-sized gf2x_mul1..gf2x_mul9 routines. It is
 * meant to be possibly included directly by applications. */

#include "gf2x/gf2x_mul1.h"
#include "gf2x/gf2x_mul2.h"
#include "gf2x/gf2x_mul3.h"
#include "gf2x/gf2x_mul4.h"

/* specialized Karatsuba, RPB 20070518 */

static void gf2x_mul5 (unsigned long *c, const unsigned long *a, const unsigned long *b)
{
#ifdef HAVE_NATIVE_gf2x_MUL2
  /* adapted from gf2x_mul6 below, with a[5]=b[5]=c[10]=c[11]=0 */
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
#else
  unsigned long aa[3], bb[3], ab[6], ab3, ab4, ab5;
  gf2x_mul2 (c+6, a+3, b+3);
  gf2x_mul3 (c, a, b);
  aa[0] = a[0] ^ a[3];
  aa[1] = a[1] ^ a[4];
  aa[2] = a[2];
  bb[0] = b[0] ^ b[3];
  bb[1] = b[1] ^ b[4];
  bb[2] = b[2];
  gf2x_mul3 (ab, aa, bb);
  ab3 = ab[3] ^ c[3];
  ab4 = ab[4] ^ c[4];
  ab5 = ab[5] ^ c[5];
  c[3] ^= ab[0] ^ c[0] ^ c[6];
  c[4] ^= ab[1] ^ c[1] ^ c[7];
  c[5] ^= ab[2] ^ c[2] ^ c[8];
  c[6] ^= ab3 ^ c[9];
  c[7] ^= ab4;
  c[8] ^= ab5;
#endif
}

/* specialized Karatsuba, RPB 20070518 */
static void gf2x_mul6 (unsigned long *c, const unsigned long *a, const unsigned long *b)
{
#ifdef HAVE_NATIVE_gf2x_MUL2
  /* This code uses the K3 formula from Weimerskirch and Paar,
     http://weimerskirch.org/papers/Weimerskirch_Karatsuba.pdf,
     which performs only 6 calls to gf2x_mul2. */
  unsigned long d01[4], d1[4], d12[4], aa[2], bb[2];
  gf2x_mul2 (c,         a,     b); /* D0 */
  gf2x_mul2 (d1,    a + 2, b + 2); /* D1 */
  gf2x_mul2 (c + 8, a + 4, b + 4); /* D2 */
  aa[0] = a[0] ^ a[2]; aa[1] = a[1] ^ a[3];
  bb[0] = b[0] ^ b[2]; bb[1] = b[1] ^ b[3];
  gf2x_mul2 (d01, aa, bb);         /* D01 */
  aa[0] = a[0] ^ a[4]; aa[1] = a[1] ^ a[5];
  bb[0] = b[0] ^ b[4]; bb[1] = b[1] ^ b[5];
  gf2x_mul2 (c + 4, aa, bb);       /* D02 */
  aa[0] = a[2] ^ a[4]; aa[1] = a[3] ^ a[5];
  bb[0] = b[2] ^ b[4]; bb[1] = b[3] ^ b[5];
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
  c[6] ^= d12[0] ^ c[10]; c[7] ^= d12[1] ^ c[11];
  c[8] ^= d12[2] ^ c[10]; c[9] ^= d12[3] ^ c[11];
#else
  unsigned long aa[3], bb[3], ab[6], ab3, ab4, ab5;
  gf2x_mul3 (c+6, a+3, b+3);
  gf2x_mul3 (c, a, b);
  aa[0] = a[0] ^ a[3];
  aa[1] = a[1] ^ a[4];
  aa[2] = a[2] ^ a[5];
  bb[0] = b[0] ^ b[3];
  bb[1] = b[1] ^ b[4];
  bb[2] = b[2] ^ b[5];
  gf2x_mul3 (ab, aa, bb);
  ab3 = ab[3] ^ c[3];
  ab4 = ab[4] ^ c[4];
  ab5 = ab[5] ^ c[5];
  c[3] ^= ab[0] ^ c[0] ^ c[6];
  c[4] ^= ab[1] ^ c[1] ^ c[7];
  c[5] ^= ab[2] ^ c[2] ^ c[8];
  c[6] ^= ab3 ^ c[9];
  c[7] ^= ab4 ^ c[10];
  c[8] ^= ab5 ^ c[11];
#endif
}

/* specialized Karatsuba, RPB 20070518 */
/* slightly faster on bogong than version with loops */

static void gf2x_mul7 (unsigned long *c, const unsigned long *a, const unsigned long *b)
{
  unsigned long aa[4], bb[4], ab[8], ab4, ab5, ab6, ab7;
  gf2x_mul3 (c+8, a+4, b+4);
  gf2x_mul4 (c, a, b);
  aa[0] = a[0] ^ a[4];
  aa[1] = a[1] ^ a[5];
  aa[2] = a[2] ^ a[6];
  aa[3] = a[3];
  bb[0] = b[0] ^ b[4];
  bb[1] = b[1] ^ b[5];
  bb[2] = b[2] ^ b[6];
  bb[3] = b[3];
  gf2x_mul4 (ab, aa, bb);
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

/* specialized Karatsuba, RPB 20070518 */
/* slightly faster on bogong than version with loops */
/* this version uses minimal temporary storage (12 = 3*n/2 words) */

static void gf2x_mul8 (unsigned long *c, const unsigned long *a, const unsigned long *b)
{
  unsigned long aa[4], bb[4], cc[4];
  gf2x_mul4 (c+8, a+4, b+4);
  gf2x_mul4 (c, a, b);
  cc[0] = c[4] ^ c[8];
  cc[1] = c[5] ^ c[9];
  cc[2] = c[6] ^ c[10];
  cc[3] = c[7] ^ c[11];
  aa[0] = a[0] ^ a[4];
  aa[1] = a[1] ^ a[5];
  aa[2] = a[2] ^ a[6];
  aa[3] = a[3] ^ a[7];
  bb[0] = b[0] ^ b[4];
  bb[1] = b[1] ^ b[5];
  bb[2] = b[2] ^ b[6];
  bb[3] = b[3] ^ b[7];
  gf2x_mul4 (c+4, aa, bb);
  c[4]  ^= c[0]  ^ cc[0];
  c[5]  ^= c[1]  ^ cc[1];
  c[6]  ^= c[2]  ^ cc[2];
  c[7]  ^= c[3]  ^ cc[3];
  c[8]  ^= c[12] ^ cc[0];
  c[9]  ^= c[13] ^ cc[1];
  c[10] ^= c[14] ^ cc[2];
  c[11] ^= c[15] ^ cc[3];
}

/* specialized Karatsuba, RPB 20070520 */
/* slightly faster on bogong than version with loops */

static void gf2x_mul9 (unsigned long *c, const unsigned long *a, const unsigned long *b)
{
  unsigned long aa[5], bb[5], ab[10], ab5, ab6, ab7, ab8, ab9;
  gf2x_mul4 (c+10, a+5, b+5);
  gf2x_mul5 (c, a, b);
  aa[0] = a[0] ^ a[5];
  aa[1] = a[1] ^ a[6];
  aa[2] = a[2] ^ a[7];
  aa[3] = a[3] ^ a[8];
  aa[4] = a[4];
  bb[0] = b[0] ^ b[5];
  bb[1] = b[1] ^ b[6];
  bb[2] = b[2] ^ b[7];
  bb[3] = b[3] ^ b[8];
  bb[4] = b[4];
  gf2x_mul5 (ab, aa, bb);
  ab5 = ab[5] ^ c[5];
  ab6 = ab[6] ^ c[6];
  ab7 = ab[7] ^ c[7];
  ab8 = ab[8] ^ c[8];
  ab9 = ab[9] ^ c[9];
  c[5] ^= ab[0] ^ c[0] ^ c[10];
  c[6] ^= ab[1] ^ c[1] ^ c[11];
  c[7] ^= ab[2] ^ c[2] ^ c[12];
  c[8] ^= ab[3] ^ c[3] ^ c[13];
  c[9] ^= ab[4] ^ c[4] ^ c[14];
  c[10] ^= ab5 ^ c[15];
  c[11] ^= ab6 ^ c[16];
  c[12] ^= ab7 ^ c[17];
  c[13] ^= ab8;
  c[14] ^= ab9;
}

#endif  /* GF2X_SMALL_H_ */
