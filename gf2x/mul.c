/* Main file for Karatsuba and Toom-Cook multiplication over GF(2)[x].

  Copyright 2007 Richard Brent & Paul Zimmermann.

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

// #define BOGONG			// Selects routines that are good on
					// bogong, a 2.2 Ghz AMD Opteron 275
#ifndef MUL_KARA_THRESHOLD  
#define MUL_KARA_THRESHOLD 10
#endif

#ifndef MUL_TOOM_THRESHOLD   
#define MUL_TOOM_THRESHOLD 17
#endif

#ifndef MUL_TOOMW_THRESHOLD   
#define MUL_TOOMW_THRESHOLD 10
#endif

#ifndef MUL_TOOM4_THRESHOLD
#define MUL_TOOM4_THRESHOLD 30
#endif

#if (NTL_BITS_PER_LONG == 64)

/* use SSE2-based mul2() */
#define HAVE_NATIVE_MUL2
#ifndef TUNE_MUL1
#ifdef BOGONG
#include "mul2rpb.c"
#else
#include "mul2t.c"
#endif /* BOGONG */
#else
#define TUNING
#endif /* TUNE_MUL1 */

#else

/* this version of mul2 is faster than mul2t on bogong,
   but slower than mul2rpb. Thus not used at present. RPB, 20070527 */
static void mul2(_ntl_ulong *c, const _ntl_ulong *a, const _ntl_ulong *b)
{
   _ntl_ulong t;
   _ntl_ulong u[2];
   
   mul1 (c, a[0], b[0]);
   mul1 (c+2, a[1], b[1]);
   t    = c[1]^c[2];
   mul1 (u, a[0]^a[1], b[0]^b[1]);
   c[1] = c[0]^u[0]^t;
   c[2] = c[3]^u[1]^t;
}
#endif /* NTL_BITS_PER_LONG == 64 */

/* specialized Karatsuba */

static void mul3kara (_ntl_ulong *c, const _ntl_ulong *a, const _ntl_ulong *b)
{
  _ntl_ulong aa[2], bb[2], ab[4], c24, c35;
  mul1 (c + 4, a[2], b[2]);
  mul2 (c, a, b);
  aa[0] = a[0] ^ a[2];
  aa[1] = a[1];
  bb[0] = b[0] ^ b[2];
  bb[1] = b[1];
  c24 = c[2] ^ c[4];
  c35 = c[3] ^ c[5];
  mul2 (ab, aa, bb);
  c[2] = ab[0] ^ c[0] ^ c24;
  c[3] = ab[1] ^ c[1] ^ c35;
  c[4] = ab[2] ^ c24;
  c[5] = ab[3] ^ c35;
}

/* specialized Karatsuba */

static void mul4kara (_ntl_ulong *c, const _ntl_ulong *a, const _ntl_ulong *b)
{
  _ntl_ulong aa[2], bb[2], ab[4], c24, c35;
  mul2 (c, a, b);
  mul2 (c + 4, a + 2, b + 2);
  aa[0] = a[0] ^ a[2];
  aa[1] = a[1] ^ a[3];
  bb[0] = b[0] ^ b[2];
  bb[1] = b[1] ^ b[3];
  c24 = c[2] ^ c[4];
  c35 = c[3] ^ c[5];
  mul2 (ab, aa, bb);
  c[2] = ab[0] ^ c[0] ^ c24;
  c[3] = ab[1] ^ c[1] ^ c35;
  c[4] = ab[2] ^ c[6] ^ c24;
  c[5] = ab[3] ^ c[7] ^ c35;
}

/* specialized Karatsuba, RPB 20070518 */

static void mul5kara (_ntl_ulong *c, const _ntl_ulong *a, const _ntl_ulong *b)
{
#ifdef HAVE_NATIVE_MUL2
  /* adapted from mul6kara below, with a[5]=b[5]=c[10]=c[11]=0 */
  _ntl_ulong d01[4], d1[4], d12[4], aa[2], bb[2];
  mul2 (c,         a,     b); /* D0 */
  mul2 (d1,    a + 2, b + 2); /* D1 */
  mul1 (c + 8, a[4], b[4]);   /* D2 has only two words */
  aa[0] = a[0] ^ a[2]; aa[1] = a[1] ^ a[3];
  bb[0] = b[0] ^ b[2]; bb[1] = b[1] ^ b[3];
  mul2 (d01, aa, bb);         /* D01 */
  aa[0] = a[0] ^ a[4]; aa[1] = a[1];
  bb[0] = b[0] ^ b[4]; bb[1] = b[1];
  mul2 (c + 4, aa, bb);       /* D02 */
  aa[0] = a[2] ^ a[4]; aa[1] = a[3];
  bb[0] = b[2] ^ b[4]; bb[1] = b[3];
  mul2 (d12, aa, bb);         /* D12 */
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
  _ntl_ulong aa[3], bb[3], ab[6], ab3, ab4, ab5;
  mul2 (c+6, a+3, b+3);
  mul3kara (c, a, b);
  aa[0] = a[0] ^ a[3];
  aa[1] = a[1] ^ a[4];
  aa[2] = a[2];
  bb[0] = b[0] ^ b[3];
  bb[1] = b[1] ^ b[4];
  bb[2] = b[2];
  mul3kara (ab, aa, bb);
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
static void mul6kara (_ntl_ulong *c, const _ntl_ulong *a, const _ntl_ulong *b)
{
#ifdef HAVE_NATIVE_MUL2
  /* This code uses the K3 formula from Weimerskirch and Paar,
     http://weimerskirch.org/papers/Weimerskirch_Karatsuba.pdf,
     which performs only 6 calls to mul2. */
  _ntl_ulong d01[4], d1[4], d12[4], aa[2], bb[2];
  mul2 (c,         a,     b); /* D0 */
  mul2 (d1,    a + 2, b + 2); /* D1 */
  mul2 (c + 8, a + 4, b + 4); /* D2 */
  aa[0] = a[0] ^ a[2]; aa[1] = a[1] ^ a[3];
  bb[0] = b[0] ^ b[2]; bb[1] = b[1] ^ b[3];
  mul2 (d01, aa, bb);         /* D01 */
  aa[0] = a[0] ^ a[4]; aa[1] = a[1] ^ a[5];
  bb[0] = b[0] ^ b[4]; bb[1] = b[1] ^ b[5];
  mul2 (c + 4, aa, bb);       /* D02 */
  aa[0] = a[2] ^ a[4]; aa[1] = a[3] ^ a[5];
  bb[0] = b[2] ^ b[4]; bb[1] = b[3] ^ b[5];
  mul2 (d12, aa, bb);         /* D12 */
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
  _ntl_ulong aa[3], bb[3], ab[6], ab3, ab4, ab5;
  mul3kara (c+6, a+3, b+3);
  mul3kara (c, a, b);
  aa[0] = a[0] ^ a[3];
  aa[1] = a[1] ^ a[4];
  aa[2] = a[2] ^ a[5];
  bb[0] = b[0] ^ b[3];
  bb[1] = b[1] ^ b[4];
  bb[2] = b[2] ^ b[5];
  mul3kara (ab, aa, bb);
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

static void mul7kara (_ntl_ulong *c, const _ntl_ulong *a, const _ntl_ulong *b)
{
  _ntl_ulong aa[4], bb[4], ab[8], ab4, ab5, ab6, ab7;
  mul3kara (c+8, a+4, b+4);
  mul4kara (c, a, b);
  aa[0] = a[0] ^ a[4];
  aa[1] = a[1] ^ a[5];
  aa[2] = a[2] ^ a[6];
  aa[3] = a[3];
  bb[0] = b[0] ^ b[4];
  bb[1] = b[1] ^ b[5];
  bb[2] = b[2] ^ b[6];
  bb[3] = b[3];
  mul4kara (ab, aa, bb);
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

static void mul8kara (_ntl_ulong *c, const _ntl_ulong *a, const _ntl_ulong *b)
{
  _ntl_ulong aa[4], bb[4], cc[4];
  long j;
  mul4kara (c+8, a+4, b+4);
  mul4kara (c, a, b);
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
  mul4kara (c+4, aa, bb);
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

static void mul9kara (_ntl_ulong *c, const _ntl_ulong *a, const _ntl_ulong *b)
{
  _ntl_ulong aa[5], bb[5], ab[10], ab5, ab6, ab7, ab8, ab9;
  mul4kara (c+10, a+5, b+5);
  mul5kara (c, a, b);
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
  mul5kara (ab, aa, bb);
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

/* an = bn = n > 0 */

static inline void mul_basecase_n (_ntl_ulong *c, const _ntl_ulong *a, 
				    const _ntl_ulong *b, long n)
{
  long i;

  switch (n)
    {
    case 0:				// This can occur in call from KarMul
      return;
    case 1:
      return mul1 (c, a[0], b[0]);
#ifndef TUNE_MUL1
    case 2:
      return mul2 (c, a, b);
    case 3:
      return mul3kara (c, a, b);
    case 4:
      return mul4kara (c, a, b);
    case 5:
      return mul5kara (c, a, b);
    case 6:
      return mul6kara (c, a, b);
    case 7:
      return mul7kara (c, a, b);
    case 8:
      return mul8kara (c, a, b);
    case 9:
      return mul9kara (c, a, b);
#endif
    default:
      printf ("mul.c: ran off end of switch, n %ld\n", n);
      if (MUL_KARA_THRESHOLD > 10)
        printf ("       MUL_KARA_THRESHOLD %ld > 10, reduce to 10\n",
                MUL_KARA_THRESHOLD);
      exit(1);
      
//    Old code for default case, should no longer be necessary
//                
//    c[n] = Mul1 (c, a, n, b[0]);
//    for (i = 1; i < n; i++)
//      c[n + i] = AddMul1 (c + i, c + i, a, n, b[i]);
//    return;

    }
}

/* Version of Karatsuba multiplication with minimal temporary storage
   sp(n) = 3*ceil(n/2) + sp(ceil(n/2)) = 3n + O(log n) words. 
   RPB, 20070522 */
   
static void KarMul (_ntl_ulong *c, const _ntl_ulong *a, const _ntl_ulong *b, 
            	      long n, _ntl_ulong *stk)
{
  _ntl_ulong t;
  _ntl_ulong *aa, *bb, *cc;
  long j, d, n2;

#if 0
  if (n <= 0)	/* if turned on this test shows that calls with n == 0 */
    		/* do occur (e.g from tunefft, FFT(19683)), but don't  */
                /* seem to be harmful if mul_basecase_n just returns.  */
    {            
    printf ("\nWarning: n %ld in call to KarMul\n", n);
    fflush (stdout);
    }
#endif    

  if (n < MUL_KARA_THRESHOLD)
    return mul_basecase_n (c, a, b, n);

  n2 = (n+1)/2;		/* ceil(n/2) */
  d  = n&1;		/* 2*n2 - n = 1 if n odd, 0 if n even */
  aa = stk;		/* Size n2   */
  bb = aa + n2;		/* Size n2   */
  cc = bb + n2;		/* Size n2   */

  stk = cc + n2;	/* sp(n) = 3*ceil(n/2)) + sp(ceil(n/2)) */

  const _ntl_ulong *a1 = a + n2;	/* a[n2] */
  const _ntl_ulong *b1 = b + n2;	/* b[n2] */
  	_ntl_ulong *c1 = c + n2;	/* c[n2]   */
  	_ntl_ulong *c2 = c1 + n2;	/* c[2*n2] */
  	_ntl_ulong *c3 = c2 + n2;	/* c[3*n2] */

  KarMul (c, a, b, n2, stk);		/* Low */

  KarMul (c2, a1, b1, n2-d, stk);	/* High */

  for (j = 0; j < n2-d; j++)
    {
    aa[j] = a[j]  ^ a1[j];
    bb[j] = b[j]  ^ b1[j];
    cc[j] = c1[j] ^ c2[j];
    }
  for (; j < n2; j++)			/* Only when n odd */
    {
    aa[j] = a[j];
    bb[j] = b[j];
    cc[j] = c1[j] ^ c2[j];
    }

  KarMul (c1, aa, bb, n2, stk);		/* Middle */

  for (j = 0; j < n2-2*d; j++)
    {
    t = cc[j];
    c1[j] ^= t ^ c[j];
    c2[j] ^= t ^ c3[j];
    }
  for (; j < n2; j++)			/* Only when n odd */
    {
    c1[j] ^= cc[j] ^ c[j];
    c2[j] ^= cc[j];
    }
}

#ifndef TUNING
#include "besttoom.h"
#endif

#ifndef TUNE_MUL1
#include "BestToom.c" 
#include "ToomSpace.c"

#include "Toom.c"
#include "mul-tc3.c"  
#include "mul-tc3w.c"
#include "mul-tc3u.c" 
#include "mul-tc4.c"  
#endif
