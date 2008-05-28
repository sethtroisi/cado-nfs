/* Multiplication over GF(2)[x]

  Copyright 2007 Richard P. Brent, Pierrick Gaudry, Paul Zimmermann, Emmanuel Thome'

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
#ifndef GF2X_H_
#define GF2X_H_

#ifdef	__i386
#define	WORDSIZE	32
#endif
#ifdef	__x86_64
#define	WORDSIZE	64
#endif
#ifndef	WORDSIZE
#error "define	WORDSIZE"
#endif

/**********************************************************************/
#include <assert.h>
#ifndef ASSERT
#define ASSERT(x)	assert(x)
#endif

/*********************************************************************/
/* Helper macros */
#ifndef	MAYBE_UNUSED
#if defined(__GNUC__)
#define MAYBE_UNUSED __attribute__ ((unused))
#else
#define MAYBE_UNUSED
#endif
#endif

#include <stdlib.h>

#include "thresholds.h"

typedef unsigned long ulong;

#ifdef __cplusplus
extern "C" {
#endif

/* Declare prototypes for basic routines */
static inline void mul1 (ulong *c, ulong a, ulong b) MAYBE_UNUSED;
/*
static ulong mul_1_n (ulong *cp, const ulong *bp, long sb, ulong a);
static ulong addmul_1_n (ulong *dp, const ulong *cp, const ulong* bp, long sb,ulong a);
*/

static void mul2(ulong *c, const ulong *a, const ulong *b) MAYBE_UNUSED;

/* Declare prototypes that build directly on top of mul1 and mul2 */
static void mul3 (ulong *c, const ulong *a, const ulong *b) MAYBE_UNUSED;
static void mul4 (ulong *c, const ulong *a, const ulong *b) MAYBE_UNUSED;
static void mul5 (ulong *c, const ulong *a, const ulong *b) MAYBE_UNUSED;
static void mul6 (ulong *c, const ulong *a, const ulong *b) MAYBE_UNUSED;
static void mul7 (ulong *c, const ulong *a, const ulong *b) MAYBE_UNUSED;
static void mul8 (ulong *c, const ulong *a, const ulong *b) MAYBE_UNUSED;
static void mul9 (ulong *c, const ulong *a, const ulong *b) MAYBE_UNUSED;

static void mul_basecase(ulong * c, const ulong * a,
			 long na, const ulong * b, long nb) MAYBE_UNUSED;


void mul_toom(unsigned long *c, const unsigned long *a,
			const unsigned long *b, long n, unsigned long *stk);
void mul_kara(ulong *c, const ulong *a, const ulong *b,
			long n, ulong *stk);
void mul_tc3(ulong *c, const ulong *a, const ulong *b,
		 	long n, ulong *stk);
void mul_tc3w(ulong *c, const ulong *a, const ulong *b,
		        long n, ulong *stk);
void mul_tc4(ulong *c, const ulong *a, const ulong *b,
			long n, ulong *stk);
short BestToom(unsigned int);
long toomspace(long);

void mul_tc3u(ulong * c, const ulong * a, long sa,
	      const ulong * b, ulong * stk);
short BestuToom(unsigned int);
long toomuspace(long);


void FFTMul0(unsigned long *c, const unsigned long *a, long an,
	     const unsigned long *b, long bn, long K, long M);
void FFTMul(unsigned long *c, const unsigned long *a, long an,
		            const unsigned long *b, long bn, long K);
void FFTMul1(unsigned long *c, long cn,
	     const unsigned long *aa, long an,
	     const unsigned long *bb, long bn, long K, long M);
void FFTMul2(unsigned long *c, const unsigned long *a, long an,
	     const unsigned long *b, long bn, long K);


struct mul_gf2x_pool_s {
	ulong * stk;
	size_t stk_size;
};
typedef struct mul_gf2x_pool_s mul_gf2x_pool_t[1];
void mul_gf2x_pool_init(mul_gf2x_pool_t);
void mul_gf2x_pool_clear(mul_gf2x_pool_t);

/* This is the toplevel multiplication routine, the only one that really
 * matters.
 *
 * The destination must not overlap with the inputs.
 *
 * c must have enough room to hold an+bn words.
 */
void mul_gf2x(unsigned long *c,
		const unsigned long *aa, unsigned int an,
		const unsigned long *bb, unsigned int bn);

/* This one is reentrant */
void mul_gf2x_r(unsigned long *c,
		const unsigned long *aa, unsigned int an,
		const unsigned long *bb, unsigned int bn, mul_gf2x_pool_t);

#ifdef __cplusplus
}
#endif

/* mul1 -- the INLINES_FILE trick is used by the tuning code. */
#ifdef	INLINES_FILE
#include INLINES_FILE
#else
#include "mul-inlines.c"
#endif

/* mul2 to mul9 */
#include "mul-small.c"

/* mul_basecase dispatcher */
#include "mul-basecase.c"

/* Karatsuba, Toom, and so on are in separate .c files */
#endif	/* GF2X_H_ */
