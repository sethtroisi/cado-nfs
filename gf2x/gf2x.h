/* This file is part of the gf2x library.

   Copyright 2007, 2008, 2009, 2010, 2013, 2014, 2015
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

/* Multiplication over GF(2)[x] */

#ifndef GF2X_H_
#define GF2X_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <stddef.h>     /* size_t */

/* This is the toplevel multiplication routine, the only one that really
 * matters. a and b point to polynomials whose coefficients span an and
 * bn machine words (unsigned longs).
 *
 * c must have enough room to hold an+bn words.
 *
 * The destination pointer c may alias either a or b (that is, one may
 * have c==a or c==b), but any other kind of overlap is unsupported).
 */
extern void gf2x_mul(unsigned long *c,
		const unsigned long *a, unsigned long an,
		const unsigned long *b, unsigned long bn);

/* The second version is reentrant */
struct gf2x_mul_pool_s {
	unsigned long * stk;
	size_t stk_size;
};
typedef struct gf2x_mul_pool_s gf2x_mul_pool_t[1];
extern void gf2x_mul_pool_init(gf2x_mul_pool_t);
extern void gf2x_mul_pool_clear(gf2x_mul_pool_t);

/* If the gf2x_mul_pool_t arg is passed as NULL, a new pool is created
 * (and freed) for that multiplication */
extern void gf2x_mul_r(unsigned long *c,
		const unsigned long *a, unsigned long an,
		const unsigned long *b, unsigned long bn, gf2x_mul_pool_t);

#ifdef __cplusplus
}
#endif

#endif	/* GF2X_H_ */
