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
#include <string.h>
#include <stdint.h>

#include "gf2x.h"

#ifdef MUL_FFT_TABLE
int64_t T_FFT_TAB[][2] = MUL_FFT_TABLE;
#endif

/* This is the toplevel multiplication routine. It handles the temporary
 * storage if necessary.
 */

/* Having this static ensures initialization to zero */
static mul_gf2x_pool_t global_pool;

void mul_gf2x_pool_init(mul_gf2x_pool_t p)
{
    memset(p, 0, sizeof(p));
}

void mul_gf2x_pool_clear(mul_gf2x_pool_t p)
{
    free(p->stk);
    p->stk_size = 0;
}

void mul_gf2x(ulong * c,
	      const ulong * a, unsigned int sa,
	      const ulong * b, unsigned int sb)
{
    mul_gf2x_r(c, a, sa, b, sb, global_pool);
}

void mul_gf2x_r(ulong * c,
		const ulong * a, unsigned int sa,
		const ulong * b, unsigned int sb, mul_gf2x_pool_t pool)
{
    unsigned int sc = sa + sb;

    if (sa > sb) {
	mul_gf2x_r(c, b, sb, a, sa, pool);
	return;
    }
    // now sa <= sb (note: sa and sb are interchanged in Toom3uMul etc

    if (sa < MUL_KARA_THRESHOLD) {
        /* This calls the hand-crafted code if sa == sb */
        mul_basecase(c, a, sa, b, sb);
	return;
    }

    /* This ugly cpp block entirely disables the FFT if it has
     * not yet been tuned */
#ifdef MUL_FFT_TABLE
    long ix, K, sab = sc / 2;
    long FFT2;
    for (ix = 0; T_FFT_TAB[ix + 1][0] <= sab; ix++);
    /* now T_FFT_TAB[ix][0] <= sab < T_FFT_TAB[ix+1][0] */
    K = T_FFT_TAB[ix][1];
    if (K < 0) {
	FFT2 = 1;
	K = -K;
    } else {
	FFT2 = 0;
    }

    /* FFTMul can handle unbalanced operands if not too
     * small: return the result in {c, sa+sb} */
    if ((K >= 3) && (sc >= MUL_FFT_THRESHOLD)) {
	if (FFT2)
	    FFTMul2(c, a, sa, b, sb, K);	// Split FFT into two
	else
	    FFTMul(c, a, sa, b, sb, K);	// Don't split here

	return;
    }
#endif

    unsigned int sp1, sp2, sp3, sp;
    sp1 = toomspace(sa);	// Space for balanced TC routines
    sp2 = toomuspace(2 * sa);	// Space for unbalanced TC routines
    sp3 = 2 * sa + toomspace(sa);	// Space for unbalanced TC routines w/ lazy cut

    sp = sp1;
    if (sp < sp2)
	sp = sp2;
    if (sp < sp3)
	sp = sp3;		// Worst-case space required

    if (pool->stk_size < sp) {
	pool->stk = realloc(pool->stk, sp * sizeof(ulong));
	pool->stk_size = sp;
    }

    if (sa == sb) {
	// Avoid copy in common case
	mul_toom(c, a, b, sa, pool->stk);
    } else if ((sa == (sb + 1) / 2) && BestuToom(sb)) {
	// Another common case
	// due to GCD algorithm
	mul_tc3u(c, b, sb, a, pool->stk);
    } else {
	ulong *v = pool->stk + toomspace(sa);

	unsigned int i, j;

	memset(c, 0, sc * sizeof(ulong));

	for (;;) {
	    if (sa == 0)
		break;

	    if (sa == 1) {
		c[sb] ^= addmul_1_n(c, c, b, sb, a[0]);
		break;
	    }

	    /* TODO: Should do addmul_2_n here, that would be an easy
	     * improvement. */

	    // finally: the general case
	    for (i = 0; i + sa <= sb; i += sa) {
		// Generic (balanced) Toom-Cook mult.
		mul_toom(v, a, b + i, sa, pool->stk);
		for (j = 0; j < 2 * sa; j++)
		    c[i + j] ^= v[j];
	    }

	    {
		const ulong *t;
		unsigned int st;

		/* Swap a and b, and go for the next spin */
		t = a;
		st = sa;
		a = b + i;
		sa = sb - i;
		b = t;
		sb = st;
	    }
	    c = c + i;
	}
    }
}

/* vim: set sw=4 sta et: */
