/* Implements 128x128 -> 256 bit product using SSE2 instructions.

  Copyright 2007 Emmanuel Thome',
  with contributions from
  	Richard Brent -- interleaving of repair steps with precomputations
	Paul Zimmermann -- precomputation with simple in-order walk.

  The code in this file is extracted from the mpfq library, Copyright 2007
  Pierrick Gaudry and Emmanuel Thome'.

  This file specifically, as part of the gf2x package, is licensed as follows.

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

#include <stdint.h>

static inline
void mul2(ulong * t, ulong const * s1, ulong const * s2)
{
	typedef uint64_t v2di __attribute__ ((vector_size (16)));
	typedef union { v2di s; ulong x[2]; } v2di_proxy;
#define SHL(x,r) (v2di)__builtin_ia32_psllqi128   ((x),(r))
#define SHR(x,r) (v2di)__builtin_ia32_psrlqi128   ((x),(r))
#define SHLD(x,r) (v2di)__builtin_ia32_pslldqi128 ((x),(r))
#define SHRD(x,r) (v2di)__builtin_ia32_psrldqi128 ((x),(r))
	v2di u;
	v2di t0;
	v2di t1;
	v2di t2;

	v2di g[16];
	v2di w;
	v2di m = (v2di) { 0xeeeeeeeeeeeeeeeeUL, 0xeeeeeeeeeeeeeeeeUL, };
	/* sequence update walk */
	g[ 0] = (v2di) { 0, };
	v2di b0 = (v2di) { s2[0], s2[1], };
	g[ 1] = b0;
	v2di v1 = (v2di) { s1[0], s1[0], };
	w = -SHR(b0,63);
	v2di v2 = (v2di) { s1[1], s1[1], };
	v1 = SHR(v1 & m, 1); t1 = v1 & w;
	g[ 2] = SHL(b0, 1); g[ 3] = g[ 2] ^ b0;
	v2 = SHR(v2 & m, 1); t2 = v2 & w;
	g[ 4] = SHL(g[ 2], 1); g[ 5] = g[ 4] ^ b0;
	w = -SHR(g[ 2],63);
	g[ 6] = SHL(g[ 3], 1); g[ 7] = g[ 6] ^ b0;
	v1 = SHR(v1 & m, 1); t1 ^= v1 & w;
	g[ 8] = SHL(g[ 4], 1); g[ 9] = g[ 8] ^ b0;
	v2 = SHR(v2 & m, 1); t2 ^= v2 & w;
	g[10] = SHL(g[ 5], 1); g[11] = g[10] ^ b0;
	w = -SHR(g[4],63);
	g[12] = SHL(g[ 6], 1); g[13] = g[12] ^ b0;
	v1 = SHR(v1 & m, 1); t1 ^= v1 & w;
	g[14] = SHL(g[ 7], 1); g[15] = g[14] ^ b0;
	v2 = SHR(v2 & m, 1); t2 ^= v2 & w;



	/* round 0 */
	u = g[s1[0]       & 15]; t0  = u;
	u = g[s1[0] >>  4 & 15]; t0 ^= SHL(u,  4); t1 ^= SHR(u, 60);
	u = g[s1[0] >>  8 & 15]; t0 ^= SHL(u,  8); t1 ^= SHR(u, 56);
	u = g[s1[0] >> 12 & 15]; t0 ^= SHL(u, 12); t1 ^= SHR(u, 52);
	u = g[s1[0] >> 16 & 15]; t0 ^= SHL(u, 16); t1 ^= SHR(u, 48);
	u = g[s1[0] >> 20 & 15]; t0 ^= SHL(u, 20); t1 ^= SHR(u, 44);
	u = g[s1[0] >> 24 & 15]; t0 ^= SHL(u, 24); t1 ^= SHR(u, 40);
	u = g[s1[0] >> 28 & 15]; t0 ^= SHL(u, 28); t1 ^= SHR(u, 36);
	u = g[s1[0] >> 32 & 15]; t0 ^= SHL(u, 32); t1 ^= SHR(u, 32);
	u = g[s1[0] >> 36 & 15]; t0 ^= SHL(u, 36); t1 ^= SHR(u, 28);
	u = g[s1[0] >> 40 & 15]; t0 ^= SHL(u, 40); t1 ^= SHR(u, 24);
	u = g[s1[0] >> 44 & 15]; t0 ^= SHL(u, 44); t1 ^= SHR(u, 20);
	u = g[s1[0] >> 48 & 15]; t0 ^= SHL(u, 48); t1 ^= SHR(u, 16);
	u = g[s1[0] >> 52 & 15]; t0 ^= SHL(u, 52); t1 ^= SHR(u, 12);
	u = g[s1[0] >> 56 & 15]; t0 ^= SHL(u, 56); t1 ^= SHR(u,  8);
	u = g[s1[0] >> 60 & 15]; t0 ^= SHL(u, 60); t1 ^= SHR(u,  4);

	/* round 1 */
	u = g[s1[1]       & 15]; t1 ^= u;
	u = g[s1[1] >>  4 & 15]; t1 ^= SHL(u,  4); t2 ^= SHR(u, 60);
	u = g[s1[1] >>  8 & 15]; t1 ^= SHL(u,  8); t2 ^= SHR(u, 56);
	u = g[s1[1] >> 12 & 15]; t1 ^= SHL(u, 12); t2 ^= SHR(u, 52);
	u = g[s1[1] >> 16 & 15]; t1 ^= SHL(u, 16); t2 ^= SHR(u, 48);
	u = g[s1[1] >> 20 & 15]; t1 ^= SHL(u, 20); t2 ^= SHR(u, 44);
	u = g[s1[1] >> 24 & 15]; t1 ^= SHL(u, 24); t2 ^= SHR(u, 40);
	u = g[s1[1] >> 28 & 15]; t1 ^= SHL(u, 28); t2 ^= SHR(u, 36);
	u = g[s1[1] >> 32 & 15]; t1 ^= SHL(u, 32); t2 ^= SHR(u, 32);
	u = g[s1[1] >> 36 & 15]; t1 ^= SHL(u, 36); t2 ^= SHR(u, 28);
	u = g[s1[1] >> 40 & 15]; t1 ^= SHL(u, 40); t2 ^= SHR(u, 24);
	u = g[s1[1] >> 44 & 15]; t1 ^= SHL(u, 44); t2 ^= SHR(u, 20);
	u = g[s1[1] >> 48 & 15]; t1 ^= SHL(u, 48); t2 ^= SHR(u, 16);
	u = g[s1[1] >> 52 & 15]; t1 ^= SHL(u, 52); t2 ^= SHR(u, 12);
	u = g[s1[1] >> 56 & 15]; t1 ^= SHL(u, 56); t2 ^= SHR(u,  8);
	u = g[s1[1] >> 60 & 15]; t1 ^= SHL(u, 60); t2 ^= SHR(u,  4);
	/* end */

	/* store result */
	{
		v2di_proxy r;
		r.s = t0 ^ SHLD(t1, 64);
		t[0] = r.x[0];
		t[1] = r.x[1];
	}

	{
		v2di_proxy r;
		r.s = t2 ^ SHRD(t1, 64);
		t[2] = r.x[0];
		t[3] = r.x[1];
	}
#undef SHL
#undef SHR
#undef SHLD
#undef SHRD
}
