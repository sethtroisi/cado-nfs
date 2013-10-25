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
#ifndef GF2X_MUL9_H_
#define GF2X_MUL9_H_

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

/* TODO: if somebody comes up with a neat way to improve the interface so
 * as to remove the false dependency on pclmul, that would be nice.
 */
static inline __v2di
GF2X_FUNC(mul9clk_c_mul1) (unsigned long a, unsigned long b)
{   
    __v2di aa = (__v2di) { (long long)a, 0 };
    __v2di bb = (__v2di) { (long long)b, 0 };
    return _mm_clmulepi64_si128(aa, bb, 0);
}

/* variant with 30 multiplications */
GF2X_STORAGE_CLASS_mul9
void gf2x_mul9 (unsigned long *c, const unsigned long *a, const unsigned long *b)
{
    /* Taken from Cenk & Ozbudak 2009 */
    /* We reserve one more to follow notations of the paper */
    __v2di ab[9] = {
        (__v2di) { (long long)a[0], (long long)b[0] },
        (__v2di) { (long long)a[1], (long long)b[1] },
        (__v2di) { (long long)a[2], (long long)b[2] },
        (__v2di) { (long long)a[3], (long long)b[3] },
        (__v2di) { (long long)a[4], (long long)b[4] },
        (__v2di) { (long long)a[5], (long long)b[5] },
        (__v2di) { (long long)a[6], (long long)b[6] },
        (__v2di) { (long long)a[7], (long long)b[7] },
        (__v2di) { (long long)a[8], (long long)b[8] },
    };
    __v2di pab[30];

#if 0
    pab[ 0] = ab[0]^ab[1]^ab[2]^ab[3]^ab[4]^ab[5]^ab[6]^ab[7]^ab[8];
    pab[ 1] = ab[0]^      ab[2]^      ab[4]^      ab[6]^      ab[8];
    pab[ 2] =       ab[1]^ab[2]^ab[3]^      ab[5]^            ab[8];
    pab[ 3] = ab[0]^      ab[2]^ab[3]^      ab[5]^ab[6]^      ab[8];
    pab[ 4] = ab[0]^ab[1]^      ab[3]^ab[4]^      ab[6]^ab[7];
    pab[ 5] = ab[0]^            ab[3]^ab[4]^ab[5]^      ab[7];
    pab[ 6] =       ab[1]^ab[2]^ab[3]^            ab[6]^      ab[8];
    pab[ 7] =             ab[2]^      ab[4]^ab[5]^ab[6];
    pab[ 8] =             ab[2]^ab[3]^ab[4]^      ab[6];
    pab[ 9] =       ab[1]^      ab[3]^      ab[5]^      ab[7];
    pab[10] = ab[0]^ab[1]^            ab[4]^      ab[6]^ab[7]^ab[8];
    pab[11] = ab[0]^            ab[3]^      ab[5]^ab[6]^ab[7];
    pab[12] = ab[0]^ab[1]^            ab[4]^ab[5]^            ab[8];
    pab[13] =       ab[1]^ab[2]^      ab[4]^ab[5]^      ab[7]^ab[8];
    pab[14] = ab[0]^ab[1]^      ab[3]^            ab[6]^ab[7]^ab[8];
    pab[15] =       ab[1]^      ab[3]^ab[4]^ab[5]^            ab[8];
    pab[16] = ab[0]^      ab[2]^ab[3]^ab[4]^            ab[7];
    pab[17] =       ab[1]^            ab[4]^ab[5]^ab[6]^      ab[8];
    pab[18] = ab[0]^      ab[2]^            ab[5]^ab[6]^ab[7];
    pab[19] =             ab[2]^ab[3]^            ab[6]^ab[7];
    pab[20] =                                     ab[6]^      ab[8];
    pab[21] = ab[0]^      ab[2];
    pab[22] = ab[0]^ab[1];
    pab[23] = ab[0];
    pab[24] =       ab[1];
    pab[25] =                                           ab[7];
    pab[26] =                                           ab[7]^ab[8];
    pab[27] =                                     ab[6];
    pab[28] =                                                 ab[8];
    pab[29] =             ab[2];
#else
    /* same as above, but optimized with Maple's codegen[optimize] function
       with 'tryhard' option: 89 XORs -> 46 XORs */
    __v2di t51, t52, t53, t54, t55, t56, t57, t58, t59, t60, t61, t62,
      t63, t64, t65, t66, t67, t68, t69, t70, t71, t72, t73, t74, t75;
    t51 = ab[8];
    t55 = ab[4];
    t75 = t51^t55;
    t54 = ab[5];
    t57 = ab[2];
    t74 = t54^t57;
    t56 = ab[3];
    t58 = ab[1];
    t73 = t56^t58;
    t59 = ab[0];
    t72 = t59^t57;
    t71 = t58^t75;
    t52 = ab[7];
    t70 = t52^t56^t59;
    t53 = ab[6];
    t69 = t53^t56^t57;
    t68 = t53^t72;
    t67 = t54^t71;
    t66 = t59^t71;
    t65 = t51^t57^t73;
    t64 = t53^t70;
    t63 = t55^t70;
    t62 = t54^t68;
    t61 = t58^t64;
    t60 = t55^t61;
    pab[0] = t51^t60^t74;
    pab[1] = t68^t75;
    pab[2] = t54^t65;
    pab[3] = t56^t51^t62;
    pab[4] = t60;
    pab[5] = t54^t63;
    pab[6] = t53^t65;
    pab[7] = t55^t53^t74;
    pab[8] = t55^t69;
    pab[9] = t54^t52^t73;
    pab[10] = t53^t52^t66;
    pab[11] = t54^t64;
    pab[12] = t54^t66;
    pab[13] = t57^t52^t67;
    pab[14] = t51^t61;
    pab[15] = t56^t67;
    pab[16] = t57^t63;
    pab[17] = t53^t67;
    pab[18] = t52^t62;
    pab[19] = t52^t69;
    pab[20] = t53^t51;
    pab[21] = t72;
    pab[22] = t59^t58;
    pab[23] = t59;
    pab[24] = t58;
    pab[25] = t52;
    pab[26] = t52^t51;
    pab[27] = t53;
    pab[28] = t51;
    pab[29] = t57;
#endif

    int i;
    for (i = 0; i < 30; ++i)
        pab[i] = _mm_clmulepi64_si128(pab[i], pab[i], 0x10);

    __v2di cc[17];

#if 0
    cc[0 ] = pab[23];
    cc[1 ] = pab[22]^pab[23]^pab[24];
    cc[2 ] = pab[21]^pab[23]^pab[24]^pab[29];
    cc[3 ] = pab[28]^pab[17]^pab[2]^pab[15]^pab[7]^pab[6]^pab[5]^pab[29]^pab[21]^pab[22]^pab[12]^pab[19]^pab[9]^pab[13]^pab[11]^pab[3]^pab[26]^pab[20]^pab[27];
    cc[4 ] = pab[4]^pab[3]^pab[10]^pab[11]^pab[6]^pab[2]^pab[8]^pab[14]^pab[9]^pab[22]^pab[23]^pab[24]^pab[1]^pab[20]^pab[27]^pab[28]^pab[25];
    cc[5 ] = pab[26]^pab[25]^pab[28]^pab[0]^pab[9]^pab[21]^pab[23]^pab[29]^pab[24]^pab[1]^pab[3]^pab[13]^pab[14]^pab[5]^pab[18]^pab[16]^pab[11]^pab[15];
    cc[6 ] = pab[26]^pab[12]^pab[19]^pab[21]^pab[23]^pab[29]^pab[4]^pab[3]^pab[14]^pab[5]^pab[18]^pab[22]^pab[1]^pab[20]^pab[27];
    cc[7 ] = pab[20]^pab[27]^pab[28]^pab[25]^pab[23]^pab[0]^pab[15]^pab[7]^pab[11]^pab[6]^pab[14]^pab[5]^pab[18];
    cc[8 ] = pab[0]^pab[23]^pab[24]^pab[10]^pab[15]^pab[7]^pab[2]^pab[18]^pab[14]^pab[17]^pab[22]^pab[26]^pab[25]^pab[28];
    cc[9 ] = pab[21]^pab[23]^pab[29]^pab[24]^pab[0]^pab[16]^pab[11]^pab[7]^pab[10]^pab[2]^pab[8]^pab[18]^pab[5]^pab[28];
    cc[10] = pab[12]^pab[0]^pab[19]^pab[9]^pab[21]^pab[29]^pab[22]^pab[3]^pab[13]^pab[16]^pab[11]^pab[7]^pab[10]^pab[20]^pab[27]^pab[28]^pab[26];
    cc[11] = pab[16]^pab[11]^pab[7]^pab[10]^pab[17]^pab[5]^pab[2]^pab[0]^pab[9]^pab[4]^pab[3]^pab[22]^pab[23]^pab[24]^pab[1]^pab[20]^pab[27]^pab[28]^pab[25];
    cc[12] = pab[26]^pab[25]^pab[28]^pab[8]^pab[14]^pab[5]^pab[17]^pab[10]^pab[6]^pab[16]^pab[15]^pab[3]^pab[13]^pab[1]^pab[9]^pab[21]^pab[23]^pab[29]^pab[24];
    cc[13] = pab[8]^pab[18]^pab[2]^pab[15]^pab[16]^pab[5]^pab[29]^pab[21]^pab[23]^pab[22]^pab[0]^pab[12]^pab[19]^pab[1]^pab[11]^pab[4]^pab[3]^pab[26]^pab[20]^pab[27];
    cc[14] = pab[20]^pab[27]^pab[28]^pab[25];
    cc[15] = pab[25]^pab[26]^pab[28];
    cc[16] = pab[28];
#else
    /* same as above, optimized with codegen[optimize] with 'tryhard' */
    __v2di t100, t101, t102, t103, t104, t105, t106, t107, t108, t109, t110,
      t111, t112, t113, t114, t115, t116, t117, t118, t119, t120, t121, t122,
      t123, t124, t125, t126, t127, t128, t129, t130, t77, t79, t80, t82, t83,
      t87, t88, t89, t90, t91, t92, t94, t95, t96, t97, t98, t99;
    t82 = pab[23];
    t87 = pab[18];
    t130 = t82^t87;
    t77 = pab[28];
    t98 = pab[7];
    t129 = t77^t98;
    t79 = pab[26];
    t83 = pab[22];
    t128 = t79^t83;
    t90 = pab[15];
    t91 = pab[14];
    t127 = t90^t91;
    t97 = pab[8];
    t99 = pab[6];
    t126 = t97^t99;
    t100 = pab[5];
    t125 = t100^t90;
    t117 = pab[27]^pab[20];
    t80 = pab[25];
    t118 = t77^t80;
    t112 = t117^t118;
    t94 = pab[11];
    t124 = t112^t94;
    t103 = pab[2];
    t105 = pab[0];
    t123 = t103^t105;
    t89 = pab[16];
    t122 = t89^t94^t97;
    t121 = t100^t105^t98;
    t102 = pab[3];
    t104 = pab[1];
    t96 = pab[9];
    t120 = t102^t104^t96;
    t119 = pab[29]^pab[21];
    t116 = pab[24]^t82;
    t115 = t79^t118;
    t114 = t83^t116;
    t113 = t116^t119;
    t95 = pab[10];
    t111 = t87^t95^t116^t123^t129;
    t110 = t102^pab[19]^pab[12]^t117^t119^t128;
    t92 = pab[13];
    t109 = t92^t94^t96^t110^t129;
    t101 = pab[4];
    t108 = t100^t101^t104^t110^t130;
    t107 = t101^t103^t95^t114^t120^t124;
    t106 = t89^t91^t92^t113^t115^t120^t125;
    t88 = pab[17];
    cc[0] = t82;
    cc[1] = t114;
    cc[2] = t113;
    cc[3] = t88^t99^t103^t109^t125;
    cc[4] = t91^t107^t126;
    cc[5] = t87^t94^t105^t106;
    cc[6] = t91^t108;
    cc[7] = t99^t121^t124^t127^t130;
    cc[8] = t88^t80^t111^t127^t128;
    cc[9] = t100^t111^t119^t122;
    cc[10] = t89^t95^t105^t109;
    cc[11] = t88^t89^t107^t121;
    cc[12] = t88^t95^t106^t126;
    cc[13] = t90^t108^t122^t123;
    cc[14] = t112;
    cc[15] = t115;
    cc[16] = t77;
#endif

    _mm_storeu_si128((__v2di*)(c),    cc[0]                           ^ _mm_slli_si128(cc[1], 8));
    _mm_storeu_si128((__v2di*)(c+2),  cc[2]  ^ _mm_srli_si128(cc[1], 8) ^ _mm_slli_si128(cc[3], 8));
    _mm_storeu_si128((__v2di*)(c+4),  cc[4]  ^ _mm_srli_si128(cc[3], 8) ^ _mm_slli_si128(cc[5], 8));
    _mm_storeu_si128((__v2di*)(c+6),  cc[6]  ^ _mm_srli_si128(cc[5], 8) ^ _mm_slli_si128(cc[7], 8));
    _mm_storeu_si128((__v2di*)(c+8),  cc[8]  ^ _mm_srli_si128(cc[7], 8) ^ _mm_slli_si128(cc[9], 8));
    _mm_storeu_si128((__v2di*)(c+10), cc[10] ^ _mm_srli_si128(cc[9], 8) ^ _mm_slli_si128(cc[11], 8));
    _mm_storeu_si128((__v2di*)(c+12), cc[12] ^ _mm_srli_si128(cc[11], 8) ^ _mm_slli_si128(cc[13], 8));
    _mm_storeu_si128((__v2di*)(c+14), cc[14] ^ _mm_srli_si128(cc[13], 8) ^ _mm_slli_si128(cc[15], 8));
    _mm_storeu_si128((__v2di*)(c+16), cc[16] ^ _mm_srli_si128(cc[15], 8));
}

#endif  /* GF2X_MUL9_H_ */
