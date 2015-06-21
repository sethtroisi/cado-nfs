/* This file is part of the gf2x library.

   Copyright 2010, 2012, 2013, 2015
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

#ifndef GF2X_MUL9_H_
#define GF2X_MUL9_H_

#include "gf2x.h"
/* All gf2x source files for lowlevel functions must include gf2x-small.h
 * This is mandatory for the tuning mechanism. */
#include "gf2x/gf2x-small.h"

#if GF2X_WORDSIZE != 64
#error "This code is for 64-bit only"
#endif

#ifndef GF2X_HAVE_PCLMUL_SUPPORT
#error "This code needs pclmul support"
#endif

#define PXOR(lop, rop) _mm_xor_si128((lop), (rop))
#define PXOR3(op1, op2, op3) PXOR(op1, PXOR(op2, op3))
#define PXOR4(op1, op2, op3, op4) PXOR(op1, PXOR3(op2, op3, op4))
#define PXOR5(op1, op2, op3, op4, op5) PXOR(op1, PXOR4(op2, op3, op4, op5))
#define PXOR6(op1, op2, op3, op4, op5, op6) PXOR(op1, PXOR5(op2, op3, op4, op5, op6))
#define PXOR7(op1, op2, op3, op4, op5, op6, op7) PXOR(op1, PXOR6(op2, op3, op4, op5, op6, op7))
#define PZERO    _mm_setzero_si128()

/* variant with 30 multiplications */
GF2X_STORAGE_CLASS_mul9
void gf2x_mul9 (unsigned long *c, const unsigned long *a, const unsigned long *b)
{
    /* Taken from Cenk & Ozbudak 2009 */
    /* We reserve one more to follow notations of the paper */
    __m128i ab[9] = {
        _gf2x_mm_setr_epi64(a[0], b[0]),
        _gf2x_mm_setr_epi64(a[1], b[1]),
        _gf2x_mm_setr_epi64(a[2], b[2]),
        _gf2x_mm_setr_epi64(a[3], b[3]),
        _gf2x_mm_setr_epi64(a[4], b[4]),
        _gf2x_mm_setr_epi64(a[5], b[5]),
        _gf2x_mm_setr_epi64(a[6], b[6]),
        _gf2x_mm_setr_epi64(a[7], b[7]),
        _gf2x_mm_setr_epi64(a[8], b[8]),
    };
    __m128i pab[30];

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
    __m128i t51, t52, t53, t54, t55, t56, t57, t58, t59, t60, t61, t62,
      t63, t64, t65, t66, t67, t68, t69, t70, t71, t72, t73, t74, t75;
    t51 = ab[8];
    t55 = ab[4];
    t75 = PXOR(t51, t55);
    t54 = ab[5];
    t57 = ab[2];
    t74 = PXOR(t54, t57);
    t56 = ab[3];
    t58 = ab[1];
    t73 = PXOR(t56, t58);
    t59 = ab[0];
    t72 = PXOR(t59, t57);
    t71 = PXOR(t58, t75);
    t52 = ab[7];
    t70 = PXOR3(t52, t56, t59);
    t53 = ab[6];
    t69 = PXOR3(t53, t56, t57);
    t68 = PXOR(t53, t72);
    t67 = PXOR(t54, t71);
    t66 = PXOR(t59, t71);
    t65 = PXOR3(t51, t57, t73);
    t64 = PXOR(t53, t70);
    t63 = PXOR(t55, t70);
    t62 = PXOR(t54, t68);
    t61 = PXOR(t58, t64);
    t60 = PXOR(t55, t61);
    pab[0] = PXOR3(t51, t60, t74);
    pab[1] = PXOR(t68, t75);
    pab[2] = PXOR(t54, t65);
    pab[3] = PXOR3(t56, t51, t62);
    pab[4] = t60;
    pab[5] = PXOR(t54, t63);
    pab[6] = PXOR(t53, t65);
    pab[7] = PXOR3(t55, t53, t74);
    pab[8] = PXOR(t55, t69);
    pab[9] = PXOR3(t54, t52, t73);
    pab[10] = PXOR3(t53, t52, t66);
    pab[11] = PXOR(t54, t64);
    pab[12] = PXOR(t54, t66);
    pab[13] = PXOR3(t57, t52, t67);
    pab[14] = PXOR(t51, t61);
    pab[15] = PXOR(t56, t67);
    pab[16] = PXOR(t57, t63);
    pab[17] = PXOR(t53, t67);
    pab[18] = PXOR(t52, t62);
    pab[19] = PXOR(t52, t69);
    pab[20] = PXOR(t53, t51);
    pab[21] = t72;
    pab[22] = PXOR(t59, t58);
    pab[23] = t59;
    pab[24] = t58;
    pab[25] = t52;
    pab[26] = PXOR(t52, t51);
    pab[27] = t53;
    pab[28] = t51;
    pab[29] = t57;
#endif

    int i;
    for (i = 0; i < 30; ++i)
        pab[i] = _mm_clmulepi64_si128(pab[i], pab[i], 0x10);

    __m128i cc[17];

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
    __m128i t100, t101, t102, t103, t104, t105, t106, t107, t108, t109, t110,
      t111, t112, t113, t114, t115, t116, t117, t118, t119, t120, t121, t122,
      t123, t124, t125, t126, t127, t128, t129, t130, t77, t79, t80, t82, t83,
      t87, t88, t89, t90, t91, t92, t94, t95, t96, t97, t98, t99;
    t82 = pab[23];
    t87 = pab[18];
    t130 = PXOR(t82, t87);
    t77 = pab[28];
    t98 = pab[7];
    t129 = PXOR(t77, t98);
    t79 = pab[26];
    t83 = pab[22];
    t128 = PXOR(t79, t83);
    t90 = pab[15];
    t91 = pab[14];
    t127 = PXOR(t90, t91);
    t97 = pab[8];
    t99 = pab[6];
    t126 = PXOR(t97, t99);
    t100 = pab[5];
    t125 = PXOR(t100, t90);
    t117 = PXOR(pab[27], pab[20]);
    t80 = pab[25];
    t118 = PXOR(t77, t80);
    t112 = PXOR(t117, t118);
    t94 = pab[11];
    t124 = PXOR(t112, t94);
    t103 = pab[2];
    t105 = pab[0];
    t123 = PXOR(t103, t105);
    t89 = pab[16];
    t122 = PXOR3(t89, t94, t97);
    t121 = PXOR3(t100, t105, t98);
    t102 = pab[3];
    t104 = pab[1];
    t96 = pab[9];
    t120 = PXOR3(t102, t104, t96);
    t119 = PXOR(pab[29], pab[21]);
    t116 = PXOR(pab[24], t82);
    t115 = PXOR(t79, t118);
    t114 = PXOR(t83, t116);
    t113 = PXOR(t116, t119);
    t95 = pab[10];
    t111 = PXOR5(t87, t95, t116, t123, t129);
    t110 = PXOR6(t102, pab[19], pab[12], t117, t119, t128);
    t92 = pab[13];
    t109 = PXOR5(t92, t94, t96, t110, t129);
    t101 = pab[4];
    t108 = PXOR5(t100, t101, t104, t110, t130);
    t107 = PXOR6(t101, t103, t95, t114, t120, t124);
    t106 = PXOR7(t89, t91, t92, t113, t115, t120, t125);
    t88 = pab[17];
    cc[0] = t82;
    cc[1] = t114;
    cc[2] = t113;
    cc[3] = PXOR5(t88, t99, t103, t109, t125);
    cc[4] = PXOR3(t91, t107, t126);
    cc[5] = PXOR4(t87, t94, t105, t106);
    cc[6] = PXOR(t91, t108);
    cc[7] = PXOR5(t99, t121, t124, t127, t130);
    cc[8] = PXOR5(t88, t80, t111, t127, t128);
    cc[9] = PXOR4(t100, t111, t119, t122);
    cc[10] = PXOR4(t89, t95, t105, t109);
    cc[11] = PXOR4(t88, t89, t107, t121);
    cc[12] = PXOR4(t88, t95, t106, t126);
    cc[13] = PXOR4(t90, t108, t122, t123);
    cc[14] = t112;
    cc[15] = t115;
    cc[16] = t77;
#endif

    _mm_storeu_si128((__m128i*)(c),    PXOR(cc[0],                           _mm_slli_si128(cc[1], 8)));
    _mm_storeu_si128((__m128i*)(c+2),  PXOR3(cc[2] , _mm_srli_si128(cc[1], 8), _mm_slli_si128(cc[3], 8)));
    _mm_storeu_si128((__m128i*)(c+4),  PXOR3(cc[4] , _mm_srli_si128(cc[3], 8), _mm_slli_si128(cc[5], 8)));
    _mm_storeu_si128((__m128i*)(c+6),  PXOR3(cc[6] , _mm_srli_si128(cc[5], 8), _mm_slli_si128(cc[7], 8)));
    _mm_storeu_si128((__m128i*)(c+8),  PXOR3(cc[8] , _mm_srli_si128(cc[7], 8), _mm_slli_si128(cc[9], 8)));
    _mm_storeu_si128((__m128i*)(c+10), PXOR3(cc[10], _mm_srli_si128(cc[9], 8), _mm_slli_si128(cc[11], 8)));
    _mm_storeu_si128((__m128i*)(c+12), PXOR3(cc[12], _mm_srli_si128(cc[11], 8), _mm_slli_si128(cc[13], 8)));
    _mm_storeu_si128((__m128i*)(c+14), PXOR3(cc[14], _mm_srli_si128(cc[13], 8), _mm_slli_si128(cc[15], 8)));
    _mm_storeu_si128((__m128i*)(c+16), PXOR(cc[16], _mm_srli_si128(cc[15], 8)));
#undef PXOR
#undef PXOR3
#undef PXOR4
#undef PXOR5
#undef PXOR6
#undef PXOR7
#undef PZERO
}

#endif  /* GF2X_MUL9_H_ */
