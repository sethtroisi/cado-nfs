/* Various matrix operations.
 *  
 * Author: Emmanuel Thomé
 * 
 * Copyright 2010 INRIA
 * 
 * This file is part of the CADO project.
 * 
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation; either version 2.1 of the License, or (at
 * your option) any later version.

 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 * License for more details.

 * You should have received a copy of the GNU Lesser General Public
 * License along with CADO-NFS; see the file COPYING.  If not, write to
 * the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA 02110-1301, USA.
 */

#include <stdio.h>
#include <stdint.h>
#include <time.h>
#include <string.h>
#include <assert.h>
#include <emmintrin.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <inttypes.h>

/* Define this to check these functions against the m4ri library */
#define xxxHAVE_M4RI
#ifdef  HAVE_M4RI
#include "m4ri/m4ri.h"
#endif  /* HAVE_M4RI */

/* This is **only** for 64 * 64 matrices */

#define WBITS   64
typedef uint64_t mat64[64];
typedef uint64_t * mat64_ptr;
typedef const uint64_t * mat64_srcptr;

#if 0
static int urandom_fd = -1;

static inline uint64_t rand64()
{
    if (urandom_fd == -1)
        urandom_fd = open("/dev/urandom", O_RDONLY);
    if (urandom_fd < 0) abort();
    uint64_t r;
    size_t s = read(urandom_fd, &r, sizeof(uint64_t));
    if (s != sizeof(uint64_t)) abort();
    return r;
}

static inline void rand64_mem(uint64_t * ptr, size_t z)
{
    if (urandom_fd == -1)
        urandom_fd = open("/dev/urandom", O_RDONLY);
    if (urandom_fd < 0) abort();
    size_t s = read(urandom_fd, ptr, z * sizeof(uint64_t));
    if (s != z) abort();
}
#else
static inline uint64_t rand64()
{
    uint64_t r0 = rand();
    uint64_t r1 = rand(); r1 <<= 22;
    uint64_t r2 = rand(); r2 <<= 44;
    return r0 ^r1 ^ r2;
}
static inline void rand64_mem(uint64_t * ptr, size_t z)
{
    for( ; z-- ; ptr[z] = rand64());
}
#endif

static inline uint64_t bitrev(uint64_t a)
{
    a = (a >> 32) ^ (a << 32);
    uint64_t m;
    m = 0x0000ffff0000ffffUL;
    a = ((a >> 16) & m) ^ ((a << 16) & ~m);
    m = 0x00ff00ff00ff00ffUL;
    a = ((a >> 8) & m) ^ ((a << 8) & ~m);
    m = 0x0f0f0f0f0f0f0f0fUL;
    a = ((a >> 4) & m) ^ ((a << 4) & ~m);
    m = 0x3333333333333333UL;
    a = ((a >> 2) & m) ^ ((a << 2) & ~m);
    m = 0x5555555555555555UL;
    a = ((a >> 1) & m) ^ ((a << 1) & ~m);
    return a;
}
/* like bitrev, but keep nibbles intact */
static inline uint64_t nibrev(uint64_t a)
{
    a = (a >> 32) ^ (a << 32);
    uint64_t m;
    m = 0x0000ffff0000ffffUL;
    a = ((a >> 16) & m) ^ ((a << 16) & ~m);
    m = 0x00ff00ff00ff00ffUL;
    a = ((a >> 8) & m) ^ ((a << 8) & ~m);
    m = 0x0f0f0f0f0f0f0f0fUL;
    a = ((a >> 4) & m) ^ ((a << 4) & ~m);
    return a;
}
#define mul_6464_6464 mul_6464_6464_sse
#define add_6464_6464 add_6464_6464_C
#define mul_N64_6464 mul_N64_6464_sse
#define mul_N64_T6464 mul_N64_6464_vec
#define addmul_To64_o64 addmul_To64_o64_lsb_sse_v1
#define mul_o64_6464 mul_o64_6464_C_lsb
#define mul_o64_T6464 mul_o64_T6464_C_parity

/* level 1 */

static inline void mul_6464_6464_sse(mat64_ptr C, mat64_srcptr A, mat64_srcptr B)
{
    int i;
    memset(C, 0, sizeof(mat64));
 
    /* As per emmintrin.h, only the __m128i type is declared with
     * attribute __may_alias__, meaning that accesses through this
     * pointer may alias other types (e.g. uint64_t's, in our cases). For
     * this reason, we _must_ use this pointer type, and not __v2di, for
     * our pointer types. Reading from a __m128 * into a __v2di, or the
     * converse, are legal operations.
     */
    __m128i *Cw = (__m128i *) C;
    __m128i *Aw = (__m128i *) A;

    for (int j = 0; j < 64; j += 2) {
	__v2di c = { 0, 0 };
	__v2di a = *Aw++;

	__v2di one = { 1, 1, };
#define SHR(x,r) _mm_srli_epi64((x),(r))
	for (i = 0; i < 64; i++) {
	    __v2di bw = { B[i], B[i], };

	    c ^= (bw & -(a & one));
	    a = SHR(a, 1);
	}
#undef  SHR
	*Cw++ = c;
    }
}

static inline void mul_6464_6464_v2(mat64_ptr C, mat64_srcptr A, mat64_srcptr B)
{
    memset(C, 0, sizeof(mat64));
 
    for (int i = 0; i < 64; i++) {
    for (int j = 0; j < 64; j++) {
        if ((A[i]>>j)&1)
            C[i]^=B[j];
    }
    }
}

static inline void add_6464_6464_C(mat64_ptr C, mat64_srcptr A, mat64_srcptr B)
{
    for (int j = 0; j < 64; j++) {
        C[j] = A[j] ^ B[j];
    }
}


static inline void addmul_To64_o64_lsb(uint64_t * r, uint64_t a, uint64_t w)
{
    /* Dans un sens */
    for (unsigned int i = 0; i < 64; i++) {
	*r++ ^= w & -(a & 1);
	a >>= 1;
    }
}

static inline void addmul_To64_o64_msb(uint64_t * r, uint64_t a, uint64_t w)
{
    /* Dans l'autre -- va un poil plus vite. */
    for (unsigned int i = 0; i < 64; i++) {
	*r++ ^= w & (((int64_t) a) >> 63);
	a <<= 1;
    }
}

static inline void addmul_To64_o64_lsb_sse_v0(uint64_t * r, uint64_t a, uint64_t w)
{
    /* À peu près comme la méthode 1, mais pas mieux */
    typedef uint64_t mvec_t[2];
    mvec_t mb[4] = {
	{0, 0}, {w, 0}, {0, w}, {w, w},
    };
    for (int i = 0; i < 64; i += 2) {
	const uint64_t *y = mb[a & 3];
	*r++ ^= y[0];
	*r++ ^= y[1];
	a >>= 2;
    }
}

static inline void addmul_To64_o64_lsb_sse_v1(uint64_t * r, uint64_t a, uint64_t w)
{
    /* Avec des sse-2 */
    __v2di mb[4] = {
	(__v2di) {0, 0},
	(__v2di) {w, 0},
	(__v2di) {0, w},
	(__v2di) {w, w},
    };
    __m128i *sr = (__m128i *) r;
    for (int i = 0; i < 64; i += 2) {
	*sr++ ^= mb[a & 3];
	a >>= 2;
    }
}

static inline void mul_o64_6464_C_lsb(uint64_t * r, uint64_t a, mat64_srcptr w)
{
    unsigned long c = 0;
    for (unsigned int i = 0; i < WBITS; i++) {
	c ^= (w[i] & -(a & 1UL));
	a >>= 1UL;
    }
    *r = c;
}

static inline void mul_o64_6464_C_msb(uint64_t *r,
                   uint64_t a,
                   mat64_srcptr w)
{
    uint64_t c = 0UL;
    for (int i = 64 - 1; i >= 0; i--) {
        c ^= (w[i] & (((long) a) >> (64 - 1)));
        a <<= 1UL;
    }
    *r = c;
}


static inline void mul_o64_T6464_C_parity(uint64_t * w, uint64_t a, mat64_srcptr b)
{
    // Uses unoptimized __builtin_parityl function -- maybe better with gcc 4.3
    // note that popcnt is faster in asm than the more restricted parity
    // functions. So if it's available, it should be tested.
    unsigned long c = 0;
    for (unsigned int i = 0; i < WBITS; i++) {
        uint64_t p = __builtin_parityl(a & b[i]);
	c ^= p << i;
    }
    *w = c;
}

/* This is stolen from code by D. Harvey. (GPL, thus can't stay here) */
#define XMIX32(b, a) (((((a) << 32) ^ (a)) >> 32) + \
                     ((((b) >> 32) ^ (b)) << 32))
#define XMIX16(b, a) (((((a) >> 16) ^ (a)) & 0x0000FFFF0000FFFFll) +     \
                     ((((b) << 16) ^ (b)) & 0xFFFF0000FFFF0000ll));
#define XMIX8(b, a) (((((a) >> 8) ^ (a)) & 0x00FF00FF00FF00FFll) + \
                    ((((b) << 8) ^ (b)) & 0xFF00FF00FF00FF00ll));
#define XMIX4(b, a) (((((a) >> 4) ^ (a)) & 0x0F0F0F0F0F0F0F0Fll) + \
                    ((((b) << 4) ^ (b)) & 0xF0F0F0F0F0F0F0F0ll));
#define XMIX2(b, a) (((((a) >> 2) ^ (a)) & 0x3333333333333333ll) + \
                    ((((b) << 2) ^ (b)) & 0xCCCCCCCCCCCCCCCCll));
#define XMIX1(b, a) (((((a) >> 1) ^ (a)) & 0x5555555555555555ll) + \
                    ((((b) << 1) ^ (b)) & 0xAAAAAAAAAAAAAAAAll));
static inline uint64_t _parity64_helper2(const uint64_t* buf, uint64_t a)
{
   uint64_t a0, a1, b0, b1, c0, c1;
   a0 = XMIX32(buf[0x20] & a, buf[0x00] & a);
   a1 = XMIX32(buf[0x30] & a, buf[0x10] & a);
   b0 = XMIX16(a1, a0);
   a0 = XMIX32(buf[0x28] & a, buf[0x08] & a);
   a1 = XMIX32(buf[0x38] & a, buf[0x18] & a);
   b1 = XMIX16(a1, a0);
   c0 = XMIX8(b1, b0);
   a0 = XMIX32(buf[0x24] & a, buf[0x04] & a);
   a1 = XMIX32(buf[0x34] & a, buf[0x14] & a);
   b0 = XMIX16(a1, a0);
   a0 = XMIX32(buf[0x2C] & a, buf[0x0C] & a);
   a1 = XMIX32(buf[0x3C] & a, buf[0x1C] & a);
   b1 = XMIX16(a1, a0);
   c1 = XMIX8(b1, b0);
   return XMIX4(c1, c0);
}

static inline void mul_o64_T6464_C_parity3(uint64_t * w, uint64_t a, mat64_srcptr b)
{
   uint64_t d0, d1, e0, e1;

   d0 = _parity64_helper2(b, a);
   d1 = _parity64_helper2(b + 2, a);
   e0 = XMIX2(d1, d0);

   d0 = _parity64_helper2(b + 1, a);
   d1 = _parity64_helper2(b + 3, a);
   e1 = XMIX2(d1, d0);

   *w = XMIX1(e1, e0);
}


static inline void transp_6464(mat64_ptr dst, mat64_srcptr src)
{
    int i, j;
    for (i = 0; i < 64; i++) {
	dst[i] = 0;
	for (j = 0; j < 64; j++) {
	    dst[i] ^= ((src[j] >> i) & 1UL) << j;
	}
    }
}

static inline void copy_6464(mat64_ptr dst, mat64_srcptr src)
{
    memcpy(dst, src, sizeof(mat64));
}


/* level 2 */

static inline void copy_N64(uint64_t * dst, const uint64_t * src, unsigned long m)
{
    memcpy(dst, src, m * sizeof(uint64_t));
}

static inline int cmp_N64(const uint64_t * dst, const uint64_t * src, unsigned long m)
{
    return memcmp(dst, src, m * sizeof(uint64_t));
}


static inline void mul_N64_6464_lookup4(uint64_t *C,
                   const uint64_t *A,
                   const uint64_t *B, unsigned long m)
{
    uint64_t Bx[16][16];
    for(int j = 0 ; j < 16 ; j++) {
        const uint64_t * bb = B + 4 * j;
        uint64_t w = 0;
        Bx[j][0]  = w; w ^= bb[0];
        Bx[j][1]  = w; w ^= bb[1];
        Bx[j][3]  = w; w ^= bb[0];
        Bx[j][2]  = w; w ^= bb[2];
        Bx[j][6]  = w; w ^= bb[0];
        Bx[j][7]  = w; w ^= bb[1];
        Bx[j][5]  = w; w ^= bb[0];
        Bx[j][4]  = w; w ^= bb[3];
        Bx[j][12] = w; w ^= bb[0];
        Bx[j][13] = w; w ^= bb[1];
        Bx[j][15] = w; w ^= bb[0];
        Bx[j][14] = w; w ^= bb[2];
        Bx[j][10] = w; w ^= bb[0];
        Bx[j][11] = w; w ^= bb[1];
        Bx[j][9]  = w; w ^= bb[0];
        Bx[j][8]  = w;
    }
    memset(C, 0, m * sizeof(unsigned long));
    for (unsigned long i = 0; i < m; i++) {
        uint64_t aa = A[i];
        C[i] = Bx[0][aa & 15]; aa>>=4;
        C[i]^= Bx[1][aa & 15]; aa>>=4;
        C[i]^= Bx[2][aa & 15]; aa>>=4;
        C[i]^= Bx[3][aa & 15]; aa>>=4;
        C[i]^= Bx[4][aa & 15]; aa>>=4;
        C[i]^= Bx[5][aa & 15]; aa>>=4;
        C[i]^= Bx[6][aa & 15]; aa>>=4;
        C[i]^= Bx[7][aa & 15]; aa>>=4;
        C[i]^= Bx[8][aa & 15]; aa>>=4;
        C[i]^= Bx[9][aa & 15]; aa>>=4;
        C[i]^= Bx[10][aa & 15]; aa>>=4;
        C[i]^= Bx[11][aa & 15]; aa>>=4;
        C[i]^= Bx[12][aa & 15]; aa>>=4;
        C[i]^= Bx[13][aa & 15]; aa>>=4;
        C[i]^= Bx[14][aa & 15]; aa>>=4;
        C[i]^= Bx[15][aa];
    }
}

static inline void mul_N64_6464_lookup8(uint64_t *C,
                   const uint64_t *A,
                   const uint64_t *B, unsigned long m)
{
    uint64_t Bx[8][256];
    for(int j = 0 ; j < 8 ; j++) {
        const uint64_t * bb = B + 8 * j;
        uint64_t w = 0;
        Bx[j][0] = w; w ^= bb[0];
        Bx[j][1] = w; w ^= bb[1];
        Bx[j][3] = w; w ^= bb[0];
        Bx[j][2] = w; w ^= bb[2];
        Bx[j][6] = w; w ^= bb[0];
        Bx[j][7] = w; w ^= bb[1];
        Bx[j][5] = w; w ^= bb[0];
        Bx[j][4] = w; w ^= bb[3];
        Bx[j][12] = w; w ^= bb[0];
        Bx[j][13] = w; w ^= bb[1];
        Bx[j][15] = w; w ^= bb[0];
        Bx[j][14] = w; w ^= bb[2];
        Bx[j][10] = w; w ^= bb[0];
        Bx[j][11] = w; w ^= bb[1];
        Bx[j][9] = w; w ^= bb[0];
        Bx[j][8] = w; w ^= bb[4];
        Bx[j][24] = w; w ^= bb[0];
        Bx[j][25] = w; w ^= bb[1];
        Bx[j][27] = w; w ^= bb[0];
        Bx[j][26] = w; w ^= bb[2];
        Bx[j][30] = w; w ^= bb[0];
        Bx[j][31] = w; w ^= bb[1];
        Bx[j][29] = w; w ^= bb[0];
        Bx[j][28] = w; w ^= bb[3];
        Bx[j][20] = w; w ^= bb[0];
        Bx[j][21] = w; w ^= bb[1];
        Bx[j][23] = w; w ^= bb[0];
        Bx[j][22] = w; w ^= bb[2];
        Bx[j][18] = w; w ^= bb[0];
        Bx[j][19] = w; w ^= bb[1];
        Bx[j][17] = w; w ^= bb[0];
        Bx[j][16] = w; w ^= bb[5];
        Bx[j][48] = w; w ^= bb[0];
        Bx[j][49] = w; w ^= bb[1];
        Bx[j][51] = w; w ^= bb[0];
        Bx[j][50] = w; w ^= bb[2];
        Bx[j][54] = w; w ^= bb[0];
        Bx[j][55] = w; w ^= bb[1];
        Bx[j][53] = w; w ^= bb[0];
        Bx[j][52] = w; w ^= bb[3];
        Bx[j][60] = w; w ^= bb[0];
        Bx[j][61] = w; w ^= bb[1];
        Bx[j][63] = w; w ^= bb[0];
        Bx[j][62] = w; w ^= bb[2];
        Bx[j][58] = w; w ^= bb[0];
        Bx[j][59] = w; w ^= bb[1];
        Bx[j][57] = w; w ^= bb[0];
        Bx[j][56] = w; w ^= bb[4];
        Bx[j][40] = w; w ^= bb[0];
        Bx[j][41] = w; w ^= bb[1];
        Bx[j][43] = w; w ^= bb[0];
        Bx[j][42] = w; w ^= bb[2];
        Bx[j][46] = w; w ^= bb[0];
        Bx[j][47] = w; w ^= bb[1];
        Bx[j][45] = w; w ^= bb[0];
        Bx[j][44] = w; w ^= bb[3];
        Bx[j][36] = w; w ^= bb[0];
        Bx[j][37] = w; w ^= bb[1];
        Bx[j][39] = w; w ^= bb[0];
        Bx[j][38] = w; w ^= bb[2];
        Bx[j][34] = w; w ^= bb[0];
        Bx[j][35] = w; w ^= bb[1];
        Bx[j][33] = w; w ^= bb[0];
        Bx[j][32] = w; w ^= bb[6];
        Bx[j][96] = w; w ^= bb[0];
        Bx[j][97] = w; w ^= bb[1];
        Bx[j][99] = w; w ^= bb[0];
        Bx[j][98] = w; w ^= bb[2];
        Bx[j][102] = w; w ^= bb[0];
        Bx[j][103] = w; w ^= bb[1];
        Bx[j][101] = w; w ^= bb[0];
        Bx[j][100] = w; w ^= bb[3];
        Bx[j][108] = w; w ^= bb[0];
        Bx[j][109] = w; w ^= bb[1];
        Bx[j][111] = w; w ^= bb[0];
        Bx[j][110] = w; w ^= bb[2];
        Bx[j][106] = w; w ^= bb[0];
        Bx[j][107] = w; w ^= bb[1];
        Bx[j][105] = w; w ^= bb[0];
        Bx[j][104] = w; w ^= bb[4];
        Bx[j][120] = w; w ^= bb[0];
        Bx[j][121] = w; w ^= bb[1];
        Bx[j][123] = w; w ^= bb[0];
        Bx[j][122] = w; w ^= bb[2];
        Bx[j][126] = w; w ^= bb[0];
        Bx[j][127] = w; w ^= bb[1];
        Bx[j][125] = w; w ^= bb[0];
        Bx[j][124] = w; w ^= bb[3];
        Bx[j][116] = w; w ^= bb[0];
        Bx[j][117] = w; w ^= bb[1];
        Bx[j][119] = w; w ^= bb[0];
        Bx[j][118] = w; w ^= bb[2];
        Bx[j][114] = w; w ^= bb[0];
        Bx[j][115] = w; w ^= bb[1];
        Bx[j][113] = w; w ^= bb[0];
        Bx[j][112] = w; w ^= bb[5];
        Bx[j][80] = w; w ^= bb[0];
        Bx[j][81] = w; w ^= bb[1];
        Bx[j][83] = w; w ^= bb[0];
        Bx[j][82] = w; w ^= bb[2];
        Bx[j][86] = w; w ^= bb[0];
        Bx[j][87] = w; w ^= bb[1];
        Bx[j][85] = w; w ^= bb[0];
        Bx[j][84] = w; w ^= bb[3];
        Bx[j][92] = w; w ^= bb[0];
        Bx[j][93] = w; w ^= bb[1];
        Bx[j][95] = w; w ^= bb[0];
        Bx[j][94] = w; w ^= bb[2];
        Bx[j][90] = w; w ^= bb[0];
        Bx[j][91] = w; w ^= bb[1];
        Bx[j][89] = w; w ^= bb[0];
        Bx[j][88] = w; w ^= bb[4];
        Bx[j][72] = w; w ^= bb[0];
        Bx[j][73] = w; w ^= bb[1];
        Bx[j][75] = w; w ^= bb[0];
        Bx[j][74] = w; w ^= bb[2];
        Bx[j][78] = w; w ^= bb[0];
        Bx[j][79] = w; w ^= bb[1];
        Bx[j][77] = w; w ^= bb[0];
        Bx[j][76] = w; w ^= bb[3];
        Bx[j][68] = w; w ^= bb[0];
        Bx[j][69] = w; w ^= bb[1];
        Bx[j][71] = w; w ^= bb[0];
        Bx[j][70] = w; w ^= bb[2];
        Bx[j][66] = w; w ^= bb[0];
        Bx[j][67] = w; w ^= bb[1];
        Bx[j][65] = w; w ^= bb[0];
        Bx[j][64] = w; w ^= bb[7];
        Bx[j][192] = w; w ^= bb[0];
        Bx[j][193] = w; w ^= bb[1];
        Bx[j][195] = w; w ^= bb[0];
        Bx[j][194] = w; w ^= bb[2];
        Bx[j][198] = w; w ^= bb[0];
        Bx[j][199] = w; w ^= bb[1];
        Bx[j][197] = w; w ^= bb[0];
        Bx[j][196] = w; w ^= bb[3];
        Bx[j][204] = w; w ^= bb[0];
        Bx[j][205] = w; w ^= bb[1];
        Bx[j][207] = w; w ^= bb[0];
        Bx[j][206] = w; w ^= bb[2];
        Bx[j][202] = w; w ^= bb[0];
        Bx[j][203] = w; w ^= bb[1];
        Bx[j][201] = w; w ^= bb[0];
        Bx[j][200] = w; w ^= bb[4];
        Bx[j][216] = w; w ^= bb[0];
        Bx[j][217] = w; w ^= bb[1];
        Bx[j][219] = w; w ^= bb[0];
        Bx[j][218] = w; w ^= bb[2];
        Bx[j][222] = w; w ^= bb[0];
        Bx[j][223] = w; w ^= bb[1];
        Bx[j][221] = w; w ^= bb[0];
        Bx[j][220] = w; w ^= bb[3];
        Bx[j][212] = w; w ^= bb[0];
        Bx[j][213] = w; w ^= bb[1];
        Bx[j][215] = w; w ^= bb[0];
        Bx[j][214] = w; w ^= bb[2];
        Bx[j][210] = w; w ^= bb[0];
        Bx[j][211] = w; w ^= bb[1];
        Bx[j][209] = w; w ^= bb[0];
        Bx[j][208] = w; w ^= bb[5];
        Bx[j][240] = w; w ^= bb[0];
        Bx[j][241] = w; w ^= bb[1];
        Bx[j][243] = w; w ^= bb[0];
        Bx[j][242] = w; w ^= bb[2];
        Bx[j][246] = w; w ^= bb[0];
        Bx[j][247] = w; w ^= bb[1];
        Bx[j][245] = w; w ^= bb[0];
        Bx[j][244] = w; w ^= bb[3];
        Bx[j][252] = w; w ^= bb[0];
        Bx[j][253] = w; w ^= bb[1];
        Bx[j][255] = w; w ^= bb[0];
        Bx[j][254] = w; w ^= bb[2];
        Bx[j][250] = w; w ^= bb[0];
        Bx[j][251] = w; w ^= bb[1];
        Bx[j][249] = w; w ^= bb[0];
        Bx[j][248] = w; w ^= bb[4];
        Bx[j][232] = w; w ^= bb[0];
        Bx[j][233] = w; w ^= bb[1];
        Bx[j][235] = w; w ^= bb[0];
        Bx[j][234] = w; w ^= bb[2];
        Bx[j][238] = w; w ^= bb[0];
        Bx[j][239] = w; w ^= bb[1];
        Bx[j][237] = w; w ^= bb[0];
        Bx[j][236] = w; w ^= bb[3];
        Bx[j][228] = w; w ^= bb[0];
        Bx[j][229] = w; w ^= bb[1];
        Bx[j][231] = w; w ^= bb[0];
        Bx[j][230] = w; w ^= bb[2];
        Bx[j][226] = w; w ^= bb[0];
        Bx[j][227] = w; w ^= bb[1];
        Bx[j][225] = w; w ^= bb[0];
        Bx[j][224] = w; w ^= bb[6];
        Bx[j][160] = w; w ^= bb[0];
        Bx[j][161] = w; w ^= bb[1];
        Bx[j][163] = w; w ^= bb[0];
        Bx[j][162] = w; w ^= bb[2];
        Bx[j][166] = w; w ^= bb[0];
        Bx[j][167] = w; w ^= bb[1];
        Bx[j][165] = w; w ^= bb[0];
        Bx[j][164] = w; w ^= bb[3];
        Bx[j][172] = w; w ^= bb[0];
        Bx[j][173] = w; w ^= bb[1];
        Bx[j][175] = w; w ^= bb[0];
        Bx[j][174] = w; w ^= bb[2];
        Bx[j][170] = w; w ^= bb[0];
        Bx[j][171] = w; w ^= bb[1];
        Bx[j][169] = w; w ^= bb[0];
        Bx[j][168] = w; w ^= bb[4];
        Bx[j][184] = w; w ^= bb[0];
        Bx[j][185] = w; w ^= bb[1];
        Bx[j][187] = w; w ^= bb[0];
        Bx[j][186] = w; w ^= bb[2];
        Bx[j][190] = w; w ^= bb[0];
        Bx[j][191] = w; w ^= bb[1];
        Bx[j][189] = w; w ^= bb[0];
        Bx[j][188] = w; w ^= bb[3];
        Bx[j][180] = w; w ^= bb[0];
        Bx[j][181] = w; w ^= bb[1];
        Bx[j][183] = w; w ^= bb[0];
        Bx[j][182] = w; w ^= bb[2];
        Bx[j][178] = w; w ^= bb[0];
        Bx[j][179] = w; w ^= bb[1];
        Bx[j][177] = w; w ^= bb[0];
        Bx[j][176] = w; w ^= bb[5];
        Bx[j][144] = w; w ^= bb[0];
        Bx[j][145] = w; w ^= bb[1];
        Bx[j][147] = w; w ^= bb[0];
        Bx[j][146] = w; w ^= bb[2];
        Bx[j][150] = w; w ^= bb[0];
        Bx[j][151] = w; w ^= bb[1];
        Bx[j][149] = w; w ^= bb[0];
        Bx[j][148] = w; w ^= bb[3];
        Bx[j][156] = w; w ^= bb[0];
        Bx[j][157] = w; w ^= bb[1];
        Bx[j][159] = w; w ^= bb[0];
        Bx[j][158] = w; w ^= bb[2];
        Bx[j][154] = w; w ^= bb[0];
        Bx[j][155] = w; w ^= bb[1];
        Bx[j][153] = w; w ^= bb[0];
        Bx[j][152] = w; w ^= bb[4];
        Bx[j][136] = w; w ^= bb[0];
        Bx[j][137] = w; w ^= bb[1];
        Bx[j][139] = w; w ^= bb[0];
        Bx[j][138] = w; w ^= bb[2];
        Bx[j][142] = w; w ^= bb[0];
        Bx[j][143] = w; w ^= bb[1];
        Bx[j][141] = w; w ^= bb[0];
        Bx[j][140] = w; w ^= bb[3];
        Bx[j][132] = w; w ^= bb[0];
        Bx[j][133] = w; w ^= bb[1];
        Bx[j][135] = w; w ^= bb[0];
        Bx[j][134] = w; w ^= bb[2];
        Bx[j][130] = w; w ^= bb[0];
        Bx[j][131] = w; w ^= bb[1];
        Bx[j][129] = w; w ^= bb[0];
        Bx[j][128] = w;
    }
    memset(C, 0, m * sizeof(unsigned long));
    for (unsigned long i = 0; i < m; i++) {
        uint64_t aa = A[i];
        C[i] = Bx[0][aa & 255]; aa>>=8;
        C[i]^= Bx[1][aa & 255]; aa>>=8;
        C[i]^= Bx[2][aa & 255]; aa>>=8;
        C[i]^= Bx[3][aa & 255]; aa>>=8;
        C[i]^= Bx[4][aa & 255]; aa>>=8;
        C[i]^= Bx[5][aa & 255]; aa>>=8;
        C[i]^= Bx[6][aa & 255]; aa>>=8;
        C[i]^= Bx[7][aa & 255];
    }
}

static inline void mul_N64_6464_vec(uint64_t *C,
                   const uint64_t *A,
                   const uint64_t *B, unsigned long m)
{

    memset(C, 0, m * sizeof(unsigned long));
    for (unsigned long i = 0; i < m; i++) {
        mul_o64_6464(C++, *A++, B);
    }
}

static inline void mul_N64_T6464_vec(uint64_t *C,
                   const uint64_t *A,
                   const uint64_t *B, unsigned long m)
{

    memset(C, 0, m * sizeof(unsigned long));
    for (unsigned long i = 0; i < m; i++) {
        mul_o64_T6464(C++, *A++, B);
    }
}

static inline void mul_N64_6464_transB(uint64_t *C,
                   const uint64_t *A,
                   const uint64_t *B, unsigned long m)
{
    uint64_t *tb = malloc(64 * sizeof(uint64_t));
    transp_6464(tb, B);
    mul_N64_T6464(C, A, B, m);
    free(tb);
}


static inline void mul_N64_6464_sse(uint64_t *C,
		 const uint64_t *A,
		 const uint64_t *B, unsigned long m)
{
    unsigned long j;
    memset(C, 0, m * sizeof(uint64_t));

    __m128i *Cw = (__m128i *) C;
    __m128i *Aw = (__m128i *) A;

    for (j = 0; j < m; j += 2) {
	__v2di c = { 0, 0 };
	__v2di a = *Aw++;

	__v2di one = { 1, 1, };
#define SHR(x,r) _mm_srli_epi64((x),(r))
	for (int i = 0; i < 64; i++) {
	    __v2di bw = { B[i], B[i], };

	    c ^= (bw & -(a & one));
	    a = SHR(a, 1);
	}
#undef  SHR
	*Cw++ = c;
    }
    C += j;
    A += j;
    for (; j < m; j++) {
	uint64_t c = 0UL;
	uint64_t a = *A++;
	for (int i = 0; i < 64; i++) {
	    c ^= (B[i] & -(a & 1UL));
	    a >>= 1UL;
	}
	*C++ = c;
    }
}

static inline void mul_64N_N64_addmul(uint64_t *r, uint64_t *a, uint64_t *w, unsigned long n)
{
    memset(r, 0, 64 * sizeof(uint64_t));
    for (unsigned long i = 0; i < n; i++) {
        addmul_To64_o64(r, a[i], w[i]);
    }
}


static inline void mul_TN32_N64_C(uint64_t * b, uint32_t * A, uint64_t * x, unsigned int ncol)
{
    uint32_t idx, i, rA;
    uint64_t rx;

    memset(b, 0, 32 * sizeof(uint64_t));
    for(idx = 0; idx < ncol; idx++) {
        rA = A[idx];
        rx = x[idx];
        for(i = 0; i < 32; i++) {
            b[i] ^= rx & -(rA & 1);
            rA >>= 1;
        }
    }
}

static inline void mul_TN64_N64_C(uint64_t * b, uint64_t * A, uint64_t * x, unsigned int ncol)
{
    uint64_t idx, i, rA;
    uint64_t rx;

    memset(b, 0, 64 * sizeof(uint64_t));
    for(idx = 0; idx < ncol; idx++) {
        rA = A[idx];
        rx = x[idx];
        for(i = 0; i < 64; i++) {
            b[i] ^= rx & -(rA & 1);
            rA >>= 1;
        }
    }
}

#if 0   /* haven't checked yet what the funny-named functions actually do... *//*{{{*/
static inline void TVUBit_v2(unsigned long m,
	       unsigned long n,
	       const unsigned long *A, const unsigned long *B,
	       unsigned long *C)
{
    unsigned long i, P, k;
    memset(C, 0, WBITS * sizeof(unsigned long));
    for (i = 0; i < m; i++) {
	P = *A++;
	for (k = 0; k < WBITS; k++) {
	    //if (P & 1UL) C[k] ^= B[i];
	    C[k] ^= (B[i] & -(P & 1UL));
	    P >>= 1UL;
	}
    }
}
#if 1
static inline void VUBit_v2(unsigned long m,
	      unsigned long n,
	      const unsigned long *A, const unsigned long *B,
	      unsigned long *C)
{
    unsigned long i, P, k;
    memset(C, 0, m * sizeof(unsigned long));
    for (i = 0; i < m; i++) {
	P = *A++;
	for (k = 0; k < n; k++) {
	    C[i] ^= (B[k] & -(P & 1UL));
	    P >>= 1UL;
	}
    }
}
#endif
#if 0
static inline void VUBit(unsigned long m,
	      unsigned long n,
	      const unsigned long *A, const unsigned long *B,
	      unsigned long *C)
{
    unsigned long i, P, k, j;
    memset(C, 0, m * sizeof(unsigned long));
    for (i = 0; i < m; i++) {
	P = *A++;
	j = C[i];
	for (k = 0; k < n; k++) {
	    j ^= (B[k] & -(P & 1UL));
	    P >>= 1UL;
	}
	C[i] = j;
    }
}
#endif
#endif/*}}}*/


/* polynomials */

/* lengths of a1 and a2 are n1 and n2 */
typedef uint64_t (*m64pol_ptr)[64];
typedef uint64_t (*const m64pol_srcptr)[64];

void m64pol_addmul(m64pol_ptr r, m64pol_srcptr a1, m64pol_srcptr a2, unsigned int n1, unsigned int n2)
{
    assert(r != a1 && r != a2);
    for(unsigned int i = 0 ; i < n1 ; i++) {
        for(unsigned int j = 0 ; j < n2 ; j++) {
            mat64 x;
            mul_6464_6464(x, a1[i], a2[j]);
            add_6464_6464(r[i+j], r[i+j], x);
        }
    }
}

void m64pol_add(m64pol_ptr r, m64pol_srcptr a1, m64pol_srcptr a2, unsigned int n)
{
    for(unsigned int i = 0 ; i < n ; i++) {
        add_6464_6464(r[i], a1[i], a2[i]);
    }
}

void m64pol_mul(m64pol_ptr r, m64pol_srcptr a1, m64pol_srcptr a2, unsigned int n1, unsigned int n2)
{
    memset(r, 0, (n1 + n2 - 1) * sizeof(mat64));
    m64pol_addmul(r, a1, a2, n1, n2);
}

void m64pol_mul_kara(m64pol_ptr r, m64pol_srcptr a1, m64pol_srcptr a2, unsigned int n1, unsigned int n2)
{
    assert(r != a1 && r != a2);
    assert(n1 == n2);
    assert((n1 & (n1 - 1)) == 0);
    if (n1 == 1) {
        m64pol_mul(r, a1, a2, n1, n2);
        return;
    }
    unsigned int h = n1 >> 1;
    memset(r, 0, (n1 + n2 - 1) * sizeof(mat64));

    m64pol_add(r, a1, a1 + h, h);
    m64pol_add(r + 2 * h, a2, a2 + h, h);

    mat64 * t = (mat64 *) malloc((2*h-1) * sizeof(mat64));
    m64pol_mul_kara(t, r, r + 2 * h, h, h);

    m64pol_mul_kara(r, a1, a2, h, h);
    m64pol_mul_kara(r + 2 * h, a1 + h, a2 + h, h, h);
    m64pol_add(t, t, r, 2 * h - 1);
    m64pol_add(t, t, r + 2 * h, 2 * h - 1);
    m64pol_add(r + h, r + h, t, 2 * h - 1);
    free(t);
}

void m64pol_addmul_kara(m64pol_ptr r, m64pol_srcptr a1, m64pol_srcptr a2, unsigned int n1, unsigned int n2)
{
    mat64 * t = (mat64 *) malloc((n1 + n2 -1) * sizeof(mat64));
    m64pol_mul_kara(t, a1, a2, n1, n2);
    m64pol_add(r, r, t, n1 + n2 - 1);
    free(t);
}


// 2.11 -- 512 512 based on basecase @ 512
// 1.49 -- 512 512 based on basecase @ 256
// 1.07 -- 512 512 based on basecase @ 128
// 0.76 -- 512 512 based on basecase @ 64
// 0.56 -- 512 512 based on basecase @ 32
// 0.41 -- 512 512 based on basecase @ 16
// 0.31 -- 512 512 based on basecase @ 8
// 0.24 -- 512 512 based on basecase @ 4
// 0.19 -- 512 512 based on basecase @ 2
// 0.14 -- 512 512 based on basecase @ 1

#define TIME1(maxtime, what, args) do {			        	\
    clock_t measuring_time = maxtime * CLOCKS_PER_SEC;			\
    clock_t t0, t1;							\
    int j;								\
    t0 = clock();							\
    for (j = 0; ; j++) {						\
        what args;							\
        t1 = clock() - t0;						\
        if (j && t1 > measuring_time)					\
            break;							\
    }									\
    double t = t1;							\
    t /= CLOCKS_PER_SEC;						\
    t /= j;								\
    char * unit = "s";							\
    if (t < 1.0e-7) { unit = "ns"; t *= 1.0e9;			\
    } else if (t < 1.0e-4) { unit = "micros"; t *= 1.0e6;		\
    } else if (t < 1.0e-1) { unit = "ms"; t *= 1.0e3; }		\
    printf(#what " \t%d times in %.4f %s each\n",       		\
            j, t, unit);		                        	\
} while (0)

/* Same spirit, but treat multiplication of 64K by 64K matrices (of
 * polynomials).
 * 
 * We assume that there is no pointer aliasing, and that matrices are
 * stored row-major, with all polynomials contiguous (and of fixed
 * lengths n1, resp n2).
 */
void m64polblock_mul(m64pol_ptr r, m64pol_srcptr a1, m64pol_srcptr a2, unsigned int n1, unsigned int n2, unsigned int K)
{
    assert(r != a1 && r != a2);
    memset(r, 0, (n1 + n2 - 1) * K * K * sizeof(mat64));
    for(unsigned int i = 0 ; i < K ; i++) {
        m64pol_srcptr ra1 = a1 + i * n1 * K;
        m64pol_srcptr rr = r + i * (n1 + n2 - 1) * K;
    for(unsigned int j = 0 ; j < K ; j++) {
        m64pol_srcptr ca2 = a2 + j * n2;
        m64pol_srcptr pr = rr + j * (n1 + n2 - 1);
    for(unsigned int k = 0 ; k < K ; k++) {
        m64pol_srcptr pa1 = ra1 + k * n1;
        m64pol_srcptr pa2 = ca2 + k * n2 * K;
        m64pol_addmul(pr, pa1, pa2, n1, n2);
    }
    }
    }
}

void m64polblock_mul_kara(m64pol_ptr r, m64pol_srcptr a1, m64pol_srcptr a2, unsigned int n1, unsigned int n2, unsigned int K)
{
    assert(r != a1 && r != a2);
    memset(r, 0, (n1 + n2 - 1) * K * K * sizeof(mat64));
    for(unsigned int i = 0 ; i < K ; i++) {
        m64pol_srcptr ra1 = a1 + i * n1 * K;
        m64pol_srcptr rr = r + i * (n1 + n2 - 1) * K;
    for(unsigned int j = 0 ; j < K ; j++) {
        m64pol_srcptr ca2 = a2 + j * n2;
        m64pol_srcptr pr = rr + j * (n1 + n2 - 1);
    for(unsigned int k = 0 ; k < K ; k++) {
        m64pol_srcptr pa1 = ra1 + k * n1;
        m64pol_srcptr pa2 = ca2 + k * n2 * K;
        m64pol_addmul_kara(pr, pa1, pa2, n1, n2);
    }
    }
    }
}

#ifdef  HAVE_M4RI
void my_mzd_randomize(mzd_t * A)
{
    for(size_t i = 0 ; i < A->nrows ; i++) {
        uint64_t * ptr = (uint64_t *) A->rows[i];
        for(size_t j = 0 ; j < A->width ; j++) {
            ptr[j] = rand64();
        }
    }
}

void mzd_set_mem(mzd_t * M, const uint64_t * s, unsigned int n)
{
    assert(M->width == 1);
    assert(M->nrows == n);
    for(size_t i = 0 ; i < M->nrows ; i++) {
        uint64_t * ptr = (uint64_t *) M->rows[i];
        ptr[0] = bitrev(s[i]);
    }
}

void mzd_set_memT(mzd_t * M, const uint64_t * s, unsigned int n)
{
    mzd_set_mem(M, s, n);
    mzd_transpose(M, M);
}

void mzd_check_mem(mzd_t * M, uint64_t * s, unsigned int n)
{
    assert(M->width == 1);
    assert(M->nrows == n);
    for(size_t i = 0 ; i < M->nrows ; i++) {
        uint64_t * ptr = (uint64_t *) M->rows[i];
        if (ptr[0] != bitrev(s[i])) {
            fprintf(stderr, "Rows %zu differ: %016"PRIx64" != %016"PRIx64"\n",
                    i, bitrev(ptr[0]), s[i]);
            abort();
        }
    }
}
void mzd_check_memT(mzd_t * M, uint64_t * s, unsigned int n)
{
    mzd_t * tmp = mzd_transpose(NULL, M);
    mzd_check_mem(M, s, n);
    mzd_free(tmp);
}
#endif  /* HAVE_M4RI */


struct l1_data_s {
    uint64_t * xr;
    uint64_t * r;
    uint64_t * a;
    uint64_t * w;
    uint64_t * wt;
#ifdef  HAVE_M4RI
    mzd_t *R;
    mzd_t *A;
    mzd_t *W;
    mzd_t *WT;
#endif  /* HAVE_M4RI */
    unsigned int n;
};

typedef struct l1_data_s l1_data[1];
typedef struct l1_data_s * l1_data_ptr;
typedef const struct l1_data_s * l1_data_srcptr;

void l1_data_init_set(l1_data_ptr D, unsigned int n)
{
    D->n = n;
    D->xr = (uint64_t *) malloc(n * sizeof(uint64_t));
    D->r = (uint64_t *) malloc(n * sizeof(uint64_t));
    D->a = (uint64_t *) malloc(n * sizeof(uint64_t));
    D->w = (uint64_t *) malloc(64 * sizeof(uint64_t));
    D->wt = (uint64_t *) malloc(64 * sizeof(uint64_t));
#ifdef  HAVE_M4RI
    D->R = mzd_init(n, 64);
    D->A = mzd_init(n, 64);
    D->W = mzd_init(64, 64);
    D->WT = mzd_init(64, 64);
#endif  /* HAVE_M4RI */

    rand64_mem(D->a, n);
    rand64_mem(D->w, 64);
    transp_6464(D->wt, D->w);
#ifdef  HAVE_M4RI
    mzd_set_mem(D->A, D->a, n);
    mzd_set_mem(D->W, D->w, 64);
    mzd_set_mem(D->WT, D->wt, 64);
#endif  /* HAVE_M4RI */
}

void l1_data_clear(l1_data_ptr D)
{
    free(D->xr);
    free(D->r);
    free(D->a);
    free(D->w);
    free(D->wt);
#ifdef  HAVE_M4RI
    mzd_free(D->R);
    mzd_free(D->A);
    mzd_free(D->W);
    mzd_free(D->WT);
#endif  /* HAVE_M4RI */
}

#define SCOPE_L1_DATA_MEMBERS_STANDARD(D)				\
    unsigned int n __attribute__((unused)) = D->n;			\
    uint64_t * xr __attribute__((unused)) = D->xr;      		\
    uint64_t * r __attribute__((unused)) = D->r;			\
    const uint64_t * a __attribute__((unused)) = D->a;			\
    const uint64_t * w __attribute__((unused)) = D->w;			\
    const uint64_t * wt __attribute__((unused)) = D->wt

#define SCOPE_L1_DATA_MEMBERS_M4RI(D)   				\
    mzd_t * R __attribute__((unused)) = D->R;				\
    mzd_t * A __attribute__((unused)) = D->A;				\
    mzd_t * W __attribute__((unused)) = D->W;				\
    mzd_t * WT __attribute__((unused)) = D->WT

#ifdef  HAVE_M4RI
#define SCOPE_L1_DATA_MEMBERS(D)                \
        SCOPE_L1_DATA_MEMBERS_STANDARD(D);      \
        SCOPE_L1_DATA_MEMBERS_M4RI(D)
#else  /* HAVE_M4RI */
#define SCOPE_L1_DATA_MEMBERS(D)                \
        SCOPE_L1_DATA_MEMBERS_STANDARD(D)
#endif  /* HAVE_M4RI */

void level1_basic_tests()
{
    l1_data D;
    l1_data_init_set(D, 64);
    SCOPE_L1_DATA_MEMBERS(D);

    /* transposition */
    transp_6464(r, a); memcpy(xr, r, 64 * sizeof(uint64_t));

    transp_6464(r, a);
    if (memcmp(xr, r, 64 * sizeof(uint64_t))) abort();
    TIME1(1, transp_6464, (r,a));

#ifdef  HAVE_M4RI
    mzd_transpose(R, A);
    mzd_check_mem(R, xr, 64);
    TIME1(1, mzd_transpose, (R, A));
#endif  /* HAVE_M4RI */

    /* copy */
    copy_6464(r, a); memcpy(xr, r, 64 * sizeof(uint64_t));

    copy_6464(r, a);
    if (memcmp(xr, r, 64 * sizeof(uint64_t))) abort();
    TIME1(1, copy_6464, (r,a));

#ifdef  HAVE_M4RI
    mzd_copy(R, A);
    mzd_check_mem(R, xr, 64);
    TIME1(1, mzd_copy, (R, A));
#endif  /* HAVE_M4RI */

    /* add */
    add_6464_6464(r, a, w); memcpy(xr, r, 64 * sizeof(uint64_t));

    add_6464_6464_C(r, a, w);
    if (memcmp(xr, r, 64 * sizeof(uint64_t))) abort();
    TIME1(1, add_6464_6464_C, (r,a,w));

#ifdef  HAVE_M4RI
    mzd_add(R, A, W);
    mzd_check_mem(R, xr, 64);
    TIME1(1, mzd_add, (R, A, W));
#endif  /* HAVE_M4RI */

    /*******/

    l1_data_clear(D);
}

void level1_mul_tests_N_list(l1_data_ptr D)
{
    SCOPE_L1_DATA_MEMBERS(D);

    mul_N64_6464_vec(r, a, w, n);
    if (memcmp(xr, r, n * sizeof(uint64_t))) abort();
    TIME1(2, mul_N64_6464_vec, (r,a,w,n));

    mul_N64_6464_transB(r, a, w, n);
    if (memcmp(xr, r, n * sizeof(uint64_t))) abort();
    TIME1(2, mul_N64_6464_transB, (r,a,w,n));

    mul_N64_6464_sse(r, a, w, n);
    if (memcmp(xr, r, n * sizeof(uint64_t))) abort();
    TIME1(2, mul_N64_6464_sse, (r,a,w,n));

    mul_N64_6464_lookup4(r, a, w, n);
    if (memcmp(xr, r, n * sizeof(uint64_t))) abort();
    TIME1(2, mul_N64_6464_lookup4, (r,a,w,n));

    mul_N64_6464_lookup8(r, a, w, n);
    if (memcmp(xr, r, n * sizeof(uint64_t))) abort();
    TIME1(2, mul_N64_6464_lookup8, (r,a,w,n));

#ifdef  HAVE_M4RI
    mzd_mul_naive(R, A, W);
    mzd_check_mem(R, xr, n);
    TIME1(1, mzd_mul_naive, (R, A, W));

    _mzd_mul_naive(R, A, WT, 1);
    mzd_check_mem(R, xr, n);
    TIME1(1, _mzd_mul_naive, (R, A, WT, 1));

    mzd_mul_m4rm(R, A, W, 0);
    mzd_check_mem(R, xr, n);
    TIME1(1, mzd_mul_m4rm, (R, A, W, 0));

    _mzd_mul_m4rm(R, A, W, 0, 1);
    mzd_check_mem(R, xr, n);
    TIME1(1, _mzd_mul_m4rm, (R, A, W, 0, 1));
#endif  /* HAVE_M4RI */
}


void level1_mul_tests_N(unsigned int n)
{
    l1_data D;
    l1_data_init_set(D, n);

    /* multiplication of a by a matrix */
    mul_N64_6464_vec(D->r, D->a, D->w, D->n);
    memcpy(D->xr, D->r, D->n * sizeof(uint64_t));

    level1_mul_tests_N_list(D);

    l1_data_clear(D);
}

void level1_mul_tests_1_list(l1_data_ptr D)
{
    SCOPE_L1_DATA_MEMBERS(D);
    assert(n == 1);

    mul_o64_6464_C_lsb(r, *a, w);
    if (memcmp(xr, r, n * sizeof(uint64_t))) abort();
    TIME1(1, mul_o64_6464_C_lsb, (r, *a, w));

    mul_o64_6464_C_msb(r, *a, w);
    if (memcmp(xr, r, n * sizeof(uint64_t))) abort();
    TIME1(1, mul_o64_6464_C_msb, (r, *a, w));

    mul_o64_T6464_C_parity(r, *a, wt);
    if (memcmp(xr, r, n * sizeof(uint64_t))) abort();
    TIME1(1, mul_o64_T6464_C_parity, (r, *a, wt));

    mul_o64_T6464_C_parity3(r, *a, wt);
    if (memcmp(xr, r, n * sizeof(uint64_t))) abort();
    TIME1(1, mul_o64_T6464_C_parity3, (r, *a, wt));
}

void level1_mul_tests_1()
{
    unsigned int n = 1;
    l1_data D;
    l1_data_init_set(D, n);

    /* multiplication of a by a matrix */
    mul_o64_6464(D->r, *D->a, D->w);
    memcpy(D->xr, D->r, D->n * sizeof(uint64_t));

    level1_mul_tests_1_list(D);
    /* Functions which can do any n can also do n=1 */
    level1_mul_tests_N_list(D);

    l1_data_clear(D);
}
void level1_mul_tests_64_list(l1_data_ptr D)
{
    SCOPE_L1_DATA_MEMBERS(D);
    assert(n == 64);

    mul_6464_6464_sse(r, a, w);
    if (memcmp(xr, r, n * sizeof(uint64_t))) abort();
    TIME1(1, mul_6464_6464_sse, (r, a, w));

    mul_6464_6464_v2(r, a, w);
    if (memcmp(xr, r, n * sizeof(uint64_t))) abort();
    TIME1(1, mul_6464_6464_v2, (r, a, w));
}

void level1_mul_tests_64()
{
    unsigned int n = 64;
    l1_data D;
    l1_data_init_set(D, n);

    /* multiplicate of two 64x64 matrices */
    mul_6464_6464(D->r, D->a, D->w);
    memcpy(D->xr, D->r, D->n * sizeof(uint64_t));
    
    level1_mul_tests_64_list(D);
    /* Functions which can do any n can also do n=64 */
    level1_mul_tests_N_list(D);

    l1_data_clear(D);
}

int main()
{
  unsigned int n = 2 * 1000 * 1000;

    if (1) {
        uint64_t * r = (uint64_t *) malloc(64 * sizeof(uint64_t));
        uint64_t a = rand64();
        uint64_t w = rand64();

        printf("-- level-1 benches --\n");
        level1_basic_tests();
        level1_mul_tests_1();
        printf("-- level-1 benches, n=64 --\n");
        level1_mul_tests_64();
        printf("-- level-1 benches, misc --\n");
        TIME1(1, addmul_To64_o64_lsb, (r, a, w));
        TIME1(1, addmul_To64_o64_msb, (r, a, w));
        TIME1(1, addmul_To64_o64_lsb_sse_v0, (r, a, w));
        TIME1(1, addmul_To64_o64_lsb_sse_v1, (r, a, w));
        TIME1(1, addmul_To64_o64, (r, a, w));
        free(r);
    }

    if (1) {
        uint64_t * r = (uint64_t *) malloc(64 * sizeof(uint64_t));
        uint64_t * a = (uint64_t *) malloc(n * sizeof(uint64_t));
        uint64_t * w = (uint64_t *) malloc(n * sizeof(uint64_t));
        rand64_mem(a, n);
        rand64_mem(w, n);

        printf("-- level-2 benches (N=%u) --\n", n);
        level1_mul_tests_N(n);
        TIME1(1, mul_64N_N64_addmul, (r,a,w,n));
        TIME1(5, mul_TN32_N64_C, (r,(uint32_t*)a,w,n));
        TIME1(5, mul_TN64_N64_C, (r,a,w,n));

        free(r); free(a); free(w);
    }

    if (0) {
        size_t n = 512;
        mat64 * A = (mat64 *) malloc(n * sizeof(mat64));
        mat64 * B = (mat64 *) malloc(n * sizeof(mat64));
        mat64 * C = (mat64 *) malloc(2 *n * sizeof(mat64));
        printf("-- polynomials (N=%zu) --\n", n);
        TIME1(5, m64pol_mul, (C,A,B,n,n));
        TIME1(5, m64pol_mul_kara, (C,A,B,n,n));
        free(A);
        free(B);
        free(C);
    }

    {
        size_t n = 512;
        mat64 * A = (mat64 *) malloc(2 * 2 * n * sizeof(mat64));
        mat64 * B = (mat64 *) malloc(2 * 2 * n * sizeof(mat64));
        mat64 * C = (mat64 *) malloc(2 * 2 * 2 *n * sizeof(mat64));
        printf("-- polynomials, larger matrices (N=%zu) --\n", n);
        TIME1(5, m64polblock_mul, (C,A,B,n,n,2));
        TIME1(5, m64polblock_mul_kara, (C,A,B,n,n,2));
        free(A);
        free(B);
        free(C);
    }

    return 0;
}
