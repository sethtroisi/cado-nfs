/* An implementation of Cantor's algorithm for multiplication of
  polynomials over GF(2).

  Copyright 2007 Pierrick Gaudry.

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

/* 20080122 -- shortened the input file somewhat. Old disabled code has
 * been deleted in this operation.
 *
 *
 * XXX The inverse transform operation should be unified between the
 * truncated and non-truncated case.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#include "cantor128.h"
#include "mpfq_2_128.h"

/* The following flags affect the behaviour of the program */

#define xxxCANTOR_GM            /* Use Gao-Mateer recursion */
#define xxxCANTOR_GM_TRUNCATE   /* Use truncated variant */
#define xxxCOUNT_MULTS          /* Count multiplications */

#define Kelt            mpfq_2_128_elt
#define Kdst_elt        mpfq_2_128_dst_elt
#define Ksrc_elt        mpfq_2_128_src_elt

int mulcount=0;
Kelt * fbase;
long int * findex;

#define Kadd(a0,a1,a2)  mpfq_2_128_add (NULL,a0,a1,a2)
#if 0
#define Kmul(a0,a1,a2)  do { \
    mpfq_2_128_mul (NULL,a0,a1,a2);     \
    mpfq_2_128_mul (NULL,a0,a1,a2);     \
} while (0)
#endif

#ifdef  COUNT_MULTS
#define Kmul(a0,a1,a2)  do { mpfq_2_128_mul (NULL,a0,a1,a2); mulcount++; } while (0)
#else
#define Kmul(a0,a1,a2)  mpfq_2_128_mul (NULL,a0,a1,a2)
#endif

#define Kset_ui(a0,a1)  mpfq_2_128_set_ui (NULL,a0,a1)
#define Kcmp(a0,a1)  mpfq_2_128_cmp (NULL,a0,a1)
#define Kset(a0,a1)  mpfq_2_128_set (NULL,a0,a1)
#define Ksqr(a0,a1)     mpfq_2_128_sqr (NULL,a0,a1)

// Some constants related to Cantor's algorithm
// Beta_i are such that Beta_{i-1} = Beta_i^2 + Beta_i
// and { Beta_i : i in [1..128] } is a GF(2)-basis of GF(2^128).
// The following is (the begining of) a solution.
// NB: the defining polynomial for GF(2^128) is x^128 + x^7 + x^2 + x + 1

#if GMP_LIMB_BITS == 32
#define BETA(x,y,z,t) { x, y, z, t, }
#else
#define BETA(x,y,z,t) { y << 32 | x, t << 32 | z }
#endif
static const Kelt Betai[32] = {
    BETA(1UL,0UL,0UL,0UL),
    BETA(2763131656UL, 1735044255UL, 4101184249UL, 693813425UL),
    BETA(465473314UL, 4279278019UL, 3443545289UL, 2719171618UL),
    BETA(45213768UL, 1404590948UL, 3790498114UL, 2035930892UL),
    BETA(1766741572UL, 3401373924UL, 3912612066UL, 150385427UL),
    BETA(1642619116UL, 1296609638UL, 891822488UL, 1237922326UL),
    BETA(1432407564UL, 2705320108UL, 2845284892UL, 87917909UL),
    BETA(3169845528UL, 1555107770UL, 977112198UL, 1297233226UL),
    BETA(3766302226UL, 3784070584UL, 2372619472UL, 751494893UL),
    BETA(3319219312UL, 1485529668UL, 1919364198UL, 554536894UL),
    BETA(2221231524UL, 3494370010UL, 2246133939UL, 3255900204UL),
    BETA(4004748792UL, 4275287155UL, 1207038226UL, 273485854UL),
    BETA(4134924972UL, 3206127284UL, 3660212267UL, 3465895584UL),
    BETA(1810844824UL, 432364429UL, 3275102095UL, 3685546794UL),
    BETA(2239846096UL, 4115864360UL, 3774983060UL, 1484196804UL),
    BETA(572310756UL, 2341975643UL, 1205907074UL, 1233095387UL),
    BETA(374654328UL, 2701982826UL, 4048512218UL, 95057245UL),
    BETA(1069419604UL, 3135871566UL, 603196539UL, 2195404101UL),
    BETA(2819376738UL, 2390789826UL, 3206934123UL, 3126567229UL),
    BETA(2667030844UL, 1546355395UL, 894286960UL, 806061570UL),
    BETA(3589162292UL, 2284726595UL, 3987361705UL, 3261030005UL),
    BETA(488495394UL, 2513641547UL, 3724856370UL, 956902965UL),
    BETA(894705622UL, 3412990576UL, 2985357841UL, 2327012336UL),
    BETA(691430644UL, 2704377464UL, 1647445035UL, 3189688948UL),
    BETA(3013060094UL, 3867471049UL, 66167789UL, 3014841169UL),
    BETA(1101590000UL, 1839657742UL, 3478054432UL, 884001900UL),
    BETA(494046998UL, 2899707219UL, 2156601667UL, 2796513881UL),
    BETA(3398486548UL, 57111011UL, 1428951442UL, 894796186UL),
    BETA(3950556970UL, 3334574894UL, 3366886047UL, 3938361620UL),
    BETA(999564750UL, 4233655517UL, 706258300UL, 428364188UL),
    BETA(1913007172UL, 1894536580UL, 871798750UL, 1621883580UL),
    BETA(2252056148UL, 2850320680UL, 3455532848UL, 234262497UL),
};

// The S_i polynomials of Cantor's algorithm: linearized polynomials such
// that S_i | S_{i+1}.
static const long ind_S[31][15] = {
    {},                         // S0
    {1},                        // S1
    {1},                        // S2
    {4, 2, 1},                  // S3
    {1},                        // S4
    {16, 2, 1},                 // S5
    {16, 4, 1},                 // S6
    {64, 32, 16, 8, 4, 2, 1},   // S7
    {1},                        // S8
    {256, 2, 1},                // S9
    {256, 4, 1},                // S10
    {1024, 512, 256, 8, 4, 2, 1},       // S11
    {256, 16, 1},               // S12
    {4096, 512, 256, 32, 16, 2, 1},     // S13
    {4096, 1024, 256, 64, 16, 4, 1},    // S14
    {16384, 8192, 4096, 2048, 1024, 512, 256, 128, 64, 32, 16, 8, 4, 2, 1},     // S15
    {1},                        // S16
    {65536, 2, 1},              // S17
    {65536, 4, 1},              // S18
    {262144, 131072, 65536, 8, 4, 2, 1},        // S19
    {65536, 16, 1},             // S20
    {1048576, 131072, 65536, 32, 16, 2, 1},     // S21
    {1048576, 262144, 65536, 64, 16, 4, 1},     // S22
    {4194304, 2097152, 1048576, 524288, 262144, 131072, 65536, 128, 64, 32, 16, 8, 4, 2, 1},    // S23
    {65536, 256, 1},            // S24
    {16777216, 131072, 65536, 512, 256, 2, 1},  // S25
    {16777216, 262144, 65536, 1024, 256, 4, 1}, // S26
    {67108864, 33554432, 16777216, 524288, 262144, 131072, 65536, 2048, 1024, 512, 256, 8, 4, 2, 1},    // S27
    {16777216, 1048576, 65536, 4096, 256, 16, 1},       // S28
    {268435456, 33554432, 16777216, 2097152, 1048576, 131072, 65536, 8192, 4096, 512, 256, 32, 16, 2, 1},       // S29
    {268435456, 67108864, 16777216, 4194304, 1048576, 262144, 65536, 16384, 4096, 1024, 256, 64, 16, 4, 1},     // S30
};

// The number of terms in S_i.
static const int ind_number[31] =
    { 0, 1, 1, 3, 1, 3, 3, 7, 1, 3, 3, 7, 3, 7, 7, 15,
    1, 3, 3, 7, 3, 7, 7, 15, 3, 7, 7, 15, 7, 15, 15
};

// At index i=i0+2*i1+4*i2+8*i3..., one computes 
//     allBetai[i] = i0*Betai[1] + i1*Betai[2] + i2*Betai[3] + ...
static inline void allBetai(Kdst_elt x, int i)
{
    int j;
    Kset_ui(x,0);
    j = 0;
    for( ; i ; i>>=1, j++) {
        if (!(i & 1)) continue;
        Kadd(x, x, Betai[j]);
    }
}

#ifdef  CANTOR_GM
#if 0
/* These versions would make sense if Kelt were defined to be _v2di *,
 * which it isn't */
static inline void allBetai1(Kdst_elt x, int i)
{
    int j;
    Kset_ui(x,0);
    j = 0;
    __v2di m[2] = { (__v2di) {0,0}, (__v2di) { -1L, -1L } };
    for( ; i ; i>>=1, j++) {
            * (__v2di*)  x ^= m[i&1] & * (__v2di*) Betai[j];
    }
}

// compute omega_i and omega_{2i}
static inline void allBetai2(Kdst_elt x, Kdst_elt y, int i)
{
    int j;
    Kset_ui(x,0);
    Kset_ui(y,0);
    j = 0;
    __v2di m[2] = { (__v2di) {0,0}, (__v2di) { -1L, -1L } };
    for( ; i ; i>>=1, j++) {
            * (__v2di*)  x ^= m[i&1] & * (__v2di*) Betai[j];
            * (__v2di*)  y ^= m[i&1] & * (__v2di*) Betai[j+1];
    }
}
static inline void allBetai4(Kdst_elt x, Kdst_elt y, Kdst_elt z, Kdst_elt t, int i)
{
    int j;
    Kset_ui(x,0);
    Kset_ui(y,0);
    Kset_ui(z,0);
    Kset_ui(t,0);
    j = 0;
    __v2di m[2] = { (__v2di) {0,0}, (__v2di) { -1L, -1L } };
    for( ; i ; i>>=1, j++) {
            * (__v2di*) x ^= m[i&1] & * (__v2di*) Betai[j];
            * (__v2di*) y ^= m[i&1] & * (__v2di*) Betai[j+1];
            * (__v2di*) z ^= m[i&1] & * (__v2di*) Betai[j+2];
            * (__v2di*) t ^= m[i&1] & * (__v2di*) Betai[j+3];
    }
}
#else
static inline void allBetai1(Kdst_elt x, int i)
{
    int j;
    Kset_ui(x,0);
    j = 0;
    for( ; i ; i>>=1, j++) {
        if (!(i & 1)) continue;
        Kadd(x, x, Betai[j]);
    }
}
static inline void allBetai2(Kdst_elt x, Kdst_elt y, int i)
{
    int j;
    Kset_ui(x,0);
    Kset_ui(y,0);
    j = 0;
    for( ; i ; i>>=1, j++) {
        if (!(i & 1)) continue;
        Kadd(x, x, Betai[j]);
        Kadd(y, y, Betai[j+1]);
    }
}
static inline void allBetai4(Kdst_elt x, Kdst_elt y, Kdst_elt z, Kdst_elt t, int i)
{
    int j;
    Kset_ui(x,0);
    Kset_ui(y,0);
    Kset_ui(z,0);
    Kset_ui(t,0);
    j = 0;
    for( ; i ; i>>=1, j++) {
        if (!(i & 1)) continue;
        Kadd(x, x, Betai[j]);
        Kadd(y, y, Betai[j+1]);
        Kadd(z, z, Betai[j+2]);
        Kadd(t, t, Betai[j+3]);
    }
}
#endif  /* v2di or not v2di */
#endif  /* CANTOR_GM */

#ifdef  CANTOR_GM
static inline void expand(Kelt * f, int t, int t0)
{
    // f has at most 2^t coefficients. Split in pieces of 2^t0
    // coefficients.
    if (t == t0) return;
    long K = 1L << (t-1);
    long K0 = 1L << (t-1-t0);
    for (long i = 2 * K - 1; i >= K; --i) {
        Kadd(f[i+K0-K], f[i+K0-K], f[i]); 
    }
    expand(f, t-1, t0);
    expand(f+K, t-1, t0);
}
static inline void expand_trunc(Kelt * f, int t, int t0, long n)
{
    // f has at most n coefficients. Split in pieces of 2^t0
    // coefficients.
    long K;
    for( ; n <= (K = (1L << (t-1))) ; t--);
    if (t <= t0) return;
    long K0 = 1L << (t-1-t0);
    for (long i = n - 1; i >= K; --i) {
        Kadd(f[i+K0-K], f[i+K0-K], f[i]); 
    }
    expand(f, t-1, t0);
    expand_trunc(f+K, t-1, t0, n-K);
}

#if 0
static inline void expand_unroll(Kelt * f, int t, int t0)
{
    // f has at most 2^t coefficients. Split in pieces of 2^t0
    // coefficients.
    if (t == t0) return;
    long K = 1L << (t-1);
    long K0 = 1L << (t-1-t0);
    for (long i = 2 * K - 1; i >= K; --i) {
        Kadd(f[i+K0-K], f[i+K0-K], f[i]); 
        fprintf(stderr, "Kadd(f[%ld],f[%ld],f[%ld]);\n",
                findex[i+K0-K+f-fbase],
                findex[i+K0-K+f-fbase],
                findex[i+f-fbase]);
    }
    expand_unroll(f, t-1, t0);
    expand_unroll(f+K, t-1, t0);
}
#endif

void transpose_inplace(Kelt * f, int t)
{
#define Kswap(a,b) do { Kelt x; Kset(x,a); Kset(a,b); Kset(b,x); } while (0)
    // transpose
    for(long i = 0 ; i < (1L<<t) ; i++) {
        for(long j = i+1 ; j < (1L<<t) ; j++) {
            Kswap(f[(i<<t)+j],f[(j<<t)+i]);
        }
    }
}

// assuming f has 1L<<t1 rows of 1L<<t2 values, transpose into g. Set striding
// in g to be 1L<<t1 values, for 1L<<t2 rows.
//
// doing this in place would be somewhat tricky, but maybe reachable
// assuming that we stick to the case of powers of 2 for the dimensions.
void transpose_outofplace(Kelt * g, Kelt * f, int t1, int t2)
{
    for(long j = 0 ; j < (1L<<t2) ; j++) {
        for(long i = 0 ; i < (1L<<t1) ; i++) {
            Kset(g[i],f[i<<t2]);
        }
        g += 1L << t1;
        f++;
    }
}

#if 0
void gm_trick_unroll(int two_t, Kelt * f, int j, int m, int d)
{
    assert((two_t & (two_t-1)) == 0);
    if (two_t == 1) {
        if (j) {
            Kelt z;
            allBetai(z,2*j);
            Kmul(z,z,f[1]);
            Kadd(f[0],z,f[0]);
            fprintf(stderr, "allBetai(z,%d*j+%d);\n",2*m,2*d);
            fprintf(stderr, "Kmul(z,z,f[%ld]);\n",findex[1 + f - fbase]);
            fprintf(stderr, "Kadd(f[%ld],z,f[%ld]);\n",findex[f - fbase], findex[f-fbase]);
        }
        Kadd(f[1],f[1],f[0]);
        fprintf(stderr, "Kadd(f[%ld],f[%ld],f[%ld]);\n",findex[f-fbase+1],findex[f-fbase+1], findex[f-fbase]);
        return;
    }

    long t = two_t >> 1;
    long tau = 1L << t;
    expand_unroll(f,2*t,t);

    for(long i = 0 ; i < tau ; i++) {
        for(long j = i+1 ; j < tau ; j++) {
            Kswap(f[(i<<t)+j],f[(j<<t)+i]);
#if 1
            long idx = findex[(i<<t)+j+f-fbase];
            findex[(i<<t)+j+f-fbase] = findex[(j<<t)+i+f-fbase];
            findex[(j<<t)+i+f-fbase] = idx;
#else
            fprintf(stderr, "Kswap(f[%ld],f[%ld]);\n",
                    (i<<t)+j+f-fbase,
                    (j<<t)+i+f-fbase);
#endif
        }
    }
    
    // evaluate ; we can use f as a buffer, now.
    for(long l = 0 ; l < tau ; l++) {
        gm_trick_unroll(t, f + (l << t), j, m, d);
    }

    for(long i = 0 ; i < tau ; i++) {
        for(long j = i+1 ; j < tau ; j++) {
            Kswap(f[(i<<t)+j],f[(j<<t)+i]);
#if 1
            long idx = findex[(i<<t)+j+f-fbase];
            findex[(i<<t)+j+f-fbase] = findex[(j<<t)+i+f-fbase];
            findex[(j<<t)+i+f-fbase] = idx;
#else
            fprintf(stderr, "Kswap(f[%ld],f[%ld]);\n",
                    (i<<t)+j+f-fbase,
                    (j<<t)+i+f-fbase);
#endif
        }
    }
    
    // evaluate again
    for(long i = 0 ; i < tau ; i++) {
        gm_trick_unroll(t, f + (i << t), (j << t) + i, m << t, (d << t) + i);
    }
}
#endif

static inline void gm_trick2(Kdst_elt f0, Kdst_elt f1, Kdst_elt f2, Kdst_elt f3, Ksrc_elt w2, Ksrc_elt w4)
{
    Kelt z;
    // Kelt w2,w4;
    Kadd(f2,f2,f3);
    Kadd(f1,f1,f2);
    // allBetai2(w2,w4,2*j);
    Kmul(z,w2,f2); Kadd(f0,z,f0); Kadd(f2,f2,f0);
    Kmul(z,w2,f3); Kadd(f1,z,f1); Kadd(f3,f3,f1);
    Kmul(z,w4,f1); Kadd(f0,z,f0); Kadd(f1,f1,f0);
    Kadd(z,w4,Betai[1]);Kmul(z,z,f3); Kadd(f2,z,f2); Kadd(f3,f3,f2);
}
static inline void gm_trick2z(Kdst_elt f0, Kdst_elt f1, Kdst_elt f2, Kdst_elt f3)
{
    Kelt z;
    Kadd(z,f2,f3);                      
    Kadd(f3,f1,f2);                     
    Kadd(f2,f0,z);                      
    Kadd(f1,f1,f2);                     
    Kmul(z,Betai[1],f3);
    Kadd(f2,z,f2);
    Kadd(f3,f3,f2);
}

void gm_trick4(Kelt * f, long j)
{
    Kelt f0; Kelt f1; Kelt f2; Kelt f3; 
    Kelt f4; Kelt f5; Kelt f6; Kelt f7; 
    Kelt f8; Kelt f9; Kelt fa; Kelt fb;
    Kelt fc; Kelt fd; Kelt fe; Kelt ff;

    Kset(f0,f[0]);	Kset(f1,f[1]);	Kset(f2,f[2]);	Kset(f3,f[3]);
    Kset(f4,f[4]);	Kset(f5,f[5]);	Kset(f6,f[6]);	Kset(f7,f[7]);
    Kset(f8,f[8]);	Kset(f9,f[9]);	Kset(fa,f[10]); Kset(fb,f[11]);
    Kset(fc,f[12]);	Kset(fd,f[13]); Kset(fe,f[14]); Kset(ff,f[15]);

    Kadd(f9,f9,ff);	Kadd(f8,f8,fe); Kadd(f7,f7,fd);
    Kadd(f6,f6,fc);	Kadd(f5,f5,fb); Kadd(f4,f4,fa);
    Kadd(f3,f3,f9);	Kadd(f2,f2,f8);	Kadd(f4,f4,f7);
    Kadd(f3,f3,f6);	Kadd(f2,f2,f5);	Kadd(f1,f1,f4);
    Kadd(fc,fc,ff);	Kadd(fb,fb,fe);	Kadd(fa,fa,fd);
    Kadd(f9,f9,fc);

    Kelt w2;
    Kelt w4;
    Kelt w8;
    Kelt w16;
    allBetai4(w2,w4,w8,w16,j<<1);

    gm_trick2(f0,f4,f8,fc,w2,w4);
    gm_trick2(f1,f5,f9,fd,w2,w4);
    gm_trick2(f2,f6,fa,fe,w2,w4);
    gm_trick2(f3,f7,fb,ff,w2,w4);
    gm_trick2(f0,f1,f2,f3,w8,w16);
    Kadd(w2,w8,Betai[1]);
    Kadd(w4,w16,Betai[2]);
    gm_trick2(f4,f5,f6,f7,w2,w4);
    Kadd(w8,w8,Betai[2]);
    Kadd(w16,w16,Betai[3]);
    gm_trick2(f8,f9,fa,fb,w8,w16);
    Kadd(w2,w2,Betai[2]);
    Kadd(w4,w4,Betai[3]);
    gm_trick2(fc,fd,fe,ff,w2,w4);
    Kset(f[0],f0);	Kset(f[1],f1);	Kset(f[2],f2);	Kset(f[3],f3);
    Kset(f[4],f4);	Kset(f[5],f5);	Kset(f[6],f6);	Kset(f[7],f7);
    Kset(f[8],f8);	Kset(f[9],f9);	Kset(f[10],fa); Kset(f[11],fb);
    Kset(f[12],fc);	Kset(f[13],fd); Kset(f[14],fe); Kset(f[15],ff);
}
void gm_trick4z(Kelt * f)
{
    Kelt f0; Kelt f1; Kelt f2; Kelt f3; 
    Kelt f4; Kelt f5; Kelt f6; Kelt f7; 
    Kelt f8; Kelt f9; Kelt fa; Kelt fb;
    Kelt fc; Kelt fd; Kelt fe; Kelt ff;

    Kset(f0,f[0]);	Kset(f1,f[1]);	Kset(f2,f[2]);	Kset(f3,f[3]);
    Kset(f4,f[4]);	Kset(f5,f[5]);	Kset(f6,f[6]);	Kset(f7,f[7]);
    Kset(f8,f[8]);	Kset(f9,f[9]);	Kset(fa,f[10]); Kset(fb,f[11]);
    Kset(fc,f[12]);	Kset(fd,f[13]); Kset(fe,f[14]); Kset(ff,f[15]);

    Kadd(f9,f9,ff);	Kadd(f8,f8,fe); Kadd(f7,f7,fd);
    Kadd(f6,f6,fc);	Kadd(f5,f5,fb); Kadd(f4,f4,fa);
    Kadd(f3,f3,f9);	Kadd(f2,f2,f8);	Kadd(f4,f4,f7);
    Kadd(f3,f3,f6);	Kadd(f2,f2,f5);	Kadd(f1,f1,f4);
    Kadd(fc,fc,ff);	Kadd(fb,fb,fe);	Kadd(fa,fa,fd);
    Kadd(f9,f9,fc);

    Kelt w2;
    Kelt w4;

    gm_trick2z(f0,f4,f8,fc);
    gm_trick2z(f1,f5,f9,fd);
    gm_trick2z(f2,f6,fa,fe);
    gm_trick2z(f3,f7,fb,ff);
    gm_trick2z(f0,f1,f2,f3);
    gm_trick2(f4,f5,f6,f7,Betai[1],Betai[2]);
    gm_trick2(f8,f9,fa,fb,Betai[2],Betai[3]);
    Kadd(w2,Betai[1],Betai[2]);
    Kadd(w4,Betai[2],Betai[3]);
    gm_trick2(fc,fd,fe,ff,w2,w4);

    Kset(f[0],f0);	Kset(f[1],f1);	Kset(f[2],f2);	Kset(f[3],f3);
    Kset(f[4],f4);	Kset(f[5],f5);	Kset(f[6],f6);	Kset(f[7],f7);
    Kset(f[8],f8);	Kset(f[9],f9);	Kset(f[10],fa); Kset(f[11],fb);
    Kset(f[12],fc);	Kset(f[13],fd); Kset(f[14],fe); Kset(f[15],ff);
}


void gm_trick(int two_t, Kelt * f, long j)
{
    assert((two_t & (two_t-1)) == 0);
    if (two_t == 4) {
        if (j) gm_trick4(f,j); else gm_trick4z(f);
        return;
    } else if (two_t == 2) {
        if (j) {
            Kelt w2,w4;
            allBetai2(w2,w4,2*j);
            gm_trick2(f[0],f[1],f[2],f[3],w2,w4);
        } else {
            gm_trick2z(f[0],f[1],f[2],f[3]);
        }
        return;
    } else if (two_t == 1) {
        if (j) {
            Kelt z;
            allBetai(z,2*j);
            Kmul(z,z,f[1]);
            Kadd(f[0],z,f[0]);
        }
        Kadd(f[1],f[1],f[0]);
        return;
    }
#if 0
    if (two_t == 4) {
        fprintf(stderr, "// Unrolling gm_trick(%d,f,%d)\n",two_t, j);
        fbase = f;
        findex = malloc(sizeof(long int) << two_t);
        for(int i = 0 ; i < (1L << two_t) ; i++) {
            findex[i] = i;
        }
        gm_trick_unroll(two_t, f, j, 1, 0);
        free(findex);
        findex = NULL;
        fbase = NULL;
        fprintf(stderr, "// Done unrolling\n\n");
        return;
    }
#endif

    int t = two_t >> 1;
    long tau = 1L << t;

    expand(f,2*t,t);

    transpose_inplace(f,t);
    
    // evaluate ; we can use f as a buffer, now.
    for(long l = 0 ; l < tau ; l++) {
        gm_trick(t, f + (l << t), j);
    }
    transpose_inplace(f,t);

    // evaluate again
    for(long i = 0 ; i < tau ; i++) {
        gm_trick(t, f + (i << t), (j << t) + i);
    }
}
#endif

#ifdef  CANTOR_GM_TRUNCATE
// compute only n' = 1L<<k values -- f must not have more than n coefficients !
//
// n' is larger than or equal to n.
//
// f and buf must hold space for 1L << k coefficients.
//
// The trailing values in f must be equal to zero.
void gm_trick_trunc(int two_t, Kelt * f, Kelt * buf, long j, int k, long n)
{
    assert((two_t & (two_t-1)) == 0);

    if (n <= 1)
        return;
    if (two_t == 1) {
        if (j) {
            Kelt z;
            allBetai(z,2*j);
            Kmul(z,z,f[1]);
            Kadd(f[0],z,f[0]);
        }
        Kadd(f[1],f[1],f[0]);
        return;
    }
    assert(n <= (1L << k));
    if (k <= (two_t>>1)) {
        gm_trick_trunc(two_t>>1, f, buf, j, k, n);
        return;
    }
    assert(k > (two_t>>1));
    assert(k <= two_t);
    assert(n > (1L << (two_t>>1)));
    assert(n <= (1L << two_t));

    int t = two_t >> 1;
    long tau = 1L << t;
    // long cn = (n+tau-1)>>t;      // Ceiling(n/tau);
    long cn = 1L << (k-t);

    expand_trunc(f,2*t,t,n);
    memset(f + n, 0, ((1L<<k)-n) * sizeof(Kelt));

    // need 1L<<k in buf
    transpose_outofplace(buf,f,k-t,t);
    
    // evaluate ; we can use f as a buffer, now.
    Kelt * h = buf;
    for(long l = 0 ; l < tau ; l++) {
        // we're going to need 1L << (k-t) in f, no big deal.
        gm_trick_trunc(t, h, f, j, k-t, cn);
        h += cn;
    }
    transpose_outofplace(f,buf,t,k-t);

    // evaluate again. Here we cannot truncate.
    for(long i = 0 ; i < cn ; i++) {
        gm_trick(t, f + (i << t), (j << t) + i);
    }
}
#endif


/* 20080122 -- hard-coded versions of reduceSi for small sizes have been
 * tried, but do not seem to pay off. Look for the svn data to
 * investigate again the potential benefit.
 */

// The generic reduction function.
static inline void reduceSi(int k, Kelt * f, Kelt beta)
{
    long i, j, K;
    Kelt coeff;
    Ksrc_elt fi;

    K = 1L << k;
    // put rest mod Sk+beta in low part of f and quo in high part.
    for (i = 2 * K - 1; i >= K; --i) {
        fi = f[i];
        for (j = 0; j < ind_number[k]; ++j) {
            int index = ind_S[k][j] + i - K;
            Kadd(f[index], f[index], fi);
        }
        Kmul(coeff, fi, beta);
        Kadd(f[i - K], f[i - K], coeff);
    }

    // The other rest is just the sum of low and hi parts
    for (i = 0; i < K; ++i)
        Kadd(f[i + K], f[i + K], f[i]);
}

#ifdef  CANTOR_GM
// same as multieval, but only to chunks of size kappa
void reduce_top_cantor(int k, int kappa, Kelt * f, long j)
{
    if (k == kappa)
        return;
    // fprintf(stderr, "top-reduce %d %d\n",k,kappa);
    j *= 2;
    Kelt x;
    allBetai(x, j);
    reduceSi(k-1, f, x);
    reduce_top_cantor(k-1,kappa,f,j);
    reduce_top_cantor(k-1,kappa,f+(1L<<(k-1)),j+1);
}
#endif


// The generic reduction function.
// truncated version
static inline void reduceSiNobeta_trunc(int k, Kelt * f, long L)
{
    long j;
    long i, K;
    Ksrc_elt fi;

    K = 1L << k;
    // put rest mod Sk+beta in low part of f and quo in high part.
    for (i = L - 1; i >= K; --i) {
        fi = f[i];
        for (j = 0; j < ind_number[k]; ++j) {
            long index = ind_S[k][j] + i - K;
            Kadd(f[index], f[index], fi);
        }
    }

    // The other rest is just the sum of low and hi parts
    for (i = 0; i < K; ++i)
        Kadd(f[i + K], f[i + K], f[i]);
}

// A generic interpolation function.
// TODO: unroll small sizes, like for reduction.
static inline void interpolateSi(int k, Kelt * f, Kelt beta)
{
    long j;
    long i, K = 1L << k;

    for (i = 0; i < K; ++i)
        Kadd(f[i + K], f[i + K], f[i]);

    for (i = K; i < K << 1; ++i) {
        Kelt coeff;
        Ksrc_elt fi;
        fi = f[i];
        for (j = 0; j < ind_number[k]; ++j) {
            long index = ind_S[k][j] + i - K;
            Kadd(f[index], f[index], fi);
        }
        Kmul(coeff, fi, beta);
        Kadd(f[i - K], f[i - K], coeff);
    }
}

static inline void interpolateSiNobeta(int k, Kelt * f)
{
    long j;
    long i, K = 1L << k;

    for (i = 0; i < K; ++i)
        Kadd(f[i + K], f[i + K], f[i]);

    for (i = K; i < K << 1; ++i) {
        Ksrc_elt fi;
        fi = f[i];
        for (j = 0; j < ind_number[k]; ++j) {
            long index = ind_S[k][j] + i - K;
            Kadd(f[index], f[index], fi);
        }
    }
}


void multievaluateKrec(Kelt * f, int i, long rep_beta)
{
    Kelt beta;
    rep_beta <<= 1;
    allBetai(beta, rep_beta);
    reduceSi(i - 1, f, beta);

    if (i == 2) {
        rep_beta <<= 1;
        allBetai(beta, rep_beta);
        reduceSi(0, f, beta);
        Kadd(beta, beta, Betai[1]);
        reduceSi(0, f + 2, beta);
        return;
    }
    if (i == 3) {
        rep_beta <<= 1;
        allBetai(beta, rep_beta);
        reduceSi(1, f, beta);
        Kadd(beta, beta, Betai[1]);
        reduceSi(1, f + 4, beta);

        rep_beta <<= 1;
        allBetai(beta, rep_beta);
        reduceSi(0, f, beta);
        Kadd(beta, beta, Betai[1]);
        reduceSi(0, f + 2, beta);
        Kadd(beta, beta, Betai[2]);
        reduceSi(0, f + 6, beta);
        Kadd(beta, beta, Betai[1]);
        reduceSi(0, f + 4, beta);
        return;
    }
    if (i > 1) {
        multievaluateKrec(f, i - 1, rep_beta);
        multievaluateKrec(f + (1L << (i - 1)), i - 1, rep_beta + 1);
    }
}

#ifdef  CANTOR_GM
// f must have 2^k coeffs exactly
void multievaluateGM(Kelt * f, int k, long length MAYBE_UNUSED)
{
    int t = 1;
    int two_t;
    for( ; (two_t = t << 1) < k ; t = two_t) ;
    assert(t < k && k <= two_t);
#ifdef CANTOR_GM_TRUNCATE
    Kelt * buf = malloc(sizeof(Kelt) << k);
    gm_trick_trunc(two_t,f,buf,0,k,length);
    free(buf);
#else
    reduce_top_cantor(k,t,f,0);
    long K = 1L << k;
    long tau = 1L << t;
    long j = 0;
    for(long u = 0 ; u < K ; u += tau) {
        gm_trick(t, f + u, j);
        j++;
    }
#endif
}
#endif

void multievaluateKrec_trunc(Kelt * f, int i, long rep_beta, long length)
{
    Kelt beta;
    rep_beta <<= 1;
    allBetai(beta, rep_beta);
    reduceSi(i - 1, f, beta);

    if (i == 2) {
        rep_beta <<= 1;
        allBetai(beta, rep_beta);
        reduceSi(0, f, beta);
        Kadd(beta, beta, Betai[1]);
        reduceSi(0, f + 2, beta);
        return;
    }
    if (i == 3) {
        rep_beta <<= 1;
        allBetai(beta, rep_beta);
        reduceSi(1, f, beta);
        Kadd(beta, beta, Betai[1]);
        reduceSi(1, f + 4, beta);

        rep_beta <<= 1;
        allBetai(beta, rep_beta);
        reduceSi(0, f, beta);
        Kadd(beta, beta, Betai[1]);
        reduceSi(0, f + 2, beta);
        Kadd(beta, beta, Betai[2]);
        reduceSi(0, f + 6, beta);
        Kadd(beta, beta, Betai[1]);
        reduceSi(0, f + 4, beta);
        return;
    }
    if (i > 1) {
        if (length < (1L << (i - 1))) {
            multievaluateKrec_trunc(f, i - 1, rep_beta, length);
        } else if (length > (1L << (i - 1))) {
            multievaluateKrec(f, i - 1, rep_beta);
            multievaluateKrec_trunc(f + (1L << (i - 1)), i - 1, rep_beta + 1,
                                    length - (1L << (i - 1)));
        } else { // length == 1L<<(i-1)
            multievaluateKrec(f, i - 1, rep_beta);
        }
    }
}


void multievaluateKnew_trunc(Kelt * f, int k, long length)
{
    reduceSiNobeta_trunc(k - 1, f, length);

    multievaluateKrec(f, k - 1, 0);
    multievaluateKrec_trunc(f + (1L << (k - 1)), k - 1, 1,
                            length - (1L << (k - 1)));
}

// Interpolation with subproduct tree.
void interpolateK(Kelt * f, int k)
{
    long i, j;
    Kelt beta;

    for (i = 0; i < k - 1; ++i) {
        for (j = 0; j < (1L << (k - 1 - i)); ++j) {
            long index = j * (1L << (i + 1));
            allBetai(beta, 2 * j);
            interpolateSi(i, f + index, beta);
        }
    }
    interpolateSiNobeta(k - 1, f);
}


// Quotrem Si. f is possibly larger than 2 times 1L<<k.
// Degree of f is given by L (with L > 1L<<k)
// Inplace computation: low part of f receive rem, high part gets quo.
static inline void quotremSi(int k, Kelt * f, long L, Kelt beta)
{
    long j;
    long i, K;
    Kelt coeff;
    Ksrc_elt fi;

    K = 1L << k;
    // put rest mod Sk+beta in low part of f and quo in high part.
    for (i = L - 1; i >= K; --i) {
        fi = f[i];
        for (j = 0; j < ind_number[k]; ++j) {
            long index = ind_S[k][j] + i - K;
            Kadd(f[index], f[index], fi);
        }
        Kmul(coeff, fi, beta);
        Kadd(f[i - K], f[i - K], coeff);
    }
}

// Mul Si. f is possibly smaller than half of 1L<<k.
// Degree of f is given by L (with L<= (1L<<k))
// Inplace computation: Top part of the result is stored in f, and low
// part of the result gets added to what is at the right of f.
static inline void mulSi(int k, Kelt * f, long L, Kelt beta)
{
    long j;
    long i, K = 1L << k;

    for (i = 0; i < L; ++i) {
        Kelt coeff;
        Ksrc_elt fi;
        fi = f[i];
        for (j = 0; j < ind_number[k]; ++j) {
            long index = ind_S[k][j] + i - K;
            Kadd(f[index], f[index], fi);
        }
        Kmul(coeff, fi, beta);
        Kadd(f[i - K], f[i - K], coeff);
    }
}

// Reduce f (with 1L<<k coeffs) modulo the product of s_i corresponding
// to length.
static inline void reduceModTrunc(Kelt * f, int k, long length)
{
    int i;
    long len, ii;
    len = length;
    Kelt *pf = f;
    Kelt beta;
    long offset_f;

    // go down, computing quotrems
    offset_f = 0;
    for (i = k - 1; i >= 0; --i) {
        if (len >> i) {
            allBetai(beta, offset_f >> i);
            quotremSi(i, f + offset_f, (1L << k) - offset_f, beta);
            offset_f += 1L << i;
            len -= 1L << i;
        }
    }

    // go up, reconstructing general rem
    len = length;
    i = __builtin_ctz(len);
    len >>= i;
    pf = f + length;
    ii = 1L << i;
    i++;
    len >>= 1;
    for (; i <= k - 1; ++i) {
        if (len & 1) {
            allBetai(beta, (length - (1L << i) - ii) >> i);
            mulSi(i, pf - ii, ii, beta);
            ii += 1L << i;
        }
        len >>= 1;
    }
}


// Interpolation with subproduct tree.
// This is a truncated version, where we reconstruct length coefficients,
// from length values. (length <= 1L<<k)
// Assume that non-computed values are set to zero.
void interpolateK_trunc(Kelt * f, int k, long length)
{
    int i;
    long j;
    Kelt beta;

    for (i = 0; i < k - 1; ++i) {
        long maxj = 1L << (k - 1 - i);
        for (j = 0; j < maxj; ++j) {
            long index = j * (1L << (i + 1));
            if (length <= (index + (1L << i))) {
                break;
            }
            allBetai(beta, 2 * j);
            interpolateSi(i, f + index, beta);
        }
    }
    if (length > (1L << (k - 1)))
        interpolateSiNobeta(k - 1, f);

    reduceModTrunc(f, k, length);
}

#if (GMP_LIMB_BITS == 64)
/* Assume that ulongs are 64 bits */
#define ULONG_BITS 64
#elif (GMP_LIMB_BITS == 32)
#define ULONG_BITS 32
#else
#error "define GMP_LIMB_BITS"
#endif

#if (ULONG_BITS == 64)
void decomposeK(Kelt * f, unsigned long * F, long Fl, int k)
{
    long i;
    for (i = 0; i < Fl; ++i) {
        f[i][0] = F[i];
        f[i][1] = 0;
    }
    memset(f + i, 0, ((1L << k) - Fl) * sizeof(Kelt));
}

void recomposeK(unsigned long * F, Kelt * f, long Fl, int k MAYBE_UNUSED)
{
    long i;

    F[0] = f[0][0];
    for (i = 1; i < Fl ; ++i)
        F[i] = f[i][0] ^ f[i - 1][1];
}
#elif (ULONG_BITS == 32)

/* untested */

void decomposeK(Kelt * f, unsigned long * F, long Fl, int k)
{
    long i;
    for (i = 0; i < Fl / 2 ; ++i) {
        f[i][0] = F[2*i];
        f[i][1] = F[2*i + 1];
        f[i][2] = 0;
        f[i][3] = 0;
    }
    if (Fl & 1) {
        f[i][0] = F[2*i];
        f[i][1] = 0;
        f[i][2] = 0;
        f[i][3] = 0;
        i++;
    }
    memset(f + i, 0, ((1L << k) - i) * sizeof(Kelt));
}

void recomposeK(unsigned long * F, Kelt * f, long Fl, int k MAYBE_UNUSED)
{
    long i;

    F[0] = f[0][0];
    F[1] = f[0][1];
    for (i = 2; i < Fl; i += 2) {
        F[i] = f[i/2][0] ^ f[i/2 - 1][2];
        F[i+1] = f[i/2][1] ^ f[i/2 - 1][3];
    }
    if (Fl & 1) {
            F[i] = f[i/2][0] ^ f[i/2 - 1][2];
    }
}
#else
#error "define ULONG_BITS"
#endif

/* nF is a number of coefficients */
void c128_setup(c128_info_t p, long nF, long nG)
{
    int k;
    long Hl;
    long n;

    /* Since internally we're working with 64-bit data, then it's really
     * a hard 64 here, not ULONG_BITS : We're just deciding on the order
     * of things.
     */
    long Fl = (nF + 63) / 64;
    long Gl = (nG + 63) / 64;

    Hl = Fl + Gl;               // nb of uint64_t of the result
    n = Hl;                     // nb of Kelt of the result.
    k = 1;
    while ((1L << k) < n)
        ++k;
    if (k < 2) {
        fprintf(stderr, "this code is for k >= 2, sorry\n");
        exit(1);
    }
    p->k = k;
    p->n = n;

}

/* nF is a number of coefficients */
void c128_dft(const c128_info_t p, c128_t x, unsigned long * F, long nF)
{
    long Fl = (nF + ULONG_BITS - 1) / ULONG_BITS;
    decomposeK(x,F,Fl,p->k);
#ifdef  CANTOR_GM
    multievaluateGM(x, p->k, p->n);
#else
    multievaluateKnew_trunc(x, p->k, p->n);
#endif
}


void c128_compose(const c128_info_t p, c128_t y, c128_src_t x1, c128_src_t x2)
{
    long j;
    for (j = 0; j < p->n; ++j) {
        Kmul(y[j], x1[j], x2[j]);
    }
}
void c128_add(const c128_info_t p, c128_t y, c128_src_t x1, c128_src_t x2)
{
    long j;
    for (j = 0; j < p->n; ++j) {
        Kadd(y[j], x1[j], x2[j]);
    }
}

void c128_cpy(const c128_info_t p, c128_t y, c128_src_t x)
{
    memcpy(y, x, (p->n)*sizeof(Kelt));
}

int c128_size(const c128_info_t p)
{
    return (p->n);
}

/* nH is a number of coefficients */
void c128_ift(
        const c128_info_t p,
        unsigned long * H,
        long nH,
        c128_t h)
{
    long Hl = (nH + ULONG_BITS - 1) / ULONG_BITS;

    // fill in with zeros to facilitate interpolation
    memset(h + p->n, 0, ((1L << p->k) - p->n) * sizeof(Kelt));
    if (p->n & (p->n - 1)) {
        /* n is not a power of 2 */
        interpolateK_trunc(h, p->k, p->n);
    } else {
        interpolateK(h, p->k);
    }
    recomposeK(H, h, Hl, p->k);
}


// Main function:
// Takes two arrays of 64 bit limbs and compute their
// product as polynomials over GF(2). 
// H must have room for the result (Fl+Gl limbs)
void mulCantor128(unsigned long * H, unsigned long * F, long Fl, unsigned long * G, long Gl)
{
    c128_info_t order;
    c128_t f,g;

    long nF = Fl * ULONG_BITS;
    long nG = Gl * ULONG_BITS;

    c128_setup(order, nF, nG);

    f = c128_alloc(order, 1);
    c128_dft(order, f, F, nF);

    g = c128_alloc(order, 1);
    c128_dft(order, g, G, nG);

    c128_compose(order, f, f, g);

    c128_free(order, g, 1);
    c128_ift(order, H, nF + nG -1, f);

    c128_free(order, f, 1);
}

c128_t c128_alloc(const c128_info_t p, int n)
{
    return (c128_t) malloc((n << p->k) * sizeof(Kelt));
}
void c128_free(
        const c128_info_t p MAYBE_UNUSED,
        c128_t x,
        int n MAYBE_UNUSED)
{
    free(x);
}
c128_t c128_get(const c128_info_t p, c128_t x, int k)
{
    return x + (k << p->k);
}
void c128_zero(const c128_info_t p, c128_t x, int n)
{
	memset(x, 0, (n << p->k) * sizeof(Kelt));
}
/* vim: set sw=4 sta et: */
