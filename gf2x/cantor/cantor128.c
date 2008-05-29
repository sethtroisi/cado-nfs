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

#define Kelt            mpfq_2_128_elt
#define Kdst_elt        mpfq_2_128_dst_elt
#define Ksrc_elt        mpfq_2_128_src_elt
#define Kadd(a0,a1,a2)  mpfq_2_128_add (NULL,a0,a1,a2)
#define Kmul(a0,a1,a2)  mpfq_2_128_mul (NULL,a0,a1,a2)
#define Kset_ui(a0,a1)  mpfq_2_128_set_ui (NULL,a0,a1)
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
static const int ind_S[31][15] = {
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
static inline void allBetai(Kelt x, int i)
{
    int j;
    Kset_ui(x,0);
    j = 0;
    while (i != 0) {
        if (i & 1) {
            Kadd(x, x, Betai[j]);
        }
        i >>= 1;
        j++;
    }
}

/* 20080122 -- hard-coded versions of reduceSi for small sizes have been
 * tried, but do not seem to pay off. Look for the svn data to
 * investigate again the potential benefit.
 */

// The generic reduction function.
static inline void reduceSi(int k, Kelt * f, Kelt beta)
{
    int i, j, K;
    Kelt coeff;
    Ksrc_elt fi;

    K = 1 << k;
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

// The generic reduction function.
// truncated version
static inline void reduceSiNobeta_trunc(int k, Kelt * f, int L)
{
    int i, j, K;
    Ksrc_elt fi;

    K = 1 << k;
    // put rest mod Sk+beta in low part of f and quo in high part.
    for (i = L - 1; i >= K; --i) {
        fi = f[i];
        for (j = 0; j < ind_number[k]; ++j) {
            int index = ind_S[k][j] + i - K;
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
    int i, j;
    int K = 1 << k;

    for (i = 0; i < K; ++i)
        Kadd(f[i + K], f[i + K], f[i]);

    for (i = K; i < K << 1; ++i) {
        Kelt coeff;
        Ksrc_elt fi;
        fi = f[i];
        for (j = 0; j < ind_number[k]; ++j) {
            int index = ind_S[k][j] + i - K;
            Kadd(f[index], f[index], fi);
        }
        Kmul(coeff, fi, beta);
        Kadd(f[i - K], f[i - K], coeff);
    }
}

static inline void interpolateSiNobeta(int k, Kelt * f)
{
    int i, j;
    int K = 1 << k;

    for (i = 0; i < K; ++i)
        Kadd(f[i + K], f[i + K], f[i]);

    for (i = K; i < K << 1; ++i) {
        Ksrc_elt fi;
        fi = f[i];
        for (j = 0; j < ind_number[k]; ++j) {
            int index = ind_S[k][j] + i - K;
            Kadd(f[index], f[index], fi);
        }
    }
}


void multievaluateKrec(Kelt * f, int i, int rep_beta)
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
        multievaluateKrec(f + (1 << (i - 1)), i - 1, rep_beta + 1);
    }
}

void multievaluateKrec_trunc(Kelt * f, int i, int rep_beta, int length)
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
        if (length < (1 << (i - 1))) {
            multievaluateKrec_trunc(f, i - 1, rep_beta, length);
        } else if (length > (1 << (i - 1))) {
            multievaluateKrec(f, i - 1, rep_beta);
            multievaluateKrec_trunc(f + (1 << (i - 1)), i - 1, rep_beta + 1,
                                    length - (1 << (i - 1)));
        } else { // length == 1<<(i-1)
            multievaluateKrec(f, i - 1, rep_beta);
        }
    }
}


void multievaluateKnew_trunc(Kelt * f, int k, int length)
{
    reduceSiNobeta_trunc(k - 1, f, length);

    multievaluateKrec(f, k - 1, 0);
    multievaluateKrec_trunc(f + (1 << (k - 1)), k - 1, 1,
                            length - (1 << (k - 1)));
}

// Interpolation with subproduct tree.
void interpolateK(Kelt * f, int k)
{
    int i, j;
    Kelt beta;

    for (i = 0; i < k - 1; ++i) {
        for (j = 0; j < (1 << (k - 1 - i)); ++j) {
            int index = j * (1 << (i + 1));
            allBetai(beta, 2 * j);
            interpolateSi(i, f + index, beta);
        }
    }
    interpolateSiNobeta(k - 1, f);
}


// Quotrem Si. f is possibly larger than 2 times 1<<k.
// Degree of f is given by L (with L > 1<<k)
// Inplace computation: low part of f receive rem, high part gets quo.
static inline void quotremSi(int k, Kelt * f, int L, Kelt beta)
{
    int i, j, K;
    Kelt coeff;
    Ksrc_elt fi;

    K = 1 << k;
    // put rest mod Sk+beta in low part of f and quo in high part.
    for (i = L - 1; i >= K; --i) {
        fi = f[i];
        for (j = 0; j < ind_number[k]; ++j) {
            int index = ind_S[k][j] + i - K;
            Kadd(f[index], f[index], fi);
        }
        Kmul(coeff, fi, beta);
        Kadd(f[i - K], f[i - K], coeff);
    }
}

// Mul Si. f is possibly smaller than half of 1<<k.
// Degree of f is given by L (with L<= (1<<k))
// Inplace computation: Top part of the result is stored in f, and low
// part of the result gets added to what is at the right of f.
static inline void mulSi(int k, Kelt * f, int L, Kelt beta)
{
    int i, j;
    int K = 1 << k;

    for (i = 0; i < L; ++i) {
        Kelt coeff;
        Ksrc_elt fi;
        fi = f[i];
        for (j = 0; j < ind_number[k]; ++j) {
            int index = ind_S[k][j] + i - K;
            Kadd(f[index], f[index], fi);
        }
        Kmul(coeff, fi, beta);
        Kadd(f[i - K], f[i - K], coeff);
    }
}

// Reduce f (with 1<<k coeffs) modulo the product of s_i corresponding
// to length.
static inline void reduceModTrunc(Kelt * f, int k, int length)
{
    int i, len, ii;
    len = length;
    Kelt *pf = f;
    Kelt beta;
    int offset_f;

    // go down, computing quotrems
    offset_f = 0;
    for (i = k - 1; i >= 0; --i) {
        if (len >> i) {
            allBetai(beta, offset_f >> i);
            quotremSi(i, f + offset_f, (1 << k) - offset_f, beta);
            offset_f += 1 << i;
            len -= 1 << i;
        }
    }

    // go up, reconstructing general rem
    len = length;
    i = __builtin_ctz(len);
    len >>= i;
    pf = f + length;
    ii = 1 << i;
    i++;
    len >>= 1;
    for (; i <= k - 1; ++i) {
        if (len & 1) {
            allBetai(beta, (length - (1 << i) - ii) >> i);
            mulSi(i, pf - ii, ii, beta);
            ii += 1 << i;
        }
        len >>= 1;
    }
}


// Interpolation with subproduct tree.
// This is a truncated version, where we reconstruct length coefficients,
// from length values. (length <= 1<<k)
// Assume that non-computed values are set to zero.
void interpolateK_trunc(Kelt * f, int k, int length)
{
    int i, j;
    Kelt beta;

    for (i = 0; i < k - 1; ++i) {
        int maxj = 1 << (k - 1 - i);
        for (j = 0; j < maxj; ++j) {
            int index = j * (1 << (i + 1));
            if (length <= (index + (1 << i))) {
                break;
            }
            allBetai(beta, 2 * j);
            interpolateSi(i, f + index, beta);
        }
    }
    if (length > (1 << (k - 1)))
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
void decomposeK(Kelt * f, unsigned long * F, int Fl, int k)
{
    int i;
    for (i = 0; i < Fl; ++i) {
        f[i][0] = F[i];
        f[i][1] = 0;
    }
    memset(f + i, 0, ((1 << k) - Fl) * sizeof(Kelt));
}

void recomposeK(unsigned long * F, Kelt * f, int Fl, int k MAYBE_UNUSED)
{
    int i;

    F[0] = f[0][0];
    for (i = 1; i < Fl ; ++i)
        F[i] = f[i][0] ^ f[i - 1][1];
}
#elif (ULONG_BITS == 32)

/* untested */

void decomposeK(Kelt * f, unsigned long * F, int Fl, int k)
{
    int i;
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
    memset(f + i, 0, ((1 << k) - i) * sizeof(Kelt));
}

void recomposeK(unsigned long * F, Kelt * f, int Fl, int k MAYBE_UNUSED)
{
    int i;

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

void c128_setup(c128_info_t p, int dF, int dG)
{
    int k;
    int Hl;
    int n;

    int Fl = (dF + 1 + 63) / 64;
    int Gl = (dG + 1 + 63) / 64;

    Hl = Fl + Gl;               // nb of uint64_t of the result
    n = Hl;                     // nb of Kelt of the result.
    k = 1;
    while ((1 << k) < n)
        ++k;
    if (k < 2) {
        fprintf(stderr, "this code is for k >= 2, sorry\n");
        exit(1);
    }
    p->k = k;
    p->n = n;
}

void c128_dft(const c128_info_t p, c128_t x, unsigned long * F, int dF)
{
    int Fl = 1 + dF / ULONG_BITS;
    decomposeK(x,F,Fl,p->k);
    multievaluateKnew_trunc(x, p->k, p->n);
}
void c128_compose(const c128_info_t p, c128_t y, c128_src_t x1, c128_src_t x2)
{
    int j;
    for (j = 0; j < p->n; ++j) {
        Kmul(y[j], x1[j], x2[j]);
    }
}
void c128_add(const c128_info_t p, c128_t y, c128_src_t x1, c128_src_t x2)
{
    int j;
    for (j = 0; j < p->n; ++j) {
        Kadd(y[j], x1[j], x2[j]);
    }
}
void c128_ift(
        const c128_info_t p,
        unsigned long * H,
        int dH,
        c128_t h)
{
    int Hl = 1 + dH / ULONG_BITS;

    // fill in with zeros to facilitate interpolation
    memset(h + p->n, 0, ((1 << p->k) - p->n) * sizeof(Kelt));
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
void mulCantor128(unsigned long * H, unsigned long * F, int Fl, unsigned long * G, int Gl)
{
    c128_info_t order;
    c128_t f,g;

    int dF = Fl * ULONG_BITS - 1;
    int dG = Gl * ULONG_BITS - 1;

    c128_setup(order, dF, dG);

    f = c128_alloc(order, 1);
    c128_dft(order, f, F, dF);

    g = c128_alloc(order, 1);
    c128_dft(order, g, G, dG);

    c128_compose(order, f, f, g);

    c128_free(order, g, 1);
    c128_ift(order, H, dF + dG, f);

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
c128_t c128_get(const c128_info_t p, c128_t x, int k) {
	return x + (k << p->k);
}
void c128_zero(const c128_info_t p, c128_t x, int n)
{
	memset(x, 0, n * (1 << p->k) * sizeof(Kelt));
}
/* vim: set sw=4 sta et: */
