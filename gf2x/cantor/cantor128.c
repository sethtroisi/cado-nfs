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

#define Kelt    MPFQ_(elt)
#define Kdst_elt        MPFQ_(dst_elt)
#define Ksrc_elt        MPFQ_(src_elt)
#define Kadd(a0,a1,a2)  MPFQ_(add) (NULL,a0,a1,a2)
#define Kmul(a0,a1,a2)  MPFQ_(mul) (NULL,a0,a1,a2)
#define Ksqr(a0,a1)     MPFQ_(sqr) (NULL,a0,a1)

// Some constants related to Cantor's algorithm
// Beta_i are such that Beta_{i-1} = Beta_i^2 + Beta_i
// and { Beta_i : i in [1..128] } is a GF(2)-basis of GF(2^128).
// The following is (the begining of) a solution.
// NB: the defining polynomial for GF(2^128) is x^128 + x^7 + x^2 + x + 1
static const Kelt Betai[32] = {
    {UINT64_C(1)},
    {UINT64_C(7451958335100816136), UINT64_C(2979905974001933049)},
    {UINT64_C(18379359142562139938), UINT64_C(11678753174964950217)},
    {UINT64_C(6032672185962850376), UINT64_C(8744256601846606146)},
    {UINT64_C(14608789766813931076), UINT64_C(645900494672607458)},
    {UINT64_C(5568895992531017964), UINT64_C(5316835906050072984)},
    {UINT64_C(11619261390503595532), UINT64_C(377604546732988956)},
    {UINT64_C(6679137017075335448), UINT64_C(5571574281931689094)},
    {UINT64_C(16252459407801923090), UINT64_C(3227645990918638800)},
    {UINT64_C(6380301344616957040), UINT64_C(2381717826074782822)},
    {UINT64_C(15008204915294424484), UINT64_C(13983984897465862323)},
    {UINT64_C(18362218515738631672), UINT64_C(1174612800055669010)},
    {UINT64_C(13770211835728229036), UINT64_C(14885908188291033131)},
    {UINT64_C(1856991084319558808), UINT64_C(15829302951382751119)},
    {UINT64_C(17677502823211816656), UINT64_C(6374576737782705044)},
    {UINT64_C(10058708795285882084), UINT64_C(5296104361219370626)},
    {UINT64_C(11604927872398312824), UINT64_C(408267762571371738)},
    {UINT64_C(13468465821495725140), UINT64_C(9429188815902477435)},
    {UINT64_C(10268364117098907234), UINT64_C(13428504000507276907)},
    {UINT64_C(6641545852185192764), UINT64_C(3462008082606701680)},
    {UINT64_C(9812826009415599412), UINT64_C(14006017226737078185)},
    {UINT64_C(10796008238720342306), UINT64_C(4109866943845289010)},
    {UINT64_C(14658682906370908118), UINT64_C(9994441883493921297)},
    {UINT64_C(11615212764610847988), UINT64_C(13699609717720089643)},
    {UINT64_C(16610661676694873598), UINT64_C(12948644223555576813)},
    {UINT64_C(7901269838824795632), UINT64_C(3796759253579916832)},
    {UINT64_C(12454147674074156822), UINT64_C(12010935663861637443)},
    {UINT64_C(245289927884982804), UINT64_C(3843120356884484498)},
    {UINT64_C(14321890119743223594), UINT64_C(16915134361088465567)},
    {UINT64_C(18183411989044536782), UINT64_C(1839810178943853948)},
    {UINT64_C(8136972654088694852), UINT64_C(6965936934891198430)},
    {UINT64_C(12242034105964537428), UINT64_C(1006149766749830960)}
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
    x[0] = 0;
    x[1] = 0;
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
    int Fll = Fl + (Fl & 1);
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
        memset(f + i, 0, ((1 << k) - Fl/2 - 1) * sizeof(Kelt));
    } else {
        memset(f + i, 0, ((1 << k) - Fl/2) * sizeof(Kelt));
    }
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
    if (k < 7) {
        fprintf(stderr, "this code is for k >= 7, sorry\n");
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
/* vim: set sw=4 sta et: */
