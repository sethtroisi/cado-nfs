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
 

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <stdint.h>
#include <gmp.h>

#include "mpfq_2_128.h"

// Hum... hopefully this one is not already declared.
typedef uint64_t uint128_t[2];

// Wrapper around mpfq's mul function. 
// Assume something about internals, here, I'm afraid...

#ifdef COUNT_MUL
static long cpt;
#endif

static inline 
void mul_f128(uint128_t z, uint128_t x, uint128_t y) {
  mpfq_2_128_mul(NULL, (mpfq_2_128_dst_elt) z, (mpfq_2_128_src_elt) x, (mpfq_2_128_src_elt) y);
#ifdef COUNT_MUL
  cpt++;
#endif
}

static inline
void add_f128(uint128_t z, const uint128_t x, const uint128_t y) {
#if 1  // SSE seems to be slower, here
  z[0] = x[0]^y[0];
  z[1] = x[1]^y[1];
#else
  typedef uint64_t v2di __attribute__ ((vector_size (16)));
  v2di xx = (v2di) {x[0], x[1]};
  v2di yy = (v2di) {y[0], y[1]};
  xx ^= yy;
  * (v2di*) (z) = xx;
#endif
}

// Some constants related to Cantor's algorithm
// Beta_i are such that Beta_{i-1} = Beta_i^2 + Beta_i
// and { Beta_i : i in [1..128] } is a GF(2)-basis of GF(2^128).
// The following is (the begining of) a solution.
// NB: the defining polynomial for GF(2^128) is x^128 + x^7 + x^2 + x + 1
static const uint128_t Betai[32] = {
    { UINT64_C(1) },
    { UINT64_C(7451958335100816136), UINT64_C(2979905974001933049) },
    { UINT64_C(18379359142562139938), UINT64_C(11678753174964950217) },
    { UINT64_C(6032672185962850376), UINT64_C(8744256601846606146) },
    { UINT64_C(14608789766813931076), UINT64_C(645900494672607458) },
    { UINT64_C(5568895992531017964), UINT64_C(5316835906050072984) },
    { UINT64_C(11619261390503595532), UINT64_C(377604546732988956) },
    { UINT64_C(6679137017075335448), UINT64_C(5571574281931689094) },
    { UINT64_C(16252459407801923090), UINT64_C(3227645990918638800) },
    { UINT64_C(6380301344616957040), UINT64_C(2381717826074782822) },
    { UINT64_C(15008204915294424484), UINT64_C(13983984897465862323) },
    { UINT64_C(18362218515738631672), UINT64_C(1174612800055669010) },
    { UINT64_C(13770211835728229036), UINT64_C(14885908188291033131) },
    { UINT64_C(1856991084319558808), UINT64_C(15829302951382751119) },
    { UINT64_C(17677502823211816656), UINT64_C(6374576737782705044) },
    { UINT64_C(10058708795285882084), UINT64_C(5296104361219370626) },
    { UINT64_C(11604927872398312824), UINT64_C(408267762571371738) },
    { UINT64_C(13468465821495725140), UINT64_C(9429188815902477435) },
    { UINT64_C(10268364117098907234), UINT64_C(13428504000507276907) },
    { UINT64_C(6641545852185192764), UINT64_C(3462008082606701680) },
    { UINT64_C(9812826009415599412), UINT64_C(14006017226737078185) },
    { UINT64_C(10796008238720342306), UINT64_C(4109866943845289010) },
    { UINT64_C(14658682906370908118), UINT64_C(9994441883493921297) },
    { UINT64_C(11615212764610847988), UINT64_C(13699609717720089643) },
    { UINT64_C(16610661676694873598), UINT64_C(12948644223555576813) },
    { UINT64_C(7901269838824795632), UINT64_C(3796759253579916832) },
    { UINT64_C(12454147674074156822), UINT64_C(12010935663861637443) },
    { UINT64_C(245289927884982804), UINT64_C(3843120356884484498) },
    { UINT64_C(14321890119743223594), UINT64_C(16915134361088465567) },
    { UINT64_C(18183411989044536782), UINT64_C(1839810178943853948) },
    { UINT64_C(8136972654088694852), UINT64_C(6965936934891198430) },
    { UINT64_C(12242034105964537428), UINT64_C(1006149766749830960) }
};

// The S_i polynomials of Cantor's algorithm: linearized polynomials such
// that S_i | S_{i+1}.
static const int ind_S[31][15] = {
  {}, 	// S0
  {1},	// S1
  {1},	// S2
  {4,2,1},	// S3
  {1},	// S4
  {16,2,1},	// S5
  {16,4,1},	// S6
  {64,32,16,8,4,2,1}, 	// S7
  {1},			// S8
  {256,2,1},			// S9
  {256,4,1},			// S10
  {1024, 512,256,8,4,2,1},	// S11
  {256,16,1},			// S12
  {4096,512,256,32,16,2,1},	// S13
  {4096,1024,256,64,16,4,1},	// S14
  {16384,8192,4096,2048,1024,512,256,128,64,32,16,8,4,2,1},	// S15
  {1},          // S16
  {65536,2,1},  // S17
  {65536,4,1},  // S18
  {262144,131072,65536,8,4,2,1},        // S19
  {65536,16,1}, // S20
  {1048576,131072,65536,32,16,2,1},     // S21
  {1048576,262144,65536,64,16,4,1},     // S22
  {4194304,2097152,1048576,524288,262144,131072,65536,128,64,32,16,8,4,2,1},    // S23
  {65536,256,1},        // S24
  {16777216,131072,65536,512,256,2,1},  // S25
  {16777216,262144,65536,1024,256,4,1}, // S26
  {67108864,33554432,16777216,524288,262144,131072,65536,2048,1024,512,256,8,4,2,1},    // S27
  {16777216,1048576,65536,4096,256,16,1},       // S28
  {268435456,33554432,16777216,2097152,1048576,131072,65536,8192,4096,512,256,32,16,2,1}, // S29
  {268435456,67108864,16777216,4194304,1048576,262144,65536,16384,4096,1024,256,64,16,4,1},     // S30
};

// The number of terms in S_i.
static const int ind_number[31] = {0,1,1,3,1,3,3,7,1,3,3,7,3,7,7,15,
  1,3,3,7,3,7,7,15,3,7,7,15,7,15,15};

static uint128_t *tableBeta;

// At index i=i0+2*i1+4*i2+8*i3..., one computes 
//     allBetai[i] = i0*Betai[1] + i1*Betai[2] + i2*Betai[3] + ...
static inline
void allBetai(uint128_t x, int i) {
  int j;
  x[0] = 0;
  x[1] = 0;
  j = 0;
  while (i != 0) {
    if ((i%2) == 1)
      add_f128(x, x, Betai[j]);
    i >>= 1;
    j++;
  }
}

static
void initTableBeta(int k) {
  int i, j;
  tableBeta = (uint128_t *)malloc((1<<k)*sizeof(uint128_t));
  for (i = 0; i < k; ++i) {
    for (j = 0; j < (1<<k); j+= (1<<i)) {
      add_f128(tableBeta[j], tableBeta[j], Betai[i]);
    }
  }
}

// Reduction functions.
// This one must be seen as a macro that will be used for small sizes.
static inline 
void reduceSk(uint128_t* f, uint128_t beta, int k, int *ind_Sk, int ind_number) {
  int i, j, K;
  uint128_t coeff;
  uint64_t *fi;

  K = 1<<k;
  // put rest mod Sk+beta in low part of f and quo in high part.
  for (i = 2*K-1; i >= K; --i) {
    fi = f[i];
    for (j = 0; j < ind_number; ++j) {
      int index = ind_Sk[j]+i-K;
      add_f128(f[index], f[index], fi);
    }
    mul_f128(coeff, fi, beta);
    add_f128(f[i-K], f[i-K], coeff);
  }

  for (i = 0; i < K; ++i)
    add_f128(f[i+K],f[i+K], f[i]);
}

// The smallest reductions are written specifically
// (don't trust too much the compiler...)

static inline
void reduceS0(uint128_t* f, uint128_t beta) {
  uint128_t tmp;
  mul_f128(tmp, beta, f[1]);
  add_f128(f[0], f[0], tmp);
  add_f128(f[1], f[1], f[0]);
}

static inline
void reduceS1(uint128_t* f, uint128_t beta) {
  int ind_S1[1] = {1};
  reduceSk(f, beta, 1, ind_S1, 1);
}

static inline
void reduceS2(uint128_t* f, uint128_t beta) {
  int ind_S2[1] = {1};
  reduceSk(f, beta, 2, ind_S2, 1);
}

static inline
void reduceS3(uint128_t* f, uint128_t beta) {
  int ind_S3[3] = {4, 2, 1};
  reduceSk(f, beta, 3, ind_S3, 3);
}

static inline
void reduceS4(uint128_t* f, uint128_t beta) {
  int ind_S4[1] = {1};
  reduceSk(f, beta, 4, ind_S4, 1);
}

static inline
void reduceS5(uint128_t* f, uint128_t beta) {
  int ind_S5[3] = {16, 2, 1};
  reduceSk(f, beta, 5, ind_S5, 3);
}

static inline
void reduceS6(uint128_t* f, uint128_t beta) {
  int ind_S6[3] = {16, 4, 1};
  reduceSk(f, beta, 6, ind_S6, 3);
}

static inline
void reduceS7(uint128_t* f, uint128_t beta) {
  int ind_S7[7] = {64, 32,16,8,4,2,1};
  reduceSk(f, beta, 7, ind_S7, 7);
}

static inline
void reduceS8(uint128_t* f, uint128_t beta) {
  int ind_S8[1] = {1};
  reduceSk(f, beta, 8, ind_S8, 1);
}

// The generic reduction function.
static inline 
void reduceSi(int k, uint128_t* f, uint128_t beta) {
  int i, j, K;
  uint128_t coeff;
  uint64_t *fi;

  K = 1<<k;
  // put rest mod Sk+beta in low part of f and quo in high part.
  for (i = 2*K-1; i >= K; --i) {
    fi = f[i];
    for (j = 0; j < ind_number[k]; ++j) {
      int index = ind_S[k][j]+i-K;
      add_f128(f[index], f[index], fi);
    }
    mul_f128(coeff, fi, beta);
    add_f128(f[i-K],f[i-K], coeff);
  }

  // The other rest is just the sum of low and hi parts
  for (i = 0; i < K; ++i)
    add_f128(f[i+K],f[i+K],f[i]);
}

// The generic reduction function.
static inline 
void reduceSiNobeta(int k, uint128_t* f) {
  int i, j, K;
  uint128_t coeff;
  uint64_t *fi;

  K = 1<<k;
  // put rest mod Sk+beta in low part of f and quo in high part.
  for (i = 2*K-1; i >= K; --i) {
    fi = f[i];
    for (j = 0; j < ind_number[k]; ++j) {
      int index = ind_S[k][j]+i-K;
      add_f128(f[index], f[index], fi);
    }
  }

  // The other rest is just the sum of low and hi parts
  for (i = 0; i < K; ++i)
    add_f128(f[i+K],f[i+K],f[i]);
}

// The generic reduction function.
// truncated version
static inline 
void reduceSiNobeta_trunc(int k, uint128_t* f, int L) {
  int i, j, K;
  uint128_t coeff;
  uint64_t *fi;

  K = 1<<k;
  // put rest mod Sk+beta in low part of f and quo in high part.
  for (i = L; i >= K; --i) {
    fi = f[i];
    for (j = 0; j < ind_number[k]; ++j) {
      int index = ind_S[k][j]+i-K;
      add_f128(f[index], f[index], fi);
    }
  }

  // The other rest is just the sum of low and hi parts
  for (i = 0; i < K; ++i)
    add_f128(f[i+K],f[i+K],f[i]);
}




// A generic interpolation function.
// TODO: unroll small sizes, like for reduction.
static inline
void interpolateSi(int k, uint128_t* f, uint128_t beta) {
  int i, j;
  int K = 1<<k;

  for (i = 0; i < K; ++i)
    add_f128(f[i+K],f[i+K],f[i]);

  for (i = K; i <= 2*K-1; ++i) {
    uint128_t coeff;
    uint64_t *fi;
    fi = f[i];
    for (j = 0; j < ind_number[k]; ++j) {
      int index = ind_S[k][j]+i-K;
      add_f128(f[index], f[index], fi);
    }
    mul_f128(coeff, fi, beta);
    add_f128(f[i-K],f[i-K],coeff);
  }
}

static inline
void interpolateSiNobeta(int k, uint128_t* f) {
  int i, j;
  int K = 1<<k;

  for (i = 0; i < K; ++i)
    add_f128(f[i+K],f[i+K],f[i]);

  for (i = K; i <= 2*K-1; ++i) {
    uint64_t *fi;
    fi = f[i];
    for (j = 0; j < ind_number[k]; ++j) {
      int index = ind_S[k][j]+i-K;
      add_f128(f[index], f[index], fi);
    }
  }
}


void multievaluateKrec(uint128_t* f, int i, int rep_beta) {
  uint128_t beta;
  rep_beta <<= 1;
  allBetai(beta, rep_beta);
  reduceSi(i-1, f, beta);

  if (i == 2) {
    rep_beta <<= 1;
    allBetai(beta, rep_beta);
    reduceSi(0, f, beta);
    add_f128(beta, beta, Betai[1]);
    reduceSi(0, f + 2, beta);
    return;
  } 
  if (i == 3) {
    rep_beta <<= 1;
    allBetai(beta, rep_beta);
    reduceSi(1, f, beta);
    add_f128(beta, beta, Betai[1]);
    reduceSi(1, f + 4, beta);
   
    rep_beta <<= 1;
    allBetai(beta, rep_beta);
    reduceSi(0, f, beta);
    add_f128(beta, beta, Betai[1]);
    reduceSi(0, f + 2, beta);
    add_f128(beta, beta, Betai[2]);
    reduceSi(0, f + 6, beta);
    add_f128(beta, beta, Betai[1]);
    reduceSi(0, f + 4, beta);
    return;
  }
  if (i > 1) {
    multievaluateKrec(f, i-1, rep_beta);
    multievaluateKrec(f + (1<<(i-1)), i-1, rep_beta+1);
  }
}

void multievaluateKrec_trunc(uint128_t* f, int i, int rep_beta, int length) {
  uint128_t beta;
  rep_beta <<= 1;
  allBetai(beta, rep_beta);
  reduceSi(i-1, f, beta);

  if (i == 2) {
    rep_beta <<= 1;
    allBetai(beta, rep_beta);
    reduceSi(0, f, beta);
    add_f128(beta, beta, Betai[1]);
    reduceSi(0, f + 2, beta);
    return;
  } 
  if (i == 3) {
    rep_beta <<= 1;
    allBetai(beta, rep_beta);
    reduceSi(1, f, beta);
    add_f128(beta, beta, Betai[1]);
    reduceSi(1, f + 4, beta);
   
    rep_beta <<= 1;
    allBetai(beta, rep_beta);
    reduceSi(0, f, beta);
    add_f128(beta, beta, Betai[1]);
    reduceSi(0, f + 2, beta);
    add_f128(beta, beta, Betai[2]);
    reduceSi(0, f + 6, beta);
    add_f128(beta, beta, Betai[1]);
    reduceSi(0, f + 4, beta);
    return;
  }
  if (i > 1) {
    if (length < (1<<(i-1)))
      multievaluateKrec_trunc(f, i-1, rep_beta, length);
    else if (length > (1<<(i-1))) {
      multievaluateKrec(f, i-1, rep_beta);
      multievaluateKrec_trunc(f + (1<<(i-1)), i-1, rep_beta+1, length - (1<<(i-1)));
    } else // length == 1<<(i-1)
      multievaluateKrec(f, i-1, rep_beta);
  }
}


void multievaluateKnew(uint128_t* f, int k) {
  reduceSiNobeta(k-1, f);

  multievaluateKrec(f, k-1, 0);
  multievaluateKrec(f + (1<<(k-1)), k-1, 1);
}

void multievaluateKnew_trunc(uint128_t* f, int k, int length) {
  reduceSiNobeta_trunc(k-1, f, length);
  //reduceSiNobeta(k-1, f);

  multievaluateKrec(f, k-1, 0);
  multievaluateKrec_trunc(f + (1<<(k-1)), k-1, 1, length - (1<<(k-1)));
}



// Multipoint evaluation at K=1<<k points
// Subproduct tree.
void multievaluateK(uint128_t* f, int k) {
  int i, j;
  uint128_t beta;
  
  reduceSiNobeta(k-1, f);

  for (i = k-2; i >=0; --i) {
    for (j = 0; j < (1<<(k-1-i)); ++j) {
      int index = j*(1<<(i+1));
      allBetai(beta, 2*j);
      reduceSi(i, f + index, beta);
    }
  }
#if 0
  for (j = 0; j < (1<<(k-6)); ++j) {
    allBetai(beta, 2*j);
    reduceS5(f+64*j, beta);
  }
  for (j = 0; j < (1<<(k-5)); ++j) {
    allBetai(beta, 2*j);
    reduceS4(f+32*j, beta);
  }
  for (j = 0; j < (1<<(k-4)); ++j) {
    allBetai(beta, 2*j);
    reduceS3(f+16*j, beta);
  }
  for (j = 0; j < (1<<(k-3)); ++j) {
    allBetai(beta, 2*j);
    reduceS2(f+8*j, beta);
  }
  for (j = 0; j < (1<<(k-2)); ++j) {
    allBetai(beta, 2*j);
    reduceS1(f+4*j, beta);
  }
  for (j = 0; j < (1<<(k-1)); ++j) {
    allBetai(beta, 2*j);
    reduceS0(f+2*j, beta);
  }
#endif
}

// Interpolation with subproduct tree.
void interpolateK(uint128_t* f, int k) {
  int i, j;
  uint128_t beta;

  for (i = 0; i < k-1; ++i) {
    for (j = 0; j < (1<<(k-1-i)); ++j) {
      int index = j*(1<<(i+1));  
      allBetai(beta, 2*j);
      interpolateSi(i, f + index, beta);
    }
  }
  interpolateSiNobeta(k-1, f);
}


// Quotrem Si. f is possibly larger than 2 times 1<<k.
// Degree of f is given by L (with L > 1<<k)
// Inplace computation: low part of f receive rem, high part gets quo.
static inline 
void quotremSi(int k, uint128_t* f, int L, uint128_t beta) {
  int i, j, K;
  uint128_t coeff;
  uint64_t *fi;

  K = 1<<k;
  // put rest mod Sk+beta in low part of f and quo in high part.
  for (i = L-1; i >= K; --i) {
    fi = f[i];
    for (j = 0; j < ind_number[k]; ++j) {
      int index = ind_S[k][j]+i-K;
      add_f128(f[index], f[index], fi);
    }
    mul_f128(coeff, fi, beta);
    add_f128(f[i-K],f[i-K], coeff);
  }
}

// Mul Si. f is possibly smaller than half of 1<<k.
// Degree of f is given by L (with L<= (1<<k))
// Inplace computation: Top part of the result is stored in f, and low
// part of the result gets added to what is at the right of f.
static inline 
void mulSi(int k, uint128_t* f, int L, uint128_t beta) {
  int i, j;
  int K = 1<<k;

  for (i = 0; i <= L-1; ++i) {
    uint128_t coeff;
    uint64_t *fi;
    fi = f[i];
    for (j = 0; j < ind_number[k]; ++j) {
      int index = ind_S[k][j]+i-K;
      add_f128(f[index], f[index], fi);
    }
    mul_f128(coeff, fi, beta);
    add_f128(f[i-K],f[i-K],coeff);
  }
}

void printPol(uint128_t* f, int fl) {
  int i;
  printf("[ ");
  for(i = 0; i < fl-1; ++i)
    printf("%lu, %lu, ", f[i][0], f[i][1]);
  printf("%lu, %lu ]", f[fl-1][0], f[fl-1][1]);
}

// Reduce f (with 1<<k coeffs) modulo the product of s_i corresponding
// to length.
static inline
void reduceModTrunc(uint128_t* f, int k, int length) {
  int i, j, len, ii;
  len = length;
  uint128_t *pf = f;
  uint128_t beta;
  int offset_f;

  // go down, computing quotrems
  offset_f = 0;
  for (i = k-1; i >= 0; --i) {
    if (len >= (1<<i)) {
      allBetai(beta, offset_f>>i);
      quotremSi(i, f + offset_f, (1<<k) - offset_f, beta);
      offset_f += 1<<i;
      len -= 1<<i;
    }
  }

  // go up, reconstructing general rem
  len = length;
  i = 0;
  while ((len%2) == 0) {
    i++;
    len /= 2;
  }
  pf = f + length;
  ii = 1<<i;
  i++;
  len /= 2;
  for (; i <= k-1; ++i) {
    if ((len%2) == 1) {
      allBetai(beta, (length-(1<<i) - ii)>>i);
//      fprintf(stderr, "i=%d, ii=%d\n", i, ii);
      mulSi(i, pf - ii, ii, beta);
      ii += 1<<i;
    }
    len /= 2;
  }
  //  { printf("fgi := "); printPol(f, 1<<k); printf(";\n");}
}


// Interpolation with subproduct tree.
// This is a truncated version, where we reconstruct length coefficients,
// from length values. (length <= 1<<k)
// Assume that non-computed values are set to zero.
void interpolateK_trunc(uint128_t* f, int k, int length) {
  int i, j;
  uint128_t beta;

  for (i = 0; i < k-1; ++i) {
    int maxj = 1<<(k-1-i);
    for (j = 0; j < maxj; ++j) {
      int index = j*(1<<(i+1));  
      if (length <= (index + (1<<i))) {
	break;
      }
      allBetai(beta, 2*j);
      interpolateSi(i, f + index, beta);
    }
  }
  if (length > (1<<(k-1)))
    interpolateSiNobeta(k-1, f);

  reduceModTrunc(f, k, length);
}

// Multipoint evaluation at K=1<<k points
// Subproduct tree.
// Truncated version: evaluate at length points.
void multievaluateK_trunc(uint128_t* f, int k, int length) {
  int i, j;
  uint128_t beta;
  
  reduceSiNobeta(k-1, f);

  for (i = k-2; i >=0; --i) {
    for (j = 0; j < (1<<(k-1-i)); ++j) {
      int index = j*(1<<(i+1));
      if (length <= index) {
	break;
      }
      allBetai(beta, 2*j);
      reduceSi(i, f + index, beta);
    }
  }
}

// Takes two poynomials of degree <  K=1<<k over GF(2^128)
// Compute their product modulo S[k].
// This is a ``truncated'' version, where only length coefficients of the
// results are wanted. Of course, length should be <= K
void mulModSk_trunc(uint128_t* h, uint128_t* f, uint128_t* g, int k, int length) {
  int j;
  multievaluateKnew_trunc(f, k, length);
  multievaluateKnew_trunc(g, k, length);
  //  printf("fi := "); printPol(f, 1<<k); printf(";\n");
  //  printf("gi := "); printPol(g, 1<<k); printf(";\n");
  
  for (j = 0; j < length; ++j)
    mul_f128(h[j], f[j], g[j]);
  
  // fill in with zeros to facilitate interpolation
  for (j = length; j < (1<<k); ++j) {
    h[j][0] = 0;
    h[j][1] = 0;
  }
  interpolateK_trunc(h, k, length);

}


// Takes two poynomials of degree <  K=1<<k over GF(2^128)
// Compute their product modulo S[k].
void mulModSk(uint128_t* h, uint128_t* f, uint128_t* g, int k) {
  int j;
  multievaluateKnew(f, k);
  multievaluateKnew(g, k);
  for (j = 0; j < (1<<k); ++j) {
    mul_f128(h[j], f[j], g[j]);
  }
  interpolateK(h, k);
}


void decomposeK(uint128_t *f, uint64_t *F, int Fl, int k) {
  int i;
  uint64_t *p;
  for (i = 0; i < Fl; ++i) {
    f[i][0] = F[i];
    f[i][1] = 0;
  }
  for (i = Fl; i < (1<<k); ++i) {
    f[i][0] = 0;
    f[i][1] = 0;
  }
}

void recomposeK(uint64_t *F, uint128_t *f, int Fl, int k) {
  int i;

  F[0] = f[0][0];
  for (i = 1; i < Fl; ++i)
    F[i] = f[i][0] ^ f[i-1][1];
}


// Main function:
// Takes two arrays of 64 bit limbs and compute their
// product as polynomials over GF(2). 
// H must have room for the result (Fl+Gl limbs)
void mulCantor128(uint64_t *H, uint64_t *F, int Fl, uint64_t *G, int Gl) {
  uint128_t *f, *g;
  int k;
  int Hl;
  int n;
  
  Hl = Fl + Gl; // nb of uint64_t of the result
  n = Hl;  // nb of uint128_t of the result.
  k = 1;
  while ((1<<k) < n)
    ++k;
//  fprintf(stderr, "Fl = %d, Gl= %d, k = %d, n = %d\n", Fl ,Gl, k, n);
  if (k < 7) {
    fprintf(stderr, "this code is for k >= 7, sorry\n");
    exit(1);
  }
#ifdef COUNT_MUL
  cpt = 0;
#endif
  
  f = (uint128_t *)malloc((1<<k)*sizeof(uint128_t));
  g = (uint128_t *)malloc((1<<k)*sizeof(uint128_t));
  
  decomposeK(f, F, Fl, k);
  decomposeK(g, G, Gl, k);
  if (n == (1<<k))
    mulModSk(f, f, g, k);
  else
    mulModSk_trunc(f, f, g, k, n);
  recomposeK(H, f, Hl, k);

#ifdef COUNT_MUL
  fprintf(stderr, "nb of muls: %ld\n", cpt);
#endif
  free(g);
  free(f);
}

// Main function:
// Takes two arrays of 64 bit limbs and compute their
// product as polynomials over GF(2). 
// H must have room for the result (Fl+Gl limbs)
void mulCantor128Plain(uint64_t *H, uint64_t *F, int Fl, uint64_t *G, int Gl) {
  uint128_t *f, *g;
  int k;
  int Hl;
  int n;
  
  Hl = Fl + Gl; // nb of uint64_t of the result
  n = Hl;  // nb of uint128_t of the result.
  k = 1;
  while ((1<<k) < n)
    ++k;
  //fprintf(stderr, "Fl = %d, Gl= %d, k = %d\n", Fl ,Gl, k);
  if (k < 7) {
    fprintf(stderr, "this code is for k >= 7, sorry\n");
    exit(1);
  }
#ifdef COUNT_MUL
  cpt = 0;
#endif
  
  f = (uint128_t *)malloc((1<<k)*sizeof(uint128_t));
  g = (uint128_t *)malloc((1<<k)*sizeof(uint128_t));
  
  decomposeK(f, F, Fl, k);
  decomposeK(g, G, Gl, k);
    mulModSk(f, f, g, k);
  recomposeK(H, f, Hl, k);

#ifdef COUNT_MUL
  fprintf(stderr, "nb of muls: %ld\n", cpt);
#endif
  free(g);
  free(f);
}



#ifdef WANT_MAIN

int main(int argc, char **argv) {
  int i, N;
  uint64_t *f, *g, *h;

  if (argc != 2) {
    fprintf(stderr, "usage: %s N\n", argv[0]);
    fprintf(stderr, "  where N is the number of 64-bit limbs of operands\n");
    exit(1);
  } 
  N = atoi(argv[1]);
  
  f = (uint64_t *)malloc(N*sizeof(uint64_t));
  g = (uint64_t *)malloc(N*sizeof(uint64_t));
  h = (uint64_t *)malloc(2*N*sizeof(uint64_t));

  mpn_random((mp_limb_t *)f, (sizeof(uint64_t)/sizeof(mp_limb_t))*N);
  mpn_random((mp_limb_t *)g, (sizeof(uint64_t)/sizeof(mp_limb_t))*N);

#if GMP_LIMB_BITS == 32
#define UINT64_FORMAT "%llu"
#else
#define UINT64_FORMAT "%lu"
#endif

  printf("f := [\n");
  for ( i = 0; i < N-1; ++i)
    printf(UINT64_FORMAT ", ", f[i]);
  printf(UINT64_FORMAT "\n];\n", f[N-1]);
  printf("\n");
  printf("g := [\n");
  for ( i = 0; i < N-1; ++i)
    printf(UINT64_FORMAT ", ", g[i]);
  printf(UINT64_FORMAT "\n];\n", g[N-1]);
  printf("\n");
  
  //for (i = 0; i < 10; ++i) 
      mulCantor128(h, f, N, g, N);

  printf("fg := [\n");
  for ( i = 0; i < 2*N-1; ++i)
    printf(UINT64_FORMAT ", ", h[i]);
  printf(UINT64_FORMAT "\n];\n",h[2*N-1]);

  free(h);
  free(g);
  free(f);
  return 0;
}
#endif

/*

Magma code that checks output:

w := 64;

load "/tmp/toto";
PP := PolynomialRing(GF(2));
F := PP!Intseq(Seqint(f, 2^w), 2);
G := PP!Intseq(Seqint(g, 2^w), 2);
FG := PP!Intseq(Seqint(fg, 2^w), 2);

FG eq F*G;


w := 64;
load "/tmp/toto";
PP := PolynomialRing(GF(2));
x := PP.1;
F128 := ext<GF(2) | x^128 + x^7 + x^2 + x + 1>;
P128<x> := PolynomialRing(F128);

function convertF128(F128, x0, x1)
  X0 := PP!Intseq(x0, 2);
  X1 := PP!Intseq(x1, 2);
  return F128!(X0 + PP.1^w*X1);
end function;

function convertP128(P128, F128, f)
  return P128! [ convertF128(F128, f[2*i-1], f[2*i]) : i in [1..(#f div 2)] ];
end function;

ffi := convertP128(P128,F128,fi);
ggi := convertP128(P128,F128,gi);
h := convertP128(P128,F128,fgi);
fg := ffi*ggi;


beta0 := F128!1;
beta1 := F128!Intseq((7451958335100816136 +2^64*2979905974001933049 ), 2);
beta2 := F128!Intseq((18379359142562139938+2^64*11678753174964950217), 2);
beta3 := F128!Intseq((6032672185962850376 +2^64*8744256601846606146 ), 2);
beta4 := F128!Intseq((14608789766813931076+2^64*645900494672607458  ), 2);
beta5 := F128!Intseq((5568895992531017964 +2^64*5316835906050072984 ), 2);
beta6 := F128!Intseq((11619261390503595532+2^64*377604546732988956  ), 2);
beta7 := F128!Intseq((6679137017075335448 +2^64*5571574281931689094 ), 2);

s4 := x^16+x;
s5 := x^32+x^16+x^2+x;
s6 := x^64+x^16+x^4+x;
s7 := x^128+x^64+x^32+x^16+x^8+x^4+x^2+x;


///////////////////////////////////////////////////////////////////////////

procedure quotremSi(i, pf, ii)
  printf "QuotremSi(%o, %o, %o)\n", i, pf, ii;
end procedure;

procedure mulSi(i, pf, ii)
  printf "mulSi(%o, %o, %o)\n", i, pf, ii;
end procedure;



procedure reduceModTrunc(k, length)
  pf := 0;
  len := length;

  // go down, computing quotrems
  printf "Going down...\n";
  ii := 2^k;
  for i := k-1 to 0 by -1 do
    if (len ge (2^i)) then
      quotremSi(i, pf, ii);
      ii -:= 2^i;
      pf +:= 2^i;
      len -:= 2^i;
    end if;
  end for;

  printf "Going up...\n";
  len := length;
  i := 0;
  while IsEven(len) do
    len div:= 2;
    i +:= 1;
  end while;
  pf := length;
  ii := 2^i;

  i +:= 1;
  len div:= 2;
  while i le k-1 do
    if ( (len mod 2) eq 1) then
      mulSi(i, pf - ii, ii);
      ii +:= 2^i;
    end if;
    len div:= 2;
    i +:= 1;
  end while;
end procedure;




*/
