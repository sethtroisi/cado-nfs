#include <stdlib.h>
#include <sys/types.h>
#include <unistd.h>
#include <stdio.h>
#include <stdint.h>
#include "sieve/las-norms.c"

#define NB_TESTS 100  /* NB_TESTS for a given degree & given logi */
#define MIN_LOGI 10 
#define MAX_LOGI 16   /* MIN_LOGI <= logi <= MAX_LOGI && MAX_LOGI <= LOG_BUCKET_REGION */

/* RAND_MAX is at least 15 bits, so 5 calls to rand() are need to have a random 64 bits */
static int64_t rand_int64_t () {
  return ((int64_t) rand () ^ ((int64_t) rand () << 15) ^ ((int64_t) rand () << 30) ^ ((int64_t) rand () << 45) ^ ((int64_t) rand () << 60));
}

int main() {
  unsigned int logI, k, l;
  uint32_t J, eJ;
  double coeff[2], cexp2[257], scale, iter, step, log2max, di, g, h;
  unsigned char *S1 = malloc_pagealigned((1U<<LOG_BUCKET_REGION) + MEMSET_MIN),
    *S2 = malloc_pagealigned((1U<<LOG_BUCKET_REGION) + MEMSET_MIN), *cS2;
  double_poly_t poly;
  
  ASSERT_ALWAYS (MIN_LOGI <= MAX_LOGI);
  ASSERT_ALWAYS (MAX_LOGI <= LOG_BUCKET_REGION);
  ASSERT_ALWAYS (S1);
  ASSERT_ALWAYS (S2);
  poly->coeff = coeff;
  poly->deg = 1;
  srand (getpid());
  for (logI = MIN_LOGI; logI <= MAX_LOGI; logI++) {
    const uint32_t I = 1U << logI;
    const double didiv2 = (double) (I >> 1);
    for (k = NB_TESTS; k--;) {
      J = rand () ^ (rand() << 15); /* 0 <= J < 2^30 at least */
      J &= ((I >> 1) - 1);          /* Original J; 0 <= J < I/2 */
      poly->coeff[0] = rand_int64_t () * (double) rand_int64_t (); /* [(-2^63)^2, (2^63-1)^2] */
      poly->coeff[1] = rand_int64_t () * (double) rand_int64_t (); /* [(-2^63)^2, (2^63-1)^2] */
      if (poly->coeff[1] == 0.) poly->coeff[1] = 1.;               /* Need by init_rat_... */
      eJ = (LOG_BUCKET_REGION - logI);
      J >>= eJ;                 /* In order to be able to compute a complete region */
      log2max = log2(fabs(get_maxnorm_alg (poly, didiv2, didiv2) + 1.));
      scale = (UCHAR_MAX - 1 - GUARD) / log2max; /* In order to have [GUARD, UCHAR_MAX-1] */
      step = 1. / scale;
      iter = -step * GUARD;
      for (l = 0; l < 257; iter += step) cexp2[l++] = exp2 (iter);
      
      /* Phase 1: compute S1 */
      memset (S1, 0, 1U<<LOG_BUCKET_REGION); /* To be sure S1 will be computed in next line */
      init_rat_norms_bucket_region_internal (S1, J, I, scale, poly->coeff[0], poly->coeff[1], cexp2);
      
      /* Phase 2: compute S2 */
      J <<= eJ;
      eJ = (1U << eJ) + J;
      memset (S2, 128, 1U<<LOG_BUCKET_REGION); /* To be sure S2 will be computed */
      cS2 = S2;
      do {
	g = poly->coeff[0] * J;
	for (di = -didiv2; di < didiv2; di++) {
	  /* u1 must be multiplied by di here, don't use next_h = h + u1! */
	  h = log2(fabs(g + poly->coeff[1] * di) + 1.) * scale + GUARD; /* log2(|x|+1.) >= 0 */
	  ASSERT_ALWAYS (h < (double) UCHAR_MAX);
	  *cS2++ = (unsigned char) h;
	}
      } while (++J < eJ);
      
      /* Phase 3: compare S1 & S2 */
      /* It's possible the difference is 1, because the log2 is approximative */
      for (l = 0; l < (1U<<LOG_BUCKET_REGION); l++) ASSERT_ALWAYS (abs(S1[l] - S2[l]) <= 1);
    }
  }
  free(S1); free(S2);
}
