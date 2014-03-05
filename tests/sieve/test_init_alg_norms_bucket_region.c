#include <stdlib.h>
#include <sys/types.h>
#include <unistd.h>
#include <stdio.h>
#include <stdint.h>
#include "sieve/las-norms.c"

#define NB_TESTS 100  /* NB_TESTS for a given degree & given logi */
#define MIN_LOGI 10  
#define MAX_LOGI 16   /* MIN_LOGI <= logi <= MAX_LOGI && MAX_LOGI <= LOG_BUCKET_REGION */
#define MAX_DEGREE 10 /* All degree between 0 and MAX_DEGREE included are tested */

/* RAND_MAX is at least 15 bits, so 5 calls to rand() are need to have a random 64 bits */
static int64_t rand_int64_t () {
  return ((int64_t) rand () ^ ((int64_t) rand () << 15) ^ ((int64_t) rand () << 30) ^ ((int64_t) rand () << 45) ^ ((int64_t) rand () << 60));
}

int main() {
  unsigned int d, logI, k, l, p, ih;
  uint32_t J, eJ;
  double coeff[MAX_DEGREE + 1], u[MAX_DEGREE + 1], scale, log2max, h, g;
  unsigned char *S1 = malloc_pagealigned((1U<<LOG_BUCKET_REGION) + MEMSET_MIN), 
    *S2 = malloc_pagealigned((1U<<LOG_BUCKET_REGION) + MEMSET_MIN), *cS2, *lineS2;
  double_poly_t poly;
  
  ASSERT_ALWAYS (MIN_LOGI <= MAX_LOGI);
  ASSERT_ALWAYS (MAX_LOGI <= LOG_BUCKET_REGION);
  ASSERT_ALWAYS (S1);
  ASSERT_ALWAYS (S2);
  ASSERT_ALWAYS (!(VERT_NORM_STRIDE & (VERT_NORM_STRIDE - 1))); /* VERT_NORM_STRIDE must be = 2^x */
  poly->coeff = coeff;
  srand (getpid());
  for (logI = MIN_LOGI; logI <= MAX_LOGI; logI++) {
    const uint32_t I = 1U << logI;
    const double didiv2 = (double) (I >> 1);
    for (d = MAX_DEGREE + 1; d--;) {
      poly->deg = d;
      for (k = NB_TESTS; k--;) {
	J = rand () ^ (rand() << 15); /* 0 <= J < 2^30 at least */
	J &= (I >> 1) - 1;  /* Original J; 0 <= J < I/2 */
	for (l = poly->deg + 1; l--;) /* (-2^63)^2 <= coeff <= (2^63-1)^2 */
	  poly->coeff[l] = rand_int64_t () * (double) rand_int64_t ();
	eJ = (LOG_BUCKET_REGION - logI);
	J >>= eJ;                 /* In order to be able to compute a complete region */
	log2max = log2(fabs(get_maxnorm_alg (poly, didiv2, didiv2) + 1.));
	scale = (UCHAR_MAX - 1 - GUARD) / log2max; /* In order to have [GUARD, UCHAR_MAX-1] */
	
	/* Phase 1: compute S1 */
	memset (S1, 0, 1U<<LOG_BUCKET_REGION); /* To be sure S1 will be computed in next line */
	init_alg_norms_bucket_region_internal (S1, J, I, poly->deg, scale, poly->coeff);
	
	/* Phase 2: compute S2 */
	J <<= eJ;
	eJ = 1U << eJ;
	p = MIN (eJ, VERT_NORM_STRIDE);
	eJ += J;
	cS2 = S2;
	memset (S2, 128, 1U<<LOG_BUCKET_REGION); /* To be sure S2 will be computed */
	do {
	  poly_scale(u, poly->coeff, poly->deg, (double) (J + (p >> 1)));
	  h = 3. - didiv2;
	  for (ih = 0; ih < I; ih += 8, h += 8.) {
	    for (g = u[poly->deg], l = poly->deg; l--; g = g * h + u[l]);
	    g = log2(fabs(g) + 1.) * scale + GUARD; /* log2(|x|+1.) >= 0 */
	    ASSERT_ALWAYS (g < (double) UCHAR_MAX);
	    memset (cS2 + ih, (unsigned char) g, 8);
	  }
	  lineS2 = cS2;
	  cS2 += I;
	  for (l = p; --l; cS2 += I) memcpy (cS2, lineS2, I);
	  J += p;
	} while (J < eJ);

	/* Phase 3: compare S1 & S2 */
	/* It's possible the difference is 1, because the log2 is approximative */
	for (l = 0; l < (1U<<LOG_BUCKET_REGION); l++) ASSERT_ALWAYS (abs(S1[l] - S2[l]) <= 1);
      }
    }
  }
  free(S1); free(S2);
}
