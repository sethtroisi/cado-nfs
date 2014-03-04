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
  unsigned int d, logi, i, idiv2, j,  k, l, ej, p;
  double coeff[MAX_DEGREE + 1], scale, log2max;
  unsigned char *S1 = malloc_pagealigned(1U<<LOG_BUCKET_REGION);
  unsigned char *S2 = malloc_pagealigned(1U<<LOG_BUCKET_REGION), *cS2;
  double_poly_t poly;
  
  ASSERT_ALWAYS (MIN_LOGI <= MAX_LOGI);
  ASSERT_ALWAYS (MAX_LOGI <= LOG_BUCKET_REGION);
  ASSERT_ALWAYS (S1);
  ASSERT_ALWAYS (S2);
  ASSERT_ALWAYS (!(VERT_NORM_STRIDE & (VERT_NORM_STRIDE - 1))); /* VERT_NORM_STRIDE must be = 2^x */
  poly->coeff = coeff;
  srand (getpid());
  for (logi = MIN_LOGI; logi <= MAX_LOGI; logi++) {
    i = 1U << logi;
    idiv2 = i >> 1;
    for (d = MAX_DEGREE + 1; d--;) {
      poly->deg = d;
      for (k = NB_TESTS; k--;) {
	j = rand() & (idiv2 - 1); /* Original j; 0 <= j < i/2 */
	for (l = poly->deg + 1; l--;)
	  poly->coeff[l] = (double) rand_int64_t () * rand_int64_t (); /* [(-2^63)^2, (2^63-1)^2] */
	ej = (LOG_BUCKET_REGION - logi);
	j >>= ej;                 /* In order to be able to compute a complete region */
	log2max = log2(fabs(get_maxnorm_alg (poly, (double) idiv2, (double) idiv2) + 1.));
	scale = (UCHAR_MAX - 1 - GUARD) / log2max; /* In order to have [GUARD, UCHAR_MAX-1] */
	
	/* Phase 1: compute S1 */
	memset (S1, 0, 1U<<LOG_BUCKET_REGION); /* To be sure S1 will be computed in next line */
	init_alg_norms_bucket_region_internal (S1, j, i, poly->deg, scale, poly->coeff);
	
	/* Phase 2: compute S2 */
	/* Next 4 lines are the beginning of init_alg_norms_bucket_region_internal */
	ej = 1U << ej;
	p = MIN (ej, VERT_NORM_STRIDE);
	j *= ej;
	ej += j;
	cS2 = S2;
	memset (S2, 128, 1U<<LOG_BUCKET_REGION); /* To be sure S2 will be computed */
	do {
	  double h, u[poly->deg + 1];
	  poly_scale(u, poly->coeff, poly->deg, (double) (j + (p >> 1)));
	  h = 3 - (int) idiv2;
	  for (unsigned int ih = 0; ih < i; ih += 8, h += 8.) {
	    double g = u[poly->deg];
	    for (l = poly->deg; l--; g = g * h + u[l]);
	    g = log2(fabs(g) + 1.) * scale + GUARD; /* log2(|x|+1.) >= 0 */
	    memset (cS2 + ih, (unsigned char) g, 8);
	  }
	  unsigned char *lineS2 = cS2;
	  cS2 += i;
	  for (l = p; --l; cS2 += i) memcpy (cS2, lineS2, i);
	  j += p;
	} while (j < ej);

	/* Phase 3: compare S1 & S2 */
	/* It's possible the difference is 1, because the log2 is approximative */
	for (l = 0; l < (1U<<LOG_BUCKET_REGION); l++) ASSERT_ALWAYS (abs(S1[l] - S2[l]) <= 1);
      }
    }
  }
  free(S1); free(S2);
}
