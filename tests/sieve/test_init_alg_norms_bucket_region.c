#define NB_TESTS 10   /* NB_TESTS for a given degree & given logi */
#define MIN_LOGI 10  
#define MAX_LOGI 16   /* MIN_LOGI <= logi <= MAX_LOGI && MAX_LOGI <= LOG_BUCKET_REGION */
#define MAX_DEGREE 10 /* All degree between 0 and MAX_DEGREE included are tested */
#define MAX_ERR 2     /* The error must be < MAX_ERR */

#include <stdlib.h>
#include <sys/types.h>
#include <unistd.h>
#include <stdio.h>
#include <stdint.h>

#include "sieve/las-config.h"
#include "sieve/las-debug.h"
#include "sieve/las-norms.h"
#include "utils/utils.h"
#include "portability.h"

/* RAND_MAX is at least 15 bits, so 5 calls to rand() are need to have a random 64 bits */
static int64_t rand_int64_t () {
  return ((int64_t) rand () ^ ((int64_t) rand () << 15) ^ ((int64_t) rand () << 30) ^ ((int64_t) rand () << 45) ^ ((int64_t) rand () << 60));
}

static inline void
my_poly_scale (double *u, const double *t, unsigned int d, double h)
{
  double hpow;
  u[d] = t[d];
  for (hpow = h; d--; hpow *= h) u[d] = t[d] * hpow;
}

int main() {
  unsigned int d, logI, k, l, p, ih;
  uint32_t original_J, J, eJ;
  double coeff[MAX_DEGREE + 1], u[MAX_DEGREE + 1], scale, log2max, h, g,
    *S0 = malloc_pagealigned((1U<<LOG_BUCKET_REGION) * sizeof (*S0)), *cS0;
  unsigned char *S1 = malloc_pagealigned((1U<<LOG_BUCKET_REGION) + MEMSET_MIN), 
    *S2 = malloc_pagealigned((1U<<LOG_BUCKET_REGION) + MEMSET_MIN), *cS2, *lineS2;
  double_poly_t poly;
  unsigned long count, err[MAX_ERR], real_err[UCHAR_MAX + 1 - GUARD];
  
  memset(err, 0, sizeof(err));
  memset(real_err, 0, sizeof(real_err));
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
	original_J = J;
	eJ = 1U << eJ;
	p = MIN (eJ, VERT_NORM_STRIDE);
	eJ += J;
	cS2 = S2;
	memset (S2, 128, 1U<<LOG_BUCKET_REGION); /* To be sure S2 will be computed */
	do {
	  my_poly_scale (u, poly->coeff, d, (double) (J + (p >> 1)));
	  h = 3. - didiv2;
	  for (ih = 0; ih < I; ih += 8, h += 8.) {
	    for (g = u[poly->deg], l = poly->deg; l--; g = g * h + u[l]);
	    g = fabs(g) + 1.;
	    g = log2(g) * scale + GUARD; /* log2(|x|+1.) >= 0 */
	    ASSERT_ALWAYS (g < (double) UCHAR_MAX);
	    memset (cS2 + ih, (unsigned char) g, 8);
	  }
	  lineS2 = cS2;
	  cS2 += I;
	  for (l = p; --l; cS2 += I) memcpy (cS2, lineS2, I);
	  J += p;
	} while (J < eJ);

	/* Phase 2 bis: compute real S0, in order to have an idea of the
	   difference of 8*VERT_NORM_STRIDE decimal real values and the integer
	   value used in S1 to approximate them */
	cS0 = S0;
	for (J = original_J; J < eJ; J++) {
	  my_poly_scale (u, poly->coeff, d, (double) (J + (p >> 1)));
	  for (h = -didiv2; h < didiv2; h++) {
	    for (g = u[poly->deg], l = poly->deg; l--; g = g * h + u[l]);
	    g = fabs(g) + 1.;
	    g = log2(g) * scale + GUARD; /* log2(|x|+1.) >= 0 */
	    ASSERT_ALWAYS (g < (double) UCHAR_MAX);
	    *cS0++ = g;
	  }
	}

	/* Phase 3: compare S1 & S2 */
	/* It's possible the difference is 1, because the log2 is approximative */
	for (l = 0; l < (1U<<LOG_BUCKET_REGION); l++) {
	  unsigned char c = abs(S1[l] - S2[l]);
	  ASSERT_ALWAYS (c < MAX_ERR);
	  err[c]++;
	}
	
	/* Phase 3 bis: compare S0 & S1 */
	/* The difference here could be very important, because the init algorithm
	   approximates 8*VERT_NORM_STRIDE values by only one value (computed in the
	   middle, 4 and VERT_NORM_STRIDE>>1). So, in the neighborhood of a polynome
	   root, the difference could be as large as UCHAR_MAX - GUARD - 1. */
	for (l = 0; l < (1U<<LOG_BUCKET_REGION); l++) {
	  unsigned int c = abs((unsigned int) S0[l] - (unsigned int) S1[l]);
	  ASSERT_ALWAYS (c < (UCHAR_MAX - GUARD));
	  real_err[c]++;
	}
      }
    }
  }
  free(S0); free(S1); free(S2);
  count = (MAX_LOGI - MIN_LOGI + 1) * (MAX_DEGREE + 1) * (NB_TESTS << LOG_BUCKET_REGION);
  fprintf (stderr, "The tests are done with %lu values. All the values are between %u and %u included.\n", count, GUARD, UCHAR_MAX - 1);
  double lc = 100. / count;
  fprintf (stderr, "\n1. Smart (fast, rough) algorithm for both. The difference is the approximate log2 versus exact log2.\nInteger differences between integer values:\n");
  for (l = 0; l < MAX_ERR; l++) if (err[l]) fprintf (stderr, "%u: %4.1f%%\n", l, err[l] * lc);
  fprintf (stderr, "\n2. Smart algorithm + approximate log2 versus complete exact algorithm:\nInteger differences between exact integer and approximate integer values:\n");
  count = 0.;
  for (l = 0; l < MAX_ERR; l++) if (real_err[l]) {
      fprintf (stderr, "%3u: %6.3f%%\n", l, real_err[l] * lc);
      count += real_err[l];
    }
  fprintf (stderr, "===> The percentage of the differences smaller than the tolerate error (%u) is: %6.3f%%\n", MAX_ERR, count * lc);
  for (; l < UCHAR_MAX + 1 - GUARD; l++)
    if (real_err[l]) fprintf (stderr, "%3u :%10.7f%%, %6lu values\n", l, real_err[l] * lc, real_err[l]);
}
