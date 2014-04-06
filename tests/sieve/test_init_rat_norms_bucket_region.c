#define NB_TESTS 100  /* NB_TESTS for a given degree & given logi */
#define MIN_LOGI 10 
#define MAX_LOGI 16   /* MIN_LOGI <= logi <= MAX_LOGI && MAX_LOGI <= LOG_BUCKET_REGION */
#define MAX_ERR  2    /* The error must be < MAX_ERR */

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

int main() {
  unsigned int logI, k, l;
  uint32_t J, eJ;
  double coeff[2], cexp2[257], scale, iter, step, log2max, di, g, h;
  unsigned char *S1 = malloc_pagealigned((1U<<LOG_BUCKET_REGION) + MEMSET_MIN),
    *S2 = malloc_pagealigned((1U<<LOG_BUCKET_REGION) + MEMSET_MIN), *cS2, *cS1;
  double_poly_t poly;
  unsigned long count, err[MAX_ERR], precision_err[MAX_ERR * 10];

  memset(err, 0, sizeof(err));
  memset(precision_err, 0, sizeof(precision_err));
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
      
      /* Phase 2: compute S2 & compare S1 & S2 */
      J <<= eJ;
      eJ = (1U << eJ) + J;
      memset (S2, 128, 1U<<LOG_BUCKET_REGION); /* To be sure S2 will be computed */
      cS2 = S2;
      cS1 = S1;
      do {
	g = poly->coeff[0] * J;
	for (di = -didiv2; di < didiv2; di++) {
	  /* u1 must be multiplied by di here, don't use next_h = h + u1! */
	  h = log2(fabs(g + poly->coeff[1] * di) + 1.) * scale + GUARD; /* log2(|x|+1.) >= 0 */
	  ASSERT_ALWAYS (h < (double) UCHAR_MAX);
	  unsigned char my_cS1 = *cS1++, my_cS2 = (unsigned char) h, c = abs(my_cS1 - my_cS2);
	  unsigned int y = abs ((int) ((h - (double) my_cS1) * 10.));
	  ASSERT_ALWAYS (c < MAX_ERR);
	  *cS2++ = my_cS2;
	  err[c]++;
	  precision_err[y]++;
	}
      } while (++J < eJ);
    }
  }
  free(S1); free(S2);
  count = (MAX_LOGI - MIN_LOGI + 1) * (NB_TESTS <<LOG_BUCKET_REGION);
  double lc = 100. / count;
  fprintf (stderr, "The tests are done with %lu values.\n", count);
  fprintf (stderr, "All the values are between %u and %u included.\n", GUARD, UCHAR_MAX - 1);
  fprintf (stderr, "The maximal difference between the values computed by exact and smart algorithms is less than %u.\n", MAX_ERR);
  fprintf (stderr, "The integer differences between real integer and roughly integer values are:\n");
  for (l = 0; l < MAX_ERR; l++)
    if (err[l]) fprintf (stderr, "%u: %4.1f%%\n", l, err[l] * lc);
  fprintf (stderr, "The decimal differences between real decimal and roughly integer values are:\n");
  count = 0;
  for (l = 0; l < 10; l++) if (precision_err[l]) {
      fprintf (stderr, "%.1f <= difference < %.1f: %6.3f%%\n",
	       (double) l / 10., (double) (l + 1) / 10., (double) precision_err[l] * lc);
      count += precision_err[l];
    }
  fprintf (stderr, "===> The percentage of the differences smaller than 1 is: %6.3f%%\n", count * lc);
  for (; l < 10 * MAX_ERR; l++)
    if (precision_err[l]) fprintf (stderr, "%.1f <= difference < %.1f: %6.3f%%\n",
			  (double) l / 10., (double) (l + 1) / 10., (double) precision_err[l] * lc);
}
