#include "cado.h"
#include <math.h>       /* rint etc */
#include "las-types.h"
#include "las-qlattice.h"
#include "las-arith.h"

/* check that the double x fits into an int32_t */
#define fits_int32_t(x) \
  ((double) INT32_MIN <= (x)) && ((x) <= (double) INT32_MAX)

/* return non-zero when the reduced lattice has entries that do not
   fit into int32_t, otherwise return 0 */
int
SkewGauss (sieve_info_ptr si, double skewness)
{
  double a[2], b[2], q, maxab0, maxab1;

  a[0] = (double) si->q;
  ASSERT_ALWAYS(a[0] < 9007199254740992.0); /* si->q should be less than 2^53
                                               so that a[0] is exact */
  b[0] = 0.0;
  a[1] = (double) si->rho;
  b[1] = skewness;
  ASSERT(b[1] != 0);
  while (1)
    {
      /* reduce vector (a[0], b[0]) with respect to (a[1], b[1]) */
      q = (a[0] * a[1] + b[0] * b[1]) / (a[1] * a[1] + b[1] * b[1]);
      q = rint (q);
      if (q == 0.0)
        break;
      a[0] -= q * a[1];
      b[0] -= q * b[1];

      /* reduce vector (a[1], b[1]) with respect to (a[0], b[0]) */
      q = (a[0] * a[1] + b[0] * b[1]) / (a[0] * a[0] + b[0] * b[0]);
      q = rint (q);
      if (q == 0.0)
        break;
      a[1] -= q * a[0];
      b[1] -= q * b[0];
    }
  if (!(fits_int32_t(a[0]) && fits_int32_t(b[0] / skewness) &&
	fits_int32_t(a[1]) && fits_int32_t(b[1] / skewness)))
    return 1;
  /* now b[0], b[1] should be of the form i*skewness, but this might not be
     exact due to rounding errors, thus we round them to the nearest integer */
  maxab0 = fabs (a[0]) > fabs (b[0]) ? fabs (a[0]) : fabs (b[0]);
  maxab1 = fabs (a[1]) > fabs (b[1]) ? fabs (a[1]) : fabs (b[1]);
  if (maxab0 <= maxab1)
    {
      si->a0 = (int32_t) a[0];
      si->b0 = (int32_t) rint (b[0] / skewness);
      si->a1 = (int32_t) a[1];
      si->b1 = (int32_t) rint (b[1] / skewness);
    }
  else /* swap (a0,b0) and (a1,b1) */
    {
      si->a1 = (int32_t) a[0];
      si->b1 = (int32_t) rint (b[0] / skewness);
      si->a0 = (int32_t) a[1];
      si->b0 = (int32_t) rint (b[1] / skewness);
      maxab1 = maxab0;
    }

  /* make sure J does not exceed I/2 */
  if (maxab1 >= si->B)
    si->J = (uint32_t) (si->B * skewness / maxab1 * (double) (si->I >> 1));
  else
    si->J = si->I >> 1;

  /* Make sure the bucket region size divides the sieve region size, 
     partly covered bucket regions may lead to problems when 
     reconstructing p from half-empty buckets. */
  /* Compute number of i-lines per bucket region, must be integer */
  ASSERT_ALWAYS(LOG_BUCKET_REGION >= si->logI);
  uint32_t i = 1U << (LOG_BUCKET_REGION - si->logI);
  i *= si->nb_threads;  /* ensures nb of bucket regions divisible by nb_threads */
  si->J = ((si->J - 1U) / i + 1U) * i; /* Round up to multiple of i */

  return 0;
}

/* return max(|a0|,|a1|)/min(|a0|,|a1|) */
double
skewness (sieve_info_ptr si)
{
  double a0 = fabs ((double) si->a0);
  double a1 = fabs ((double) si->a1);

  return (a0 > a1) ? a0 / a1 : a1 / a0;
}

