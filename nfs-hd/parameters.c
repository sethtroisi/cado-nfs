#include <stdint.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include "cado.h"
#include "utils.h"
#include "sieving_bound.h"
#include "parameters.h"
#include "polyselect/rho.h"
#include "polynomials.h"

/* ---------- */

/* From sieve/strategies/gen_decomp.c (07/05/2015)*/
/* number of primes <= 2^n */
/* All values below fit in 53-bit mantissa, so that makes well-defined
 * floating point literals.
 */
double A7053[256] = {0,1,2,4,6,11,18,31,54,97,172,309,564,1028,1900,
  3512,6542,12251,23000,43390,82025,155611,295947,
  564163,1077871,2063689,3957809,7603553,14630843,
  28192750,54400028,105097565,203280221,393615806,
  762939111,1480206279,2874398515,5586502348,
  10866266172,21151907950,41203088796,80316571436,
  156661034233,305761713237,597116381732,
  1166746786182,2280998753949,4461632979717,
  8731188863470,17094432576778,33483379603407,
  65612899915304,128625503610475,0,};

/* return prime_pi(2^i) */
static double prime_pi (unsigned long i)
{
  if (A7053[i] != 0)
    return (double) A7053[i];
  else
    {
      double x = ldexp (1.0, i);
      return x / log (x);
    }
}

/* ---------- */

uint64_t size_sieving_region(sieving_bound_srcptr H)
{
  uint64_t res = 1;
  for (unsigned int i = 0; i < H->t - 1; i++) {
    res = res * 2 * (uint64_t) H->h[i];
  }
  res = res * (uint64_t) H->h[H->t - 1];

  ASSERT(res > 0);

  return res;
}

void sieving_region_adapted(sieving_bound_ptr H, uint64_t number_element)
{
  double tmp = (double)number_element / pow(2, (double)(H->t - 1));
  unsigned int H0 = (unsigned int) nearbyint(pow(tmp, 1/(double)H->t));
  for (unsigned int i = 0; i < H->t; i++) {
    sieving_bound_set_hi(H, i, H0);
  }
}

void sieving_region_classical(sieving_bound_ptr H, mpz_srcptr p, unsigned int n,
    uint64_t number_element)
{
  double tmp = (double)number_element / pow(2, (double)(H->t - 1));
  double p_d = mpz_get_d(p);
  double power = ( (double) (H->t - H->t * H->t) ) / (2.0 * (double) n);
  tmp = tmp / pow(p_d, power);
  unsigned int H0 = (unsigned int) nearbyint(pow(tmp, 1/(double)H->t));
 
  sieving_bound_set_hi(H, 0, H0);
  for (unsigned int i = 1; i < H->t - 1; i++) {
    tmp = H0 * pow(p_d, - (double) i / (double) n);
    sieving_bound_set_hi(H, i, (unsigned int) nearbyint(tmp));
  }
  tmp = H0 * pow(p_d, - (double) (H->t - 1) / (double) n);
  if ((unsigned int) nearbyint(tmp) == 0) {
    sieving_bound_set_hi(H, H->t - 1, 1);
  } else {
    sieving_bound_set_hi(H, H->t - 1, (unsigned int) nearbyint(tmp));
  }
}

void rand_mpz_poly(mpz_poly_ptr a, sieving_bound_srcptr H)
{
  for (int i = 0; i < (int)H->t - 1; i++) {
    mpz_poly_setcoeff_si(a, i, rand() % (2 * (int)H->h[i]) - (int)H->h[i]);
  }
  mpz_poly_setcoeff_si(a, (int)H->t - 1, rand() % ((int)H->h[H->t - 1]));

#ifndef NDEBUG
  mpz_t tmp;
  mpz_init(tmp);
  for (int i = 0; i < a->deg; i++) {
    mpz_poly_getcoeff(tmp, i, a);
    mpz_abs(tmp, tmp);
    ASSERT(mpz_cmp_ui(tmp, H->h[i]) <= 0);
  }
  if (a->deg == (int)H->t - 1) {
    mpz_poly_getcoeff(tmp, H->t - 1, a);
    ASSERT(mpz_cmp_ui(tmp, H->h[H->t - 1]) <= 0);
    ASSERT(mpz_cmp_ui(tmp, 0) >= 0);
    mpz_clear(tmp);
  }
#endif
}

void mean_approx_number(mpz_t * mean, mpz_poly_t * f, unsigned int nb_fields,
    uint64_t number_a, sieving_bound_srcptr H)
{
  for (unsigned int i = 0; i < nb_fields; i++) {
    mpz_set_ui(mean[i], 0);
  }

  srand(time(NULL));
  mpz_poly_t a;
  mpz_poly_init(a, (int)H->t);

  mpz_t res;
  mpz_init(res);
  for (uint64_t i = 0; i < number_a; i++) {
    rand_mpz_poly(a, H);
    for(unsigned int j = 0; j < nb_fields; j++) {
      mpz_poly_resultant(res, f[j], a);
      mpz_abs(res, res);
      mpz_add(mean[j], mean[j], res);
    }
  }

  for(unsigned int i = 0; i < nb_fields; i++) {
    mpz_cdiv_q_ui(mean[i], mean[i], number_a);
  }

  mpz_clear(res);
  mpz_poly_clear(a);
}

static unsigned int test(mpz_t * mean, unsigned int i, mpz_poly_t * f,
    unsigned int nb_fields, uint64_t number_a, sieving_bound_srcptr H,
    unsigned int lpb)
{
  ASSERT(nb_fields == 2);

  mean_approx_number(mean, f, nb_fields, number_a, H); 
  double goal = ((double)i); 
  double val = log2( 2 * prime_pi((unsigned long)i) );
  double smoothness = dickman_rho(log2(mpz_get_d(mean[0])) / (double)lpb);
  smoothness = smoothness * dickman_rho((log2(mpz_get_d(mean[1])) -
        (double) lpb - 1) / (double)lpb);
  val = val - log2( smoothness );
  val = ceil(val);
  if (val <= goal) {
    return 1;
  }
  return 0;
} 

void find_parameters_adapted(mpz_srcptr p, unsigned int n, uint64_t number_a,
    unsigned int lpb_min, unsigned int lpb_max, unsigned int t_min,
    unsigned int t_max, uint64_t size_start, mpz_poly_srcptr h, int coeff0,
    int coeff1, unsigned int nb_times, double weight_0, double weight_1)
{
  sieving_bound_t H;
  mpz_t a;
  mpz_t b;
  mpz_t c;
  mpz_init(a);
  mpz_init(b);
  mpz_init(c);
  mpz_poly_t g;
  mpz_poly_init(g, n);
  mpz_poly_t f0;
  mpz_poly_init(f0, n);
  mpz_poly_t f1;
  mpz_poly_init(f1, n);
  mpz_poly_t * f = malloc(2 * sizeof(mpz_poly_t));
  mpz_t * mean = malloc(2 * sizeof(mpz_t));
  for (unsigned int k = 0; k < 2; k++) {
    mpz_poly_init(f[k], n);
    mpz_init(mean[k]);
  }

  for (unsigned int lpb = lpb_min; lpb < lpb_max; lpb++) {
    for (unsigned int t = t_min; t < t_max; t++) {
      sieving_bound_init(H, t);
      uint64_t size_current = size_start;

      sieving_region_adapted(H, size_current);
      function_special_q(f0, f1, g, a, b, c, p, h, coeff0, coeff1, lpb - 1, t,
          number_a, weight_0, weight_1);
      mpz_poly_set(f[0], f0);
      mpz_poly_set(f[1], f1);
      
      test(mean, size_current, f, 2, nb_times, H, lpb);
      sieving_bound_clear(H);
    }
  }
  mpz_clear(a);
  mpz_clear(b);
  mpz_clear(c);
  mpz_poly_clear(g);
  mpz_poly_clear(f0);
  mpz_poly_clear(f1);
  for (unsigned int k = 0; k < 2; k++) {
    mpz_poly_clear(f[k]);
    mpz_clear(mean[k]);
  }
  free(f);
  free(mean);
}

