#include <stdint.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include "cado.h"
#include "utils.h"
#include "utils.h"
#include "sieving_bound.h"
#include "parameters.h"
#include "polyselect/rho.h"

void sieving_region_adapted(sieving_bound_ptr H, uint64_t number_element)
{
  double H_t = (double)number_element / pow(2, (double)(H->t - 1));
  unsigned int H0 = (unsigned int) nearbyint(pow(H_t, 1/(double)H->t));
  for (unsigned int i = 0; i < H->t; i++) {
    sieving_bound_set_hi(H, i, H0);
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
#ifndef NDEBUG
  for (unsigned int i = 0; i < nb_fields; i++) {
    ASSERT(mpz_cmp_ui(mean[i], 0) == 0);
  }
#endif

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

/* unsigned int test(unsigned int i, mpz_poly_t * f, unsigned int nb_fields, */
/*                   uint64_t number_a, sieving_bound_srcptr H) */
/* { */
/*   mpz_t  * mean = malloc(sizeof(mpz_t) * nb_fields); */
/*   mean_approx_number(mean, f, nb_fields, number_a, H); */
/*   double goal = ((double)i); */
/*   double val = 2 * mpz_get_d(); */
/*   free(mean); */
/* } */
