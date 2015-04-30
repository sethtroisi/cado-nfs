#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <stdint.h>
#include "sieving_bound.h"

void sieving_region_adapted(sieving_bound_ptr H, uint64_t number_element);

void rand_mpz_poly(mpz_poly_ptr a, sieving_bound_srcptr H);

void mean_approx_number(mpz_t * mean, mpz_poly_t * f, unsigned int nb_fields,
                        uint64_t number_a, sieving_bound_srcptr H);

#endif // PARAMETERS_H
