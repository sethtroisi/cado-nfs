#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <stdint.h>
#include "sieving_bound.h"

void sieving_region_adapted(sieving_bound_ptr H, uint64_t number_element);

void sieving_region_classical(sieving_bound_ptr H, mpz_srcptr p, unsigned int n,
    uint64_t number_element);

void rand_mpz_poly(mpz_poly_ptr a, sieving_bound_srcptr H);

void mean_approx_number(mpz_t * mean, mpz_poly_t * f, unsigned int nb_fields,
                        uint64_t number_a, sieving_bound_srcptr H);

void find_parameters_adapted(mpz_srcptr p, unsigned int n, uint64_t number_a,
    unsigned int lpb_min, unsigned int lpb_max, unsigned int t_min,
    unsigned int t_max, uint64_t size_start, mpz_poly_srcptr h, int coeff0,
    int coeff1, unsigned int nb_times, double weight_0, double weight_1);

#endif // PARAMETERS_H
