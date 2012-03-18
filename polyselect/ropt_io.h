#ifndef ROPT_IO_H
#define ROPT_IO_H

#include "utils.h"
#include "murphyE.h"
#include "ropt_param.h"

/* -- declarations -- */
void read_ggnfs ( mpz_t N,
                  mpz_t *f,
                  mpz_t *g,
                  mpz_t M );

#if SKIP_ROOTSIEVE_M
double print_poly_info_short ( mpz_t *f,
                               mpz_t *g,
                               int d,
                               mpz_t N );

#endif

#endif
