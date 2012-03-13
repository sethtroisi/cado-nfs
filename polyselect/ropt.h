#ifndef ROPT_H
#define ROPT_H

#include "ropt_stage1.h"
#include "ropt_stage2.h"

/* -- declarations -- */
void ropt_main ( rsstr_t rs,
                 bestpoly_t bestpoly,
                 param_t param,
                 int verbose );

void ropt_polyselect ( mpz_t *f,
                       int d,
                       mpz_t m,
                       mpz_t l,
                       mpz_t N,
                       int max_k,
                       int verbose );

#endif /* ROPT_H */
