#ifndef ROPT_H
#define ROPT_H

#include "ropt_linear.h"
#include "ropt_quadratic.h"

/* -- declarations -- */

int cachesize_cpuid(int verbose);

int cachesize_guess(int verbose);

void ropt ( ropt_poly_t rs,
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

void ropt_return_bestpoly ( ropt_poly_t rs,
                            MurphyE_pq *global_E_pqueue,
                            bestpoly_t bestpoly );


#endif /* ROPT_H */
