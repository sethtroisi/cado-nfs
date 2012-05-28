#ifndef ROPT_LINEAR_H
#define ROPT_LINEAR_H

#include "ropt.h"
#include "ropt_stage1.h"
#include "ropt_stage2.h"

/* -- declarations -- */

void ropt_linear ( ropt_poly_t poly,
                   ropt_bestpoly_t bestpoly,
                   ropt_param_t param,
                   ropt_info_t info );


void ropt_linear_sieve ( ropt_poly_t poly,
                         ropt_bound_t bound,
                         ropt_s1param_t s1param,
                         ropt_param_t param,
                         ropt_info_t info,
                         alpha_pq *alpha_pqueue,
                         MurphyE_pq *global_E_pqueue );

#endif /* ROPT_LINEAR_H */
