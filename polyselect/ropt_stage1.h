#ifndef ROPT_STAGE1_H
#define ROPT_STAGE1_H

#include "ropt_param.h"
#include "ropt_tree.h"
#include "ropt_arith.h"

int ropt_stage1 ( ropt_poly_t rs,
                  ropt_param_t rsparam,
                  sub_alpha_pq *alpha_pqueue,
                  int verbose,
                  int w );


// tune parameters for sublattices?
#if TUNE_FIND_SUBLATTICE
#include "ropt_stage2.h"
double
ropt_param_tune_findlat ( ropt_poly_t rs,
                       ropt_param_t rsparam,
                       param_t param,
                       int num_trials,
                       int w,
                       int verbose );
#endif


// rank found sublattices by test sieving?
#if TUNE_RANK_SUBLATTICE
#include "ropt_stage2.h"
void
ropt_param_tune_ranklat ( ropt_poly_t rs,
                       ropt_param_t rsparam,
                       param_t param,
                       sub_alpha_pq *alpha_pqueue,
                       int used,
                       sub_alpha_pq *tsieve_MurphyE_pqueue,
                       int verbose );
#endif


#endif /* ROPT_STAGE1_H */
