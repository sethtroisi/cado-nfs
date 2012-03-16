#ifndef ROPT_STAGE1_H
#define ROPT_STAGE1_H

#include "ropt_param.h"
#include "ropt_tree.h"
#include "ropt_arith.h"

int ropt_stage1 ( rsstr_t rs,
                  rsparam_t rsparam,
                  sub_alpha_pq *alpha_pqueue,
                  int verbose,
                  int w );


// tune parameters for sublattices?
#if TUNE_FIND_SUBLATTICE
#include "ropt_stage2.h"
double
rsparam_tune_findlat ( rsstr_t rs,
                       rsparam_t rsparam,
                       param_t param,
                       int num_trials,
                       int w,
                       int verbose );
#endif


// rank found sublattices by test sieving?
#if TUNE_RANK_SUBLATTICE
#include "ropt_stage2.h"
void
rsparam_tune_ranklat ( rsstr_t rs,
                       rsparam_t rsparam,
                       param_t param,
                       sub_alpha_pq *alpha_pqueue,
                       int used,
                       sub_alpha_pq *tsieve_MurphyE_pqueue,
                       int verbose );
#endif


#endif /* ROPT_STAGE1_H */
