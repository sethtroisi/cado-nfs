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

#if WANT_TUNE // declare this if WANT_TUNE==1.
#include "ropt_stage2.h"
double
rsparam_tune ( rsstr_t rs,
               rsparam_t rsparam,
               param_t param,
               int num_trials,
               int w,
               int verbose );
#endif

#endif /* ROPT_STAGE1_H */
