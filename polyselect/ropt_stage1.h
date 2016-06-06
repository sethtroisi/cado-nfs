#ifndef ROPT_STAGE1_H
#define ROPT_STAGE1_H

#include "ropt_param.h"
#include "ropt_tree.h"
#include "ropt_arith.h"


int ropt_stage1 ( ropt_poly_t poly,
                  ropt_bound_t bound,
                  ropt_s1param_t s1param,
                  ropt_param_t param,
                  alpha_pq *alpha_pqueue,
                  int current_w );


#endif /* ROPT_STAGE1_H */
