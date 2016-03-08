#ifndef ROPT_QUADRATIC_H
#define ROPT_QUADRATIC_H

#include "ropt.h"
#include "ropt_stage1.h"
#include "ropt_stage2.h"

/* -- declarations -- */

void ropt_quadratic ( ropt_poly_t rs,
                      ropt_bestpoly_t bestpoly,
                      ropt_param_t param,
                      ropt_info_t info );


#endif /* ROPT_QUADRATIC_H */
