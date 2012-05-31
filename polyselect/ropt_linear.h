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

#endif /* ROPT_LINEAR_H */
