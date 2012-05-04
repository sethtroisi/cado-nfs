#ifndef ROPT_QUADRATIC_H
#define ROPT_QUADRATIC_H

#include "ropt.h"
#include "ropt_stage1.h"
#include "ropt_stage2.h"

/* -- declarations -- */
double ropt_quadratic ( ropt_poly_t rs,
                        bestpoly_t bestpoly,
                        param_t param,
                        int verbose );



#endif /* ROPT_QUADRATIC_H */
