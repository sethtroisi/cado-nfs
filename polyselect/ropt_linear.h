#ifndef ROPT_LINEAR_H
#define ROPT_LINEAR_H

#include "ropt.h"
#include "ropt_stage1.h"
#include "ropt_stage2.h"

/* -- declarations -- */
double ropt_linear ( ropt_poly_t rs,
                     bestpoly_t bestpoly,
                     param_t param,
                     int verbose );


#endif /* ROPT_LINEAR_H */
