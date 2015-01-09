#ifndef ROPT_H
#define ROPT_H

#include "ropt_linear.h"
#include "ropt_quadratic.h"
#include "ropt_param.h"
#include "ropt_arith.h"
#include "ropt_io.h"
#include "ropt_str.h"


/* -- declarations -- */


void ropt ( ropt_poly_t poly,
            ropt_bestpoly_t bestpoly,
            ropt_param_t param,
            ropt_info_t info );

void ropt_get_bestpoly ( ropt_poly_t poly,
                         MurphyE_pq *global_E_pqueue,
                         ropt_bestpoly_t bestpoly );

void ropt_polyselect (cado_poly_ptr output_poly, cado_poly_ptr input_poly,
                      const int effort, const int verbose);

#endif /* ROPT_H */
