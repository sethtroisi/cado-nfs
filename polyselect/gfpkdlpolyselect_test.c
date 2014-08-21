/**
 * \file cado-nfs/polyselect/polyselect_test.c
 * 
 * \date 21/08/2014
 * \author Aurore Guillevic
 * \email guillevic@lix.polytechnique.fr
 * \brief test gfpkdlpolyselect functions
 *
 * \test ask for polynomials for p of 60, 80, 100... dd for Fp2.
 * (TODO)
 */

#include "cado.h"
#include "auxiliary.h"
#include "area.h"
#include "utils.h"
#include "portability.h"
#include "murphyE.h"
#include <ctype.h>
#include <stdlib.h>
//#include <time.h> 
#include "gfpkdlpolyselect.h"

// for MurphyE value
double area=AREA;
double bound_f=BOUND_F;
double bound_g=BOUND_G;


// maybe: in a separate file.
// gfpkdlpolyselect_utils.c --> with all the auxiliary functions
// gfpkdlpolyselect.c       --> main function (interface with Python ?)
// gfpkdlpolyselect_test.c  --> with also a main test function
int main(int argc, char* argv[])
{


}
