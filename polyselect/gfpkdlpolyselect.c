/**
 * \file cado-nfs/polyselect/gfpkdlpolyselect.c
 * 
 * \date 08/08/2014
 * \author Aurore Guillevic
 * \email guillevic@lix.polytechnique.fr
 * \brief Compute two polynomials f, g suitable for discrete logarithm in 
 *        extension fields, with the conjugation method.
 *
 * \test TODO
 */


#include "cado.h"
#include "auxiliary.h"
#include "area.h"
#include "utils.h"
#include "portability.h"
#include "murphyE.h"
#include <ctype.h>
#include <stdlib.h>
#include <stdbool.h>
//#include <time.h> 

// for MurphyE value
double area=AREA;
double bound_f=BOUND_F;
double bound_g=BOUND_G;

/**
 * \brief select suitable polynomial f 
 * 
 * \param[in] p multiprecision integer, prime
 * \param[in] k integer greater than 1 ( k = 2 or 3 at the moment)
 * \param[out] f polynomial (the one of higher degree)

 */
// note: mpz_poly_t is a 1-dim 1-element array. So this is a pointer to 
// an initialized thing.
void polygen_CONJ_f ( mpz_t p, unsigned int k, mpz_poly_t f )
{
  /*which f table to choose ? */
  int** TAB = polygen_CONJ_get_tab_f(k);

  while ()// nothing found{
    
    }
}

/**
 * \brief select suitable polynomial g according to f, p, k 
 * 
 * \param[in] p multiprecision integer, prime
 * \param[in] k integer greater than 1 ( k = 2 or 3 at the moment)
 * \param[in] f polynomial (the one of higher degree)
 * \param[out] g polynomial (the one of lower degree)
 *
 */
void polygen_CONJ_g ( mpz_t p, unsigned int k, mpz_poly_t f, mpz_poly_t g )
{

}

