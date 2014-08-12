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

#include "cado_poly.h"

#define DEG_PY 2

typedef struct {
  int t;     // parameter t
  int PY[DEG_PY + 1]; // no --> use one of the poly stuct of cado-nfs!
  int f[MAXDEGREE + 1];  // polynomial f of degree at most MAXDEGREE 
                         // set to 10 at the moment in cado_poly.h
}table_f_poly_t;

// tables containing polynomials f
// how to encode varphi form ?

#include "table_t_Py_f_deg4_type0_h1_t-200--200.c"


// for MurphyE value
double area=AREA;
double bound_f=BOUND_F;
double bound_g=BOUND_G;

/**
 * \brief return a pointer to a table [{t, PY, f}] with f of degree k.
 * 
 * \param[in] 
 * \param[in] deg_f integer greater than 1 ( deg_f = 4 or 6 implemented for now)
 * \param[out] table_f (pointer to a) table of polynomials f (the one of higher degree)
 *                     of degree deg_f, NULL if there is no such table for deg_f.
 * \param[out] table_f_size* size of the returned table 
 * 
 */
bool polygen_CONJ_get_tab_f(unsigned int deg_f, \
			    table_f_poly_t* table_f, \
			    unsigned int* table_f_size)
{
  switch (deg_f){
  case 4:
    table_f = TAB_f4_CYCLIC;
    *table_f_size = TAB_f4_CYCLIC_SIZE;
    return true;
  case 6:
    table_f = TAB_f6_CYCLIC;
    *table_f_size = TAB_f6_CYCLIC_SIZE;
    return true;
  default:
    table_f = NULL;
    *table_f_size = O;
    return false;
  }
}



/**
 * \brief test if varphi is a suitable candidate for giving g or not
 * 
 * \param[in] p multiprecision integer, prime
 * \param[in] k integer greater than 1 ( k = 2 or 3 at the moment)
 * \param[in] varphi polynomial with coefficients mod p (so very large)
 *
 * \return bool true if varphi is of degree k and irreducible mod p
 * \return false otherwise
 */
bool is_good_varphi(mpz_poly_t varphi, unsigned int k, mpz_t p)
{
  //    return (Degree(varphi_p) eq k) and IsIrreducible(varphi_p);


}

/**
 * \brief 
 * 
 * \param[in] 
 * \param[in] 
 * \param[out] 
 * \param[out] 
 * 
 */
bool is_good_f_PY()
{


}

/**
 * \brief select suitable polynomial f 
 * 
 * \param[in] p multiprecision integer, prime
 * \param[in] k integer greater than 1 ( k = 2 or 3 at the moment)
 * \param[out] f polynomial (the one of higher degree)
 * 
 */
// note: mpz_poly_t is a 1-dim 1-element array. So this is a pointer to 
// an initialized thing.
void polygen_CONJ_f ( mpz_t p, unsigned int k, mpz_poly_t f )
{
  /*which f table to choose ? */
  table_f_poly_t* table_f;
  unsigned int table_f_size;
  polygen_CONJ_get_tab_f(k, &table_f, &table_f_size);
  
  while ()// nothing found
    // i.e. the i-th polynomial in f is not irreducible mod p,
    // or is but PY has no root.
    {

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

