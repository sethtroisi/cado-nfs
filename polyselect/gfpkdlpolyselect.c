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
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
//#include <time.h> 

#include "cado_poly.h"
#include "gfpkdlpolyselect.h"

#include "table_t_Py_f_deg4_type0_h1_t-200--200.c"
/*extern const row_f_poly_t tab_t_Py_f_deg4_type0_h1[26];
extern const unsigned int table_f4_size;
extern const table_f_poly_t table_f4;*/

/**
 * \brief Return appropriate degrees of polynomials f and g for CONJUGATION
 * method according to extension degree field k > 1.
 * 
 * \param[in] k integer greater than 1 ( k = 2 or 3 at the moment)
 * \param[out] deg_f degree of polynomial f, this is 2*k
 * \param[out] deg_g degree of polynomial g, this is k
 * 
 * note that deg_f >= deg_g by general convention in all this cado-nfs library.
 */
void get_degree_CONJ_f_g(unsigned int k, unsigned int *deg_f, unsigned int *deg_g){
  *deg_f = 2*k;
  *deg_g = k;
}

/**
 * \brief evaluate varphi at (u,v) and outputs a polynomial g, 
 *        assuming that Py is of degree 2.
 * 
 * \param[in] table_f structure that contains varphi and its degree
 * \param[in] u  mpz_t integer, of size half the size of p
 * \param[in] v  mpz_t integer, of size half the size of p
 *               u/v = y mod p with y a root of PY mod p.
 * \param[out] g polynomial g (already initialized)
 * 
 * NOTE: maybe input varphi and deg_varphi directly ? Or a poly structure ?
 */
// works only if PY is of degree 2
void eval_varphi_mpz(mpz_poly_t g, mpz_t** varphi_coeff, unsigned int deg_varphi, mpz_t u, mpz_t v)
{
  unsigned int i;
  for (i=0; i <= deg_varphi; i++){
    mpz_mul(g->coeff[i], v, varphi_coeff[i][0]);    // gi <- varphi_i0 * v
    mpz_addmul(g->coeff[i], u, varphi_coeff[i][1]); // gi <- gi + varphi_i1 * u
  }
}
    // t[0] + t[1]*X + ... + t[deg_varphi]*X^deg_varphi
    // with t[i] = t[i][0] + t[i][1]*Y + ... t[i][deg_Py - 1]*Y^(deg_Py-1)
    // t[i][j] are int
    // u, v are mpz_t (multi precision integers)
    // table_f.varphi[i][0] * v +
    // table_f.varphi[i][1] * u
    // varphi_i = varphi_i0 + varphi_i1 * Y
    // the function does not use modular arithmetic but exact integer arithmetic
void eval_varphi_si(mpz_poly_t g, long int** varphi_coeff, unsigned int deg_varphi, mpz_t u, mpz_t v)
{
  unsigned int i;
  for (i=0; i <= deg_varphi; i++){
    // if varphi has coefficients of type (signed) long int
    mpz_mul_si(g->coeff[i], v, varphi_coeff[i][0]);    // gi <- varphi_i0 * v
    mpz_addmul_si(g->coeff[i], u, varphi_coeff[i][1]); // gi <- gi + varphi_i1 * u
  }
}


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
			    table_f_poly_t** table_f, \
			    unsigned int* table_f_size)
{
  switch (deg_f){
  case 4:
    *table_f = &(table_f4);
    *table_f_size = table_f4_size;
    return true;
    /*
  case 6:
    table_f = TAB_f6_CYCLIC;
    *table_f_size = TAB_f6_CYCLIC_SIZE;
    return true;*/
  default:
    table_f = NULL;
    *table_f_size = 0;
    return false;
  }
}

/**
 * \brief test irreducibility over integers of polynomial varphi
 * 
 * \param[in] varphi: polynomial
 * \return true if the polynomial is irreducible
 * \return false if the polynomial is not irreducible
 * 
 * 1st version : for Degree(varphi) <= 2.
 */
bool is_irreducible_ZZ(mpz_poly_t varphi)
{
  /* if deg <= 1, return true
     if deg == 2, test discriminant 
        --> if D < 0, varphi is irreducible
        --> if D = 0, varphi is a square and is not irreducible
        --> if D > 0 is a square, then varphi has two degree 2 factors
	             otherwise varphi is irreducible
     if deg >= 3: maybe use the method of finding a factorization modulo 
     small primes, then conclude according to that
  */

  int degree, sign_D, exact_sqrt;
  bool is_irreducible = false;
  mpz_t D, tmp;
  degree = varphi->deg;
  if (degree <= 1){
    is_irreducible = true;
  }
  else{
    if (degree == 2){
      //compute discriminant, varphi = a*x^2 + b*x + c
      mpz_init(D);
      mpz_init(tmp);
      mpz_mul(D, varphi->coeff[1], varphi->coeff[1]); //   D <- b^2
      mpz_mul(tmp, varphi->coeff[0], varphi->coeff[2]); // tmp <- a*c
      mpz_mul_ui(tmp, tmp, 4);      // tmp <- 4*a*c
      sign_D = mpz_cmp(D, tmp);     // sign_D = 1 if b^2 > 4ac, 0 if b^2 = 4ac, -1 if b^2 < 4ac
      if (sign_D == 0){
	// very easy, D = 0 and varphi is reducible
	is_irreducible = true;
      }else{
	if (sign_D < 0){
	  is_irreducible = false;
	}else{
	  mpz_sub(D, D, tmp); //   D <- b^2 - 4*a*c and is > 0
	  // is_square ?
	  exact_sqrt = mpz_perfect_square_p(D);
	  if (exact_sqrt == 0){
	    is_irreducible = false;
	  }else{
	    is_irreducible = true;
	  }
	}// sign_D < 0
      }// sign_D == 0
    }// degree == 2
  }// degree <= 1
  return is_irreducible;
}

/**
 * \brief test irreducibility modulo a prime p of polynomial varphi
 * 
 * \param[in] varphi: polynomial
 * \param[in] p: (large) prime number
 * \return true if the polynomial is irreducible
 * \return false if the polynomial is not irreducible
 * 
 * 1st version : for Degree(varphi) <= 2.
 */
bool is_irreducible_mod_p(mpz_poly_t varphi, mpz_t p)
{
  /* if deg <= 1, return true
     if deg == 2, test discriminant 
        --> if D is a square mod p, (compute Legendre symbol)
            then varphi has two degree 2 factors
	    otherwise varphi is irreducible
     if deg >= 3: use MPFQ library ?
  */
  int degree, sign_D, Legendre_D_p;
  bool is_irreducible = false;
  mpz_t D, tmp;
  degree = varphi->deg;
  if (degree <= 1){
    is_irreducible = true;
  }
  else{
    if (degree == 2){
      //compute discriminant, varphi = a*x^2 + b*x + c
      mpz_init(D);
      mpz_init(tmp);
      mpz_mul(D, varphi->coeff[1], varphi->coeff[1]); //   D <- b^2
      mpz_mul(tmp, varphi->coeff[0], varphi->coeff[2]); // tmp <- a*c
      mpz_mul_ui(tmp, tmp, 4);      // tmp <- 4*a*c
      sign_D = mpz_cmp(D, tmp);     // sign_D = 1 if b^2 > 4ac, 0 if b^2 = 4ac, -1 if b^2 < 4ac
      if (sign_D == 0){
	// very easy, D = 0 and varphi is reducible
	is_irreducible = true;
      }else{
	mpz_sub(D, D, tmp); //   D <- b^2 - 4*a*c and is > 0
	// is_square mod p ?
	// compute Legendre symbol
	Legendre_D_p = mpz_legendre(D, p);
	  if (Legendre_D_p >= 0){ // indeed, (D/p) can be 0 if D is a multiple of p.
	    is_irreducible = false;
	  }else{ // (D/p) = -1
	    is_irreducible = true;
	  }
      }// sign_D == 0
    }else{
      // degree > 2, use a function that tests irreducibility. But the function does not exists yet.
      is_irreducible = false;
    }// degree == 2
  }// degree <= 1
  return is_irreducible;
}

/**
 * \brief test wether varphi is a suitable candidate for computing g from varphi
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
  // Now, find an appropriate poly f in that table.
  // start with the first one, etc because the polynomials are sorted in 
  // decreasing order of interest (decreasing order of Murphy E value
  
 
  /*
  while ()// nothing found
    // i.e. the i-th polynomial in f is not irreducible mod p,
    // or is but PY has no root.
    {

    }
  */
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
void polygen_CONJ_g ( mpz_t p, unsigned int k, mpz_poly_t f, mpz_poly_t g ){

}


/**
 * \brief The main function that computes f and g and print them in the .poly file.
 * 
 * \param[in] p multiprecision integer, prime
 * \param[in] k integer greater than 1 ( k = 2 or 3 at the moment)
 * \param[in] label tring to label the output .poly file
 * 
 */

void gfpkdlpolyselect( mpz_t p, unsigned int k, char* label){

  int dd_p, deg_f = 2*k, deg_g = 2*k;
  // take the largest possibility as default init for deg_f and deg_g.
  mpz_poly_t f, g;
  
  get_degree_CONJ_f_g(k, &deg_f, &deg_g);

  mpz_poly_init(f, deg_f);
  mpz_poly_init(g, deg_g);

  // is there a mpz_t function to compute the size in decimal digits of a number ?
  // look at Miscellanous Functions on Integers in GMP doc
  dd_p = mpz_sizeinbase(p, 10);

  char* filename = (char*)malloc(strlen(label) + 16);

  sprintf(filename, "p%ddd%d%s.poly", k, dd_p, label);
  FILE* outputpoly = fopen(filename, "w");
  if (outputpoly != NULL){
    // do stuff

    // is there a printing function for f and g ?
    // dlpolyselect.c: NO, there is 
    // print_nonlinear_poly_info(f, g, deg_f, deg_g, 1);
    // but that does NOT do what I want.
    fclose(outputpoly);
  }else{
    fprintf(stderr, "Error output file: %s .\n", filename);
  }

  mpz_poly_clear(f);
  mpz_poly_clear(g);

}

