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
 * \version 1: dedicated for GF(p^2).
 */

#include "cado.h"
#include <stdio.h>
#include "portability.h"
#include "utils.h"

#include "gfpkdlpolyselect.h"

#include "table_f_Py_phi__f_deg4_s02_C4V4_h1_Py_s20_f01.c"
// contains: fPyphi_poly_t list_f4

/**
 * \brief Return appropriate degrees of polynomials f and g for CONJUGATION
 * method according to extension degree field n > 1.
 * 
 * \param[in] n integer greater than 1 ( n = 2 at the moment)
 * \param[out] deg_f degree of polynomial f, this is 2*n
 * \param[out] deg_g degree of polynomial g, this is n
 * 
 * note that deg_f >= deg_g by general convention in all this cado-nfs library.
 */
void get_degree_CONJ_f_g(unsigned int n, unsigned int *deg_f, unsigned int *deg_g){
    *deg_f = 2*n;
    *deg_g = n;
}

/**
 * \brief evaluate phi at (u,v) and outputs a polynomial g, 
 *        assuming that Py is of degree 2.
 * 
 * \param[in] phi_coeff double entry array. The poly phi is expressed as (ZZ[y])[X]
 *                         i.e. sum_{i=0}^{deg phi} ((a_{i0} + a_{i1}*y)*X^i)
 * \param[in] deg_phi degree of polynomial phi, for the moment this must be less than 2.
 * \param[in] u  mpz_t integer, of size half the size of p
 * \param[in] v  mpz_t integer, of size half the size of p
 *               u/v = y mod p with y a root of PY mod p.
 * \param[out] g polynomial g of same degree as phi (already initialized)
 * 
 * NOTE: maybe input phi and deg_phi directly ? Or a poly structure ?
 */
// works only if PY is of degree 2
void eval_mpz_phi_mpz_uv(mpz_poly_t g, mpz_t** phi_coeff, unsigned int deg_phi, mpz_t u, mpz_t v)
{
  unsigned int i;
  for (i=0; i <= deg_phi; i++){
    mpz_mul(g->coeff[i], v, phi_coeff[i][0]);    // gi <- phi_i0 * v
    mpz_addmul(g->coeff[i], u, phi_coeff[i][1]); // gi <- gi + phi_i1 * u
  }
}
    // t[0] + t[1]*X + ... + t[deg_phi]*X^deg_phi
    // with t[i] = t[i][0] + t[i][1]*Y + ... t[i][deg_Py - 1]*Y^(deg_Py-1)
    // t[i][j] are int
    // u, v are mpz_t (multi precision integers)
    // list_f.phi[i][0] * v +
    // list_f.phi[i][1] * u
    // phi_i = phi_i0 + phi_i1 * Y
    // the function does not use modular arithmetic but exact integer arithmetic
void eval_si_phi_mpz_uv(mpz_poly_t g, long int** phi_coeff, unsigned int deg_phi, mpz_t u, mpz_t v)
{
  unsigned int i;
  for (i=0; i <= deg_phi; i++){
    // if phi has coefficients of type (signed) long int
    mpz_mul_si(g->coeff[i], v, phi_coeff[i][0]);    // gi <- phi_i0 * v
    mpz_addmul_si(g->coeff[i], u, phi_coeff[i][1]); // gi <- gi + phi_i1 * u
  }
}

void eval_si_phi_mpz_y(mpz_poly_t g, long int** phi_coeff, unsigned int deg_phi, mpz_t y)
{
  unsigned int i;
  for (i=0; i <= deg_phi; i++){
    // if phi has coefficients of type (signed) long int
    mpz_set_si(g->coeff[i], phi_coeff[i][0]);    // gi <- phi_i0
    mpz_addmul_si(g->coeff[i], y, phi_coeff[i][1]); // gi <- phi_i0 + phi_i1 * y
  }
}

/**
 * \brief return a pointer to a table [{t, PY, f}] with f of degree n.
 * 
 * \param[in] 
 * \param[in] deg_f integer greater than 1 ( deg_f = 4 or 6 implemented for now)
 * \param[out] list_f (pointer to a) table of polynomials f (the one of higher degree)
 *                     of degree deg_f, NULL if there is no such table for deg_f.
 * \param[out] list_f_size* size of the returned table 
 * 
 */
bool polygen_CONJ_get_tab_f(unsigned int deg_f, \
			    fPyphi_poly_t** list_f)
{
  switch (deg_f){
  case 4:
    *list_f = &(list_f4);
    return true;
    /*
  case 6:
    list_f = TAB_f6_CYCLIC;
    return true;*/
  default:
    list_f = NULL;
    return false;
  }
}

// [almost] same as is_good_poly in gfpk/magma/polyselect_utils.mag
// return: error_code, see above
//int is_good_poly_light(pp_t params, ppf_t poly_params, mpz_poly_t f){}

  /* check: 
     + deg f
     + f irreducible over Q
     + f has a degree params->n factor mod params->p
     o [ Galois group of f: assume that the table is already checked.]
     o [ h=1 ]
     + signature(f) : how to do that in C ??? [look at the Cohen, then ask Pierrick and Paul.]
     o signature subfield: how to find the subfield ? [id]
     o badideals: call the magma function...
     o smexp: does this already exists in C ? call magma maybe.
  */


// same as is_good_poly_check_all in gfpk/magma/polyselect_utils.mag
// return: error_code, see above
//int is_good_poly(pp_t params, ppf_t poly_params, mpz_poly_t f){}

  /* check: 
     - deg f
     - f irreducible over Q
     - f has a degree params->n factor mod params->p
       for CONJ: check 
     - [ Galois group of f: assume that the table is already checked.]
     - [ h=1 ]
     - signature(f) : how to do that in C ??? [look at the Cohen, then ask Pierrick and Paul.]
     - signature subfield: how to find the subfield ? [id]
     - badideals: call the magma function...
     - smexp: does this already exists in C ? call magma maybe.
  */


/**
 * \brief 
 * 
 * \param[in] list_f a table of suitable polynomials f
 * \param[in] list_f_size the size of the table, so ::index should be strictly less than that
 *
 * \return index the index of the next suitable f in the table, if not found, then 
 *               index = ::list_f_size
 * 
 * This function returns the index of the first suitable f, starting at index = list_f_size
 * A suitable f is such that
 * f is irreducible over the integers ZZ
 * f has an irreducible degree n factor modulo p
 * By design, the table is made of polynomial pairs (f, Py), so that
 * if Py (of degree 2) has a root mod p, then f factors and has a degree n factor phi.
 *
 */
//unsigned int get_index_next_poly_in_tab_f(ppf_t f, pp_t pp, fPyphi_poly_t * list_f){}


/**
 * \brief test irreducibility over integers of polynomial phi
 * 
 * \param[in] phi: polynomial
 * \return true if the polynomial is irreducible
 * \return false if the polynomial is not irreducible
 * 
 * 1st version : for Degree(phi) <= 2.
 */
bool is_irreducible_ZZ(mpz_poly_t phi)
{
  /* if deg <= 1, return true
     if deg == 2, test discriminant 
        --> if D < 0, phi is irreducible
        --> if D = 0, phi is a square and is not irreducible
        --> if D > 0 is a square, then phi has two degree 2 factors
	             otherwise phi is irreducible
     if deg >= 3: maybe use the method of finding a factorization modulo 
     small primes, then conclude according to that
  */

  int degree, sign_D, exact_sqrt;
  bool is_irreducible = false;
  mpz_t D, tmp;
  degree = phi->deg;
  if (degree <= 1){
    is_irreducible = true;
  }
  else{
    if (degree == 2){
      //compute discriminant, phi = a*x^2 + b*x + c
      mpz_init(D);
      mpz_init(tmp);
      mpz_mul(D, phi->coeff[1], phi->coeff[1]); //   D <- b^2
      mpz_mul(tmp, phi->coeff[0], phi->coeff[2]); // tmp <- a*c
      mpz_mul_ui(tmp, tmp, 4);      // tmp <- 4*a*c
      sign_D = mpz_cmp(D, tmp);     // sign_D = 1 if b^2 > 4ac, 0 if b^2 = 4ac, -1 if b^2 < 4ac
      if (sign_D == 0){
	// very easy, D = 0 and phi is reducible
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
 * \brief test irreducibility modulo a prime p of polynomial phi
 * 
 * \param[in] phi: polynomial with large coefficients
 * \param[in] p: (large) prime number
 * \return true if the polynomial is irreducible
 * \return false if the polynomial is not irreducible
 * 
 * 1st version : for Degree(phi) <= 2.
 */
bool is_irreducible_mod_p(mpz_poly_t phi, mpz_t p)
{
  /* if deg <= 1, return true
     if deg == 2, test discriminant 
        --> if D is a square mod p, (compute Legendre symbol)
            then phi has two degree 2 factors
	    otherwise phi is irreducible
     if deg >= 3: use MPFQ library ?
  */
  int degree, sign_D, Legendre_D_p;
  bool is_irreducible = false;
  mpz_t D, tmp;
  degree = phi->deg;
  if (degree <= 1){
    is_irreducible = true;
  }
  else{
    if (degree == 2){
      //compute discriminant, phi = a*x^2 + b*x + c
      mpz_init(D);
      mpz_init(tmp);
      mpz_mul(D, phi->coeff[1], phi->coeff[1]); //   D <- b^2
      mpz_mul(tmp, phi->coeff[0], phi->coeff[2]); // tmp <- a*c
      mpz_mul_ui(tmp, tmp, 4);      // tmp <- 4*a*c
      sign_D = mpz_cmp(D, tmp);     // sign_D = 1 if b^2 > 4ac, 0 if b^2 = 4ac, -1 if b^2 < 4ac
      if (sign_D == 0){
	// very easy, D = 0 and phi is reducible
	is_irreducible = false;
      }else{
	mpz_sub(D, D, tmp); //   D <- b^2 - 4*a*c and is != 0
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

// same as above but with a poly with "small" coeffs given in a tab, with size <-> deg.
bool is_irreducible_mod_p_si(long int * f, int deg_f, mpz_t p){
  int sign_D, Legendre_D_p;
  bool is_irreducible = false;
  mpz_t D, tmp;
  mpz_t a, b, c;
  if (deg_f <= 1){
    is_irreducible = true;
  }
  else{
    if (deg_f == 2){
      //compute discriminant, phi = a*x^2 + b*x + c
      mpz_init(D);
      mpz_init(tmp);
      mpz_init_set_si(a, f[0]);
      mpz_init_set_si(b, f[1]);
      mpz_init_set_si(c, f[2]);
      mpz_mul(D, b, b); //   D <- b^2
      mpz_mul(tmp, a, c); // tmp <- a*c
      mpz_mul_ui(tmp, tmp, 4);      // tmp <- 4*a*c
      sign_D = mpz_cmp(D, tmp);     // sign_D = 1 if b^2 > 4ac, 0 if b^2 = 4ac, -1 if b^2 < 4ac
      if (sign_D == 0){
	// very easy, D = 0 and phi is reducible
	is_irreducible = false;
      }else{
	mpz_sub(D, D, tmp); //   D <- b^2 - 4*a*c and is != 0
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
 * \brief test wether phi is a suitable candidate for computing g from phi
 * 
 * \param[in] p multiprecision integer, prime
 * \param[in] n integer greater than 1 ( n = 2 at the moment)
 * \param[in] phi polynomial with coefficients mod p (so very large)
 *
 * \return bool true if phi is of degree n and irreducible mod p
 * \return false otherwise
 */
//bool is_good_phi(mpz_poly_t phi, unsigned int n, mpz_t p){}
  //    return (Degree(phi_p) eq n) and IsIrreducible(phi_p);


/**
 * \brief select suitable polynomial f 
 * 
 * \param[in] p multiprecision integer, prime
 * \param[out] phi_list a list of possible phi
 * 
 * phi is a degree n factor of f, with coefficients in y with y a root of Py.
 * \return true if a suitable polynomial f was found
 * \return false otherwise 
 */

int get_f_CONJ( mpz_t p, fPyphi_poly_t * list_f, mpz_poly_t * list_phi){
  int i = 0, j=0;// index of first good f in table
  bool found_good_f = false;
  //long int *fi = NULL;
  long int *Pyi = NULL;
  long int **phii = NULL;
  mpz_t y;


  if(list_f != NULL){
    // Now, find an appropriate poly f in that table.
    // start with the first one, etc because the polynomials are sorted in 
    // decreasing order of interest (decreasing order of Murphy E value
    while ((found_good_f != true) && (i < list_f->size)){// nothing found
      // i.e. the i-th polynomial in f is not irreducible mod p,
      // or is but PY has no root.
      //fi = ((list_f->table_f[i]).f);
      Pyi = ((list_f->table_f[i]).Py);
      phii = ((list_f->table_f[i]).phi);

      if (is_irreducible_mod_p_si(Pyi, list_f->deg_Py, p)){
	// compute the two roots of Pyi, 
	// --> compute a square root mod p
	// HOW TO DO THAT ?
	// we obtain y a root of Py mod p.

	// compute the possible phi(s)
	j = 0;
	eval_si_phi_mpz_y(list_phi[j], phii, list_f->deg_phi, y);
	// test if at least one phi is irreducible
	// is_irreducible_mod_p(phii_j, p)
	// add in table each irreducible phii_j
	
	
      }
    }// either poly f found or end of list reached.
  }
  return i;
}

/**
 * \brief The main function that computes f and g and print them in the .poly file.
 * 
 * \param[in] p multiprecision integer, prime
 * \param[in] n integer greater than 1 ( n = 2 at the moment)
 * \param[in] label tring to label the output .poly file
 * 
 */

// , mpz_t ell, unsigned int mnfs
int gfpkdlpolyselect( unsigned int n, mpz_t p, char* label){

  int dd_p;
  unsigned int deg_f = 2*n, deg_g = 2*n;
  // take the largest possibility as default init for deg_f and deg_g.
  mpz_poly_t f, g;
  fPyphi_poly_t* list_f = NULL;
  mpz_poly_t table_phi[DEG_PY];

  if (n==2){
    get_degree_CONJ_f_g(n, &deg_f, &deg_g);
    polygen_CONJ_get_tab_f(deg_f, &list_f);
    // list_f pointe sur une structure qui contient tout ce qu'il faut : une table de {f, Py, phi} 
    // plus quelques autres parametres
    if (list_f != NULL){
      // is there a mpz_t function to compute the size in decimal digits of a number ?
      // look at Miscellanous Functions on Integers in GMP doc
      dd_p = mpz_sizeinbase(p, 10);
      mpz_poly_init(f, deg_f);
      // 1. find a good f in table. If no f is found, return failed and exit 0.
      int f_id = get_f_CONJ( p, list_f, table_phi);
      

      if ((f_id >= 0) && (f_id < (list_f->size))){
	// 2. find a good g with LLL.
	// ppf_t pg = {{2, {0,1}, 0, {0,0}, 2, -1, 0}}; 
	mpz_poly_init(g, deg_g);

	
	char* filename = (char*)malloc(strlen(label) + 16);
	
	sprintf(filename, "p%ddd%d%s.poly", n, dd_p, label);
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
	mpz_poly_clear(g);
      }else{
	return 0; // f not found.
      }
      mpz_poly_clear(f);
    }else{
      return 0; // NO list_f available for this n.
    }
    
  }else{
    return 0; // n not supported (n=2 only at the moment)
  }
  return 0;
}

