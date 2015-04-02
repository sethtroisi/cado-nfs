/**
 * \file cado-nfs/polyselect/gfpkdlpolyselect.h
 * 
 * \date 21/08/2014
 * \author Aurore Guillevic
 * \email guillevic@lix.polytechnique.fr
 * \brief header file of gfpkdlpolyselect.c
 *        Compute two polynomials f, g suitable for discrete logarithm in 
 *        extension fields, with the conjugation method.
 *
 * \test TODO
 */

#ifndef GFPKDLPOLYSELECT_H
#define GFPKDLPOLYSELECT_H

#define ERROR_UNSPECIFIED 0
#define POLY_OK 1
#define ERROR_DEGREE 2
#define ERROR_NOT_IRREDUCIBLE 3
#define ERROR_NO_DEGREE_FACTOR 4
#define ERROR_SIGNATURE 5
#define ERROR_SUBFIELD 6
#define ERROR_SMEXP 7
#define ERROR_EASYBADIDEALS 8
#define ERROR_VERYBADIDEALS 9
#define ERROR_MAX 10


#define DEG_PY 2
#if DEG_PY > 2
#error "the code works only for Py of degree <= 2, sorry."
#else
//#define PHI_COEFF_INT 1 // ?

// table structure, old version.

typedef struct {
  long int t;     // parameter t
  long int PY[DEG_PY + 1]; // no --> use one of the poly stuct of cado-nfs!
  long int f[MAXDEGREE + 1];  // polynomial f of degree at most MAXDEGREE 
                         // set to 10 at the moment in utils/cado_poly.h
}row_ftPy_poly_t;

// tables containing polynomials f
// how to encode phi form ?
typedef struct {
  unsigned int deg_f;
  unsigned int deg_Py;
  unsigned int deg_phi;
  int size;
  const row_ftPy_poly_t* table_f; // la table de {f, Py, phi}
  long int phi[MAXDEGREE + 1][DEG_PY]; // poly whose coefficients are themselves poly in Y 
  //(a root of PY) modulo PY so of degree at most DEG_PY-1 --> of DEG_PY coeffs.
}ftPy_poly_t;


// table structure for CONJ with Fp2 (and JLSV1 Fp4 ?).

typedef struct {
  long int f[MAXDEGREE + 1];  // polynomial f of degree at most MAXDEGREE 
                              // set to 10 at the moment in utils/cado_poly.h
  long int Py[DEG_PY + 1];    // no --> use one of the poly stuct of cado-nfs!
  long int phi[MAXDEGREE + 1][DEG_PY]; // poly whose coefficients are themselves poly in Y 
}row_fPyphi_poly_t;

typedef struct {
  unsigned int deg_f;
  unsigned int deg_Py;
  unsigned int deg_phi;
  int size;
  const row_fPyphi_poly_t* table_f;
}fPyphi_poly_t;

// for keeping parameters of each poly (each side) along the process of polyselect

typedef struct {
  unsigned int deg;
  unsigned int sgtr[2];
  unsigned int deg_subfield;
  unsigned int sgtr_subfield[2];
  unsigned int smexp; // from 1 to deg
  int nb_max_easybadideals;
  int nb_max_verybadideals;
}polyselect_parameter_poly_t;

typedef polyselect_parameter_poly_t ppf_t[1];

typedef struct {
  unsigned int n; // also written k in earlier versions
  mpz_t p;
  mpz_t ell;
  unsigned int mnfs;
}polyselect_parameters_t;

typedef polyselect_parameters_t pp_t[1];

void get_degree_CONJ_f_g(unsigned int n, unsigned int *deg_f, unsigned int *deg_g);

// works only if PY is of degree 2
void eval_mpz_phi_mpz_uv(mpz_poly_t g, mpz_t** phi_coeff, unsigned int deg_phi, mpz_t u, mpz_t v);
void eval_si_phi_mpz_y(mpz_poly_t g, long int** phi_coeff, unsigned int deg_phi, mpz_t y);
void eval_si_phi_mpz_uv(mpz_poly_t g, long int** phi_coeff, unsigned int deg_phi, mpz_t u, mpz_t v);
// return the table of suitable polynomials f according to deg_f.
// apparently, for GF(p^2) with CONJ and GF(p^4) with JLSV1, the same table is used.
bool polygen_CONJ_get_tab_f(unsigned int deg_f, \
			    fPyphi_poly_t** table_f);
bool is_irreducible_ZZ(mpz_poly_t phi);
bool is_irreducible_mod_p(mpz_poly_t phi, mpz_t p);
bool is_irreducible_mod_p_si(long int* Py, int deg_Py, mpz_t p);

// same as is_good_poly_check_all in gfpk/magma/polyselect_utils.mag
// return: error_code, see above
int is_good_poly(pp_t pp, ppf_t ppf, mpz_poly_t f);
// I need a function that works also for polynomials of small coefficients (signed long int)
bool is_good_phi(mpz_poly_t phi, unsigned int n, mpz_t p);
bool is_good_f_PY(row_fPyphi_poly_t* row_f_Py_phi, mpz_poly_t** phi);

// same as ***<to be completed>*** function in gfpk/magma/polyselect_utils.mag
unsigned int get_index_next_poly_in_tab_f(fPyphi_poly_t * list_f);

// renvoie i l'indice de la bonne ligne du tableau de {F, Py, phi}
int get_f_CONJ( mpz_t p, fPyphi_poly_t * list_f, mpz_poly_t * list_phi);
// list_f->list_f[i] is the line with a right f.

// the function to call for generating a .poly file.
// works only for n=2 at the moment.
// , mpz_t ell, unsigned int mnfs
int gfpkdlpolyselect( unsigned int n, mpz_t p, char* label);
#endif // DEG_PY > 2 is not supported
#endif // GFPKDLPOLYSELECT_H
