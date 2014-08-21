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

#define DEG_PY 2
#if DEG_PY > 2
#error "the code works only for Py of degree <= 2, sorry."
#endif
#define VARPHI_COEFF_INT 1

typedef long int coeff_t;

typedef struct {
  coeff_t t;     // parameter t
  coeff_t PY[DEG_PY + 1]; // no --> use one of the poly stuct of cado-nfs!
  coeff_t f[MAXDEGREE + 1];  // polynomial f of degree at most MAXDEGREE 
                         // set to 10 at the moment in cado_poly.h
}row_f_poly_t;

// tables containing polynomials f
// how to encode varphi form ?
typedef struct {
  unsigned int deg_f;
  unsigned int deg_Py;
  unsigned int deg_varphi;
  row_f_poly_t* table_f;
  coeff_t varphi[MAXDEGREE + 1][DEG_PY]; // poly whose coefficients are themselves poly in Y 
  //(a root of PY) modulo PY so of degree at most DEG_PY-1 --> of DEG_PY coeffs.
}table_f_poly_t;

// works only if PY is of degree 2
void eval_varphi_mpz(mpz_poly_t g, mpz_t** varphi_coeff, unsigned int deg_varphi, mpz_t u, mpz_t v);
void eval_varphi_si(mpz_poly_t g, coeff_t** varphi_coeff, unsigned int deg_varphi, mpz_t u, mpz_t v);
bool polygen_CONJ_get_tab_f(unsigned int deg_f, \
			    table_f_poly_t* table_f, \
			    unsigned int* table_f_size);
bool is_irreducible_ZZ(mpz_poly_t varphi);
bool is_irreducible_mod_p(mpz_poly_t varphi, mpz_t p);
bool is_good_varphi(mpz_poly_t varphi, unsigned int k, mpz_t p);
bool is_good_f_PY();
void polygen_CONJ_f ( mpz_t p, unsigned int k, mpz_poly_t f );
void polygen_CONJ_g ( mpz_t p, unsigned int k, mpz_poly_t f, mpz_poly_t g );
