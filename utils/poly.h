/**
   FIXME [ET]: Someday, this should me moved to utils/ ; other polynomial
   libraries are there (admittedly different, since geared towards tiny
   moduli), and we could perhaps fint it useful to adapt the interface
   here for rootfinding code above ULONG_MAX
*/

#ifndef POLY_H_
#define POLY_H_

#include <gmp.h>
#ifdef __cplusplus
extern "C" {
#endif


/* Note, deg = -1 means P=0; otherwise, one should have coeff[deg] != 0.
   Warning: a polynomial of degree d needs d+1 allocation. */
typedef struct {
  int alloc;
  int deg;
  mpz_t *coeff;
} poly_struct_t;

typedef poly_struct_t poly_t[1];
typedef poly_struct_t * poly_ptr;
typedef const poly_struct_t * poly_srcptr;

/* Let F(x) be a (non-monic) polynomial of degree d:
   F(x) = f_d x^d + f_{d-1} x^{d-1} + .... + f_1 x + f_0
   The following type represents a polynomial modulo F(x):
   If P is such an element, it means: P = P.p / f_d^P.v */
typedef struct {
  poly_t p;
  int v;
} polymodF_struct_t;

typedef polymodF_struct_t polymodF_t[1];

/* Special structure to represent coeffs of a polynomial in base p^k. */
typedef struct {
  int deg;
  char **coeff;
} poly_base_struct_t;

typedef poly_base_struct_t poly_base_t[1];


/* Management of the structure, set and print coefficients. */
void poly_alloc(poly_t f, int d);
void poly_free(poly_t f);
void poly_print(const poly_t f);
void cleandeg(poly_t f, int deg);
void poly_set(poly_t f, mpz_t * coeffs, int d);
void poly_set_zero(poly_t f);
void poly_setcoeff(poly_t f, int i, const mpz_t z);
void poly_setcoeff_si(poly_t f, int i, int z);
void poly_setcoeff_str(poly_t f, int i,char *str, int base);
void poly_get(poly_t f, mpz_t * coeffs, int d);
void poly_getcoeff(mpz_t res, int i, const poly_t f);
void poly_copy(poly_t g, const poly_t f);

/* Polynomial arithmetic */
void poly_add(poly_t f, const poly_t g, const poly_t h);
void poly_sub(poly_t f, const poly_t g, const poly_t h);
void poly_sub_ui(poly_t f, unsigned long a);
void poly_sub_mod_mpz(poly_t f, const poly_t g, const poly_t h,
                      const mpz_t m);
void poly_mul(poly_t f, const poly_t g, const poly_t h);
void poly_mul_ui(poly_t f, const poly_t g, unsigned long a);
void poly_mul_mpz(poly_t Q, const poly_t P, const mpz_t a);
void poly_mul_mod_f_mod_mpz(poly_t Q, const poly_t P1, const poly_t P2,
                            const poly_t f, const mpz_t m,
                            const mpz_t invm);
void poly_div_r (poly_t h, const poly_t f, const mpz_t p);
void poly_div_ui(poly_t f, const poly_t g, unsigned long a);
void poly_div_ui_mod_ui(poly_t f, const poly_t g, unsigned long a,
                        const unsigned long m);
void poly_div_2_mod_mpz(poly_t f, const poly_t g, const mpz_t m);
void poly_divexact (poly_t q, poly_t h, const poly_t f, const mpz_t p);
  
void poly_eval(mpz_t res, const poly_t f, const mpz_t x);
void poly_eval_mod_mpz(mpz_t res, const poly_t f, const mpz_t x,
                       const mpz_t m);
void polymodF_mul(polymodF_t Q, const polymodF_t P1, const polymodF_t P2,
                  const poly_t F);
void poly_reduce_mod_mpz(poly_t Q, const poly_t P, const mpz_t m);
void poly_reduce_makemonic_mod_mpz(poly_t Q, const poly_t P,
                                   const mpz_t m);
void poly_sqr_mod_f_mod_mpz(poly_t Q, const poly_t P, const poly_t f,
                            const mpz_t m, const mpz_t invm);
void poly_power_mod_f_mod_ui(poly_t Q, const poly_t P, const poly_t f,
                             const mpz_t a, unsigned long p);
void poly_power_mod_f_mod_mpz (poly_t Q, const poly_t P, const poly_t f,
                               const mpz_t a, const mpz_t p);
void poly_derivative(poly_t df, const poly_t f);
void barrett_init (mpz_ptr invm, mpz_srcptr m);
void barrett_mod (mpz_ptr a, mpz_srcptr b, mpz_srcptr m,
                  mpz_srcptr invm);
poly_t* poly_base_modp_init (const poly_t P0, int p, int *K, int l);
void poly_base_modp_clear (poly_t *P);
void poly_base_modp_lift (poly_t a, poly_t *P, int k, mpz_t pk);
int poly_is_constant(const poly_t f); // for rootsieve
size_t poly_sizeinbase (poly_t f, int d, int base);
void poly_swap (poly_t f, poly_t g);
void poly_gcd_mpz (poly_t f, poly_t g, const mpz_t p);
  
int poly_cantor_zassenhaus (mpz_t *r, poly_t f, const mpz_t p, int depth);
int poly_roots_mpz (mpz_t *r, mpz_t *f, int d, const mpz_t p);
void poly_swap (poly_t f, poly_t g);


#ifdef __cplusplus
}
#endif

#endif	/* POLY_H_ */
