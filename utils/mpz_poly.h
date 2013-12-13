
#ifndef MPZ_POLY_H_
#define MPZ_POLY_H_

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
} mpz_poly_struct_t;

typedef mpz_poly_struct_t mpz_poly_t[1];
typedef mpz_poly_struct_t * mpz_poly_ptr;
typedef const mpz_poly_struct_t * mpz_poly_srcptr;

/* -------------------------------------------------------------------------- */

/* [LI] This should also be renamed to mpz_polymodF.  */
/* And may be moved outside mpz_poly.[ch]? */
 
/* Let F(x) be a (non-monic) polynomial of degree d:
   F(x) = f_d x^d + f_{d-1} x^{d-1} + .... + f_1 x + f_0
   The following type represents a polynomial modulo F(x):
   If P is such an element, it means: P = P.p / f_d^P.v */
typedef struct {
  mpz_poly_t p;
  int v;
} polymodF_struct_t;

typedef polymodF_struct_t polymodF_t[1];
/* -------------------------------------------------------------------------- */

/* Special structure to represent coeffs of a polynomial in base p^k. */
typedef struct {
  int deg;
  char **coeff;
} poly_base_struct_t;

typedef poly_base_struct_t poly_base_t[1];


/* Management of the structure, set and print coefficients. */
void mpz_poly_init(mpz_poly_t, int d);
void mpz_poly_free(mpz_poly_t f);
void mpz_poly_print(const mpz_poly_t f);
void mpz_poly_cleandeg(mpz_poly_t f, int deg);
void mpz_poly_set(mpz_poly_t f, mpz_t * coeffs, int d);
void mpz_poly_set_zero(mpz_poly_t f);
void mpz_poly_setcoeff(mpz_poly_t f, int i, const mpz_t z);
void mpz_poly_setcoeff_si(mpz_poly_t f, int i, int z);
void mpz_poly_setcoeff_int64(mpz_poly_t f, int i, int64_t z);
void mpz_poly_setcoeff_uint64(mpz_poly_t f, int i, uint64_t z);
void mpz_poly_setcoeff_str(mpz_poly_t f, int i,char *str, int base);
void mpz_poly_get(mpz_poly_t f, mpz_t * coeffs, int d);
void mpz_poly_getcoeff(mpz_t res, int i, const mpz_poly_t f);
void mpz_poly_copy(mpz_poly_t g, const mpz_poly_t f);

/* Polynomial arithmetic */
void mpz_poly_add(mpz_poly_t f, const mpz_poly_t g, const mpz_poly_t h);
void mpz_poly_sub(mpz_poly_t f, const mpz_poly_t g, const mpz_poly_t h);
void mpz_poly_sub_ui(mpz_poly_t f, unsigned long a);
void mpz_poly_sub_mod_mpz(mpz_poly_t f, const mpz_poly_t g, const mpz_poly_t h,
                      const mpz_t m);
void mpz_poly_mul(mpz_poly_t f, const mpz_poly_t g, const mpz_poly_t h);
void mpz_poly_mul_ui(mpz_poly_t f, const mpz_poly_t g, unsigned long a);
void mpz_poly_mul_mpz(mpz_poly_t Q, const mpz_poly_t P, const mpz_t a);
int mpz_poly_mod_f_mod_mpz (mpz_t *R, int d, mpz_t *f, int df, const mpz_t m,
                        const mpz_t invm);
void mpz_poly_mul_mod_f_mod_mpz(mpz_poly_t Q, const mpz_poly_t P1, const mpz_poly_t P2,
                            const mpz_poly_t f, const mpz_t m,
                            const mpz_t invm);
void mpz_poly_div_r (mpz_poly_t h, const mpz_poly_t f, const mpz_t p);
void mpz_poly_div_qr (mpz_poly_t q, mpz_poly_t r, const mpz_poly_t f, const mpz_poly_t g, const mpz_t p);
void mpz_poly_div_ui(mpz_poly_t f, const mpz_poly_t g, unsigned long a);
void mpz_poly_div_ui_mod_ui(mpz_poly_t f, const mpz_poly_t g, unsigned long a,
                        const unsigned long m);
void mpz_poly_div_2_mod_mpz(mpz_poly_t f, const mpz_poly_t g, const mpz_t m);
void mpz_poly_divexact (mpz_poly_t q, mpz_poly_t h, const mpz_poly_t f, const mpz_t p);
  
void mpz_poly_eval(mpz_t res, const mpz_poly_t f, const mpz_t x);
void mpz_poly_eval_mod_mpz(mpz_t res, const mpz_poly_t f, const mpz_t x,
                       const mpz_t m);
void polymodF_mul(polymodF_t Q, const polymodF_t P1, const polymodF_t P2,
                  const mpz_poly_t F);
void mpz_poly_reduce_mod_mpz(mpz_poly_t Q, const mpz_poly_t P, const mpz_t m);
void mpz_poly_reduce_makemonic_mod_mpz(mpz_poly_t Q, const mpz_poly_t P,
                                   const mpz_t m);
void mpz_poly_sqr_mod_f_mod_mpz(mpz_poly_t Q, const mpz_poly_t P, const mpz_poly_t f,
                            const mpz_t m, const mpz_t invm);
void mpz_poly_power_mod_f_mod_ui(mpz_poly_t Q, const mpz_poly_t P, const mpz_poly_t f,
                             const mpz_t a, unsigned long p);
void mpz_poly_power_mod_f_mod_mpz (mpz_poly_t Q, const mpz_poly_t P, const mpz_poly_t f,
                               const mpz_t a, const mpz_t p);
void mpz_poly_derivative(mpz_poly_t df, const mpz_poly_t f);
void barrett_init (mpz_ptr invm, mpz_srcptr m);
void barrett_mod (mpz_ptr a, mpz_srcptr b, mpz_srcptr m,
                  mpz_srcptr invm);
mpz_poly_t* mpz_poly_base_modp_init (const mpz_poly_t P0, int p, int *K, int l);
void mpz_poly_base_modp_clear (mpz_poly_t *P);
void mpz_poly_base_modp_lift (mpz_poly_t a, mpz_poly_t *P, int k, mpz_t pk);
int mpz_poly_is_constant(const mpz_poly_t f); // for rootsieve
size_t mpz_poly_sizeinbase (mpz_poly_t f, int d, int base);
void mpz_poly_swap (mpz_poly_t f, mpz_poly_t g);
void mpz_poly_gcd_mpz (mpz_poly_t f, mpz_poly_t g, const mpz_t p);
  void mpz_poly_xgcd_mpz(mpz_poly_t gcd, const mpz_poly_t f, const mpz_poly_t g, mpz_poly_t u, mpz_poly_t v, const mpz_t p);
int mpz_poly_cantor_zassenhaus (mpz_t *r, mpz_poly_t f, const mpz_t p, int depth);
int mpz_poly_roots_mpz (mpz_t *r, mpz_t *f, int d, const mpz_t p);
void mpz_poly_swap (mpz_poly_t f, mpz_poly_t g);


#ifdef __cplusplus
}
#endif

#endif	/* MPZ_POLY_H_ */
