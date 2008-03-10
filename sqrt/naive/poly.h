#include <gmp.h>

// deg = -1 means P=0
// otherwise, one should have coeff[deg] != 0
// Warning: a polynomial of degree d needs d+1 allocation.
typedef struct {
  int alloc;
  int deg;
  mpz_t *coeff;
} poly_struct_t;

typedef poly_struct_t poly_t[1];


// Let F(x) be a (non-monic) polynomial of degree d:
//   F(x) = f_d x^d + f_{d-1} x^{d-1} + .... + f_1 x + f_0
// The following type represents a polynomial modulo F(x):
// If P is such an element, it means:
//   P = P.p / f_d^P.v
typedef struct {
  poly_t p;
  int v;
} polymodF_struct_t;

typedef polymodF_struct_t polymodF_t[1];

/* special structure to represent the coefficients of a polynomial in base p^k
 */
typedef struct {
  int deg;
  char **coeff;
} poly_base_struct_t;

typedef poly_base_struct_t poly_base_t[1];

void poly_alloc(poly_t f, int d);
void poly_free(poly_t f);
void poly_print(const poly_t f);
void cleandeg(poly_t f, int deg);
void poly_setcoeff(poly_t f, int i, const mpz_t z);
void poly_copy(poly_t g, const poly_t f);
void poly_add(poly_t f, const poly_t g, const poly_t h);
void poly_sub(poly_t f, const poly_t g, const poly_t h);
void poly_sub_mod_mpz(poly_t f, const poly_t g, const poly_t h, const mpz_t m);
void poly_mul_ui(poly_t f, const poly_t g, unsigned long a);
// void poly_div_ui_mod_mpz(poly_t f, const poly_t g, unsigned long a, const mpz_t m);
void poly_div_2_mod_mpz(poly_t f, const poly_t g, const mpz_t m);
void poly_eval_mod_mpz(mpz_t res, const poly_t f, const mpz_t x, const mpz_t m);
void poly_mul(poly_t f, const poly_t g, const poly_t h);
void poly_reducemodF(polymodF_t P, poly_t p, const poly_t F);
void polymodF_mul(polymodF_t Q, const polymodF_t P1, const polymodF_t P2, const poly_t F);
void poly_mul_mpz(poly_t Q, const poly_t P, const mpz_t a);
void poly_reduce_mod_mpz (poly_t Q, const poly_t P, const mpz_t m);
void poly_reduce_mod_mpz_fast (poly_t Q, const mpz_t m, poly_base_t S, int p, int k);
void poly_reduce_makemonic_mod_mpz(poly_t Q, const poly_t P, const mpz_t m);
void poly_mul_mod_f_mod_mpz(poly_t Q, const poly_t P1, const poly_t P2,
                            const poly_t F, const mpz_t m, const mpz_t invm);
void poly_sqr_mod_f_mod_mpz(poly_t Q, const poly_t P, const poly_t F,
                            const mpz_t m, const mpz_t invm);
void poly_power_mod_f_mod_ui(poly_t Q, const poly_t P, const poly_t F,
        const mpz_t a, unsigned long p);
void barrett_init (mpz_t invm, const mpz_t m);
void poly_base_init_set (poly_base_t Q, poly_t P, int p);
void poly_base_clear (poly_base_t Q);
