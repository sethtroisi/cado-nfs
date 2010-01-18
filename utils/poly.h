#include <gmp.h>

/* FIXME [ET]: Someday, this should me moved to utils/ ; other polynomial
 * libraries are there (admittedly different, since geared towards tiny
 * moduli), and we could perhaps fint it useful to adapt the interface
 * here for rootfinding code above ULONG_MAX
 */

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

// Management of the structure, printing
void poly_alloc(poly_t f, int d);
void poly_free(poly_t f);
void poly_print(const poly_t f);
void cleandeg(poly_t f, int deg);
void poly_setcoeff(poly_t f, int i, const mpz_t z);
void poly_setcoeff_str(poly_t f, int i,char *str, int base);
void poly_set(poly_t f, mpz_t * coeffs, int d);
void poly_getcoeff(mpz_t res, int i, const poly_t f); // Added for rootsieve
void poly_copy(poly_t g, const poly_t f);

// Polynomial operations

int poly_is_constant(const poly_t f); // Added for rootsieve
void poly_add(poly_t f, const poly_t g, const poly_t h);  // uncommented for rootsieve
void poly_sub(poly_t f, const poly_t g, const poly_t h);  // uncommented for rootsieve
void poly_sub_mod_mpz(poly_t f, const poly_t g, const poly_t h, const mpz_t m);
void poly_mul_ui(poly_t f, const poly_t g, unsigned long a); // uncommented for rootsieve
void poly_sub_ui(poly_t f, unsigned long a);
// void poly_div_ui_mod_mpz(poly_t f, const poly_t g, unsigned long a, const mpz_t m);
void poly_div_ui_mod_ui(poly_t f, const poly_t g, unsigned long a, const unsigned long m); // done - Added for rootsieve
void poly_div_ui(poly_t f, const poly_t g, unsigned long a);
void poly_div_2_mod_mpz(poly_t f, const poly_t g, const mpz_t m);
void poly_eval(mpz_t res, const poly_t f, const mpz_t x);
void poly_eval_mod_mpz(mpz_t res, const poly_t f, const mpz_t x, const mpz_t m);
void poly_mul(poly_t f, const poly_t g, const poly_t h); // uncommented for rootsieve

// void poly_reducemodF(polymodF_t P, poly_t p, const poly_t F);
void polymodF_mul(polymodF_t Q, const polymodF_t P1, const polymodF_t P2, const poly_t F);
void poly_mul_mpz(poly_t Q, const poly_t P, const mpz_t a);
void poly_reduce_mod_mpz (poly_t Q, const poly_t P, const mpz_t m);
void poly_reduce_makemonic_mod_mpz(poly_t Q, const poly_t P, const mpz_t m);
void poly_mul_mod_f_mod_mpz(poly_t Q, const poly_t P1, const poly_t P2,
                            const poly_t f, const mpz_t m, const mpz_t invm);
void poly_sqr_mod_f_mod_mpz(poly_t Q, const poly_t P, const poly_t f,
                            const mpz_t m, const mpz_t invm);
void poly_power_mod_f_mod_ui(poly_t Q, const poly_t P, const poly_t f,
        const mpz_t a, unsigned long p);

void poly_derivative(poly_t df, const poly_t f);

void barrett_init (mpz_t invm, const mpz_t m);
void barrett_mod (mpz_t a, mpz_t b, const mpz_t m, const mpz_t invm);
poly_t* poly_base_modp_init (const poly_t P0, int p, int *K, int l);
void poly_base_modp_clear (poly_t *P);
void poly_base_modp_lift (poly_t a, poly_t *P, int k, mpz_t pk);
size_t poly_sizeinbase (poly_t f, int d, int base);

