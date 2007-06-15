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

void poly_alloc(poly_t f, int d);
void poly_free(poly_t f);
void cleandeg(poly_t f, int deg);
void poly_setcoeff(poly_t f, int i, mpz_t z);
void poly_copy(poly_t g, poly_t f);
void poly_add(poly_t f, poly_t g, poly_t h);
void poly_sub(poly_t f, poly_t g, poly_t h);
void poly_mul(poly_t f, poly_t g, poly_t h);
void poly_reducemodF(polymodF_t P, poly_t p, poly_t F);
void polymodF_mul(polymodF_t Q, polymodF_t P1, polymodF_t P2, poly_t F);
