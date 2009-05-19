#include "gmp.h"

// Keep in mind the polynomial p(x) = \sum_{i=0}^{deg} a[i]*x^i
typedef struct {
	int deg;
	mpz_t * a;
} rspoly_struct;

typedef rspoly_struct rspoly[1];

// ************************
// BASIC METHODS (ON PLACE)
// ************************
void rspoly_init(rspoly p, int deg);
void rspoly_set(rspoly p, mpz_t * coeffs);
void rspoly_init_set(rspoly p, int deg, mpz_t * coeffs);
void rspoly_set_to_identity(rspoly p);
void rspoly_clear(rspoly p);
void rspoly_copy(rspoly target, rspoly source);
void rspoly_trim(rspoly p);
int  rspoly_is_constant(rspoly p);
void rspoly_constant_coeff(mpz_t res, rspoly p);
int  rspoly_equal(rspoly p, rspoly q);
int  rspoly_equivalent(rspoly p, rspoly q);

// Algebraic operations
void rspoly_diff(rspoly res, rspoly op);
void rspoly_mod_coeffs(rspoly res, rspoly op, unsigned long m);
void rspoly_div_coeffs(rspoly res, rspoly op, unsigned long m);
void rspoly_coeff_product_si(rspoly res, rspoly op, long s);
void rspoly_sum(rspoly res, rspoly op1, rspoly op2);
void rspoly_prod(rspoly res, rspoly op1, rspoly op2);

