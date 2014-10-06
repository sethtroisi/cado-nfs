#ifndef SM_UTILS_H_
#define SM_UTILS_H_

#define MAX_LEN_RELSET 1024

typedef struct {
  mpz_poly_t num[2];
  mpz_poly_t denom[2];
} sm_relset_struct_t;

typedef sm_relset_struct_t sm_relset_t[1];
typedef sm_relset_struct_t * sm_relset_ptr;
typedef const sm_relset_struct_t * sm_relset_srcptr;

void sm_relset_init (sm_relset_t r, int *d);
void sm_relset_clear (sm_relset_t r);

// (a,b) -> a - b*x
void mpz_poly_init_set_ab (mpz_poly_ptr, int64_t, uint64_t);

// array of rows and exponents -> rational fractions corresponding to the
// combination of these rows, stored in a sm_relset structure.
void sm_build_one_relset (sm_relset_ptr rel, uint64_t *r, int64_t *e, int len,
        mpz_poly_t * abpolys, mpz_poly_ptr *F, const mpz_t ell2);

// (a,b) -> SM as a polynomial.
void sm_single_rel (mpz_poly_ptr *SM, int64_t a, uint64_t b,
        mpz_poly_ptr * F, const mpz_ptr *eps,
        const mpz_t ell, const mpz_t ell2, const mpz_t invl2);

// Taking a polynomial modulo F as input, compute the corresponding SM
// as a polynomial.
void compute_sm (mpz_poly_t, mpz_poly_t, const mpz_poly_t, const mpz_t,
                 const mpz_t, const mpz_t, const mpz_t);

// Print coeffs of the SM polynomial
void print_sm (FILE *, mpz_poly_t, int, int);

#endif /* SM_UTILS_H_ */
