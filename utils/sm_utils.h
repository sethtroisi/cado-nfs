#ifndef SM_UTILS_H_
#define SM_UTILS_H_

#define MAX_LEN_RELSET 1024

#include "mpz_poly.h"

struct sm_side_info_s {
    int nsm;
    mpz_t ell;
    mpz_t ell2;
    mpz_t invl2;        /* barrett precomputed inverse.
                           Not always used; see compute_sm_lowlevel */
    mpz_poly_srcptr f0;
    mpz_poly_t f;       /* monic */
    mpz_poly_factor_list fac;
    mpz_t exponent;
    mpz_t * exponents;
    mpz_t * matrix;
};
typedef struct sm_side_info_s sm_side_info[1];
typedef struct sm_side_info_s * sm_side_info_ptr;
typedef const struct sm_side_info_s * sm_side_info_srcptr;

void sm_side_info_init(sm_side_info_ptr sm, mpz_poly_srcptr f0, mpz_srcptr ell);
void sm_side_info_clear(sm_side_info_ptr sm);
void sm_side_info_print(FILE * out, sm_side_info_srcptr sm);

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
// If F[0] or F[1] is NULL, then no computation is done on the
// corresponding side.
void sm_build_one_relset (sm_relset_ptr rel, uint64_t *r, int64_t *e, int len,
        mpz_poly_t * abpolys, mpz_poly_ptr *F, const mpz_t ell2);

// Taking a polynomial modulo F as input, compute the corresponding SM
// as a polynomial.
void compute_sm_straightforward(mpz_poly_ptr dst, mpz_poly_srcptr u, sm_side_info_srcptr sm);

// Print coeffs of the SM polynomial
void print_sm (FILE *, mpz_poly_t, int, int);

/* This does the same as compute_sm_straightforward, except that it works piecewise on
 * the different components. It is thus noticeably faster. Results are
 * compatible, as the change of basis is precomputed within the
 * sm_side_info structure.
 */
void compute_sm_piecewise(mpz_poly_ptr dst, mpz_poly_srcptr u, sm_side_info_srcptr sm);


#endif /* SM_UTILS_H_ */
