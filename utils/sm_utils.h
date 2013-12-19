#ifndef FILTER_SM_H_
#define FILTER_SM_H_

typedef struct {
  mpz_poly_t num;
  mpz_poly_t denom;
} sm_relset_struct_t;

typedef sm_relset_struct_t sm_relset_t[1];
typedef sm_relset_struct_t * sm_relset_ptr;
typedef const sm_relset_struct_t * sm_relset_srcptr;


void mpz_poly_power_mod_f_mod_mpz_Barrett (mpz_poly_t, const mpz_poly_t, const mpz_poly_t,
                                       const mpz_t, const mpz_t, const mpz_t);
void mpz_poly_init_set_ab (mpz_poly_ptr, int64_t, uint64_t);
void mpz_poly_reduce_frac_mod_f_mod_mpz (sm_relset_ptr, const mpz_poly_t, const mpz_t,
                                     mpz_t, mpz_poly_t, mpz_poly_t, mpz_poly_t);
void compute_sm (mpz_poly_t, mpz_poly_t, const mpz_poly_t, const mpz_t, const mpz_t,
                 const mpz_t, const mpz_t);
void print_sm (FILE *, mpz_poly_t, int);
void sm_single_rel(mpz_poly_t, int64_t, uint64_t, mpz_poly_t, const mpz_t, const mpz_t,
                   const mpz_t, const mpz_t);

#endif /* FILTER_SM_H_ */
