#ifndef FILTER_SM_H_
#define FILTER_SM_H_

typedef struct {
  poly_t num;
  poly_t denom;
} sm_relset_struct_t;

typedef sm_relset_struct_t sm_relset_t[1];
typedef sm_relset_struct_t * sm_relset_ptr;
typedef const sm_relset_struct_t * sm_relset_srcptr;


void poly_power_mod_f_mod_mpz_Barrett (poly_t, const poly_t, const poly_t,
                                       const mpz_t, const mpz_t, const mpz_t);
void poly_alloc_and_set_from_ab (poly_ptr, int64_t, uint64_t);
void poly_reduce_frac_mod_f_mod_mpz (sm_relset_ptr, const poly_t, const mpz_t,
                                     mpz_t, poly_t, poly_t, poly_t);
void compute_sm (poly_t, poly_t, const poly_t, const mpz_t, const mpz_t,
                 const mpz_t, const mpz_t);
void print_sm (FILE *, poly_t, int);
void sm_single_rel(poly_t, int64_t, uint64_t, poly_t, const mpz_t, const mpz_t,
                   const mpz_t, const mpz_t);
#endif /* FILTER_SM_H_ */
