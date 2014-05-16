#ifndef SM_UTILS_H_
#define SM_UTILS_H_

#define MAX_LEN_RELSET 1024

typedef struct {
  mpz_poly_t num;
  mpz_poly_t denom;
} sm_relset_struct_t;

typedef sm_relset_struct_t sm_relset_t[1];
typedef sm_relset_struct_t * sm_relset_ptr;
typedef const sm_relset_struct_t * sm_relset_srcptr;

void sm_relset_init (sm_relset_t r, int d);
void sm_relset_clear (sm_relset_t r);
void sm_build_one_relset (sm_relset_ptr, uint64_t *, int64_t *, int,
                          mpz_poly_t *, mpz_poly_t, const mpz_t);
void mpz_poly_init_set_ab (mpz_poly_ptr, int64_t, uint64_t);
void compute_sm (mpz_poly_t, mpz_poly_t, const mpz_poly_t, const mpz_t,
                 const mpz_t, const mpz_t, const mpz_t);
void print_sm (FILE *, mpz_poly_t, int, int);
void sm_single_rel (mpz_poly_t, int64_t, uint64_t, mpz_poly_t, const mpz_t,
                    const mpz_t, const mpz_t, const mpz_t);

#endif /* SM_UTILS_H_ */
