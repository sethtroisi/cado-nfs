#ifndef LAS_NORMS_H_
#define LAS_NORMS_H_

#include <stdint.h>
#include "las-types.h"
#include "double_poly.h"

#ifdef __cplusplus
extern "C" {
#endif

/* initializing norms */
/* Knowing the norm on the rational side is bounded by 2^(2^k), compute
   lognorms approximations for k bits of exponent + NORM_BITS-k bits
   of mantissa */
void init_norms (sieve_info_ptr si, int side);

/* Initialize lognorms for the bucket_region number J. It's a wrapper.
 * For the moment, nothing clever, wrt discarding (a,b) pairs that are
 * not coprime.
 * The sieve area S must be preallocated with at least (BUCKET_REGION +
 * MEMSET_MIN) space. Indeed, the algorithm that is used might write a
 * bit beyond the meaningful area.
 */
void init_norms_bucket_region (unsigned char *S, uint32_t J, sieve_info_ptr si, unsigned int side, unsigned int smart);

/* These functions are internals. Don't use them. Use the wrapper above.
   It's need to declare them here for units & coverage tests.
 */
void init_degree_one_norms_bucket_region_internal     (unsigned char *S, uint32_t J, uint32_t I, double scale, double u0, double u1, double *cexp2);
void init_exact_degree_X_norms_bucket_region_internal (unsigned char *S, uint32_t J, uint32_t I, double scale, unsigned int d, double *fijd);
void init_smart_degree_X_norms_bucket_region_internal (unsigned char *S, uint32_t J, uint32_t I, double scale, unsigned int d, double *fijd, unsigned int nroots, double *roots);
void init_norms_roots_internal (unsigned int degree, double *coeff, double max_abs_root, unsigned int *nroots, double *roots);

double get_maxnorm_alg (double_poly_srcptr src_poly, const double X, const double Y);

/* This prepares the auxiliary data which is used by
 * init_rat_norms_bucket_region and init_alg_norms_bucket_region
 */
void sieve_info_init_norm_data(sieve_info_ptr si);

void sieve_info_clear_norm_data(sieve_info_ptr si);

void sieve_info_update_norm_data (FILE *, sieve_info_ptr, int);

int sieve_info_adjust_IJ(sieve_info_ptr si, int nb_threads);
void sieve_info_init_norm_data_sq (sieve_info_ptr si, unsigned long q);

#ifdef __cplusplus
}
#endif

#endif	/* LAS_NORMS_H_ */
