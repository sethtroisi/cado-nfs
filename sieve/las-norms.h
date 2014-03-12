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

/*  initialize norms for bucket regions */
/* Initialize lognorms on the rational side for the bucket_region
 * number J.
 * For the moment, nothing clever, wrt discarding (a,b) pairs that are
 * not coprime.
 * The sieve area S must be preallocated with at least (BUCKET_REGION +
 * MEMSET_MIN) space. Indeed, the algorithm that is used might write a
 * bit beyond the meaningful area.
 */
void init_rat_norms_bucket_region (unsigned char *S, uint32_t J, sieve_info_ptr si);

/* Initialize lognorms on the algebraic side for the bucket
 * number J.
 * Case GCD(i,j)!=1 gets 255.
 * Since we can have two algebraic sides, we need to say explicitly which
 * side we want to be initialized.
 */
void init_alg_norms_bucket_region (unsigned char *alg_S, uint32_t J, sieve_info_ptr si, int side);

/* These functions are internals. Don't use them. It's need to declare them
 * here for units & coverage tests.
 */
void init_alg_norms_bucket_region_internal (unsigned char *S, uint32_t J, uint32_t I, unsigned int d, double scale, double *fijd);

void init_rat_norms_bucket_region_internal (unsigned char *S, uint32_t J, uint32_t I, double scale, double u0, double u1, double *cexp2);

double get_maxnorm_alg (double_poly_srcptr src_poly, const double X, const double Y);

/* This prepares the auxiliary data which is used by
 * init_rat_norms_bucket_region and init_alg_norms_bucket_region
 */
void sieve_info_init_norm_data(sieve_info_ptr si);

void sieve_info_clear_norm_data(sieve_info_ptr si);

void sieve_info_update_norm_data (FILE *, sieve_info_ptr, int);

void sieve_info_init_norm_data_sq (sieve_info_ptr si, unsigned long q);

#ifdef __cplusplus
}
#endif

#endif	/* LAS_NORMS_H_ */
