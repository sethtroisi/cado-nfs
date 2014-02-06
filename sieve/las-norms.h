#ifndef LAS_NORMS_H_
#define LAS_NORMS_H_

#include <stdint.h>
#include "las-types.h"

#ifdef __cplusplus
extern "C" {
#endif

/*  initializing norms */
/* Knowing the norm on the rational side is bounded by 2^(2^k), compute
   lognorms approximations for k bits of exponent + NORM_BITS-k bits
   of mantissa */
void
init_norms (sieve_info_ptr si, int side);


/*  initialize norms for bucket regions */
/* Initialize lognorms on the rational side for the bucket_region
 * number N.
 * For the moment, nothing clever, wrt discarding (a,b) pairs that are
 * not coprime, except for the line j=0.
 */
void
init_rat_norms_bucket_region (unsigned char *S, unsigned int N, sieve_info_ptr si);

/* Initialize lognorms on the algebraic side for the bucket
 * number N.
 * Only the survivors of the rational sieve will be initialized, the
 * others are set to 255. Case GCD(i,j)!=1 also gets 255.
 * return nothing because the number of reports (= number of norm
 * initialisations) is algorithm dependent of ALG_RAT & ALG_LAZY.
 */
void
init_alg_norms_bucket_region (unsigned char *alg_S, 
                              unsigned char *rat_S,  unsigned int N, 
                              sieve_info_ptr si);

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
