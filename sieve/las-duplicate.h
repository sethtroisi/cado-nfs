#ifndef LAS_DUPLICATE_H_
#define LAS_DUPLICATE_H_

#include "las-types.h"
#include "relation.h"
#include "mpz_array.h"

#ifdef __cplusplus
extern "C" {
#endif

sieve_info_ptr
fill_in_sieve_info(const mpz_t q, const mpz_t rho,
                   const int sq_side, const uint32_t I, const uint32_t J,
                   const unsigned long limits[2], facul_strategy_t *strategy[2],
                   cado_poly_ptr cpoly, siever_config_srcptr conf);
int relation_is_duplicate(FILE *, relation_t *, int, sieve_info_srcptr);

/* FIXME: These function are defined in las.c. The prototypes should not be 
   here. The functions should be moved to an appropriate file 
   (las-qlattice.c?) and the prototypes to the corresponding header file. */
int sieve_info_adjust_IJ(sieve_info_ptr, int);
int check_leftover_norm (const mpz_t n, sieve_info_srcptr si, int side);
int factor_both_leftover_norms(mpz_t *, const mpz_t, mpz_array_t **, 
                               uint32_array_t **, sieve_info_srcptr);

#ifdef __cplusplus
}
#endif

#endif
