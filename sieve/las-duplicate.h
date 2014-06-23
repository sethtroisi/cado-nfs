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
void clear_sieve_info(sieve_info_ptr);
int relation_is_duplicate(FILE *, relation_t *, int, sieve_info_srcptr);

#ifdef __cplusplus
}
#endif

#endif
