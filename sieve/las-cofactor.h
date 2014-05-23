#ifndef LAS_COFACTOR_H_
#define LAS_COFACTOR_H_

#include "gmp.h"
#include "mpz_array.h"
#include "las-types.h"

int check_leftover_norm (const mpz_t n, sieve_info_srcptr si, int side);
int factor_both_leftover_norms(mpz_t *, const mpz_t, mpz_array_t **,
    uint32_array_t **, sieve_info_srcptr);
#endif
