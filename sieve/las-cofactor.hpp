#ifndef LAS_COFACTOR_HPP_
#define LAS_COFACTOR_HPP_

#include "gmp.h"
#include "mpz_array.h"
#include "las-types.hpp"

int check_leftover_norm (const mpz_t n, sieve_info const & si, int side);
int factor_both_leftover_norms(mpz_t *, mpz_array_t **, uint32_array_t **,
        sieve_info const &);
#endif
