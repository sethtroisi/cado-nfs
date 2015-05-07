#include "utils_mpz.h"
#include "macros.h"

int main()
{
  mpz_t z;
  mpz_init(z);

  mpz_set_str(z, "1365", 10);
  factor_t factor;
  ASSERT_ALWAYS(gmp_brute_force_factorize(factor, z) == 0);
  ASSERT_ALWAYS(factor_assert(factor, z) == 0);
  ASSERT_ALWAYS(mpz_cmp_ui(factor->factorization[3], 13) == 0);
  ASSERT_ALWAYS(factor->number == 4);
  factor_clear(factor);

  mpz_set_str(z, "21039", 10);
  ASSERT_ALWAYS(gmp_brute_force_factorize(factor, z) == 0);
  ASSERT_ALWAYS(factor_assert(factor, z) == 0);
  ASSERT_ALWAYS(mpz_cmp_ui(factor->factorization[1], 7013) == 0);
  ASSERT_ALWAYS(factor->number == 2);
  factor_clear(factor);
  
  mpz_set_str(z, "54631", 10);
  ASSERT_ALWAYS(gmp_brute_force_factorize(factor, z) == 0);
  ASSERT_ALWAYS(factor_assert(factor, z) == 0);
  factor_clear(factor);
  
  mpz_set_str(z, "87919", 10);
  ASSERT_ALWAYS(gmp_brute_force_factorize(factor, z) == 0);
  ASSERT_ALWAYS(factor_assert(factor, z) == 0);
  factor_clear(factor);
  
  mpz_set_str(z, "81097", 10);
  ASSERT_ALWAYS(gmp_brute_force_factorize(factor, z) == 0);
  ASSERT_ALWAYS(factor_assert(factor, z) == 0);
  factor_clear(factor);
  
  mpz_set_str(z, "2361", 10);
  ASSERT_ALWAYS(gmp_brute_force_factorize(factor, z) == 0);
  ASSERT_ALWAYS(factor_assert(factor, z) == 0);
  factor_clear(factor);
  
  mpz_set_str(z, "30658", 10);
  ASSERT_ALWAYS(gmp_brute_force_factorize(factor, z) == 0);
  ASSERT_ALWAYS(factor_assert(factor, z) == 0);
  factor_clear(factor);
  
  mpz_set_str(z, "64759", 10);
  ASSERT_ALWAYS(gmp_brute_force_factorize(factor, z) == 0);
  ASSERT_ALWAYS(factor_assert(factor, z) == 0);
  factor_clear(factor);
  
  mpz_set_str(z, "64285", 10);
  ASSERT_ALWAYS(gmp_brute_force_factorize(factor, z) == 0);
  ASSERT_ALWAYS(factor_assert(factor, z) == 0);
  factor_clear(factor);
  
  mpz_set_str(z, "75099", 10);
  ASSERT_ALWAYS(gmp_brute_force_factorize(factor, z) == 0);
  ASSERT_ALWAYS(factor_assert(factor, z) == 0);
  factor_clear(factor);
  
  mpz_set_str(z, "44190", 10);
  ASSERT_ALWAYS(gmp_brute_force_factorize(factor, z) == 0);
  ASSERT_ALWAYS(factor_assert(factor, z) == 0);
  ASSERT_ALWAYS(factor->number == 5);
  ASSERT_ALWAYS(mpz_cmp(factor->factorization[1],
        factor->factorization[2]) == 0);
  factor_clear(factor);
  
  mpz_set_str(z, "44190", 10);
  ASSERT_ALWAYS(gmp_factorize(factor, z) == 0);
  ASSERT_ALWAYS(factor_assert(factor, z) == 0);
  ASSERT_ALWAYS(factor->number == 5);
  ASSERT_ALWAYS(mpz_cmp(factor->factorization[1],
        factor->factorization[2]) == 0);
  factor_clear(factor);

  mpz_clear(z);
  return 0;
}
