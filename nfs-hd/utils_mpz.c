#include "utils_mpz.h"
#include <stdint.h>
#include "gmp-ecm.h"
#include <stdlib.h>
#include "macros.h"
#include "getprime.h"

void factor_init(factor_ptr factor, unsigned int alloc)
{
  ASSERT(alloc > 0);

  factor->alloc = alloc;
  factor->number = 0;
  factor->factorization = (mpz_t * ) malloc(sizeof(mpz_t) * alloc);
  for (unsigned int i = 0; i < alloc; i++) {
    mpz_init(factor->factorization[i]);
  }
}

void factor_clear(factor_ptr factor)
{
  for (unsigned int i = 0; i < factor->alloc; i++) {
    mpz_clear(factor->factorization[i]);
  }
  free(factor->factorization);
  factor->alloc = 0;
  factor->number = 0;
}

void factor_append(factor_ptr factor, mpz_srcptr z)
{
  if (factor->alloc == factor->number) {
    unsigned int newalloc = factor->alloc + 10;
    factor->factorization = (mpz_t * ) realloc(factor->factorization,
        sizeof(mpz_t) * newalloc);
    for (unsigned int i = factor->alloc; i < newalloc; ++i)
      mpz_init(factor->factorization[i]);
    factor->alloc = newalloc;
  }

  mpz_set(factor->factorization[factor->number], z);
  factor->number++;
}

void factor_fprintf(FILE * file, factor_srcptr factor)
{
  fprintf(file, "[");
  for (unsigned int i = 0; i < factor->number - 1; i++) {
    gmp_fprintf(file, "%Zd, ", factor->factorization[i]);
  }
  gmp_fprintf(file, "%Zd]\n", factor->factorization[factor->number - 1]);
}

unsigned int factor_is_smooth(factor_srcptr factor, mpz_t B, unsigned int sort)
{
  if (sort) {
    ASSERT(sort == 1);

    if (mpz_cmp(factor->factorization[factor->number - 1], B) > 0) {
      return 0;
    }
  } else {
    ASSERT(sort == 0);

    for (unsigned int i = 0; i < factor->number; i++) {
      if (mpz_cmp(factor->factorization[i], B) > 0) {
        return 0;
      }
    }
  }

  return 1;
}

unsigned int factor_assert(factor_srcptr factor, mpz_srcptr z)
{
  mpz_t tmp;
  mpz_init(tmp);
  mpz_set_ui(tmp, 1);
  for (unsigned int i = 0; i < factor->number; i++) {
    mpz_mul(tmp, tmp, factor->factorization[i]);
  }
  unsigned int assert_facto = 0;
  if (!mpz_cmp(tmp, z)) {
    ASSERT(mpz_cmp(tmp, z) == 0);

    assert_facto = 1;
  }
  mpz_clear(tmp);

  return assert_facto;
}

#if 0
BROKEN
static unsigned int factorize(factor_ptr factor, mpz_srcptr z_root,
    unsigned int * number)
{
  mpz_t z;
  mpz_init(z);
  mpz_set(z, z_root);
  int ret = 0;
  mpz_t res;
  mpz_init(res);
  for (unsigned int i = 0; i < 10; i++) {
    ret = ecm_factor(res, z, 200, NULL);
    if (ret) {
      if (mpz_cmp_ui(res, 0) != 0 && mpz_cmp_ui(res, 1) != 0) {
        if (mpz_perfect_power_p(res)){
          unsigned int j = 2;
          mpz_t root;
          mpz_t rem;
          mpz_init(root);
          mpz_init(rem);
          do {
            mpz_rootrem(root, rem, res, j);
            j++;
          } while (mpz_cmp_ui(rem, 0) != 0);
          if (mpz_probab_prime_p (root, 25)) {
            for (unsigned int k = 1; k < j; k++) {
              mpz_set(factor->factorization[* number], root);
              * number = * number + 1;
            }
          } else {
            for (unsigned int k = 1; k < j; k++) {
              factorize(factor, root, number);
            }
          }
          mpz_clear(root);
          mpz_clear(rem);
        } else if (mpz_probab_prime_p (res, 25)) {
          mpz_set(factor->factorization[* number], res);
          * number = * number + 1;
        } else {
          factorize(factor, res, number);
        }
        mpz_divexact(res, z, res);
        if (mpz_probab_prime_p (res, 25)) {
          mpz_set(factor->factorization[* number], res);
          * number = * number + 1;
        } else if (mpz_cmp_ui(res, 1) != 0) {
          factorize(factor, res, number);
        }
        break;
      }
    }
  }
  
  unsigned int factorise = 0;
  if (mpz_cmp_ui(z, 1) == 0) {
    factorise = 1;
  }
  mpz_clear(res);
  mpz_clear(z);

  return factorise;
}

static unsigned int brute_force_factorize(factor_ptr factor, mpz_ptr z,
    unsigned int * number, mpz_srcptr z_root, mpz_srcptr bound)
{
  mpz_set(z, z_root);
  mpz_t prime;
  mpz_init(prime);
  mpz_set_ui(prime, 1);
  while(mpz_cmp_ui(z, 1) != 0 && mpz_cmp(prime, bound) <= 0) {
    mpz_nextprime(prime, prime);
    mpz_t q;
    mpz_init(q);
    mpz_t r;
    mpz_init(r);
    mpz_fdiv_qr(q, r, z, prime);
    while (mpz_cmp_ui(r, 0) == 0) {
      mpz_set(z, q);
      mpz_set(factor->factorization[* number], prime);
      * number = * number + 1;
      mpz_fdiv_qr(q, r, z, prime);
    }
    mpz_clear(q);
    mpz_clear(r);
  }

  unsigned int factorise = 1;
  if (mpz_cmp_ui(z, 1) == 0) {
    factorise = 0;
  }

  mpz_clear(prime);

  return factorise;
}
#endif

int brute_force_factorize_ul(factor_ptr factor, mpz_ptr z,
    mpz_srcptr z_root, unsigned long bound)
{
  factor_init(factor, 10);
  prime_info pi;
  prime_info_init (pi);

  mpz_set(z, z_root);
  mpz_t prime_Z;
  mpz_init(prime_Z);
  unsigned long prime = 2;

  for (prime = 2; prime <= bound; prime = getprime_mt(pi)) {
    if (mpz_cmp_ui(z, 1) != 0) {
      mpz_set_ui(prime_Z, prime);
      mpz_t q;
      mpz_init(q);
      mpz_t r;
      mpz_init(r);
      mpz_fdiv_qr(q, r, z, prime_Z);
      while (mpz_cmp_ui(r, 0) == 0) {
        mpz_set(z, q);
        factor_append(factor, prime_Z);
        mpz_fdiv_qr(q, r, z, prime_Z);
      }
      mpz_clear(q);
      mpz_clear(r);
    }
  }

  mpz_clear(prime_Z);
  prime_info_clear (pi);

  if (mpz_cmp_ui(z, 1) == 0) {
    return 1;
  }
  return 0;
}

static int compare_factor(const void * p0, const void * p1)
{
  return(mpz_cmp((mpz_srcptr) p0, (mpz_srcptr) p1));
}

void sort_factor(factor_ptr factor)
{
  qsort(factor->factorization, factor->number,
      sizeof(factor->factorization[0]), compare_factor);
}

#if 0
unsigned int gmp_brute_force_factorize(factor_ptr factor, mpz_srcptr z)
{
  factor_init(factor, 10);
  unsigned int assert_facto = 1;
  if (mpz_probab_prime_p (z, 25)) {
    factor_append(factor, z);
    assert_facto = 0;
  } else {
    mpz_t tmp;
    mpz_init(tmp);
    assert_facto = brute_force_factorize(factor, tmp, &nb, z, z);
    ASSERT(mpz_cmp_ui(tmp, 1) == 0);
    mpz_clear(tmp);
  }
  factor_realloc(factor, nb);
  sort_factor(factor);
  XXXXXXXXXXXXXXXXXXXXXXXX BROKEN

  ASSERT(assert_facto == 0);
#ifndef NDEBUG
  ASSERT(factor_assert(factor, z) == 0);
#endif

  return assert_facto;
}
#endif

#if 0

XXXXXXXXXXXXXXX BROKEN
unsigned int gmp_factorize(factor_ptr factor, mpz_srcptr z)
{
  unsigned int number = mpz_sizeinbase(z, 2);
  factor_init(factor, number);
  unsigned int nb = 0;
  unsigned int assert_facto = 1; 
  if (mpz_probab_prime_p (z, 25)) {
    mpz_set(factor->factorization[nb], z);
    nb = 1;
  } else {
    mpz_t z_res;
    mpz_init(z_res);

#ifdef FACTO_BRUTE_FORCE
    mpz_t bound;
    mpz_init(bound);
    mpz_set_ui(bound, 1048576);
    assert_facto = brute_force_factorize(factor, z_res, &nb, z, bound);
#else
    mpz_set(z_res, z);
#endif

    if (assert_facto) {
      ASSERT(assert_facto == 1);

      assert_facto = factorize(factor, z_res, &nb);
    }
#ifdef FACTO_BRUTE_FORCE
    mpz_clear(bound);
#endif

    mpz_clear(z_res);
  }
  if (nb == 0) {
    mpz_set(factor->factorization[nb], z);
    nb = 1;
  }
  factor_realloc(factor, nb);
  sort_factor(factor);

#ifndef NDEBUG
  ASSERT(factor_assert(factor, z) == assert_facto);
#endif

  return assert_facto;
}
#endif

int mpz_invert_ui(mpz_ptr rop, mpz_srcptr op1, const uint64_t op2)
{
  mpz_t op3;
  mpz_init(op3);
  mpz_set_ui(op3, op2);
  int res = mpz_invert(rop, op1, op3);
  mpz_clear(op3);

  return res;
}
