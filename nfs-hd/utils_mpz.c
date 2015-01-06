#include "utils_mpz.h"
#include <stdint.h>
#include "gmp-ecm.h"
#include <stdlib.h>
#include "macros.h"

void gmp_abs(mpz_ptr a, mpz_srcptr b)
{
  if (mpz_cmp_si(b, 0) < 0) {
    mpz_mul_si(a, b, -1);
    return;
  }
  mpz_set(a, b);
}

void gmp_max(mpz_ptr max, mpz_srcptr a, mpz_srcptr b)
{
  if (mpz_cmp(a, b) < 0) {
    mpz_set(max, b);
  } else if (mpz_cmp(a, b) > 0) {
    mpz_set(max, a);
  } else {
    mpz_set(max, a);
  }
}

void gmp_factorial_int(mpz_ptr res, int a)
{
  mpz_set_si(res, 1);
  for (int i = 2; i <= a; i++) {
    mpz_mul_si(res, res, i);
  }
}

void factor_init(factor_ptr factor, unsigned int number)
{
  factor->number = number;
  factor->factorization = (mpz_t * ) malloc(sizeof(mpz_t) * number);
  for (unsigned int i = 0; i < number; i++) {
    mpz_init(factor->factorization[i]);
  }
}

void factor_clear(factor_ptr factor)
{
  for (unsigned int i = 0; i < factor->number; i++) {
    mpz_clear(factor->factorization[i]);
  }
  free(factor->factorization);
  factor->number = 0;
}

void factor_realloc(factor_ptr factor, unsigned int number)
{
  ASSERT(factor->number > number);
  for (unsigned int i = number; i < factor->number; i++) {
    mpz_clear(factor->factorization[i]);
  }
  factor->factorization = realloc(factor->factorization, sizeof(mpz_t) *
				  number);
  factor->number = number;
}

static void factorize(factor_ptr factor, mpz_t z, unsigned int * number)
{
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
  mpz_clear(res);
}

static int compare(const void * p1, const void * p2)
{
  return(mpz_cmp((mpz_srcptr) p1, (mpz_srcptr) p2));
}

void sort_factor(factor_ptr factor)
{
  qsort(factor->factorization, factor->number, sizeof(factor->factorization[0]),
	compare);
}

void gmp_factorize(factor_ptr factor, mpz_t z)
{
  unsigned int number = mpz_sizeinbase(z, 2);
  factor_init(factor, number);
  unsigned int nb = 0;
  if (mpz_probab_prime_p (z, 25)) {
    mpz_set(factor->factorization[nb], z);
    nb = 1;
  } else {
    factorize(factor, z, &nb);
  }
  factor_realloc(factor, nb);
  sort_factor(factor);
}

void factor_printf(factor_srcptr factor)
{
  printf("[");
  for (unsigned int i = 0; i < factor->number - 1; i++) {
    gmp_printf("%Zd, ", factor->factorization[i]);
  }
  gmp_printf("%Zd]\n", factor->factorization[factor->number - 1]);
}

void factor_fprintf(FILE * file, factor_srcptr factor)
{
  fprintf(file, "[");
  for (unsigned int i = 0; i < factor->number - 1; i++) {
    gmp_fprintf(file, "%Zd, ", factor->factorization[i]);
  }
  gmp_fprintf(file, "%Zd]\n", factor->factorization[factor->number - 1]);
}

int factor_is_smooth(factor_srcptr factor, mpz_t B)
{
  for (unsigned int i = 0; i < factor->number; i++) {
    if (mpz_cmp(factor->factorization[i], B) >= 0) {
      return 0;
    }
  }
  return 1;
}

int mpz_invert_ui(mpz_ptr rop, mpz_srcptr op1, const uint64_t op2)
{
  mpz_t op3;
  mpz_init(op3);
  mpz_set_ui(op3, op2);
  int res = mpz_invert(rop, op1, op3);
  mpz_clear(op3);
  return res;
}
