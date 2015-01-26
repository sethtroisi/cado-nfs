#ifndef UTILS_MPZ_H
#define UTILS_MPZ_H

#include <stdio.h>
#include <stdint.h>
#include <gmp.h>

//Define an array with all the factors.
typedef struct
{
  mpz_t * factorization;
  unsigned int number;
} s_factor_t;

typedef s_factor_t factor_t[1];
typedef s_factor_t * factor_ptr;
typedef const s_factor_t * factor_srcptr;

/*
  WARNING: redundant with mpz_abs.
  Set a = abs(b).

  b: mpz_srcptr, a number.
  a: mpz_ptr, the absolute value of b.
*/
void gmp_abs(mpz_ptr a, mpz_srcptr b);

/*
  Set the maximum between a and b in max.

  max: the maximum.
  a: a number.
  b: a number.
 */
void gmp_max(mpz_ptr max, mpz_srcptr a, mpz_srcptr b);

/*
  Compute the factorial of an integer a, and set the result in a mpz_t res.

  res: the result.
   a: the integer.
 */
void gmp_factorial_int(mpz_ptr res, int a);

/*
  Initialize an array of factor with number, the maximum number of elements.

  factor: the array.
  number: the maximum number of elements.
*/
void factor_init(factor_ptr factor, unsigned int number);

/*
  Delete an array of factors.

  factor: the array.
*/
void factor_clear(factor_ptr factor);

/*
  Realloc an array of factor. Number must be less than factor->number.

  factor: the array.
  number: new number of elements.
*/
void factor_realloc(factor_ptr factor, unsigned int number);

void sort_factor(factor_ptr factor);

/*
  Factorize a mpz and set the factors in factor. Return 0 if the factorization
   is correct, 1 else.

  factor: an array of factors.
  z: the number we want to factorize.
*/
unsigned int gmp_factorize(factor_ptr factor, mpz_t z);

/*
  Print an array of factors.

  factor: factor_srcptr, an array of factors.
*/
void factor_fprintf(FILE * file, factor_srcptr factor);

int factor_is_smooth(factor_srcptr factor, mpz_t B);

int mpz_invert_ui(mpz_ptr rop, mpz_srcptr op1, const uint64_t op2);

#endif	/* UTILS_MPZ_H */
