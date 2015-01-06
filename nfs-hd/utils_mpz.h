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
  IN:
    b: mpz_srcptr, a number.
    a: mpz_ptr, the absolute value of b.
*/
void gmp_abs(mpz_ptr a, mpz_srcptr b);

/*
  Set the maximum between a and b in max.
  IN:
    max: mpz_ptr, the maximum.
    a: mpz_srcptr, a number.
    b: mpz_srcptr, a number.
 */
void gmp_max(mpz_ptr max, mpz_srcptr a, mpz_srcptr b);

/*
  Compute the factorial of an integer a, and set the result in a mpz_t res.
  IN:
    res: gmp_ptr, the result.
    a: int, the integer.
 */
void gmp_factorial_int(mpz_ptr res, int a);

/*
  Initialize an array of factor with number, the maximum number of elements.
  IN:
    factor: factor_ptr, the array.
    number: unsigned int, the maximum number of elements.
*/
void factor_init(factor_ptr factor, unsigned int number);

/*
  Delete an array of factors.
  IN:
    factor: factor_ptr, the array.
*/
void factor_clear(factor_ptr factor);

/*
  Realloc an array of factor. Number must be less than factor->number.
  IN:
    factor: factor_ptr, the array.
    number: unsigned int, new number of elements.
*/
void factor_realloc(factor_ptr factor, unsigned int number);

void sort_factor(factor_ptr factor);

/*
  Factorize a mpz and set the factors in factor.
  IN:
    factor: factor_ptr, an array of factors.
    z: mpz_t, the number we want to factorize.
*/
void gmp_factorize(factor_ptr factor, mpz_t z);

/*
  Print an array of factors.
  IN:
    factor: factor_srcptr, an array of factors.
*/
void factor_printf(factor_srcptr factor);

void factor_fprintf(FILE * file, factor_srcptr factor);

int factor_is_smooth(factor_srcptr factor, mpz_t B);

int mpz_invert_ui(mpz_ptr rop, mpz_srcptr op1, const uint64_t op2);

#endif	/* UTILS_MPZ_H */
