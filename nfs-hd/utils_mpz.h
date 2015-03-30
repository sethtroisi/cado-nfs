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
 * Initialize an array of factors with number, the maximum number of elements.
 *
 * factor: the array of factors.
 * number: the maximum number of elements.
 */
void factor_init(factor_ptr factor, unsigned int number);

/*
 * Delete an array of factors.
 *
 * factor: the array of factors.
 */
void factor_clear(factor_ptr factor);

/*
 * Realloc an array of factors. Number must be less than factor->number.
 *
 * factor: the array of factors.
 * number: new number of elements.
 */
void factor_realloc(factor_ptr factor, unsigned int number);

/*
 * To sort by acending value the element of the factor array.
 *
 * factor: array of factors.
 */
void sort_factor(factor_ptr factor);

/*
 * Factorize a mpz and set the factors in factor. Return 0 if the factorization
 *  is correct, 1 else.
 *
 * factor: an array of factors.
 * z: the number we want to factorize.
 */
unsigned int gmp_factorize(factor_ptr factor, mpz_t z);

/*
 * Print an array of factors.
 *
 * factor: an array of factors.
 */
void factor_fprintf(FILE * file, factor_srcptr factor);

/*
 * Test if the maximum factor of the array factor is less or equal to B. Return
 *  1 if true, 0 otherwise. If factor is sorted, set sort to 1, 0 otherwise.
 *
 * factor: an array of factors.
 * B: the smoothness bound.
 * sort: 1 if factor is sorted, 0 otherwise.
 */
unsigned int factor_is_smooth(factor_srcptr factor, mpz_t B, unsigned int sort);

/*
 * Compute the inversion of op1 mod op2 and set the result to rop.
 *
 * rop: the result of the inversion.
 * op1: the number for which we want to compute the inversion.
 * op2: the modulo.
 */
int mpz_invert_ui(mpz_ptr rop, mpz_srcptr op1, const uint64_t op2);

#endif	/* UTILS_MPZ_H */
