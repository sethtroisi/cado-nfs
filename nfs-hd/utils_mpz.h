#ifndef UTILS_MPZ_H
#define UTILS_MPZ_H

#include <stdio.h>
#include <stdint.h>
#include <gmp.h>

/*
 * Compute the inversion of op1 mod op2 and set the result to rop.
 *
 * rop: the result of the inversion.
 * op1: the number for which we want to compute the inversion.
 * op2: the modulo.
 */
int mpz_invert_ui(mpz_ptr rop, mpz_srcptr op1, const uint64_t op2);

#endif	/* UTILS_MPZ_H */
