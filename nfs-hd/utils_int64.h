#ifndef UTILSINT64_T_H
#define UTILSINT64_T_H

#include <stdint.h>
#include "cado.h"
#include "utils.h"
#include "mod_ul_default.h"
#include "mod_ul.h"

/*
 * Compute d^e.
 *
 * d: a number.
 * e: a number.
 */
uint64_t pow_uint64_t(uint64_t d, int64_t e);

/*
 * Compute the modular inverse of xx mod mm.
 *
 * xx: the number for which we want the inverse.
 * mm: the modulo.
 */
uint64_t invmod_uint64(uint64_t xx, uint64_t mm);

#endif /* UTILSINT64_T_H */
