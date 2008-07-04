#ifndef CADO_UTILS_ROOTFINDER_H_
#define CADO_UTILS_ROOTFINDER_H_

#include <gmp.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* This is the entry point for the root finding routines.
 *
 * It relies on either the mod_ul inline assembly layer, or the plain C
 * layer in plain_poly.
 */

extern int poly_roots(mpz_t * r, mpz_t * f, int d, mpz_t p);
extern int poly_roots_long(long * r, mpz_t * f, int d, unsigned long p);
extern int poly_roots_ulong(unsigned long * r, mpz_t * f, int d, unsigned long p);
extern int poly_roots_uint64(uint64_t * r, mpz_t * f, int d, uint64_t p);


#ifdef __cplusplus
}
#endif

#endif	/* CADO_UTILS_ROOTFINDER_H_ */
