#ifndef LAS_DUPLICATE_H_
#define LAS_DUPLICATE_H_

#include "las-types.h"
#include "relation.h"

#ifdef __cplusplus
extern "C" {
#endif

int relation_is_duplicate(relation_t *, double, int, sieve_info_srcptr);

/* FIXME: These function is defined in las.c. The prototypes should not be 
   here. The functions should be moved to an appropriate file 
   (las-qlattice.c?) and the prototypes to the corresponding header file. */
int sieve_info_adjust_IJ(sieve_info_ptr, double, int);
int check_leftover_norm (const mpz_t n, sieve_info_srcptr si, int side);


#ifdef __cplusplus
}
#endif

#endif
