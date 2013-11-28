#ifndef LAS_DUPLICATE_H_
#define LAS_DUPLICATE_H_

#include "las-types.h"
#include "relation.h"

#ifdef __cplusplus
extern "C" {
#endif

int relation_is_duplicate(relation_t *, double, int, sieve_info_srcptr);

/* FIXME: This function is defined in las.c. The prototype should not be here.
   The function should be moved to an appropriate file (las-qlattice.c?)
   and the prototype to the corresponding header file. */
int sieve_info_adjust_IJ(sieve_info_ptr, double, int);


#ifdef __cplusplus
}
#endif

#endif
