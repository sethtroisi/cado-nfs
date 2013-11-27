#ifndef LAS_DUPLICATE_H_
#define LAS_DUPLICATE_H_

#include "las-types.h"
#include "relation.h"

#ifdef __cplusplus
extern "C" {
#endif

int relation_is_duplicate(relation_t *, double, sieve_info_srcptr);

#ifdef __cplusplus
}
#endif

#endif
