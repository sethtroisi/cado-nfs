#ifndef GPF_H_
#define GPF_H_

/* A look-up table of the largest prime factor of an integer

Before calling gpf_get(i), gpf_init(m) must have been called with m >= i.

*/

#include "macros.h"

#ifdef __cplusplus
extern "C" {
#endif

extern unsigned int *gpf;
extern unsigned int gpf_max;

/* caution: gpf_init is not MT-safe. Caller must take appropriate measures */
void gpf_init(unsigned int);
static inline unsigned int gpf_get(const unsigned long i) {
    ASSERT(i <= gpf_max);
    return gpf[i];
}

void gpf_clear();

#ifdef __cplusplus
}
#endif

#endif
