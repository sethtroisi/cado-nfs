#ifndef GPF_H_
#define GPF_H_

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
