#ifndef __FB_H__
#define __FB_H__

#include "types.h"
#include "sublat.h"

// Initialize a factor base, reading the ideals from a file.
// Return 1 if successful.
int factor_base_init(factor_base_ptr FB, const char *filename,
                     unsigned maxdeg);

// Precompute lambda for each element of the factor base.
void factor_base_precomp_lambda(factor_base_ptr FB, qlat_srcptr qlat,
        sublat_ptr sublat);

// Clean up memory.
void factor_base_clear(factor_base_ptr FB);

#endif   /* __FB_H__ */
