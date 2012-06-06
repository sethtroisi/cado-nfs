#ifndef __FB_H__
#define __FB_H__

#include "types.h"
#include "sublat.h"



// Initialize a factor base, reading the ideals from a file.
// The factor base should be sorted by increasing degree, starting at degree
// sorted_min_degp; ideals of lower degree don't have to be sorted.
// If max_degp is nonzero, the factor base is read until reaching ideals of
// degree max_degp.
// FIXME: fix max_degp so that this bound is exclusive, not inclusive!
// Return 1 if successful.
int factor_base_init(factor_base_ptr FB, const char *filename,
                     unsigned sorted_min_degp, unsigned max_degp);

int factor_base_init2(large_factor_base_ptr LFB, small_factor_base_ptr SFB,
        const char *filename, unsigned sorted_min_degp, unsigned max_degp);

// Precompute lambda for each element of the factor base.
void factor_base_precomp_lambda(factor_base_ptr FB, qlat_srcptr qlat,
        sublat_ptr sublat);

// Return the largest degree of the ideals in the factor base.
// /!\ Assume that the factor base is sorted, at least for the ideals higher
//     degree. The highest degree should thus be the degree of the last
//     ideal in the factor base.
// FIXME: should this be inclusive or exclusive?
unsigned factor_base_max_degp(factor_base_srcptr FB);

// Clean up memory.
void factor_base_clear(factor_base_ptr FB);

#endif   /* __FB_H__ */
