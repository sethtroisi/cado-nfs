#ifndef __LATSIEVE_H__
#define __LATSIEVE_H__

#include "types.h"
#include "fb.h"
#include "qlat.h"
#include "sublat.h"

// this will probably be split into 2 or more functions once the bucket
// sieving is operational.

void sieveFB(uint8_t *S, small_factor_base_ptr FB, unsigned I, unsigned J,
             ij_t j0, ijpos_t pos0, ijpos_t size, sublat_ptr sublat);


#endif   /* __LATSIEVE_H__ */
