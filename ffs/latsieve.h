#ifndef __LATSIEVE_H__
#define __LATSIEVE_H__

#include "types.h"
#include "fb.h"

// this will probably be split into 2 or more functions once the bucket
// sieving is operational.
void sieveFB(uint8_t *S, factorbase_t FB, int I, int J, qlat_t qlat);



#endif   /* __LATSIEVE_H__ */
