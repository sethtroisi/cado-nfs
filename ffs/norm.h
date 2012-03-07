#ifndef __NORM_H__
#define __NORM_H__

#ifdef __cplusplus
extern "C" {
#endif


#include "fppol.h"
#include "ffspol.h"
#include "types.h"

void ffspol_norm(fppol_t norm, ffspol_ptr ffspol, fppol_t a, fppol_t b);
void init_norms(uint8_t *S, ffspol_ptr ffspol, int I, int J, qlat_t qlat,
        int sqside);

#ifdef __cplusplus
}
#endif

#endif   /* __NORM_H__ */
