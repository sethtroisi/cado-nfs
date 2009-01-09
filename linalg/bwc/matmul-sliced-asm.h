#ifndef MATMUL_SLICED_ASM_H_
#define MATMUL_SLICED_ASM_H_

#include "abase.h"

#ifdef __cplusplus
extern "C" {
#endif

const uint16_t * matmul_sliced_asm(abobj_ptr x, abt * where, const abt * from, const uint16_t * q);

#ifdef __cplusplus
}
#endif

#endif	/* MATMUL_SLICED_ASM_H_ */

