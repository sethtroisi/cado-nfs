#ifndef MATMUL_SUB_LARGE_FBD_H_
#define MATMUL_SUB_LARGE_FBD_H_

#include "abase.h"

#ifdef __cplusplus
extern "C" {
#endif

extern void matmul_sub_large_fbd(abobj_ptr x, abt ** sb, const abt * z, const uint8_t * q, unsigned int n);

#ifdef __cplusplus
}
#endif

#endif	/* MATMUL_SUB_LARGE_FBD_H_ */
