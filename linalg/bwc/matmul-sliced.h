#ifndef MATMUL_SLICED_H_
#define MATMUL_SLICED_H_

#include "matmul.h"

#ifdef __cplusplus
extern "C" {
#endif

extern matmul_ptr matmul_sliced_build(abobj_ptr, const char * filename);
extern matmul_ptr matmul_sliced_reload_cache(abobj_ptr, const char * filename);
extern void matmul_sliced_save_cache(matmul_ptr, const char * filename);
extern void matmul_sliced_mul(matmul_ptr, abt *, abt const *, int);
extern void matmul_sliced_report(matmul_ptr);
extern void matmul_sliced_clear(matmul_ptr mm);

#ifdef __cplusplus
}
#endif

#endif	/* MATMUL_SLICED_H_ */
