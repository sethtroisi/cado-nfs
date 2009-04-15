#ifndef MATMUL_BASIC_H_
#define MATMUL_BASIC_H_

#include "matmul.h"

#ifdef __cplusplus
extern "C" {
#endif

extern matmul_ptr matmul_basic_build(abobj_ptr, const char * filename);
extern matmul_ptr matmul_basic_reload_cache(abobj_ptr, const char * filename);
extern void matmul_basic_save_cache(matmul_ptr, const char * filename);
extern void matmul_basic_mul(matmul_ptr, abt *, abt const *, int);
extern void matmul_basic_report(matmul_ptr);
extern void matmul_basic_clear(matmul_ptr mm);
extern void matmul_basic_aux(matmul_ptr mm, int op, ...);
extern void matmul_basic_auxv(matmul_ptr mm, int op, va_list ap);

#ifdef __cplusplus
}
#endif

#endif	/* MATMUL_BASIC_H_ */
