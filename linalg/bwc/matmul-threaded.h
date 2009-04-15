#ifndef MATMUL_THREADED_H_
#define MATMUL_THREADED_H_

#include "matmul.h"

#ifdef __cplusplus
extern "C" {
#endif

extern matmul_ptr matmul_threaded_build(abobj_ptr, const char * filename);
extern matmul_ptr matmul_threaded_reload_cache(abobj_ptr, const char * filename);
extern void matmul_threaded_save_cache(matmul_ptr, const char * filename);
extern void matmul_threaded_mul(matmul_ptr, abt *, abt const *, int);
extern void matmul_threaded_report(matmul_ptr);
extern void matmul_threaded_clear(matmul_ptr mm);
extern void matmul_threaded_aux(matmul_ptr mm, int op, ...);
extern void matmul_threaded_auxv(matmul_ptr mm, int op, va_list ap);

#ifdef __cplusplus
}
#endif

#endif	/* MATMUL_THREADED_H_ */
