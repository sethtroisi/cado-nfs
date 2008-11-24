#ifndef MATMUL_H_
#define MATMUL_H_

#include "abase.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef void * matmul_t;
typedef void * matmul_ptr;

extern matmul_ptr matmul_build(abobj_ptr, const char * filename);
extern matmul_ptr matmul_reload_cache(abobj_ptr, const char * filename);
extern void matmul_save_cache(matmul_ptr, const char * filename);
extern void matmul_report(matmul_ptr);
extern void matmul(matmul_ptr, abt *, abt const *);
extern void matmul_clear(matmul_ptr mm);

#ifdef __cplusplus
}
#endif

#endif	/* MATMUL_H_ */
