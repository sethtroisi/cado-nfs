#ifndef MATMUL_BASIC_H_
#define MATMUL_BASIC_H_

#include "matmul.h"

#ifdef __cplusplus
extern "C" {
#endif

struct matmul_basic_data_s;

extern struct matmul_basic_data_s * matmul_basic_build(abobj_ptr, const char * filename, param_list pl, int);
extern struct matmul_basic_data_s * matmul_basic_reload_cache(abobj_ptr, const char * filename, param_list pl, int);
extern void matmul_basic_save_cache(struct matmul_basic_data_s *, const char * filename);
extern void matmul_basic_mul(struct matmul_basic_data_s *, abt *, abt const *, int);
extern void matmul_basic_report(struct matmul_basic_data_s *);
extern void matmul_basic_clear(struct matmul_basic_data_s * mm);
extern void matmul_basic_aux(struct matmul_basic_data_s * mm, int op, ...);
extern void matmul_basic_auxv(struct matmul_basic_data_s * mm, int op, va_list ap);

#ifdef __cplusplus
}
#endif

#endif	/* MATMUL_BASIC_H_ */
