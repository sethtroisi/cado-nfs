#ifndef MATMUL_SLICED_H_
#define MATMUL_SLICED_H_

#include "matmul.h"

#ifdef __cplusplus
extern "C" {
#endif

struct matmul_sliced_data_s;

extern struct matmul_sliced_data_s * matmul_sliced_build(abobj_ptr, const char * filename, param_list pl);
extern struct matmul_sliced_data_s * matmul_sliced_reload_cache(abobj_ptr, const char * filename, param_list pl);
extern void matmul_sliced_save_cache(struct matmul_sliced_data_s *, const char * filename);
extern void matmul_sliced_mul(struct matmul_sliced_data_s *, abt *, abt const *, int);
extern void matmul_sliced_report(struct matmul_sliced_data_s *);
extern void matmul_sliced_clear(struct matmul_sliced_data_s * mm);
extern void matmul_sliced_aux(struct matmul_sliced_data_s * mm, int op, ...);
extern void matmul_sliced_auxv(struct matmul_sliced_data_s * mm, int op, va_list ap);

#ifdef __cplusplus
}
#endif

#endif	/* MATMUL_SLICED_H_ */
