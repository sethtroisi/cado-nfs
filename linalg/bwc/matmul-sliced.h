#ifndef MATMUL_SLICED_H_
#define MATMUL_SLICED_H_

#include "matmul.h"

#ifdef __cplusplus
extern "C" {
#endif

struct matmul_sliced_data_s;

extern struct matmul_sliced_data_s * matmul_sliced_init(void*, param_list pl, int);
extern void matmul_sliced_build_cache(struct matmul_sliced_data_s *, uint32_t *);
extern int matmul_sliced_reload_cache(struct matmul_sliced_data_s *);
extern void matmul_sliced_save_cache(struct matmul_sliced_data_s *);
extern void matmul_sliced_mul(struct matmul_sliced_data_s *, void *, void const *, int);
extern void matmul_sliced_report(struct matmul_sliced_data_s *, double);
extern void matmul_sliced_clear(struct matmul_sliced_data_s * mm);
extern void matmul_sliced_aux(struct matmul_sliced_data_s * mm, int op, ...);
extern void matmul_sliced_auxv(struct matmul_sliced_data_s * mm, int op, va_list ap);

#ifdef __cplusplus
}
#endif

#endif	/* MATMUL_SLICED_H_ */
