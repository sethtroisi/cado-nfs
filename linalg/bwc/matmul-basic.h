#ifndef MATMUL_BASIC_H_
#define MATMUL_BASIC_H_

#include "matmul.h"

#ifdef __cplusplus
extern "C" {
#endif

struct matmul_basic_data_s;

// The void* here is the concrete implementation from mpfq, which has to
// be passed as the v->obj, where v is the OO interface to the same
// object.
extern struct matmul_basic_data_s * matmul_basic_init(void*, param_list pl, int);
extern void matmul_basic_build_cache(struct matmul_basic_data_s *, uint32_t *);
extern int matmul_basic_reload_cache(struct matmul_basic_data_s *);
extern void matmul_basic_save_cache(struct matmul_basic_data_s *);
extern void matmul_basic_mul(struct matmul_basic_data_s *, void *, void const *, int);
extern void matmul_basic_report(struct matmul_basic_data_s *, double);
extern void matmul_basic_clear(struct matmul_basic_data_s * mm);
extern void matmul_basic_aux(struct matmul_basic_data_s * mm, int op, ...);
extern void matmul_basic_auxv(struct matmul_basic_data_s * mm, int op, va_list ap);

#ifdef __cplusplus
}
#endif

#endif	/* MATMUL_BASIC_H_ */
