#ifndef MATMUL_BUCKET_H_
#define MATMUL_BUCKET_H_

#include "matmul.h"

#ifdef __cplusplus
extern "C" {
#endif

struct matmul_bucket_data_s;

extern struct matmul_bucket_data_s * matmul_bucket_build(abobj_ptr, const char * filename, param_list pl, int);
extern struct matmul_bucket_data_s * matmul_bucket_reload_cache(abobj_ptr, const char * filename, param_list pl, int);
extern void matmul_bucket_save_cache(struct matmul_bucket_data_s *, const char * filename);
extern void matmul_bucket_mul(struct matmul_bucket_data_s *, abt *, abt const *, int);
extern void matmul_bucket_report(struct matmul_bucket_data_s *, double);
extern void matmul_bucket_clear(struct matmul_bucket_data_s * mm);
extern void matmul_bucket_aux(struct matmul_bucket_data_s * mm, int op, ...);
extern void matmul_bucket_auxv(struct matmul_bucket_data_s * mm, int op, va_list ap);

#ifdef __cplusplus
}
#endif

#endif	/* MATMUL_BUCKET_H_ */
