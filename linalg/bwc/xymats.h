#ifndef XYMATS_H_
#define XYMATS_H_

#include "abase.h"
#include "matmul_top.h"
#include "parallelizing_info.h"

#ifdef __cplusplus
extern "C" {
#endif

void reduce_generic_threadlevel(abobj_t abase, mmt_vec_ptr v, pi_wiring_ptr wr, unsigned int n);
void reduce_generic_mpilevel(abobj_t abase, mmt_vec_ptr v, pi_wiring_ptr wr, unsigned int n);
void reduce_generic(abobj_t abase, mmt_vec_ptr v, pi_wiring_ptr wr, unsigned int n);

#ifdef __cplusplus
}
#endif

#endif	/* XYMATS_H_ */
