#ifndef XYMATS_H_
#define XYMATS_H_

#include "abase.h"
#include "matmul_top.h"
#include "parallelizing_info.h"

#ifdef __cplusplus
extern "C" {
#endif

// void allreduce_generic_threadlevel(abobj_t abase, mmt_vec_ptr v, pi_wiring_ptr wr, unsigned int n);
// void allreduce_generic_mpilevel(abobj_t abase, mmt_vec_ptr v, pi_wiring_ptr wr, unsigned int n);
void allreduce_generic(abobj_t abase, mmt_vec_ptr v, pi_wiring_ptr wr, unsigned int n);

// slightly different calling interface because broadcasting is not
// abase-dependent.
void broadcast_generic(mmt_generic_vec_ptr v, pi_wiring_ptr wr, size_t siz, unsigned int j0, unsigned int t0);

#ifdef __cplusplus
}
#endif

#endif	/* XYMATS_H_ */
