#ifndef XYMATS_H_
#define XYMATS_H_

#include "matmul_top.h"
#include "parallelizing_info.h"

#ifdef __cplusplus
extern "C" {
#endif

/* This is used as a very simple reduction routine, mostly for the small
 * m*n matrices. It has nothing to do with the more involved reduce
 * operations relevant to vectors. Here, the data we're dealing with is
 * tiny.
 *
 * Clearly the function names should be changed.
 *
 * This operation serializes threads.
 */

// void allreduce_generic_threadlevel(mmt_vec_ptr v, pi_wiring_ptr wr, unsigned int n);
// void allreduce_generic_mpilevel(mmt_vec_ptr v, pi_wiring_ptr wr, unsigned int n);
void allreduce_generic(mmt_vec_ptr v, pi_wiring_ptr wr, unsigned int n);

void broadcast_generic(mmt_vec_ptr v, pi_wiring_ptr wr, unsigned int n, unsigned int j0, unsigned int t0);

#ifdef __cplusplus
}
#endif

#endif	/* XYMATS_H_ */
