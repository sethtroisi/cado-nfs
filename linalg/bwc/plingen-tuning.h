#ifndef PLINGEN_TUNING_H_
#define PLINGEN_TUNING_H_

#include "mpfq_layer.h"
#include "params.h"
#include "select_mpi.h"

#ifdef __cplusplus
extern "C" {
#endif

void plingen_tuning(abdst_field ab, unsigned int m, unsigned int n, MPI_Comm comm, param_list_ptr pl);

#ifdef __cplusplus
}
#endif

#endif	/* PLINGEN_TUNING_H_ */
