#ifndef BW_COMMON_MPI_H_
#define BW_COMMON_MPI_H_

#include "bw-common.h"

#ifdef __cplusplus
extern "C" {
#endif

extern int bw_common_init_mpi(struct bw_params * bw, param_list pl, int argc, char * argv[]);
extern int bw_common_clear_mpi(struct bw_params * bw);

#ifdef __cplusplus
}
#endif

#endif	/* BW_COMMON_MPI_H_ */
