#ifndef BALANCING_H_
#define BALANCING_H_

#include "parallelizing_info.h"
#include "select_mpi.h"

#include "utils.h"

#include "mf.h"

#ifdef __cplusplus
extern "C" {
#endif

struct matrix_u32_s {
    // input arguments.
    const char * mfile;    // matrix file name
    const char * bfile;    // balancing file name ; NULL will mean auto-detect
    int transpose;
    // output arguments.
    uint32_t * p;
    size_t size;
};
typedef struct matrix_u32_s matrix_u32[1];
typedef struct matrix_u32_s * matrix_u32_ptr;

void * get_matrix_u32(parallelizing_info_ptr pi, param_list pl, matrix_u32_ptr arg);

#ifdef __cplusplus
}
#endif

#endif	/* BALANCING_H_ */
