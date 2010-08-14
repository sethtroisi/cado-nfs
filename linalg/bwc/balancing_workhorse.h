#ifndef BALANCING_WORKHORSE_H_
#define BALANCING_WORKHORSE_H_

#include <stdint.h>

#include "parallelizing_info.h"
#include "utils.h"

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

#ifdef __cplusplus
extern "C" {
#endif

void * balancing_get_matrix_u32(parallelizing_info_ptr pi, param_list pl, matrix_u32_ptr arg);


#ifdef __cplusplus
}
#endif

#endif	/* BALANCING_WORKHORSE_H_ */
