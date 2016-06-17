#ifndef MATMUL_MF_H_
#define MATMUL_MF_H_

#include "matmul.h"
#include "raw_matrix_u32.h"

#ifdef __cplusplus
extern "C" {
#endif

void mf_prepare_matrix_u32(matmul_ptr mm, matrix_u32_ptr m, const char * file, int withcoeffs);

#ifdef __cplusplus
}
#endif

#endif	/* MATMUL_MF_H_ */
