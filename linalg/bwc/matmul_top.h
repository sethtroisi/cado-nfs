#ifndef MATMUL_TOP_H_
#define MATMUL_TOP_H_

#include <stdint.h>

#include "abase.h"
#include "select_mpi.h"
#include "intersections.hpp"
#include "parallelizing_info.h"
#include "matmul.h"

#define CONJUGATED_PERMUTATIONS

/* A wiring is the information concerning our matrix in one of its two
 * dimensions. Everything is made so that inverting the direction flag
 * should make it possible to compute with the transposed matrix easily.
 * 
 * Wiring number 0 is ``horizontal''. In detail, what this means is:
 *
 * n is the total number of rows.
 * v is the LEFT vector (result of matrix-times-vector, or input of
 *   vectora-times-matrix).
 * fences denote the row indices within the matrix where the different
 *   jobs are split.
 * x denotes the intersection data for communication across rows.
 * i0 denotes the row indices of our local matrix chunk.
 *
 * Horizontal indices i0 to i1 of the LEFT column vector are stored locally.
 * --> v is a chunk of i1-i0 data elements.
 *
 * Communication with the left vector. The intent is to obtain a right
 * vector.
 * --> reduce is done across rows. The destination for row index i is
 *   column index i on conjugation mode(*)
 * XXX TBC.
 *
 * (*) it's somewhat different in non-conjugation mode, since a
 * permutation has to be taken (in which case allreduce is probably
 * preferrable).
 */

struct mmt_wiring {
    // global
    unsigned int n;

    abt * v;
#ifdef  CONJUGATED_PERMUTATIONS
    /* This is likely to be relevant only for conjugated matrices */
    unsigned int * fences;
    struct isect_info * x;
    unsigned int xlen;
#endif
    unsigned int i0;
    unsigned int i1;
};

struct matmul_top_data_s {
    abobj_ptr abase;
    parallelizing_info_ptr pi;
    matmul_ptr mm;

    // w <- M v, or v^T <- w^T M would make sense. Of course v^T and w^T
    // are represented in exactly the same way.

    // global stuff
    unsigned int ncoeffs_total;

    // local stuff
    unsigned int ncoeffs;

    char * filename;

    struct mmt_wiring wr[2][1];
};

typedef struct matmul_top_data_s matmul_top_data[1];
typedef struct matmul_top_data_s * matmul_top_data_ptr;
typedef struct matmul_top_data_s const * matmul_top_data_srcptr;

#ifdef __cplusplus
extern "C" {
#endif

void matmul_top_init(matmul_top_data_ptr mmt,
        abobj_ptr abase,
        parallelizing_info_ptr pi,
        const char * filename);
void matmul_top_read_matrix(matmul_top_data_ptr mmt);
void matmul_top_clear(matmul_top_data_ptr mmt, abobj_ptr abase);
void matmul_top_fill_random_right(matmul_top_data_ptr mmt);
void matmul_top_fill_random_left(matmul_top_data_ptr mmt);


#ifdef __cplusplus
}
#endif

#endif	/* MATMUL_TOP_H_ */
