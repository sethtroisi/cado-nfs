#ifndef MATMUL_TOP_H_
#define MATMUL_TOP_H_

#include <stdint.h>

#include "abase.h"
#include "select_mpi.h"
#include "intersections.h"
#include "parallelizing_info.h"
#include "matmul.h"

#define CONJUGATED_PERMUTATIONS

/* A wiring is the information concerning our matrix in one of its two
 * dimensions. Everything is made so that inverting the direction flag
 * should make it possible to compute with the transposed matrix easily.
 * 
 * Wiring number 0 is ``horizontal''. In detail, what this means is:
 *
 * n is the total number of rows ; the total number of instances of this
 * ``wiring'', so to say.
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
 *
 *
 *
 * Note also that we have the pi_wiring to consider as well. Those
 * have nothing to do with the matrix, and are relevant with regard to
 * the parallelization of the work.
 */

// the ``all_v'' field collects all the pointers to the per-thread vector
// values. There are exactly pi->wr[0]->ncores such pointers in
// mmt->wr[0]. In some cases, these pointers may be equal (for source
// vectors, never used as destination), and in some cases not.

struct mmt_wiring_s {
    abt * v;
    abt * * all_v;
    unsigned int i0;
    unsigned int i1;
#ifdef  CONJUGATED_PERMUTATIONS
    /* This is likely to be relevant only for conjugated matrices */
    struct isect_info * x;
    unsigned int xlen;

    // this second intersection list is used for reduction, where the row
    // threads collectively build the result vector from their respective
    // parts -- deciding where the computed sum eventually goes is
    // determined by the intersection of i0+k/n*(i1-i0)..i0+(k+1)/n*(i1-i0)
    // with the vertical fences.
    //
    // note that the way reduction is performed is likely to cause some
    // cache invalidates unless some fences end up correctly aligned.
    struct isect_info * y;
    unsigned int ylen;
#endif
};
typedef struct mmt_wiring_s mmt_wiring[1];
typedef struct mmt_wiring_s * mmt_wiring_ptr;
typedef struct mmt_wiring_s const * mmt_wiring_srcptr;

struct matmul_top_data_s {
    abobj_ptr abase;
    parallelizing_info_ptr pi;
    matmul_ptr mm;

    // w <- M v, or v^T <- w^T M would make sense. Of course v^T and w^T
    // are represented in exactly the same way.

    // global stuff. n[] is going to be misleading in any situation. It's
    // the size of the matrix in the given direction. Meaning that n[0]
    // is the horizontal size -- the size of the rows --. This means the
    // number of columns.
    unsigned int n[2];
    unsigned int ncoeffs_total;

    // local stuff
    unsigned int ncoeffs;
    
    char * locfile;

    mmt_wiring wr[2];

#ifdef  CONJUGATED_PERMUTATIONS
    // fences[0]: horizontal fences == row indices (horizontal separating line)
    // fences[1]: vertical fences == col indices (vertical separating line)
    // fences[0] contains 1 + pi->wr[1]->totalsize values
    // fences[1] contains 1 + pi->wr[0]->totalsize values
    //
    // it turns out later on that this numbering is slightly orthogonal
    // to the one used in most of the rest of the code, so maybe it's a
    // bad decision. Fortunately, fences[] is rarely used.
    unsigned int * fences[2];
#endif

    // internal use
    int flags[2];
};

/* THREAD_MULTIPLE_VECTOR is when several threads will be writing to the
 * vector simultaneously. This implies having separate areas.
 *
 * TODO: Think about waiving this restriction in the spirit of bucket
 * sieving, maybe */
#define THREAD_MULTIPLE_VECTOR    0
#define THREAD_SHARED_VECTOR    1

typedef struct matmul_top_data_s matmul_top_data[1];
typedef struct matmul_top_data_s * matmul_top_data_ptr;
typedef struct matmul_top_data_s const * matmul_top_data_srcptr;

#ifdef __cplusplus
extern "C" {
#endif

void matmul_top_init(matmul_top_data_ptr mmt,
        abobj_ptr abase,
        /* matmul_ptr mm, */
        parallelizing_info_ptr pi,
        int const * flags,
        const char * filename);

void mmt_finish_init(matmul_top_data_ptr mmt);

void matmul_top_read_submatrix(matmul_top_data_ptr mmt);
void matmul_top_clear(matmul_top_data_ptr mmt, abobj_ptr abase);
void matmul_top_fill_random_source(matmul_top_data_ptr mmt, int d);
void matmul_top_load_vector(matmul_top_data_ptr mmt, int d, unsigned int index, unsigned int iter);
void matmul_top_save_vector(matmul_top_data_ptr mmt, int d, unsigned int index, unsigned int iter);
void matmul_top_mul(matmul_top_data_ptr mmt, int d);


#ifdef __cplusplus
}
#endif

#endif	/* MATMUL_TOP_H_ */
