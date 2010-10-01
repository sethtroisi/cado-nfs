#ifndef BLOCKMATRIX_H_
#define BLOCKMATRIX_H_

#include "bit_matrices.h"
#include "macros.h"

#ifdef __cplusplus
extern "C" {
#endif

struct blockmatrix_s {
    mat64 * mb;
    unsigned int nrblocks;   // row blocks
    unsigned int ncblocks;   // col blocks
    unsigned int nrows;      // only for convenience
    unsigned int ncols;      // only for convenience
    int sub;                 // marker indicating a submatrix.
    unsigned int stride;     // may differ from nrblocks if sub!=0
};

typedef struct blockmatrix_s * blockmatrix;

extern blockmatrix blockmatrix_alloc(unsigned int nrows, unsigned int ncols);
extern void blockmatrix_free(blockmatrix b);
extern void blockmatrix_zero(blockmatrix b);
extern void blockmatrix_mul_Ta_b(blockmatrix c,
        const blockmatrix a,
        const blockmatrix b);
extern void blockmatrix_mul_smallb(blockmatrix c,
        const blockmatrix a,
        const blockmatrix b);
extern void blockmatrix_copy_to_flat(uint64_t * tiny, unsigned int stride,
        int i0, int j0, blockmatrix m);
extern void blockmatrix_copy_transpose_to_flat(uint64_t * tiny, unsigned int stride,
        int i0, int j0, blockmatrix m);
extern void blockmatrix_copy_transpose_from_flat(blockmatrix m, uint64_t * tiny, unsigned int stride, int i0, int j0);
extern void blockmatrix_copy_from_flat(blockmatrix m, uint64_t * tiny, unsigned int stride, int i0, int j0);
extern void blockmatrix_read_from_flat_file(blockmatrix k, int i0, int j0, const char * name, unsigned int fnrows, unsigned int fncols);
extern void blockmatrix_write_to_flat_file(const char * name, blockmatrix k, int i0, int j0, unsigned int fnrows, unsigned int fncols);
extern void blockmatrix_transpose(blockmatrix b, blockmatrix a);
extern blockmatrix blockmatrix_submatrix(blockmatrix k, int i0, int j0, unsigned int nrows, unsigned int ncols);
// extern void blockmatrix_read_transpose_from_flat_file(blockmatrix k, int i0, int j0, const char * name, unsigned int fnrows, unsigned int fncols);
extern void blockmatrix_swap(blockmatrix B, blockmatrix A);

/* Use this macro to allocate flat matrix areas with proper readahead
 * padding. In some situations, it is also necessary to use zero out the
 * padding data as well, so this macro must also be used in the memset()
 * calls.
 */
#define FLAT_BYTES_WITH_READAHEAD(nr, nc) \
    iceildiv((nr), 64) * iceildiv((nc), 64) * sizeof(mat64)

#ifdef __cplusplus
}
#endif

#endif	/* BLOCKMATRIX_H_ */
