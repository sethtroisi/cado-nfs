#ifndef RAW_MATRIX_U32_H_
#define RAW_MATRIX_U32_H_

#include <stdint.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

struct matrix_u32_s {
    // input arguments.
    const char * mfile;    // matrix file name
    const char * bfile;    // balancing file name ; NULL will mean auto-detect
    int transpose;
    int withcoeffs;
    // output arguments.
    uint32_t * p;
    size_t size;
    uint32_t (*twist)[2];
    size_t ntwists;
};
typedef struct matrix_u32_s matrix_u32[1];
typedef struct matrix_u32_s * matrix_u32_ptr;

/* Structure of type matrix_u32 are created empty and zeroed out by the
 * caller functions. mfile and bfile are expected from the caller. The
 * rest is obtained by balancing_get_matrix_u32.
 *
 * There is no destructor function for this type. The global assumption
 * is that the mm layer which is fed with this structure has the right to
 * claim ownership of the data referenced from this structure. If it does
 * not, then it is this layer's responsibility to release the memory. In
 * a sense, matmul_build_cache *is* the destructor...
 *
 * To illustrate this, the layer mm_basic reuses (and in fact, also
 * modifies) the m->p data, and thus it remains in memory.  In contrast,
 * the bucket layer discards this data, as the post-processed matrix form
 * is more efficient eventually.
 */
 
extern void matrix_u32_init_from_file(matrix_u32_ptr m, const char * file, int stored_transposed);


#ifdef __cplusplus
}
#endif

#endif	/* RAW_MATRIX_U32_H_ */
