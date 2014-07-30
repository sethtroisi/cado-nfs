#ifndef MATMUL_TOP_H_
#define MATMUL_TOP_H_

#include <stddef.h>
#include <stdint.h>

#include "select_mpi.h"
#include "intersections.h"
#include "parallelizing_info.h"
#include "matmul.h"
#include "params.h"
#include "balancing.h"
#include "misc.h"
#include "mpfq/mpfq_vbase.h"

/* Don't touch this. */
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
 *   vector-times-matrix).
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

struct mmt_vec_s {
    mpfq_vbase_ptr abase;
    size_t stride;      // shorcut to this->abase->vec_elt_stride(this->abase,1)
    void * v;
    void * * all_v;
    // internal use
    int flags;
};
typedef struct mmt_vec_s mmt_vec[1];
typedef struct mmt_vec_s * mmt_vec_ptr;

/*
struct mmt_generic_vec_s {
    void * v;
    void ** all_v;
    // internal use
    int flags;
};
typedef struct mmt_generic_vec_s mmt_generic_vec[1];
typedef struct mmt_generic_vec_s * mmt_generic_vec_ptr;
*/

/* some handy inlines to access portions of mmt_vec's */
static inline void* mmt_vec_subvec(mmt_vec_ptr v, void * p, ptrdiff_t offset)
{
    return pointer_arith(p, v->abase->vec_elt_stride(v->abase, offset));
}
static inline const void* mmt_vec_subvec_const(mmt_vec_ptr v, const void * p, ptrdiff_t offset)
{
    return pointer_arith_const(p, v->abase->vec_elt_stride(v->abase, offset));
}

#define SUBVEC(v,w,offset) mmt_vec_subvec((v), (v)->w, (offset))
#define SUBVEC_const(v,w,offset) mmt_vec_subvec_const((v), (v)->w, (offset))


struct mmt_wiring_s {
    unsigned int i0;
    unsigned int i1;
    /* Note that while i0 and i1 are obviously abase-dependent, clearly
     * the two following fields are not ! */
    mmt_vec v;
    size_t rsbuf_size;          // auto-expanded on demand.
    void * rsbuf[2];             // only for USE_ALTERNATIVE_REDUCE_SCATTER
};
typedef struct mmt_wiring_s mmt_wiring[1];
typedef struct mmt_wiring_s * mmt_wiring_ptr;
typedef struct mmt_wiring_s const * mmt_wiring_srcptr;

/* TODO: It's unhandy to have abase here, since it forces some files to
 * be abase-dependent, without any real reason.
 */
struct matmul_top_data_s {
    parallelizing_info_ptr pi;

    // w <- M v, or v^T <- w^T M would make sense. Of course v^T and w^T
    // are represented in exactly the same way.

    // global stuff. n[] is going to be misleading in any situation. It's
    // the size of the matrix in the given direction. Meaning that n[0]
    // is the horizontal size -- the size of the rows --. This means the
    // number of columns.
    unsigned int n[2];
    unsigned int n0[2]; // n0: unpadded.
    // unsigned int ncoeffs_total;

    // local stuff
    // // unsigned int ncoeffs;
    
    // this really ends up within the mm field.
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
    balancing bal;

    matmul_ptr mm;

    mpfq_vbase_ptr abase;

    int io_truncate;
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

extern void matmul_top_init(matmul_top_data_ptr mmt,
        mpfq_vbase_ptr abase,
        parallelizing_info_ptr pi,
        int const * flags,
        param_list pl,
        int optimized_direction);


extern void matmul_top_clear(matmul_top_data_ptr mmt);
#if 0
extern void matmul_top_fill_random_source(matmul_top_data_ptr mmt, int d);
#endif
extern void matmul_top_load_vector(matmul_top_data_ptr mmt, const char * name, int d, unsigned int iter, unsigned int itemsondisk);
extern void matmul_top_save_vector(matmul_top_data_ptr mmt, const char * name, int d, unsigned int iter, unsigned int itemsondisk);
extern void matmul_top_mul_cpu(matmul_top_data_ptr mmt, int d);
extern void matmul_top_mul_comm(matmul_top_data_ptr mmt, int d);
static inline void matmul_top_mul(matmul_top_data_ptr mmt, int d)
{
    matmul_top_mul_cpu(mmt, d);
    matmul_top_mul_comm(mmt, d);
}

/* decoding ``apply_P_apply_S''. If bw->dir==0, vectors are row vectors.
 * Multiplication is on the right. For d==0, the result is then v*P*S (we
 * multiply by P, then by S). For d==1, the applications are transposed,
 * so the transpose of the result is S*P*trsp(v). In all cases,
 * ``applying'' means multiplying by the named matrix(ces), in order, in
 * the case d==bw->dir, and to the left or to the right depending which
 * is appropriate. The multiplications are transposed in the case of
 * d==!bw->dir
 */
extern void matmul_top_apply_P_apply_S(matmul_top_data_ptr mmt, int d);
extern void matmul_top_unapply_S_unapply_P(matmul_top_data_ptr mmt, int d);
extern void matmul_top_apply_P(matmul_top_data_ptr mmt, int d);
extern void matmul_top_unapply_P(matmul_top_data_ptr mmt, int d);
extern void matmul_top_apply_S(matmul_top_data_ptr mmt, int d);
extern void matmul_top_unapply_S(matmul_top_data_ptr mmt, int d);
extern void matmul_top_twist_vector(matmul_top_data_ptr mmt, int d);
extern void matmul_top_untwist_vector(matmul_top_data_ptr mmt, int d);
extern void indices_apply_S(matmul_top_data_ptr mmt, uint32_t * xs, unsigned int n, int d);
extern void indices_apply_P(matmul_top_data_ptr mmt, uint32_t * xs, unsigned int n, int d);

extern void matmul_top_vec_init(matmul_top_data_ptr mmt, int d, int flags);
extern void matmul_top_vec_clear(matmul_top_data_ptr mmt, int d);

/* Now some of the generic interface calls. By design, not everything is
 * possible with these calls. In particular, nothing critical is doable.
 * As a general convention, the only thing these calls need to know about
 * the abase is the stride value, which corresponds to sizeof(abelt)
 * Besides that, the generic calls take the vector pointer as argument.
 * Specifying NULL as vector argument is equivalent to taking the default
 * vector defined in the mmt data for the given direction flag. */


extern void matmul_top_vec_init_generic(matmul_top_data_ptr mmt, mpfq_vbase_ptr abase, mmt_vec_ptr v, int d, int flags);
extern void matmul_top_vec_clear_generic(matmul_top_data_ptr mmt, mmt_vec_ptr v, int d);
#if 0
extern void matmul_top_fill_random_source_generic(matmul_top_data_ptr mmt, size_t stride, mmt_vec_ptr v, int d);
#endif
extern void matmul_top_load_vector_generic(matmul_top_data_ptr mmt, mmt_vec_ptr v, const char * name, int d, unsigned int iter, unsigned int itemsondisk);
extern void matmul_top_save_vector_generic(matmul_top_data_ptr mmt, mmt_vec_ptr v, const char * name, int d, unsigned int iter, unsigned int itemsondisk);

/* These two do not really belong here, but come as a useful complement */
extern void vec_init_generic(pi_wiring_ptr, mpfq_vbase_ptr, mmt_vec_ptr, int, unsigned int);
extern void vec_clear_generic(pi_wiring_ptr, mmt_vec_ptr, unsigned int);

/* we should refrain from exposing these. At least for mksol,
 * allreduce_across is really useful, though.
 */
// extern void broadcast_down(matmul_top_data_ptr mmt, int d);
// extern void reduce_across(matmul_top_data_ptr mmt, int d);
extern void allreduce_across(matmul_top_data_ptr mmt, int d);
// extern void apply_permutation(matmul_top_data_ptr mmt, int d);
// extern void apply_identity(matmul_top_data_ptr mmt, int d);


#ifdef __cplusplus
}
#endif

#endif	/* MATMUL_TOP_H_ */
