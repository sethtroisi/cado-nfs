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

/* A communicator is the information concerning our matrix in one of its
 * two dimensions. Everything is made so that inverting the direction
 * flag should make it possible to compute with the transposed matrix
 * easily.
 * 
 * [note: communicators, both for this layer and the parallelizing_info
 * layer, used to be called "wirings", hence the variable name "wr"]
 *
 * Communicator number 0 is ``horizontal''. In detail, what this means
 * is:
 *
 * n is the total number of rows ; the total number of instances of this
 * communicator, so to say.
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
 * Note also that we have the pi_comm to consider as well. Those
 * have nothing to do with the matrix, and are relevant with regard to
 * the parallelization of the work.
 */

// the ``all_v'' field collects all the pointers to the per-thread vector
// values. There are exactly pi->wr[0]->ncores such pointers in
// mmt->wr[0]. In some cases, these pointers may be equal (for source
// vectors, never used as destination), and in some cases not.

struct mmt_vec_s {
    mpfq_vbase_ptr abase;
    pi_datatype_ptr pitype;
    // size_t stride;      // shorcut to this->abase->vec_elt_stride(this->abase,1)
    void * v;
    void * * all_v;
    // internal use
    int flags;
};
typedef struct mmt_vec_s mmt_vec[1];
typedef struct mmt_vec_s * mmt_vec_ptr;

/* some handy macros to access portions of mmt_vec's */
#define SUBVEC(v,w,offset) (v)->abase->vec_subvec(v->abase, (v)->w, offset)
#define SUBVEC_const(v,w,offset) (v)->abase->vec_subvec_const(v->abase, (v)->w, offset)

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

struct mmt_comm_s {
    unsigned int i0;
    unsigned int i1;
    /* Note that while i0 and i1 are obviously abase-dependent, clearly
     * the two following fields are not ! */
    mmt_vec v;
    size_t rsbuf_size;          // auto-expanded on demand.
    void * rsbuf[2];            // only for USE_ALTERNATIVE_REDUCE_SCATTER
};
typedef struct mmt_comm_s mmt_comm[1];
typedef struct mmt_comm_s * mmt_comm_ptr;
typedef struct mmt_comm_s const * mmt_comm_srcptr;

/* yikes. yet another list structure. Wish we had the STL. */
struct permutation_data {
    uint32_t n;
    size_t alloc;
    uint32_t (*x)[2];
};
typedef struct permutation_data permutation_data[1];
typedef struct permutation_data * permutation_data_ptr;
/* all methods are private */

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
    //
    // The number of items in the vector mmt->wr[0]->v is the number of
    // rows, so it's mmt->n0[1].
    //
    unsigned int n[2];
    unsigned int n0[2]; // n0: unpadded.

    // this really ends up within the mm field.
    char * locfile;

    mmt_comm wr[2];

    // fences[0]: horizontal fences == row indices (horizontal separating line)
    // fences[1]: vertical fences == col indices (vertical separating line)
    // fences[0] contains nh = 1 + pi->wr[1]->totalsize values
    // fences[1] contains nv = 1 + pi->wr[0]->totalsize values
    unsigned int * fences[2];

    balancing bal;

    /* These are local excerpts of the balancing permutation: arrays of
     * pairs (row index in the (sub)-matrix) ==> (row index in the
     * original matrix), but only when both coordinates of this pair are
     * in the current row and column range. This can be viewed as the set
     * of non-zero positions in the permutation matrix if it were split
     * just like the current matrix is. */
    permutation_data_ptr perm[2];       /* rowperm, colperm */

    matmul_ptr mm;

    mpfq_vbase_ptr abase;
    pi_datatype_ptr pitype;
};

/* These are flags for the distributed vectors. For the moment we have
 * only one flag */
#define THREAD_SHARED_VECTOR    1

typedef struct matmul_top_data_s matmul_top_data[1];
typedef struct matmul_top_data_s * matmul_top_data_ptr;
typedef struct matmul_top_data_s const * matmul_top_data_srcptr;

#ifdef __cplusplus
extern "C" {
#endif

extern void matmul_top_init(matmul_top_data_ptr mmt,
        mpfq_vbase_ptr abase,
        pi_datatype_ptr pitype,
        parallelizing_info_ptr pi,
        int const * flags,
        param_list pl,
        int optimized_direction);


void matmul_top_decl_usage(param_list_ptr pl);
void matmul_top_lookup_parameters(param_list_ptr pl);
extern void matmul_top_clear(matmul_top_data_ptr mmt);
#if 0
extern void matmul_top_fill_random_source(matmul_top_data_ptr mmt, int d);
#endif
extern size_t mmt_my_own_size_in_items(matmul_top_data_ptr mmt, int d);
extern size_t mmt_my_own_size_in_bytes(matmul_top_data_ptr mmt, mmt_vec_ptr v, int d);
extern size_t mmt_my_own_offset_in_items(matmul_top_data_ptr mmt, int d);
extern size_t mmt_my_own_offset_in_bytes(matmul_top_data_ptr mmt, mmt_vec_ptr v, int d);
extern void * mmt_my_own_subvec(matmul_top_data_ptr mmt, mmt_vec_ptr v, int d);
extern void mmt_vec_set_random_through_file(matmul_top_data_ptr mmt, mmt_vec_ptr v, const char * name, int d, unsigned int iter, unsigned int itemsondisk, gmp_randstate_t rstate);
/* do not use this function if you want consistency when the splitting
 * changes ! */
extern void mmt_vec_set_random_inconsistent(matmul_top_data_ptr mmt, mmt_vec_ptr v, int d, gmp_randstate_t rstate);


extern void matmul_top_mul_cpu(matmul_top_data_ptr mmt, int d);
extern void matmul_top_comm_bench(matmul_top_data_ptr mmt, int d);
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
extern void matmul_top_zero_vec_area(matmul_top_data_ptr mmt, int d);
extern void matmul_top_apply_P_apply_S(matmul_top_data_ptr mmt, int d);
extern void matmul_top_unapply_S_unapply_P(matmul_top_data_ptr mmt, int d);
extern void matmul_top_apply_P(matmul_top_data_ptr mmt, int d);
extern void matmul_top_unapply_P(matmul_top_data_ptr mmt, int d);
extern void matmul_top_apply_S(matmul_top_data_ptr mmt, int d);
extern void matmul_top_unapply_S(matmul_top_data_ptr mmt, int d);
extern void matmul_top_apply_T(matmul_top_data_ptr mmt, int d);
extern void matmul_top_unapply_T(matmul_top_data_ptr mmt, int d);
extern void matmul_top_twist_vector(matmul_top_data_ptr mmt, int d);
extern void matmul_top_untwist_vector(matmul_top_data_ptr mmt, int d);
extern void indices_apply_S(matmul_top_data_ptr mmt, uint32_t * xs, unsigned int n, int d);
extern void indices_apply_P(matmul_top_data_ptr mmt, uint32_t * xs, unsigned int n, int d);

/* Now some of the generic interface calls. By design, not everything is
 * possible with these calls. In particular, nothing critical is doable.
 * As a general convention, the only thing these calls need to know about
 * the abase is the stride value, which corresponds to sizeof(abelt)
 * Besides that, the generic calls take the vector pointer as argument.
 * Specifying NULL as vector argument is equivalent to taking the default
 * vector defined in the mmt data for the given direction flag. */


extern void mmt_vec_init(matmul_top_data_ptr mmt, mpfq_vbase_ptr abase, pi_datatype_ptr, mmt_vec_ptr v, int d, int flags);
extern void mmt_vec_clear(matmul_top_data_ptr mmt, mmt_vec_ptr v, int d);
extern void mmt_vec_set(matmul_top_data_ptr mmt, mmt_vec_ptr w, mmt_vec_ptr v, int d);
extern void mmt_vec_set_zero(matmul_top_data_ptr mmt, mmt_vec_ptr v, int d);
#if 0
extern void matmul_top_fill_random_source_generic(matmul_top_data_ptr mmt, size_t stride, mmt_vec_ptr v, int d);
#endif
extern void mmt_vec_load(matmul_top_data_ptr mmt, mmt_vec_ptr v, const char * name, int d, unsigned int iter, unsigned int itemsondisk);
extern void mmt_vec_save(matmul_top_data_ptr mmt, mmt_vec_ptr v, const char * name, int d, unsigned int iter, unsigned int itemsondisk);

/* These two do not really belong here, but come as a useful complement */
extern void vec_init_generic(pi_comm_ptr, mpfq_vbase_ptr, pi_datatype_ptr, mmt_vec_ptr, int, unsigned int);
extern void vec_clear_generic(pi_comm_ptr, mmt_vec_ptr, unsigned int);

/* we should refrain from exposing these. At least for mksol,
 * allreduce_across is really useful, though. We use it for the block
 * Lanczos iterations, too.
 */
/* FIXME: well, for BL I fear that I might need broadcast_down */
extern void broadcast_down(matmul_top_data_ptr mmt, mmt_vec_ptr v, int d);
extern void reduce_across(matmul_top_data_ptr mmt, mmt_vec_ptr w, mmt_vec_ptr v, int d);
extern void allreduce_across(matmul_top_data_ptr mmt, mmt_vec_ptr v, int d);
// extern void apply_permutation(matmul_top_data_ptr mmt, int d);
// extern void apply_identity(matmul_top_data_ptr mmt, int d);


#ifdef __cplusplus
}
#endif

#endif	/* MATMUL_TOP_H_ */
