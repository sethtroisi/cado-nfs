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

// the ``all_v'' field collects all the pointers to the per-thread vector
// values. There are exactly pi->wr[0]->ncores such pointers in
// mmt->wr[0]. In some cases, these pointers may be equal (for source
// vectors, never used as destination), and in some cases not.

struct mmt_vec_s;
typedef struct mmt_vec_s * mmt_vec_ptr;

struct mmt_vec_s {
    mpfq_vbase_ptr abase;
    parallelizing_info_ptr pi;
    int d;
    pi_datatype_ptr pitype;
    void * v;
    mmt_vec_ptr * siblings;     /* pi->wr[d]->ncores siblings ;
                                   only in case all cores in the communicator
                                   have their own data area v */
    unsigned int n;             // total size in items
    unsigned int i0;
    unsigned int i1;
    size_t rsbuf_size;          // auto-expanded on demand.
    void * rsbuf[2];            // only for RS_CHOICE == RS_CHOICE_MINE
    int consistency;    /* 0 == inconsistent ; 1 == partial ; 2 == full */
};
typedef struct mmt_vec_s mmt_vec[1];

/* some handy macros to access portions of mmt_vec's */
#define SUBVEC(v,w,offset) (v)->abase->vec_subvec(v->abase, (v)->w, offset)
#define SUBVEC_const(v,w,offset) (v)->abase->vec_subvec_const(v->abase, (v)->w, offset)

/* yikes. yet another list structure. Wish we had the STL. */
struct permutation_data {
    size_t n;
    size_t alloc;
    unsigned int (*x)[2];
};
typedef struct permutation_data permutation_data[1];
typedef struct permutation_data * permutation_data_ptr;
/* all methods are private */

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
        parallelizing_info_ptr pi,
        param_list pl,
        int optimized_direction);


void matmul_top_decl_usage(param_list_ptr pl);
void matmul_top_lookup_parameters(param_list_ptr pl);
extern void matmul_top_clear(matmul_top_data_ptr mmt);
#if 0
extern void matmul_top_fill_random_source(matmul_top_data_ptr mmt, int d);
#endif
extern size_t mmt_my_own_size_in_items(mmt_vec_ptr v);
extern size_t mmt_my_own_size_in_bytes(mmt_vec_ptr v);
extern size_t mmt_my_own_offset_in_items(mmt_vec_ptr v);
extern size_t mmt_my_own_offset_in_bytes(mmt_vec_ptr v);
extern void * mmt_my_own_subvec(mmt_vec_ptr v);

extern void mmt_vec_set_random_through_file(mmt_vec_ptr v, const char * name, unsigned int iter, unsigned int itemsondisk, gmp_randstate_t rstate);
/* do not use this function if you want consistency when the splitting
 * changes ! */
extern void mmt_vec_set_random_inconsistent(mmt_vec_ptr v, gmp_randstate_t rstate);
extern void mmt_vec_set_x_indices(mmt_vec_ptr y, uint32_t * gxvecs, int m, unsigned int nx);

extern void matmul_top_mul_cpu(matmul_top_data_ptr mmt, mmt_vec_ptr w, mmt_vec_ptr v);
extern void matmul_top_comm_bench(matmul_top_data_ptr mmt, int d);
extern void matmul_top_mul_comm(mmt_vec_ptr w, mmt_vec_ptr v);

/* v is both input and output. w is temporary */
static inline void matmul_top_mul(matmul_top_data_ptr mmt, mmt_vec_ptr w, mmt_vec_ptr v)
{
    ASSERT_ALWAYS(v->consistency == 2);
    matmul_top_mul_cpu(mmt, w, v);
    ASSERT_ALWAYS(w->consistency == 0);
    matmul_top_mul_comm(v, w);
    ASSERT_ALWAYS(v->consistency == 2);
}

/* Now some of the generic interface calls. By design, not everything is
 * possible with these calls. In particular, nothing critical is doable.
 * As a general convention, the only thing these calls need to know about
 * the abase is the stride value, which corresponds to sizeof(abelt)
 * Besides that, the generic calls take the vector pointer as argument.
 * Specifying NULL as vector argument is equivalent to taking the default
 * vector defined in the mmt data for the given direction flag. */


extern void mmt_vec_init(matmul_top_data_ptr mmt, mpfq_vbase_ptr abase, pi_datatype_ptr pitype, mmt_vec_ptr v, int d, int flags, unsigned int n);
extern void mmt_vec_clear(matmul_top_data_ptr mmt, mmt_vec_ptr v);
extern void mmt_own_vec_set(mmt_vec_ptr w, mmt_vec_ptr v);
extern void mmt_own_vec_set2(mmt_vec_ptr z, mmt_vec_ptr w, mmt_vec_ptr v);
extern void mmt_full_vec_set(mmt_vec_ptr w, mmt_vec_ptr v);
extern void mmt_full_vec_set_zero(mmt_vec_ptr v);
#if 0
extern void matmul_top_fill_random_source_generic(matmul_top_data_ptr mmt, size_t stride, mmt_vec_ptr v, int d);
#endif
extern void mmt_vec_load(mmt_vec_ptr v, const char * name, unsigned int iter, unsigned int itemsondisk);
extern void mmt_vec_save(mmt_vec_ptr v, const char * name, unsigned int iter, unsigned int itemsondisk);

extern void mmt_vec_broadcast(mmt_vec_ptr v);
extern void mmt_vec_reduce(mmt_vec_ptr w, mmt_vec_ptr v);
extern void mmt_vec_reduce_sameside(mmt_vec_ptr v);
extern void mmt_vec_allreduce(mmt_vec_ptr v);
extern void mmt_vec_twist(matmul_top_data_ptr mmt, mmt_vec_ptr y);
extern void mmt_vec_untwist(matmul_top_data_ptr mmt, mmt_vec_ptr y);
extern void mmt_vec_apply_T(matmul_top_data_ptr mmt, mmt_vec_ptr y);
extern void mmt_vec_unapply_T(matmul_top_data_ptr mmt, mmt_vec_ptr y);
extern void mmt_vec_apply_S(matmul_top_data_ptr mmt, mmt_vec_ptr y);
extern void mmt_vec_unapply_S(matmul_top_data_ptr mmt, mmt_vec_ptr y);
extern void mmt_vec_apply_P(matmul_top_data_ptr mmt, mmt_vec_ptr y);
extern void mmt_vec_unapply_P(matmul_top_data_ptr mmt, mmt_vec_ptr y);
extern void mmt_apply_identity(mmt_vec_ptr w, mmt_vec_ptr v);
extern void indices_twist(matmul_top_data_ptr mmt, uint32_t * xs, unsigned int n, int d);

#ifdef __cplusplus
}
#endif

#endif	/* MATMUL_TOP_H_ */
