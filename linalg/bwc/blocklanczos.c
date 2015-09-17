#include "cado.h"

#include <stdio.h>
#include "bwc_config.h"
#include "parallelizing_info.h"
#include "matmul_top.h"
#include "select_mpi.h"
#include "params.h"
#include "xvectors.h"
#include "portability.h"
#include "misc.h"
#include "bw-common.h"
#include "async.h"
#include "filenames.h"
#include "xdotprod.h"
#include "rolling.h"
#include "matops.h"
#include "mpfq/mpfq.h"
#include "mpfq/mpfq_vbase.h"

struct blstate {
    matmul_top_data mmt;
    mpfq_vbase A;
    pi_datatype_ptr A_pi;
    gmp_randstate_t rstate;

    /* We'll need several intermediary n*n matrices. These will be
     * allocated everywhere (we set flags to 0, so as to avoid having
     * shared vectors) */
    mmt_vec V[3];
    void * L[3];
    bit_vector D[3];

    /* Here are the semantics of the data fields above.
     *
     * For iteration n, we let n0 = n%3, n1 = (n-1)%3, n2 = (n-2)%3.
     *
     * V[n0] is the V vector which is the input to iteration n. It does
     *       not consist of independent  vectors.
     * D[n0] is the extracted set of columns, computed from V[n0] (and
     *       also from D[n1]. Together, V[n0] and D[n0] allow to compute
     *       the basis of the n-th sub vector space W, although no
     *       explicit data is reserved to its storage.
     * L[n0] is computed from V[n0] and D[n0], and is some sort of local
     *       inverse (iirc, it's W_i^{inv} in most accounts).
     *
     * of course *[n1] and *[n2] are the same for the previous steps.
     *
     * In order to compute the vector for step n + 1, the data to be used
     * is V[*], L[*], and D[n0, n1]
     */
};

/* given a 64 by 64 symmetric matrix A (given as 64 uint64_t's), and an input
 * bitmap S given in the form of one uint64_t, compute B and T so that
 * the following hold:
 *
 *    (i)    B = T * B = B * T
 *    (ii)   B * A * T = T
 *    (iii)  rank(B) = rank(A)
 *    (iv)   (1-S)*T = 1-S
 *
 * where S is identified with the diagonal matrix whose entries are given by S.
 *
 * In other words, this inverts (ii) a maximal minor (iii) in the matrix A,
 * selecting in priority (iv) for defining this minor the row indices
 * which are zero in S.
 *
 * The matrix A is clobbered.
 *
 * This routine does some allocation on the stack. Speed is not critical.
 */
uint64_t extraction_step(uint64_t * B, uint64_t * A, uint64_t S)
{
    int order[64], reorder[64];
    uint64_t B0[64];
    uint64_t T = 0;
    /* convert to a priority list, in a "save trees" style.  */
    for(int o=64,z=0,i=64;i-->0;) order[S&(UINT64_C(1)<<i)?--o:z++]=i;
    for(int i = 0 ; i < 64 ; i++) B0[i] = UINT64_C(1)<<i;
    for(int i = 0 ; i < 64 ; i++) reorder[i]=-1;
    for(int j = 0 ; j < 64 ; j++) {
        int oj = order[j];
        uint64_t mj = UINT64_C(1)<<oj;
        int p = -1;
        for(int i = 0 ; i < 64 ; i++) {
            int oi = order[i];
            uint64_t mi = UINT64_C(1)<<oi;
            if (T & mi) continue;
            if (A[oi] & mj) {
                p = i;
                break;
            }
        }
        if (p < 0) continue;

        int op = order[p];
        /* Of course it's important to use indices op and oj here ! */
        reorder[op] = oj;
        uint64_t mp = UINT64_C(1) << op;
        /* We have a pivot, great. */
        ASSERT_ALWAYS(!(T & mp));
        T |= mp;
        /* add row op to all rows except itself */
        for(int i = 0 ; i < 64 ; i++) {
            if (i == p) continue;
            int oi = order[i];
            uint64_t x = ~-!(A[oi] & mj);
            B0[oi] ^= B0[op] & x;
            A[oi] ^= A[op] & x;
        }
    }
    /* Now at this point, we have some more work to do.
     *
     * A*T is now almost the identity matrix -- its square is diagonal
     * with zeros and ones, so A*T is just an involution. The reorder[]
     * array just has this involution.
     *
     * B is such that B*original_A = new_A.
     *
     * The matrix we want to return is new_A*T*b*T. We use reorder[] to
     * actually copy this into A.
     */
    memset(B, 0, sizeof(B0));
    for(int i = 0 ; i < 64 ; i++) {
        if (reorder[i] >= 0)
            B[reorder[i]]=B0[i]&T;
    }
    return T;
}


void blstate_init(struct blstate * bl, parallelizing_info_ptr pi, param_list_ptr pl)
{
    matmul_top_data_ptr mmt = bl->mmt;
    mpfq_vbase_ptr A = bl->A;
    /* Note that THREAD_SHARED_VECTOR can't work in a block Lanczos
     * context, since both ways are used for input and output.
     *
     * Well, at least it does not seem to be as easy as it is in the
     * block Wiedemann context. Would be neat to find a way, though,
     * since otherwise this represents quite a significant memory
     * footprint in the end.
     */
    int flags[2] = {0,0};

    gmp_randinit_default(bl->rstate);

    mpfq_vbase_oo_field_init_byfeatures(A,
            MPFQ_PRIME_MPZ, bw->p,
            MPFQ_GROUPSIZE, bw->ys[1]-bw->ys[0],
            MPFQ_DONE);

    bl->A_pi = pi_alloc_mpfq_datatype(pi, A);

    matmul_top_init(mmt, A, bl->A_pi, pi, flags, pl, bw->dir);

    /* So far we have constrained A to have width exactly n, but it may
     * also conceivably be a divisor of n. */
    size_t nelts_for_nnmat = bw->n * (bw->n / A->groupsize(A));

    for(int i = 0 ; i < 3 ; i++) {
        A->vec_init(A, &(bl->L[i]), nelts_for_nnmat);
        /* We also need D_n, D_{n-1}, D_{n-2}. Those are in fact bitmaps.
         * Not clear that the bitmap type is really the one we want, though. */
        bit_vector_init(bl->D[i], bw->n);
        /* We need as well the two previous vectors. For these, distributed
         * storage will be ok. */
        mmt_vec_init(mmt, NULL, NULL, bl->V[i], bw->dir, flags[bw->dir]);
    }

}

void blstate_clear(struct blstate * bl)
{
    matmul_top_data_ptr mmt = bl->mmt;
    mpfq_vbase_ptr A = bl->A;
    parallelizing_info_ptr pi = mmt->pi;
    size_t nelts_for_nnmat = bw->n * (bw->n / A->groupsize(A));

    serialize(mmt->pi->m);
    for(int i = 0 ; i < 3 ; i++) {
        A->vec_clear(A, &(bl->L[i]), nelts_for_nnmat);
        /* We also need D_n, D_{n-1}, D_{n-2}. Those are in fact bitmaps.
         * Not clear that the bitmap type is really the one we want, though. */
        bit_vector_clear(bl->D[i]);
        /* We need as well the two previous vectors. For these, distributed
         * storage will be ok. */
        mmt_vec_clear(mmt, bl->V[i], bw->dir);
    }
    pi_free_mpfq_datatype(pi, bl->A_pi);
    matmul_top_clear(bl->mmt);
    A->oo_field_clear(bl->A);
    gmp_randclear(bl->rstate);
}

void blstate_set_start(struct blstate * bl)
{
    matmul_top_data_ptr mmt = bl->mmt;
    mpfq_vbase_ptr A = bl->A;
    size_t nelts_for_nnmat = bw->n * (bw->n / A->groupsize(A));

    int g = A->groupsize(A);
    /* D = identity, L too, and V = 0 */
    for(int i = 0 ; i < 3 ; i++) {
        bit_vector_set(bl->D[i], 1);
        /* Set L[i] to identity -- it's awkward, yes. */
        A->vec_set_zero(A, bl->L[i], nelts_for_nnmat);
        for(int j = 0 ; j < bw->n ; j++) {
            int p = j * bw->n + j;
            A->set_ui_at(A, A->vec_coeff_ptr(A, bl->L[i], p / g), p % g, 1);
        }
    }
    /* matmul_top_vec_init has already set V to zero */
    /* for bw->dir=0, mmt->n0[0] is the number of rows. */
    mmt_vec_set_random_through_file(mmt, bl->V[0], "blstart", bw->dir, 0, mmt->n0[bw->dir], bl->rstate);
    mmt_vec_set(mmt, 0, bl->V[0], bw->dir);
}

void blstate_load(struct blstate * bl, unsigned int iter)
{
    unsigned int i0 = iter % 3;
    unsigned int i1 = (iter+3-1) % 3;
    unsigned int i2 = (iter+3-2) % 3;
    matmul_top_data_ptr mmt = bl->mmt;
    mpfq_vbase_ptr A = bl->A;
    parallelizing_info_ptr pi = mmt->pi;

    size_t nelts_for_nnmat = bw->n * (bw->n / A->groupsize(A));
    /* bw->dir=0: mmt->n0[bw->dir] = number of rows */
    /* bw->dir=1: mmt->n0[bw->dir] = number of columns */
    size_t sizeondisk = A->vec_elt_stride(A, mmt->n0[bw->dir]);
    size_t mysize = mmt_my_own_size_in_bytes(bl->mmt, NULL, bw->dir);

    char * filename;
    int rc = asprintf(&filename, "blstate.%u", iter);
    ASSERT_ALWAYS(rc >= 0);
    pi_file_handle f;
    int tcan_print = bw->can_print && pi->m->trank == 0;
    if (tcan_print) { printf("Loading %s...", filename); fflush(stdout); }

    pi_file_open(f, mmt->pi, bw->dir, filename, "rb");
    if (pi->m->jrank == 0 && pi->m->trank == 0) {
        bit_vector_read_from_stream(bl->D[i1], f->f);
        size_t rc;
        rc = fread(bl->L[i1], A->vec_elt_stride(A,1), nelts_for_nnmat, f->f);
        ASSERT_ALWAYS(rc == (size_t) nelts_for_nnmat);
        rc = fread(bl->L[i2], A->vec_elt_stride(A,1), nelts_for_nnmat, f->f);
        ASSERT_ALWAYS(rc == (size_t) nelts_for_nnmat);
    }
    ssize_t s;

    void * V0 = mmt_my_own_subvec(bl->mmt, bl->V[i0], bw->dir);
    s = pi_file_read(f, V0, mysize, sizeondisk);
    ASSERT_ALWAYS(s >= 0 && (size_t) s == sizeondisk);

    void * V1 = mmt_my_own_subvec(bl->mmt, bl->V[i1], bw->dir);
    s = pi_file_read(f, V1, mysize, sizeondisk);
    ASSERT_ALWAYS(s >= 0 && (size_t) s == sizeondisk);
 
    void * V2 = mmt_my_own_subvec(bl->mmt, bl->V[i2], bw->dir);
    s = pi_file_read(f, V2, mysize, sizeondisk);
    ASSERT_ALWAYS(s >= 0 && (size_t) s == sizeondisk);

    pi_file_close(f);
    if (tcan_print) { printf("done\n"); fflush(stdout); }
    free(filename);
}

void blstate_save(struct blstate * bl, unsigned int iter)
{
    unsigned int i0 = iter % 3;
    unsigned int i1 = (iter+3-1) % 3;
    unsigned int i2 = (iter+3-2) % 3;
    matmul_top_data_ptr mmt = bl->mmt;
    mpfq_vbase_ptr A = bl->A;
    parallelizing_info_ptr pi = mmt->pi;

    size_t nelts_for_nnmat = bw->n * (bw->n / A->groupsize(A));
    /* bw->dir=0: mmt->n0[bw->dir] = number of rows */
    /* bw->dir=1: mmt->n0[bw->dir] = number of columns */
    size_t sizeondisk = A->vec_elt_stride(A, mmt->n0[bw->dir]);
    size_t mysize = mmt_my_own_size_in_bytes(bl->mmt, NULL, bw->dir);

    char * filename;
    int rc = asprintf(&filename, "blstate.%u", iter);
    ASSERT_ALWAYS(rc >= 0);
    pi_file_handle f;
    int tcan_print = bw->can_print && pi->m->trank == 0;
    if (tcan_print) { printf("Saving %s...", filename); fflush(stdout); }

    pi_file_open(f, mmt->pi, bw->dir, filename, "wb");
    if (pi->m->jrank == 0 && pi->m->trank == 0) {
        bit_vector_write_to_stream(bl->D[i1], f->f);
        size_t rc;
        rc = fwrite(bl->L[i1], A->vec_elt_stride(A,1), nelts_for_nnmat, f->f);
        ASSERT_ALWAYS(rc == (size_t) nelts_for_nnmat);
    }
    ssize_t s;
    void * V0 = mmt_my_own_subvec(bl->mmt, bl->V[i0], bw->dir);
    s = pi_file_write(f, V0, mysize, sizeondisk);
    ASSERT_ALWAYS(s >= 0 && (size_t) s == sizeondisk);

    void * V1 = mmt_my_own_subvec(bl->mmt, bl->V[i1], bw->dir);
    s = pi_file_write(f, V1, mysize, sizeondisk);
    ASSERT_ALWAYS(s >= 0 && (size_t) s == sizeondisk);
 
    void * V2 = mmt_my_own_subvec(bl->mmt, bl->V[i2], bw->dir);
    s = pi_file_write(f, V2, mysize, sizeondisk);
    ASSERT_ALWAYS(s >= 0 && (size_t) s == sizeondisk);

    pi_file_close(f);
    if (tcan_print) { printf("done\n"); fflush(stdout); }
    free(filename);
}

/* Given a *BINARY* vector block of *EXACTLY* 64 vectors (that is, we
 * force bw->n == 0), compute a 64*64 matrix such that m * v is in row
 * reduced echelon form. Also return the rank.
 *
 * The vector block has to be according to direction d. If v is NULL, the
 * built-in vector block of the mmt structure in that direction is taken.
 *
 * Note that this is a collective operation.
 *
 * Vector v is modified by this process (and put in RREF).
 */
int mmt_vec_echelon(matmul_top_data_ptr mmt, mat64_ptr m, mmt_vec_ptr v0, int d)
{
    memset(m, 0, sizeof(mat64));
    uint64_t * v = mmt_my_own_subvec(mmt, v0, d);
    size_t eblock = mmt_my_own_size_in_items(mmt, d);
    /* This is the total number of non-zero coordinates of the vector v */
    size_t n = mmt->n0[!d];
    /* In all what follows, we'll talk about v being a 64*n matrix, with
     * [v[i]&1] being "the first row", and so on.  */
    uint64_t usedrows = 0;
    int rank = 0;
    for(int i = 0 ; i < 64 ; i++) {
        uint64_t mi = UINT64_C(1) << i;
        /* Find the earliest column which has non-zero in the first row */
        unsigned int j;
        for(j = 0 ; j < eblock ; j++) {
            if (v[j] & mi) break;
        }
        if (j == eblock) j = n;
        else j += mmt_my_own_offset_in_items(mmt, d);
        unsigned int jmin;
        pi_allreduce(&jmin, &j, 1, BWC_PI_UNSIGNED, BWC_PI_MIN, mmt->pi->m);
        if (jmin == n) {
            /* zero row */
            continue;
        }
        usedrows |= mi;
        rank++;
        /* zero out the other coefficients of this column. This means
         * that we need to collect the control data from the thread which
         * owns this column, and then act accordingly */
        uint64_t control = 0;
        if (jmin == j) {
            control = v[j - mmt_my_own_offset_in_items(mmt, d)];
            ASSERT_ALWAYS(control & mi);
            control ^= mi;
        }
        /* TODO: once we require mpi-3.0, use MPI_UINT64_T instead */
        ASSERT_ALWAYS(sizeof(unsigned long long) == sizeof(uint64_t));
        pi_allreduce(NULL, &control, 1, BWC_PI_UNSIGNED_LONG, BWC_PI_MAX, mmt->pi->m);
        for(unsigned int k = 0 ; k < eblock ; k++) {
            v[k] ^= control & -!!(v[k] & mi);
        }
    }
    /* TODO: compute m... */
    /* TODO: realign, to put in RREF */
    return rank;
}



/* This saves first the vector V, and then the vector V*A. By
 * construction, the two are orthogonal. It is the responsibility of the
 * external program to check exactly where we have kernel vectors in V,
 * using the comparison with V*A.
 */
void blstate_save_result(struct blstate * bl, unsigned int iter)
{
    unsigned int i0 = iter % 3;
    matmul_top_data_ptr mmt = bl->mmt;
    mpfq_vbase_ptr A = bl->A;
    parallelizing_info_ptr pi = mmt->pi;

    /* bw->dir=0: mmt->n0[bw->dir] = number of rows */
    /* bw->dir=1: mmt->n0[bw->dir] = number of columns */

    char * filename;
    int rc = asprintf(&filename, "blsolution.bin"); // %u", iter);
    ASSERT_ALWAYS(rc >= 0);
    pi_file_handle f;
    int tcan_print = bw->can_print && pi->m->trank == 0;
    if (tcan_print) { printf("Saving %s...", filename); fflush(stdout); }

    pi_file_open(f, mmt->pi, bw->dir, filename, "wb");

    {
        size_t sizeondisk = A->vec_elt_stride(A, mmt->n0[bw->dir]);
        size_t mysize = mmt_my_own_size_in_bytes(bl->mmt, NULL, bw->dir);
        void * myvec = mmt_my_own_subvec(bl->mmt, bl->V[i0], bw->dir);

        ssize_t s = pi_file_write(f, myvec, mysize, sizeondisk);
        ASSERT_ALWAYS(s >= 0 && (size_t) s == sizeondisk);
    }

    /* Save V*M as well. Note that it's a bit peculiar, because we will
     * need to save something in the other direction.
     */
    {
        mmt_vec_set(mmt, 0, bl->V[i0], bw->dir);
        matmul_top_twist_vector(mmt, bw->dir);
        matmul_top_mul_cpu(mmt, bw->dir);
        allreduce_across(mmt, NULL, !bw->dir);
        matmul_top_untwist_vector(mmt, !bw->dir);

        mat64 m;
        int r = mmt_vec_echelon(mmt, m, NULL, bw->dir);
        printf("rank(V) == %d\n", r);
        r = mmt_vec_echelon(mmt, m, NULL, !bw->dir);
        printf("rank(V*M) == %d\n", r);

        size_t sizeondisk = A->vec_elt_stride(A, mmt->n0[!bw->dir]);
        size_t mysize = mmt_my_own_size_in_bytes(bl->mmt, NULL, !bw->dir);
        void * myvec = mmt_my_own_subvec(bl->mmt, NULL, !bw->dir);

        ssize_t s = pi_file_write(f, myvec, mysize, sizeondisk);
        ASSERT_ALWAYS(s >= 0 && (size_t) s == sizeondisk);
    }

    pi_file_close(f);
    if (tcan_print) { printf("done\n"); fflush(stdout); }
    free(filename);
}


void * bl_prog(parallelizing_info_ptr pi, param_list pl, void * arg MAYBE_UNUSED)
{
    int tcan_print = bw->can_print && pi->m->trank == 0;
    struct timing_data timing[1];

    struct blstate bl[1];
    /* shortcuts */
    matmul_top_data_ptr mmt = bl->mmt;
    mpfq_vbase_ptr A = bl->A;

    /* so many features we do not support ! */
    ASSERT_ALWAYS(bw->m == bw->n);
    ASSERT_ALWAYS(bw->ys[0] == 0);
    ASSERT_ALWAYS(bw->ys[1] == bw->n);

    /* Don't think we stand any chance with interleaving with block
     * Lanczos... */
    ASSERT_ALWAYS(!pi->interleaved);

    // int withcoeffs = mpz_cmp_ui(bw->p, 2) > 0;
    block_control_signals();

    blstate_init(bl, pi, pl);

    mpfq_vbase_tmpl AxA;
    mpfq_vbase_oo_init_templates(AxA, A, A);

    size_t nelts_for_nnmat = bw->n * (bw->n / bl->A->groupsize(bl->A));

    mmt_comm_ptr mcol = mmt->wr[bw->dir];
    mmt_comm_ptr mrow = mmt->wr[!bw->dir];

    ASSERT_ALWAYS(mmt->wr[bw->dir]->i0 == 0);
    ASSERT_ALWAYS(mmt->wr[!bw->dir]->i0 == 0);
    /* hmm, right. So here we're in effect saying that we do not support
     * threading yet...
     */
    ASSERT_ALWAYS(mmt->wr[bw->dir]->i1 == mmt->n[bw->dir]);
    ASSERT_ALWAYS(mmt->wr[!bw->dir]->i1 == mmt->n[!bw->dir]);

    serialize(pi->m);

    /* note that we don't know how to do checking for BL. Sure, we can
     * check for mutual orthogonality of vector blocks at the
     * checkpoints, but that is rather lame.
     */

    matmul_top_comm_bench(mmt, bw->dir);

    if (bw->start == 0) {
        blstate_set_start(bl);
        blstate_save(bl, 0);
    } else {
        blstate_load(bl, bw->start);
    }

    if (bw->end == 0) {
        /* Decide on an automatic ending value */
        unsigned int length;
        /* FIXME: is it really mmt->n[bw->dir] ? Looks like it should be
         * mmt->n[!bw->dir] instead (number of columns, for bw->dir==1).
         */
        length = mmt->n[bw->dir] / (bw->n - 0.76) + 10;
        /* allow some deviation */
        length += 2*integer_sqrt(length);
        /* Because bw is a global variable, we protect its use */
        if (serialize_threads(pi->m)) {
            bw->end = length;
        }
        serialize(pi->m);
    }

    if (tcan_print) {
        fprintf(stderr, "Target iteration is %u ; going to %u\n", bw->end,
                bw->interval * iceildiv(bw->end, bw->interval));
    }


    timing_init(timing, bw->start, bw->interval * iceildiv(bw->end, bw->interval));

    void * vav;
    void * vaav;

    A->vec_init(A, &vav, nelts_for_nnmat);
    A->vec_init(A, &vaav, nelts_for_nnmat);

    /* TODO: Put that in the state file. */
    int sum_Ni = 0;

    for(int s = bw->start ; s < bw->end ; s += bw->interval ) {
        /* XXX For BL, how much do we care about twisting ? */
        /* The BW iteration, combined with shuffled_product option, is
         * such that for a stored matrix M', we do the iteration with
         * P*M'. So in order to have an equivalent of multiplying by M,
         * we want P*M' = S^-1 * M * S, i.e. M' = (SP)^-1 M S.
         *
         * Now if in fact, we are doing the iteration followed by its
         * transpose, this may end up somewhat different.
         *
         * But it's incorrect to think like this. In the BL iteration,
         * shuffled_product makes no sense at all. All communications are
         * done with allreduce, so there is no inherent shuffling in the
         * parallelized matrix multiplication. So when we have a stored
         * matrix M', we do multiply by M'. Now let us say that we have
         * run the balancing process just like it would have been for BW.
         * So we have M' = (SP)^-1 M S. If we multiply by M'^T
         * afterwards, we end up multiplying by S^-1 * M^T * M * S. So
         * twisting a vector means computing S^-1*V -- that should be
         * essentially the same as for bwc.
         *
         * [description above is for solving M*w=0, while for BL we are
         * really in the factoring case where we care about w*M=0. So
         * then, twisting would mean v --> v*S*P, which is what
         * matmul_top_twist_vector does for bw->dir==0 and nullspace=left]
         *
         * Bottom line: it's harmless to do the balancing the very same
         * way we do with BW.
         */

        for(int k = 0 ; k < 3 ; k++) {
            mmt_vec_set(mmt, 0, bl->V[k], bw->dir);
            matmul_top_twist_vector(mmt, bw->dir);
            mmt_vec_set(mmt, bl->V[k], 0, bw->dir);
        }
        mmt_vec_set(mmt, 0, bl->V[s % 3], bw->dir);

        /* I *think* that this is necessary */
        broadcast_down(mmt, NULL, bw->dir);

        /* FIXME: for BL, we use 8 timers while only 4 are defined in the
         * structure.
         */
        serialize(pi->m);
        int i, i0, i1, i2;

        for(i = 0 ; i < bw->interval ; i++) {
            i0 = (s + i) % 3;
            i1 = (s + i + 3 - 1) % 3;
            i2 = (s + i + 3 - 2) % 3;

            A->vec_set_zero(A, mrow->v->v, mrow->i1 - mrow->i0);

            for(int d = 0 ; d < 2 ; d++) {
                /* The first part of this loop must be guaranteed to be free
                 * of any mpi calls */
                /* at this point, timer is [0] (CPU) */

                {
                    /* for d==0: in:mcol->v, out:mrow->v. */
                    /* for d==1: in:mrow->v, out:mcol->v. */
                    matmul_top_mul_cpu(mmt, bw->dir ^ d);

                    timing_next_timer(timing);  /* now timer is [1] (cpu-wait) */
                }
                serialize(pi->m);           /* for measuring waits only */

                timing_next_timer(timing);  /* now timer is [2] (COMM) */
                allreduce_across(mmt, NULL, 1 ^ d ^ bw->dir);

                timing_next_timer(timing);  /* now timer is [3] (comm-wait) */
                serialize(pi->m);           /* for measuring waits only */

                timing_next_timer(timing);  /* now timer is [0] (CPU) */
            }

            /* We have now:
             *  V_{n}           bl->V[0]
             *  V_{n-1}         bl->V[1]
             *  V_{n-2}         bl->V[2]
             *  A * V_n         mcol->v
             */

            AxA->dotprod(A->obj, A->obj, vav,
                    mmt_my_own_subvec(bl->mmt, bl->V[i0], bw->dir), mcol->v->v,
                    mmt_my_own_size_in_items(bl->mmt, bw->dir));

            pi_allreduce(NULL, vav,
                    nelts_for_nnmat, bl->A_pi, BWC_PI_BXOR, pi->m);

            AxA->dotprod(A->obj, A->obj, vaav,
                    mcol->v->v, mcol->v->v,
                    mmt_my_own_size_in_items(bl->mmt, bw->dir));

            pi_allreduce(NULL, vaav, nelts_for_nnmat,
                    bl->A_pi, BWC_PI_BXOR, pi->m);


            ASSERT_ALWAYS(bl->D[i0]->n == 64);



            {
                pi_comm_ptr picol = mmt->pi->wr[bw->dir];
                ASSERT_ALWAYS((mcol->i1 - mcol->i0) % picol->totalsize == 0);
                size_t eblock = (mcol->i1 - mcol->i0) /  picol->totalsize;
                ASSERT_ALWAYS(mmt->abase->vec_elt_stride(mmt->abase, 1) == sizeof(uint64_t));

                // Here are the operations we will now perform
                //
                // X := (IN/*-Dn0*/)*Vn + Dn0*VA;
                // C0 := ((IN/*-Dn0*/)*vav + Dn0*vaav) * Ln0 * Vn;
                // C1 := Dn0 * vav * Ln1 * Vn1;
                // C2 := Dn0 * vav * (IN-Dn1) * Ln2 * Vn2;
                //
                // all are performed locally.

                // we set up X in mcol->v

                uint64_t * V0 = mmt_my_own_subvec(bl->mmt, bl->V[i0], bw->dir);
                uint64_t * V1 = mmt_my_own_subvec(bl->mmt, bl->V[i1], bw->dir);
                uint64_t * V2 = mmt_my_own_subvec(bl->mmt, bl->V[i2], bw->dir);
                uint64_t * VA  = mmt_my_own_subvec(bl->mmt, NULL, bw->dir);
                uint64_t * X   = mmt_my_own_subvec(bl->mmt, NULL, bw->dir);
                uint64_t D0;
                uint64_t D1 = bl->D[i1]->p[0];
                mat64_ptr mvav = (mat64_ptr) vav;
                mat64_ptr mvaav = (mat64_ptr) vaav;
                mat64_ptr mL0 = (mat64_ptr) bl->L[i0];
                mat64_ptr mL1 = (mat64_ptr) bl->L[i1];
                mat64_ptr mL2 = (mat64_ptr) bl->L[i2];
                mat64 m0, m1, m2, t;

                /* We need to save vav for use a wee bit later in this loop. */
                memcpy(t, mvav, sizeof(mat64));

                D0 = bl->D[i0]->p[0] = extraction_step(mL0, t, D1);

                int Ni = bit_vector_popcount(bl->D[i0]);
                sum_Ni += Ni;
                // int defect = bw->n - Ni;
                //if (tcan_print) printf("step %d, dim=%d\n", s+i, Ni);
                if (Ni == 0) {
                    break;
                }

            
                for(size_t i = 0 ; i < eblock ; i++) {
                    X[i] = V0[i] ^ (VA[i] & D0);
                }

                for(int i = 0 ; i < 64 ; i++) t[i] = mvav[i] ^ (mvaav[i] & D0);
                mul_6464_6464(m0, mL0, t);
                for(int i = 0 ; i < 64 ; i++) mvav[i] &= D0;
                mul_6464_6464(m1, mL1, mvav);
                for(int i = 0 ; i < 64 ; i++) mL2[i] &= ~D1;
                mul_6464_6464(m2, mL2, mvav);

                addmul_N64_6464(X, V0, m0, eblock);
                addmul_N64_6464(X, V1, m1, eblock);
                addmul_N64_6464(X, V2, m2, eblock);

                /* So. We need to get ready for the next step, don't we ? */
                mmt_vec_set(mmt, bl->V[i2], 0, bw->dir);
            }
            timing_check(pi, timing, s+i+1, tcan_print);
        }
        serialize(pi->m);

        for(int k = 0 ; k < 3 ; k++) {
            mmt_vec_set(mmt, 0, bl->V[k], bw->dir);
            matmul_top_untwist_vector(mmt, bw->dir);
            mmt_vec_set(mmt, bl->V[k], 0, bw->dir);
        }

        serialize(pi->m);
        if (i < bw->interval) {
            if (tcan_print) printf("Finished at iteration %d\n", s+i);

            blstate_save_result(bl, s+i);
            /* We need to cheat somewhat */
            bw->end = s + i;
            timing->end_mark = s + i;
            break;
        }


        blstate_save(bl, s+bw->interval);
        serialize(pi->m);

        if (tcan_print) printf("N=%d ; sum_dim=%d=N*(%u-%.3f)\n", s+i, sum_Ni, bw->n, bw->n-(double)sum_Ni/(s+i));
        // reached s + bw->interval. Count our time on cpu, and compute the sum.
        timing_disp_collective_oneline(pi, timing, s + bw->interval, mmt->mm->ncoeffs, tcan_print, "blocklanczos");
    }

    // we can't do as we do with BW, since the process really stops in
    // the end, continuing the last block to the checkpoint mark does not
    // make sense.
    timing_final_tally(pi, timing, mmt->mm->ncoeffs, tcan_print, "blocklanczos");

    A->vec_clear(A, &vav, nelts_for_nnmat);
    A->vec_clear(A, &vaav, nelts_for_nnmat);

#if 0
    pi_log_clear(pi->m);
    pi_log_clear(pi->wr[0]);
    pi_log_clear(pi->wr[1]);
#endif

    if (tcan_print) {
        printf("Done blocklanczos.\n");
    }
    serialize(pi->m);

    blstate_clear(bl);

    timing_clear(timing);

    return NULL;
}

/* The unit test for extraction_step, which is a static function here,
 * actually includes the source file. We just want to make sure we don't
 * expose main()
 */
#ifndef BL_TESTING
int main(int argc, char * argv[])
{
    param_list pl;

    bw_common_init_new(bw, &argc, &argv);
    param_list_init(pl);

    bw_common_decl_usage(pl);
    parallelizing_info_decl_usage(pl);
    matmul_top_decl_usage(pl);
    /* declare local parameters and switches: none here (so far). */

    bw_common_parse_cmdline(bw, pl, &argc, &argv);

    bw_common_interpret_parameters(bw, pl);
    parallelizing_info_lookup_parameters(pl);
    matmul_top_lookup_parameters(pl);
    /* interpret our parameters */
    if (bw->ys[0] < 0) { fprintf(stderr, "no ys value set\n"); exit(1); }
    catch_control_signals();

    if (param_list_warn_unused(pl)) {
        param_list_print_usage(pl, bw->original_argv[0], stderr);
        exit(EXIT_FAILURE);
    }

    pi_go(bl_prog, pl, 0);

    param_list_clear(pl);
    bw_common_clear_new(bw);

    return 0;
}
#endif
