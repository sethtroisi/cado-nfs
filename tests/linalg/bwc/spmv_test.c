#include "cado.h"
#include <stdio.h>
#include "bwc_config.h"
#include "parallelizing_info.h"
#include "matmul_top.h"
#include "select_mpi.h"
#include "params.h"
#include "xvectors.h"
#include "bw-common.h"
#include "filenames.h"
#include "mpfq/mpfq.h"
#include "mpfq/mpfq_vbase.h"
#include "portability.h"

/*
extern void broadcast_down(matmul_top_data_ptr mmt, int d);
extern void reduce_across(matmul_top_data_ptr mmt, int d);
extern void allreduce_across(matmul_top_data_ptr mmt, int d);
extern void apply_permutation(matmul_top_data_ptr mmt, int d);
extern void apply_identity(matmul_top_data_ptr mmt, int d);
*/

/* This only does a multiplication */

void * tst_prog(parallelizing_info_ptr pi, param_list pl, void * arg MAYBE_UNUSED)
{
    /* Interleaving does not make sense for this program. So the second
     * block of threads just leave immediately */
    if (pi->interleaved && pi->interleaved->idx)
        return NULL;

    int tcan_print = bw->can_print && pi->m->trank == 0;
    matmul_top_data mmt;

    int flags[2];
    flags[bw->dir] = THREAD_SHARED_VECTOR;
    flags[!bw->dir] = 0;

    int withcoeffs = mpz_cmp_ui(bw->p, 2) > 0;
    int nchecks = withcoeffs ? NCHECKS_CHECK_VECTOR_GFp : NCHECKS_CHECK_VECTOR_GF2;
    mpfq_vbase A;
    mpfq_vbase_oo_field_init_byfeatures(A, 
            MPFQ_PRIME_MPZ, bw->p,
            MPFQ_GROUPSIZE, nchecks,
            MPFQ_DONE);

    pi_datatype_ptr A_pi = pi_alloc_mpfq_datatype(pi, A);

    matmul_top_init(mmt, A, A_pi, pi, flags, pl, bw->dir);
    unsigned int unpadded = MAX(mmt->n0[0], mmt->n0[1]);

    mmt_vec_set_zero(mmt, NULL, 0);
    mmt_vec_set_zero(mmt, NULL, 1);

    serialize(pi->m);

    /* saves data to wr[bw->dir]->v == mcol->v) */

    int d=bw->dir;

    if  (tcan_print)
        printf("bw->dir==%d d==%d\n", bw->dir, d);
    // known strategies:
    
    mmt_vec_load(mmt, NULL, "Y", d, 0, unpadded);
    matmul_top_twist_vector(mmt, d);
    matmul_top_mul(mmt, d);
    matmul_top_untwist_vector(mmt, d);
    mmt_vec_save(mmt, NULL, "MY", d, 0, unpadded);

#if 0
    matmul_top_apply_P(mmt, d);
    mmt_vec_save(mmt, NULL, "xA", d, 0);

    mmt_vec_load(mmt, NULL, "Y", d, 0);
    matmul_top_unapply_P(mmt, d);
    mmt_vec_save(mmt, NULL, "xB", d, 0);

    mmt_vec_load(mmt, NULL, "Y", d, 0);
    matmul_top_apply_P_apply_S(mmt, d);
    mmt_vec_save(mmt, NULL, "xC", d, 0);

    mmt_vec_load(mmt, NULL, "Y", d, 0);
    matmul_top_unapply_S_unapply_P(mmt, d);
    mmt_vec_save(mmt, NULL, "xD", d, 0);

    mmt_vec_load(mmt, NULL, "Y", d, 0);
    matmul_top_apply_S(mmt, d);
    mmt_vec_save(mmt, NULL, "xE", d, 0);

    mmt_vec_load(mmt, NULL, "Y", d, 0);
    matmul_top_unapply_S(mmt, d);
    mmt_vec_save(mmt, NULL, "xF", d, 0);
#endif

#if 0
    // load to mcol->v, apply the permutation from the balancing.
    // Communicate, save mcol->v
    // rr, r => zv eq Pr*Sr^-1*Pr^-1*zu
    // rr, m->fw_colperm[rr], => zv eq Pr*Sr*zu
    //
    // other direction: zv - Pr^-1*Sr^-1*zu
    // twist
    // bw->dir==1 d==1 Pr^-1*Sr^-1*zu
    // bw->dir==1 d==0 Pr*Sr*zu
    mmt_vec_load(mmt, NULL, "U", d, 0);
    matmul_top_twist_vector(mmt, d);
    mmt_vec_save(mmt, NULL, "V", d, 0);
#endif

#if 0
    // untwist
    // rr, m->fw_colperm[rr], => zv - Sr^-1*Pr^-1*zu;
    // other direction: Sr*Pr*zu
    // bw->dir==1 d==0 Sr^-1*Pr^-1*zu;
    // bw->dir==1 d==1 Sr*Pr*zu;
    mmt_vec_load(mmt, NULL, "U", d, 0);
    matmul_top_untwist_vector(mmt, d);
    mmt_vec_save(mmt, NULL, "V", d, 0);
#endif
    // apply_identity(mmt, d);
    // mmt_vec_save(mmt, NULL, "V", !d, 0);

    serialize(pi->m);

#if 0
    mmt_vec oldv[2];
    memcpy(oldv[0], mmt->wr[0]->v, sizeof(mmt_vec));
    memcpy(oldv[1], mmt->wr[1]->v, sizeof(mmt_vec));
    memset(mmt->wr[0]->v, 0, sizeof(mmt_vec));
    memset(mmt->wr[1]->v, 0, sizeof(mmt_vec));
    mmt_vec_init(mmt, NULL, NULL, NULL, 0, flags[1 ^ 0]);
    mmt_vec_init(mmt, NULL, NULL, NULL, 1, flags[1 ^ 1]);

    serialize_threads(mmt->pi->m);

    size_t stride = abbytes(mmt->abase, 1);
    ASSERT_ALWAYS((mrow->v->flags & THREAD_SHARED_VECTOR) == 0);
    if (pirow->trank == 0 && pirow->jrank == 0) {
        mpfq_generic_copy(stride, mmt->wr[!bw->dir]->v->v, oldv[!bw->dir]->v, mmt->wr[!bw->dir]->i1 - mmt->wr[!bw->dir]->i0);
    }
    serialize_threads(mmt->pi->m);

    matmul_top_mul_comm(mmt, bw->dir);
    mmt_vec_save(mmt, NULL, CHECK_FILE_BASE, bw->dir, bw->interval);

    mmt_vec_clear(mmt, NULL, 0);
    mmt_vec_clear(mmt, NULL, 1);
    memcpy(mmt->wr[0]->v, oldv[0], sizeof(mmt_vec));
    memcpy(mmt->wr[1]->v, oldv[1], sizeof(mmt_vec));
    memset(oldv[0], 0, sizeof(mmt_vec));
    memset(oldv[1], 0, sizeof(mmt_vec));
#endif
    matmul_top_clear(mmt);
    pi_free_mpfq_datatype(pi, A_pi);
    A->oo_field_clear(A);

    return NULL;
}


void usage()
{
    exit(1);
}

int main(int argc, char * argv[])
{
    param_list pl;

    bw_common_init_new(bw, &argc, &argv);
    param_list_init(pl);
    parallelizing_info_init();

    bw_common_decl_usage(pl);
    parallelizing_info_decl_usage(pl);
    matmul_top_decl_usage(pl);
    /* declare local parameters and switches: none here (so far). */

    bw_common_parse_cmdline(bw, pl, &argc, &argv);

    bw_common_interpret_parameters(bw, pl);
    parallelizing_info_lookup_parameters(pl);
    matmul_top_lookup_parameters(pl);

    if (param_list_warn_unused(pl)) {
        param_list_print_usage(pl, bw->original_argv[0], stderr);
        exit(EXIT_FAILURE);
    }

    pi_go(tst_prog, pl, 0);

    parallelizing_info_finish();
    param_list_clear(pl);
    bw_common_clear_new(bw);

    return 0;
}

