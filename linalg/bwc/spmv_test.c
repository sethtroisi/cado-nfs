#include "cado.h"
#include <stdio.h>
#include "bwc_config.h"
#include "parallelizing_info.h"
#include "matmul_top.h"
#include "select_mpi.h"
#include "params.h"
#include "xvectors.h"
#include "bw-common-mpi.h"
#include "filenames.h"
#include "mpfq/mpfq.h"
#include "mpfq/abase_vbase.h"
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

    int withcoeffs = param_list_lookup_string(pl, "prime") != NULL;
    int nchecks = withcoeffs ? NCHECKS_CHECK_VECTOR_GFp : NCHECKS_CHECK_VECTOR_GF2;

    mpz_t p;
    mpz_init_set_ui(p, 2);
    param_list_parse_mpz(pl, "prime", p);
    abase_vbase A;
    abase_vbase_oo_field_init_byfeatures(A, 
            MPFQ_PRIME_MPZ, p,
            MPFQ_GROUPSIZE, nchecks,
            MPFQ_DONE);
    mpz_clear(p);


    matmul_top_init(mmt, A, pi, flags, pl, bw->dir);
    unsigned int unpadded = MAX(mmt->n0[0], mmt->n0[1]);

    mmt_wiring_ptr mcol = mmt->wr[bw->dir];
    mmt_wiring_ptr mrow = mmt->wr[!bw->dir];
    // pi_wiring_ptr picol = mmt->pi->wr[bw->dir];
    // pi_wiring_ptr pirow = mmt->pi->wr[!bw->dir];

    A->vec_set_zero(A, mrow->v->v, mrow->i1 - mrow->i0);
    A->vec_set_zero(A, mcol->v->v, mcol->i1 - mcol->i0);

    serialize(pi->m);

    /* saves data to wr[bw->dir]->v == mcol->v) */

    int d=bw->dir;

    if  (tcan_print)
        printf("bw->dir==%d d==%d\n", bw->dir, d);
    // known strategies:
    
    matmul_top_load_vector(mmt, "Y", d, 0, unpadded);
    matmul_top_twist_vector(mmt, d);
    matmul_top_mul(mmt, d);
    matmul_top_untwist_vector(mmt, d);
    matmul_top_save_vector(mmt, "MY", d, 0, unpadded);

#if 0
    matmul_top_apply_P(mmt, d);
    matmul_top_save_vector(mmt, "xA", d, 0);

    matmul_top_load_vector(mmt, "Y", d, 0);
    matmul_top_unapply_P(mmt, d);
    matmul_top_save_vector(mmt, "xB", d, 0);

    matmul_top_load_vector(mmt, "Y", d, 0);
    matmul_top_apply_P_apply_S(mmt, d);
    matmul_top_save_vector(mmt, "xC", d, 0);

    matmul_top_load_vector(mmt, "Y", d, 0);
    matmul_top_unapply_S_unapply_P(mmt, d);
    matmul_top_save_vector(mmt, "xD", d, 0);

    matmul_top_load_vector(mmt, "Y", d, 0);
    matmul_top_apply_S(mmt, d);
    matmul_top_save_vector(mmt, "xE", d, 0);

    matmul_top_load_vector(mmt, "Y", d, 0);
    matmul_top_unapply_S(mmt, d);
    matmul_top_save_vector(mmt, "xF", d, 0);
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
    matmul_top_load_vector(mmt, "U", d, 0);
    matmul_top_twist_vector(mmt, d);
    matmul_top_save_vector(mmt, "V", d, 0);
#endif

#if 0
    // untwist
    // rr, m->fw_colperm[rr], => zv - Sr^-1*Pr^-1*zu;
    // other direction: Sr*Pr*zu
    // bw->dir==1 d==0 Sr^-1*Pr^-1*zu;
    // bw->dir==1 d==1 Sr*Pr*zu;
    matmul_top_load_vector(mmt, "U", d, 0);
    matmul_top_untwist_vector(mmt, d);
    matmul_top_save_vector(mmt, "V", d, 0);
#endif
    // apply_identity(mmt, d);
    // matmul_top_save_vector(mmt, "V", !d, 0);

    serialize(pi->m);

#if 0
    mmt_vec oldv[2];
    memcpy(oldv[0], mmt->wr[0]->v, sizeof(mmt_vec));
    memcpy(oldv[1], mmt->wr[1]->v, sizeof(mmt_vec));
    memset(mmt->wr[0]->v, 0, sizeof(mmt_vec));
    memset(mmt->wr[1]->v, 0, sizeof(mmt_vec));
    matmul_top_vec_init(mmt, 0, flags[1 ^ 0]);
    matmul_top_vec_init(mmt, 1, flags[1 ^ 1]);

    serialize_threads(mmt->pi->m);

    size_t stride = abbytes(mmt->abase, 1);
    ASSERT_ALWAYS((mrow->v->flags & THREAD_SHARED_VECTOR) == 0);
    if (pirow->trank == 0 && pirow->jrank == 0) {
        abase_generic_copy(stride, mmt->wr[!bw->dir]->v->v, oldv[!bw->dir]->v, mmt->wr[!bw->dir]->i1 - mmt->wr[!bw->dir]->i0);
    }
    serialize_threads(mmt->pi->m);

    matmul_top_mul_comm(mmt, bw->dir);
    matmul_top_save_vector(mmt, CHECK_FILE_BASE, bw->dir, bw->interval);

    matmul_top_vec_clear(mmt, 0);
    matmul_top_vec_clear(mmt, 1);
    memcpy(mmt->wr[0]->v, oldv[0], sizeof(mmt_vec));
    memcpy(mmt->wr[1]->v, oldv[1], sizeof(mmt_vec));
    memset(oldv[0], 0, sizeof(mmt_vec));
    memset(oldv[1], 0, sizeof(mmt_vec));
#endif
    matmul_top_clear(mmt);

    return NULL;
}


void usage()
{
    exit(1);
}

int main(int argc, char * argv[])
{
    param_list pl;
    param_list_init(pl);
    bw_common_init_mpi(bw, pl, &argc, &argv);
    if (param_list_warn_unused(pl)) usage();

    setvbuf(stdout,NULL,_IONBF,0);
    setvbuf(stderr,NULL,_IONBF,0);

    pi_go(tst_prog, pl, 0);

    param_list_clear(pl);
    bw_common_clear_mpi(bw);
    return 0;
}

