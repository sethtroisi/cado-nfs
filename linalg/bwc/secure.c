#define _POSIX_C_SOURCE 200112L
#define _GNU_SOURCE

#include <stdio.h>

#include "parallelizing_info.h"
#include "matmul_top.h"
#include "abase.h"
#include "select_mpi.h"
#include "manu.h"

#include "params.h"
#include "xvectors.h"
#include "bw-common-mpi.h"
#include "filenames.h"

abobj_t abase;
struct bw_params bw[1];

void * sec_prog(parallelizing_info_ptr pi, void * arg MAYBE_UNUSED)
{
    /* OK. Now we merely have to build up a check vector. We've got a
     * fairly easy candidate for that: the x vector.
     */
    int tcan_print = bw->can_print && pi->m->trank == 0;
    matmul_top_data mmt;

    /* XXX Here we're roking in the opposite direction compared to
     * prep/krylov/mksol ! */
    int flags[2];
    flags[!bw->dir] = THREAD_SHARED_VECTOR;
    flags[bw->dir] = 0;

    matmul_top_init(mmt, abase, pi, flags, bw->matrix_filename);

    abzero(abase, mmt->wr[!bw->dir]->v->v, mmt->wr[!bw->dir]->i1 - mmt->wr[!bw->dir]->i0);

    uint32_t * gxvecs = malloc(bw->nx * bw->m * sizeof(uint32_t));

    load_x(gxvecs, bw->m, bw->nx, pi);

    ASSERT_ALWAYS(bw->m >= NCHECKS_CHECK_VECTOR);
    for(int j = 0 ; j < NCHECKS_CHECK_VECTOR ; j++) {
        for(unsigned int k = 0 ; k < bw->nx ; k++) {
            // set bit j of entry gxvecs[j*bw->nx+k] to 1.
            uint32_t i = gxvecs[j*bw->nx+k];
            if (i < mmt->wr[!bw->dir]->i0 || i >= mmt->wr[!bw->dir]->i1)
                continue;
            abt * where;
            where = mmt->wr[!bw->dir]->v->v + aboffset(abase, i-mmt->wr[!bw->dir]->i0);
            abset_ui(abase, where, j, 1);
        }
    }

    if (tcan_print) {
        printf("Computing trsp(x)*bw->M^%d\n",bw->interval);
    }

    serialize(pi->m);

    for(int k = 0 ; k < bw->interval ; k++) {
        matmul_top_mul(mmt, !bw->dir);
        if (tcan_print) {
            putchar('.');
            fflush(stdout);
        }
    }
    if (tcan_print) {
        printf("\n");
    }
    serialize(pi->m);

    matmul_top_save_vector(mmt, CHECK_FILE_BASE, !bw->dir, bw->interval);

    matmul_top_clear(mmt, abase);

    free(gxvecs);

    return NULL;
}


void usage()
{
    fprintf(stderr, "Usage: ./secure <options>\n");
    fprintf(stderr, bw_common_usage_string());
    fprintf(stderr, "Relevant options here: wdir cfg m n mpi thr matrix interval\n");
    fprintf(stderr, "Note: data files must be found in wdir !\n");
    exit(1);
}

int main(int argc, char * argv[])
{
    param_list pl;
    param_list_init(pl);
    bw_common_init_mpi(bw, pl, argc, argv);
    if (param_list_warn_unused(pl)) usage();
    param_list_clear(pl);

    if (bw->nx == 0) { fprintf(stderr, "no nx value set\n"); exit(1); } 

    abobj_init(abase);
    abobj_set_nbys(abase, NCHECKS_CHECK_VECTOR);

    pi_go(sec_prog, bw->mpi_split[0], bw->mpi_split[1], bw->thr_split[0], bw->thr_split[1], 0);

    bw_common_clear_mpi(bw);
    return 0;
}

