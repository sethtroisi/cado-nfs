#define _POSIX_C_SOURCE 200112L
#define _GNU_SOURCE

#include <stdio.h>

#include "parallelizing_info.h"
#include "matmul_top.h"
#include "abase.h"
#include "select_mpi.h"
#include "gauss.h"
//#include "hexstring.h"
#include "manu.h"

#include "params.h"
#include "xvectors.h"
#include "bw-common.h"

abobj_t abase;

void * sec_prog(parallelizing_info_ptr pi, void * arg MAYBE_UNUSED)
{
    /* OK. Now we merely have to build up a check vector. We've got a
     * fairly easy candidate for that: the x vector.
     */
    int tcan_print = can_print && pi->m->trank == 0;
    matmul_top_data mmt;

    int flags[2];
    flags[!dir] = THREAD_SHARED_VECTOR;
    flags[dir] = 0;

    matmul_top_init(mmt, abase, pi, flags, matrix_filename);

    abzero(abase, mmt->wr[!dir]->v->v, mmt->wr[!dir]->i1 - mmt->wr[!dir]->i0);

    uint32_t * gxvecs = malloc(nx * m * sizeof(uint32_t));

    load_x(gxvecs, m, nx, pi);

    ASSERT_ALWAYS(m >= NCHECKS_CHECK_VECTOR);
    for(int j = 0 ; j < NCHECKS_CHECK_VECTOR ; j++) {
        for(unsigned int k = 0 ; k < nx ; k++) {
            // set bit j of entry gxvecs[j*nx+k] to 1.
            uint32_t i = gxvecs[j*nx+k];
            if (i < mmt->wr[!dir]->i0 || i >= mmt->wr[!dir]->i1)
                continue;
            abt * where;
            where = mmt->wr[!dir]->v->v + aboffset(abase, i-mmt->wr[!dir]->i0);
            abset_ui(abase, where, j, 1);
        }
    }

    if (tcan_print) {
        printf("Computing trsp(x)*M^%d\n",interval);
    }

    serialize(pi->m);

    for(int k = 0 ; k < interval ; k++) {
        matmul_top_mul(mmt, !dir);
        if (tcan_print) {
            putchar('.');
            fflush(stdout);
        }
    }
    if (tcan_print) {
        printf("\n");
    }
    serialize(pi->m);

    matmul_top_save_vector(mmt, "C", !dir, interval);

    serialize(pi->m);
    matmul_top_clear(mmt, abase);
    serialize(pi->m);

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
    bw_common_init_mpi(argc, argv);

    if (nx == 0) { fprintf(stderr, "no nx value set\n"); exit(1); } 

    abobj_init(abase);
    abobj_set_nbys(abase, NCHECKS_CHECK_VECTOR);

    pi_go(sec_prog, mpi_split[0], mpi_split[1], thr_split[0], thr_split[1], 0);

    bw_common_clear();
    return 0;
}

