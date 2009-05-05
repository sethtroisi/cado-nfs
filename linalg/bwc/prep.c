#define _POSIX_C_SOURCE 200112L
#define _GNU_SOURCE

#include <stdio.h>
#include "bwc_config.h"
#include "parallelizing_info.h"
#include "matmul_top.h"
#include "abase.h"
#include "select_mpi.h"
#include "random_generation.h"
#include "gauss.h"
#include "gauss.h"
#include "params.h"
#include "xvectors.h"
#include "xymats.h"
#include "bw-common-mpi.h"
#include "filenames.h"

/* Number of copies of m by n matrices to use for trying to obtain a
 * full-rank matrix (rank m).
 *
 * Note that it must be at least m/n, otherwise we stand no chance !
 */
#define NBITER  2

abobj_t abase;

void * prep_prog(parallelizing_info_ptr pi, param_list pl, void * arg MAYBE_UNUSED)
{
    // Doing the ``hello world'' test is a very good way of testing the
    // global mpi/pthreads setup. So despite its apparent irrelevance, I
    // suggest leaving it here as a cheap sanity check.
    hello(pi);

    // avoid cluttering output too much.
    int tcan_print = bw->can_print && pi->m->trank == 0;

    unsigned int my_nx = 1;

    matmul_top_data mmt;

    int flags[2];
    flags[bw->dir] = THREAD_SHARED_VECTOR;
    flags[!bw->dir] = 0;

    matmul_top_init(mmt, abase, pi, flags, pl, MATRIX_INFO_FILE, bw->dir);

    mmt_wiring_ptr mcol = mmt->wr[bw->dir];
    mmt_wiring_ptr mrow = mmt->wr[!bw->dir];

    uint32_t * xvecs = malloc(my_nx * bw->m * sizeof(uint32_t));

    mmt_vec xymats;
    size_t stride =  abbytes(abase, 1);

    /* We're cheating on the generic init routines */
    vec_init_generic(pi->m, stride, (mmt_generic_vec_ptr) xymats, 0, bw->m*NBITER);

    for (unsigned ntri = 0;; ntri++) {
        serialize_threads(pi->m);

        if (tcan_print) {
            printf("// Generating new x,y vector pair (trial # %u)\n", ntri);
        }
        if (ntri >= my_nx * 10) {
            ++my_nx;
            if (tcan_print) {
                printf("// Getting bored. Trying %u x vectors\n", my_nx);
            }
            xvecs = realloc(xvecs, my_nx * bw->m * sizeof(uint32_t));
            ASSERT_ALWAYS(xvecs != NULL);
        }

        // if we're looking for the right nullspace, then x is on the left.
        // Otherwise, it's on the right.
        setup_x_random(xvecs, bw->m, my_nx, mmt->n[bw->dir], pi);

        // Compute y.
        matmul_top_fill_random_source(mmt, bw->dir);

        // we need to save this starting vector for later use if it turns out
        // that we need to save it for real.
        matmul_top_save_vector(mmt, Y_FILE_BASE, bw->dir, 0);

        // We must compute x^T M y, x^T M^2 y, and so on.
        // XXX Note that x^Ty does not count here, because it does nto
        // take part to the sequence computed by lingen !
        matmul_top_mul(mmt, bw->dir);
        
        // we have indices mmt->wr[1]->i0..i1 available.
        abzero(abase, xymats->v, bw->m * NBITER);

        for(unsigned int k = 0 ; k < NBITER ; k++) {
            for(int j = 0 ; j < bw->m ; j++) {
                abt * where = xymats->v + aboffset(abase, j * NBITER + k);
                for(unsigned int t = 0 ; t < my_nx ; t++) {
                    uint32_t i = xvecs[j*my_nx+t];
                    if (i < mcol->i0 || i >= mcol->i1)
                        continue;
                    /* We want the data to match our bw->interval on both
                     * directions, because otherwise we'll end up
                     * computing rubbish -- recall that no broadcast_down
                     * has occurred yet.
                     */
                    if (i < mrow->i0 || i >= mrow->i1)
                        continue;
                    abadd(abase, where,
                            mcol->v->v + aboffset(abase, i - mcol->i0));
                }
            }
            matmul_top_mul(mmt, bw->dir);
        }

        /* Make sure computation is over for everyone ! */
        serialize_threads(pi->m);

        /* Now all threads and jobs must collectively reduce the zone
         * pointed to by xymats */
        allreduce_generic(abase, xymats, pi->m, bw->m * NBITER);

        /* OK -- now everybody has the same data */

        int dimk;
        int * pdimk;
        
        /* the kernel() call is not reentrant */
        if (pi->m->trank == 0) {
            dimk = kernel((mp_limb_t *) xymats->v, NULL,
                    bw->m, NBITER * abnbits(abase), 
                abbytes(abase,NBITER) / sizeof(mp_limb_t), 0);
            pdimk = &dimk;
        }
        thread_agreement(pi->m, (void **) &pdimk, 0);
        dimk = * pdimk;

        if (tcan_print)
            printf("// Dimension of kernel: %d\n", dimk);

        if (dimk == 0) {
            if (tcan_print)
                printf("// Found good x,y vector pair after %u trials\n",
                        ntri+1);
            break;
        }
    }

    if (pi->m->trank == 0) {
        bw->nx = my_nx;
    }

    save_x(xvecs, bw->m, my_nx, pi);

    matmul_top_clear(mmt, abase);

    /* clean up xy mats stuff */
    vec_clear_generic(pi->m, stride, (mmt_generic_vec_ptr) xymats, bw->m*NBITER);

    free(xvecs);
    return NULL;
}

void usage()
{
    fprintf(stderr, "Usage: ./prep <options>\n");
    fprintf(stderr, bw_common_usage_string());
    fprintf(stderr, "Relevant options here: wdir cfg m n mpi thr matrix\n");
    exit(1);
}

int main(int argc, char * argv[])
{
    param_list pl;
    param_list_init(pl);
    bw_common_init_mpi(bw, pl, &argc, &argv);
    if (param_list_warn_unused(pl)) usage();

    if (bw->seed) setup_seeding(bw->seed);

    abobj_init(abase);
    abobj_set_nbys(abase, bw->n);

    pi_go(prep_prog, pl, 0);

    // we save the parameter list once again, because the prep program
    // generates some useful info, bw->nx in particular.
    param_list_save_parameter(pl, PARAMETER_FROM_FILE, "nx", "%u", bw->nx);
    param_list_save(pl, BW_CONFIG_FILE);
    param_list_clear(pl);

    bw_common_clear_mpi(bw);
    return 0;
}

