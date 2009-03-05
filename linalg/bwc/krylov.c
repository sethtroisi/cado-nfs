#define _POSIX_C_SOURCE 200112L
#define _GNU_SOURCE

#include <stdio.h>

#include "parallelizing_info.h"
#include "matmul_top.h"
#include "abase.h"
#include "select_mpi.h"
#include "manu.h"
#include "debug.h"
#include "params.h"
#include "xvectors.h"
#include "xymats.h"
#include "bw-common-mpi.h"
#include "async.h"
// #include "rusage.h"
#include "filenames.h"
#include "xdotprod.h"

abobj_t abase;

void * krylov_prog(parallelizing_info_ptr pi, void * arg MAYBE_UNUSED)
{
    int tcan_print = bw->can_print && pi->m->trank == 0;
    matmul_top_data mmt;
    struct timing_data timing[1];

    int flags[2];
    flags[bw->dir] = THREAD_SHARED_VECTOR;
    flags[!bw->dir] = 0;

    block_control_signals();

    matmul_top_init(mmt, abase, pi, flags, bw->matrix_filename);

    mmt_wiring_ptr mcol = mmt->wr[bw->dir];
    mmt_wiring_ptr mrow = mmt->wr[!bw->dir];

    serialize(pi->m);
    
    abzero(abase, mrow->v->v, mrow->i1 - mrow->i0);

    uint32_t * gxvecs = malloc(bw->nx * bw->m * sizeof(uint32_t));

    load_x(gxvecs, bw->m, bw->nx, pi);

    char * v_name;
    int rc = asprintf(&v_name, V_FILE_BASE_PATTERN, bw->ys[0], bw->ys[1]);

    if (tcan_print) { printf("Loading %s...", v_name); fflush(stdout); }
    matmul_top_load_vector(mmt, v_name, bw->dir, bw->start);
    if (tcan_print) { printf("done\n"); }

    // We must also fetch the check vector C.
    abvobj_t abase_check;
    abvobj_set_nbys(abase_check, NCHECKS_CHECK_VECTOR);
    size_t vstride = abvbytes(abase_check, 1);
    size_t stride =  abbytes(abase, 1);

    mmt_generic_vec check_vector;
    matmul_top_vec_init_generic(mmt, abvbytes(abase_check, 1),
            check_vector, !bw->dir, THREAD_SHARED_VECTOR);
    if (tcan_print) { printf("Loading check vector..."); fflush(stdout); }
    matmul_top_load_vector_generic(mmt, vstride,
            check_vector, CHECK_FILE_BASE, !bw->dir, bw->interval);
    if (tcan_print) { printf("done\n"); }

    mmt_vec ahead;
    vec_init_generic(pi->m, stride, (mmt_generic_vec_ptr) ahead, 0, NCHECKS_CHECK_VECTOR);

    /* We'll store all xy matrices locally before doing reductions. Given
     * the small footprint of these matrices, it's rather innocuous.
     */
    mmt_vec xymats;
    if (tcan_print) {
        printf("Each thread allocates %zu kb for the A matrices\n",
                abbytes(abase, bw->m*bw->interval) >> 10);
    }
    vec_init_generic(pi->m, stride, (mmt_generic_vec_ptr) xymats, 0, bw->m*bw->interval);
    
    timing_init(timing, bw->start);

    if (bw->end == 0) {
        /* Decide on an automatic ending value */
        unsigned int length;
        length = MAX(mmt->n[0], mmt->n[1]);
        length = iceildiv(length, bw->m) + iceildiv(length, bw->n);
        length += 2 * iceildiv(bw->m, bw->n);
        length += 2 * iceildiv(bw->n, bw->m);
        length += 10;
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

    for(int s = bw->start ; s < bw->end ; s += bw->interval ) {
        // Plan ahead. The check vector is here to predict the final A matrix.
        // Note that our share of the dot product is determined by the
        // intersections of the i0..i1 intervals on both sides.
        
        unsigned int how_many;
        unsigned int offset_c;
        unsigned int offset_v;
        how_many = intersect_two_intervals(&offset_c, &offset_v,
                mrow->i0, mrow->i1,
                mcol->i0, mcol->i1);

        abzero(abase, ahead->v, NCHECKS_CHECK_VECTOR);
        if (how_many) {
            size_t bytes_c =  abvbytes(abase_check, offset_c);
            abvt * c = abase_generic_ptr_add(check_vector->v, bytes_c);
            abt * v = mcol->v->v + aboffset(abase, offset_v);
            abvdotprod(abase, abase_check, ahead->v, c, v, how_many);
        }

        /*
        debug_write(ahead->v, NCHECKS_CHECK_VECTOR * stride, "ahead.%u.j%u.t%u",
                s, mmt->pi->m->jrank, mmt->pi->m->trank);
         */

        abzero(abase, xymats->v, bw->m*bw->interval);

        serialize(pi->m);
        for(int i = 0 ; i < bw->interval ; i++) {
            /* Compute the product by x */
            x_dotprod(mmt, gxvecs,
                    xymats->v + aboffset(abase, i * bw->m), bw->m);
            matmul_top_mul(mmt, bw->dir);

            /*
            debug_write(mcol->v->v,
                    (mcol->i1-mcol->i0) * stride, "pV%u.j%u.t%u",
                    s+i+1, mmt->pi->m->jrank, mmt->pi->m->trank);
             */
            timing_check(pi, timing, s+i+1, tcan_print);
        }
        serialize(pi->m);

        /*
        debug_write(xymats->v, (bw->m * bw->interval) * stride, "xy.%u-%u.j%u.t%u",
                s, s + bw->interval, mmt->pi->m->jrank, mmt->pi->m->trank);
         */

        /* Last dot product. This must cancel ! */
        x_dotprod(mmt, gxvecs, ahead->v, NCHECKS_CHECK_VECTOR);

        /*
        debug_write(ahead->v, NCHECKS_CHECK_VECTOR * stride, "post%u.j%u.t%u",
                s, mmt->pi->m->jrank, mmt->pi->m->trank);
         */

        allreduce_generic(abase, ahead, pi->m, NCHECKS_CHECK_VECTOR);
        if (!abis_zero(abase, ahead->v, NCHECKS_CHECK_VECTOR)) {
            printf("Failed check at iteration %d\n", s + bw->interval);
            exit(1);
        }

        /* Now (and only now) collect the xy matrices */
        allreduce_generic(abase, xymats, pi->m, bw->m * bw->interval);

        if (pi->m->trank == 0 && pi->m->jrank == 0) {
            char * tmp;
            rc = asprintf(&tmp, A_FILE_PATTERN, bw->ys[0], bw->ys[1], s, s+bw->interval);
            FILE * f = fopen(tmp, "w");
            rc = fwrite(xymats->v, stride, bw->m*bw->interval, f);
            if (rc != bw->m*bw->interval) {
                fprintf(stderr, "Ayee -- short write\n");
                // make sure our input data won't be deleted -- this
                // chunk will have to be redone later, maybe the disk
                // failure is temporary (?)
            }
            fclose(f);
            free(tmp);
        }

        matmul_top_save_vector(mmt, v_name, bw->dir, s + bw->interval);

        serialize(pi->m);

        // reached s + bw->interval. Count our time on cpu, and compute the sum.
        timing_disp_collective_oneline(pi, timing, s + bw->interval, tcan_print);
    }

    if (tcan_print) {
        printf("Done.\n");
    }
    serialize(pi->m);

    matmul_top_vec_clear_generic(mmt, abvbytes(abase_check, 1), check_vector, !bw->dir);
    vec_clear_generic(pi->m, stride, (mmt_generic_vec_ptr) xymats, bw->m*bw->interval);
    vec_clear_generic(pi->m, stride, (mmt_generic_vec_ptr) ahead, NCHECKS_CHECK_VECTOR);
    free(gxvecs);
    free(v_name);

    matmul_top_clear(mmt, abase);

    timing_clear(timing);

    return NULL;
}


void usage()
{
    fprintf(stderr, "Usage: ./krylov <options>\n");
    fprintf(stderr, bw_common_usage_string());
    fprintf(stderr, "Relevant options here: wdir cfg m n mpi thr matrix interval start ys\n");
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
    if (bw->ys[0] < 0) { fprintf(stderr, "no ys value set\n"); exit(1); }

    abobj_init(abase);
    abobj_set_nbys(abase, bw->ys[1]-bw->ys[0]);

    catch_control_signals();
    pi_go(krylov_prog, bw->mpi_split[0], bw->mpi_split[1], bw->thr_split[0], bw->thr_split[1], 0);

    bw_common_clear_mpi(bw);
    return 0;
}

