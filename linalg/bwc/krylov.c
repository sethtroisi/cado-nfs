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
#include "bw-common.h"
#include "async.h"
// #include "rusage.h"
#include "filenames.h"

abobj_t abase;
struct bw_params bw[1];

void x_dotprod(matmul_top_data_ptr mmt, uint32_t * xv, abt * v)
{
    /* We're reading from the shared right vector data -- this area is
     * written to by the other threads in the column. Some of them might
     * be lingering in reduce operations, so we have to wait for them
     */
    if (mmt->wr[bw->dir]->v->flags & THREAD_SHARED_VECTOR) {
        serialize_threads(mmt->pi->wr[bw->dir]);
    } else {
        /* I presume that no locking is needed here. But it's unchecked
         */
        BUG();
    }

    for(int j = 0 ; j < bw->m ; j++) {
        abt * where = v + aboffset(abase, j);
        for(unsigned int t = 0 ; t < bw->nx ; t++) {
            uint32_t i = xv[j*bw->nx+t];
            unsigned int vi0 = mmt->wr[bw->dir]->i0;
            unsigned int vi1 = mmt->wr[bw->dir]->i1;
            unsigned int hi0 = mmt->wr[!bw->dir]->i0;
            unsigned int hi1 = mmt->wr[!bw->dir]->i1;
            if (i < vi0 || i >= vi1)
                continue;
            /* We want the data to match our interval on both
             * directions, because otherwise we'll end up
             * computing rubbish -- recall that no broadcast_down
             * has occurred yet.
             */
            if (i < hi0 || i >= hi1)
                continue;
            abadd(mmt->abase, where,
                    mmt->wr[bw->dir]->v->v + aboffset(abase, i - vi0));
        }
    }
}

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

    serialize(pi->m);
    
    abzero(abase, mmt->wr[!bw->dir]->v->v, mmt->wr[!bw->dir]->i1 - mmt->wr[!bw->dir]->i0);

    uint32_t * gxvecs = malloc(bw->nx * bw->m * sizeof(uint32_t));

    load_x(gxvecs, bw->m, bw->nx, pi);

    char * v_name;
    int rc = asprintf(&v_name, V_FILE_BASE_PATTERN, bw->ys[0], bw->ys[1]);

    if (tcan_print) { printf("Loading %s...", v_name); fflush(stdout); }
    matmul_top_load_vector(mmt, v_name, bw->dir, bw->start);
    if (tcan_print) { printf("done\n"); }

    // We must also fetch the check vector C.
    abvobj_t placeholder;
    abvobj_set_nbys(placeholder, NCHECKS_CHECK_VECTOR);
    size_t vstride = abvbytes(placeholder, 1);
    size_t stride =  abbytes(abase, 1);

    mmt_generic_vec check_vector;
    matmul_top_vec_init_generic(mmt, abvbytes(placeholder, 1),
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
        fprintf(stderr, "Target iteration is %u\n", bw->end);
    }

    for(int s = bw->start ; s < bw->end ; s += bw->interval ) {
        // Plan ahead. The check vector is here to predict the final A matrix.
        // Note that our share of the dot product is determined by the
        // intersections of the i0..i1 intervals on both sides.
        
        unsigned int how_many;
        unsigned int offset_c;
        unsigned int offset_v;
        how_many = intersect_two_intervals(&offset_c, &offset_v,
                mmt->wr[!bw->dir]->i0, mmt->wr[!bw->dir]->i1,
                mmt->wr[bw->dir]->i0, mmt->wr[bw->dir]->i1);

        abzero(abase, ahead->v, NCHECKS_CHECK_VECTOR);
        if (how_many) {
            size_t bytes_c =  abvbytes(placeholder, offset_c);
            abvt * c = abase_generic_ptr_add(check_vector->v, bytes_c);
            abt * v = mmt->wr[bw->dir]->v->v + aboffset(abase, offset_v);
            abvdotprod(abase, placeholder, ahead->v, c, v, how_many);
        }

        /*
        debug_write(abase, ahead->v, NCHECKS_CHECK_VECTOR, "ahead.%u.j%u.t%u",
                s, mmt->pi->m->jrank, mmt->pi->m->trank);
         */

        abzero(abase, xymats->v, bw->m*bw->interval);

        serialize(pi->m);
        for(int i = 0 ; i < bw->interval ; i++) {
            /* Compute the product by x */
            x_dotprod(mmt, gxvecs, xymats->v + aboffset(abase, i * bw->m));
            matmul_top_mul(mmt, bw->dir);

            /*
            debug_write(abase, mmt->wr[bw->dir]->v->v,
                    mmt->wr[bw->dir]->i1-mmt->wr[bw->dir]->i0, "pV%u.j%u.t%u",
                    s+i+1, mmt->pi->m->jrank, mmt->pi->m->trank);
             */
            timing_check(pi, timing, s+i+1, tcan_print);
        }
        serialize(pi->m);

        /*
        debug_write(abase, xymats->v, bw->m * bw->interval, "xy.%u-%u.j%u.t%u",
                s, s + bw->interval, mmt->pi->m->jrank, mmt->pi->m->trank);
         */

        /* Last dot product. This must cancel ! */
        x_dotprod(mmt, gxvecs, ahead->v);

        /*
        debug_write(abase, ahead->v, NCHECKS_CHECK_VECTOR, "post%u.j%u.t%u",
                s, mmt->pi->m->jrank, mmt->pi->m->trank);
         */

        reduce_generic(abase, ahead, pi->m, NCHECKS_CHECK_VECTOR);
        if (!abis_zero(abase, ahead->v, NCHECKS_CHECK_VECTOR)) {
            printf("Failed check at iteration %d\n", s + bw->interval);
            exit(1);
        }

        /* Now (and only now) collect the xy matrices */
        reduce_generic(abase, xymats, pi->m, bw->m * bw->interval);

        if (pi->m->trank == 0 && pi->m->jrank == 0) {
            char * tmp;
            rc = asprintf(&tmp, A_FILE_PATTERN, bw->ys[0], bw->ys[1], s, s+bw->interval);
            FILE * f = fopen(tmp, "w");
            rc = fwrite(xymats->v, stride, aboffset(abase, bw->m*bw->interval), f);
            if (rc != aboffset(abase, bw->m*bw->interval)) {
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
        printf("\n");
    }
    serialize(pi->m);

    matmul_top_vec_clear_generic(mmt, abvbytes(placeholder, 1), check_vector, !bw->dir);
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
    bw_common_init_mpi(bw, argc, argv);

    if (bw->nx == 0) { fprintf(stderr, "no nx value set\n"); exit(1); } 
    if (bw->ys[0] < 0) { fprintf(stderr, "no ys value set\n"); exit(1); }

    abobj_init(abase);
    abobj_set_nbys(abase, bw->ys[1]-bw->ys[0]);

    catch_control_signals();
    pi_go(krylov_prog, bw->mpi_split[0], bw->mpi_split[1], bw->thr_split[0], bw->thr_split[1], 0);

    bw_common_clear(bw);
    return 0;
}

