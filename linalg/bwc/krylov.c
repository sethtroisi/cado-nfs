#define _POSIX_C_SOURCE 200112L
#define _GNU_SOURCE

#include <stdio.h>

#include "parallelizing_info.h"
#include "matmul_top.h"
#include "abase.h"
#include "select_mpi.h"
#include "gauss.h"
#include "manu.h"
#include "debug.h"

#include "params.h"
#include "xvectors.h"
#include "xymats.h"

abobj_t abase;

#include "bw-common.h"

void x_dotprod(matmul_top_data_ptr mmt, abt * v, int m, uint32_t * xv, unsigned int nx)
{
    /* We're reading from the shared right vector data -- this area is
     * written to by the other threads in the column. Some of them might
     * be lingering in reduce operations, so we have to wait for them
     */
    serialize_threads(mmt->pi->wr[dir]);

    for(int j = 0 ; j < m ; j++) {
        abt * where = v + aboffset(abase, j);
        for(unsigned int t = 0 ; t < nx ; t++) {
            uint32_t i = xv[j*nx+t];
            if (i < mmt->wr[dir]->i0 || i >= mmt->wr[dir]->i1)
                continue;
            /* We want the data to match our interval on both
             * directions, because otherwise we'll end up
             * computing rubbish -- recall that no broadcast_down
             * has occurred yet.
             */
            if (i < mmt->wr[!dir]->i0 || i >= mmt->wr[!dir]->i1)
                continue;
            abadd(mmt->abase, where,
                    mmt->wr[dir]->v->v + aboffset(abase, i - mmt->wr[dir]->i0));
        }
    }
}

void * krylov_prog(parallelizing_info_ptr pi, void * arg MAYBE_UNUSED)
{
    int tcan_print = can_print && pi->m->trank == 0;
    matmul_top_data mmt;

    int flags[2];
    flags[dir] = THREAD_SHARED_VECTOR;
    flags[!dir] = 0;

    matmul_top_init(mmt, abase, pi, flags, matrix_filename);

    serialize(pi->m);
    
    abzero(abase, mmt->wr[!dir]->v->v, mmt->wr[!dir]->i1 - mmt->wr[!dir]->i0);

    uint32_t * gxvecs = malloc(nx * m * sizeof(uint32_t));

    load_x(gxvecs, m, nx, pi);

    char * v_name;
    asprintf(&v_name, "V%u-%u", ys[0], ys[1]);

    if (tcan_print) { printf("Loading %s...", v_name); fflush(stdout); }
    matmul_top_load_vector(mmt, v_name, dir, start);
    if (tcan_print) { printf("done\n"); }

    // We must also fetch the check vector C.
    abvobj_t placeholder;
    abvobj_set_nbys(placeholder, NCHECKS_CHECK_VECTOR);
    size_t vstride = abvbytes(placeholder, 1);
    size_t stride =  abbytes(abase, 1);

    mmt_generic_vec check_vector;
    matmul_top_vec_init_generic(mmt, abvbytes(placeholder, 1),
            check_vector, !dir, THREAD_SHARED_VECTOR);
    if (tcan_print) { printf("Loading check vector..."); fflush(stdout); }
    matmul_top_load_vector_generic(mmt, vstride,
            check_vector, "C", !dir, interval);
    if (tcan_print) { printf("done\n"); }

    mmt_vec ahead;
    vec_init_generic(pi->m, stride, (mmt_generic_vec_ptr) ahead, 0, NCHECKS_CHECK_VECTOR);

    /* We'll store all xy matrices locally before doing reductions. Given
     * the small footprint of these matrices, it's rather innocuous.
     */
    mmt_vec xymats;
    if (tcan_print) {
        printf("Each thread allocates %zu kb for the A matrices\n",
                abbytes(abase, m*interval) >> 10);
    }
    vec_init_generic(pi->m, stride, (mmt_generic_vec_ptr) xymats, 0, m*interval);
    
    for(int s = start ; ; s += interval ) {
        // Plan ahead. The check vector is here to predict the final A matrix.
        // Note that our share of the dot product is determined by the
        // intersections of the i0..i1 intervals on both sides.
        
        unsigned int how_many;
        unsigned int offset_c;
        unsigned int offset_v;
        how_many = intersect_two_intervals(&offset_c, &offset_v,
                mmt->wr[!dir]->i0, mmt->wr[!dir]->i1,
                mmt->wr[dir]->i0, mmt->wr[dir]->i1);

        abzero(abase, ahead->v, NCHECKS_CHECK_VECTOR);
        if (how_many) {
            size_t bytes_c =  abvbytes(placeholder, offset_c);
            abvt * c = abase_generic_ptr_add(check_vector->v, bytes_c);
            abt * v = mmt->wr[dir]->v->v + aboffset(abase, offset_v);
            abvdotprod(abase, placeholder, ahead->v, c, v, how_many);
        }

        /*
        debug_write(abase, ahead->v, NCHECKS_CHECK_VECTOR, "ahead.%u.j%u.t%u",
                s, mmt->pi->m->jrank, mmt->pi->m->trank);
         */

        abzero(abase, xymats->v, m*interval);

        serialize(pi->m);
        for(int i = 0 ; i < interval ; i++) {
            /* Compute the product by x */
            x_dotprod(mmt, xymats->v + aboffset(abase, i * m), m, gxvecs, nx);
            matmul_top_mul(mmt, dir);

            /*
            debug_write(abase, mmt->wr[dir]->v->v,
                    mmt->wr[dir]->i1-mmt->wr[dir]->i0, "pV%u.j%u.t%u",
                    s+i+1, mmt->pi->m->jrank, mmt->pi->m->trank);
             */
        }
        serialize(pi->m);

        /*
        debug_write(abase, xymats->v, m * interval, "xy.%u-%u.j%u.t%u",
                s, s + interval, mmt->pi->m->jrank, mmt->pi->m->trank);
         */

        /* Last dot product. This must cancel ! */
        x_dotprod(mmt, ahead->v, m, gxvecs, nx);

        /*
        debug_write(abase, ahead->v, NCHECKS_CHECK_VECTOR, "post%u.j%u.t%u",
                s, mmt->pi->m->jrank, mmt->pi->m->trank);
         */

        reduce_generic(abase, ahead, pi->m, NCHECKS_CHECK_VECTOR);
        if (!abis_zero(abase, ahead->v, NCHECKS_CHECK_VECTOR)) {
            printf("Failed check at iteration %d\n", s + interval);
            exit(1);
        }

        /* Now (and only now) collect the xy matrices */
        reduce_generic(abase, xymats, pi->m, m * interval);

        if (pi->m->trank == 0 && pi->m->jrank == 0) {
            char * tmp;
            asprintf(&tmp, "A%u-%u.%u-%u", ys[0], ys[1], s, s+interval);
            FILE * f = fopen(tmp, "w");
            fwrite(xymats->v, stride, aboffset(abase, m*interval), f);
            fclose(f);
        }
        matmul_top_save_vector(mmt, "Y", dir, s + interval);

        if (serialize(pi->m)) {
            printf("Reached iteration %d\n", s + interval);
        }
    }

    if (tcan_print) {
        printf("\n");
    }
    serialize(pi->m);

    matmul_top_vec_clear_generic(mmt, abvbytes(placeholder, 1), check_vector, !dir);
    vec_clear_generic(pi->m, stride, (mmt_generic_vec_ptr) xymats, m*interval);
    vec_clear_generic(pi->m, stride, (mmt_generic_vec_ptr) ahead, NCHECKS_CHECK_VECTOR);
    free(gxvecs);

    serialize(pi->m);
    matmul_top_clear(mmt, abase);
    serialize(pi->m);

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
    bw_common_init_mpi(argc, argv);

    if (nx == 0) { fprintf(stderr, "no nx value set\n"); exit(1); } 
    if (ys[0] < 0) { fprintf(stderr, "no ys value set\n"); exit(1); }

    abobj_init(abase);
    abobj_set_nbys(abase, ys[1]-ys[0]);

    pi_go(krylov_prog, mpi_split[0], mpi_split[1], thr_split[0], thr_split[1], 0);

    bw_common_clear();
    return 0;
}

