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
#include "mpfq/mpfq.h"
#include "mpfq/mpfq_vbase.h"
#include "cheating_vec_init.h"

void * krylov_prog(parallelizing_info_ptr pi, param_list pl, void * arg MAYBE_UNUSED)
{
    int fake = param_list_lookup_string(pl, "random_matrix") != NULL;
    int tcan_print = bw->can_print && pi->m->trank == 0;
    matmul_top_data mmt;
    mmt_vec y, my;
    struct timing_data timing[1];

    int ys[2] = { bw->ys[0], bw->ys[1], };
    if (pi->interleaved) {
        ASSERT_ALWAYS((bw->ys[1]-bw->ys[0]) % 2 == 0);
        ys[0] = bw->ys[0] + pi->interleaved->idx * (bw->ys[1]-bw->ys[0])/2;
        ys[1] = ys[0] + (bw->ys[1]-bw->ys[0])/2;
    }

    int withcoeffs = mpz_cmp_ui(bw->p, 2) > 0;
    int nchecks = withcoeffs ? NCHECKS_CHECK_VECTOR_GFp : NCHECKS_CHECK_VECTOR_GF2;
    mpfq_vbase A;
    mpfq_vbase_oo_field_init_byfeatures(A, 
            MPFQ_PRIME_MPZ, bw->p,
            MPFQ_GROUPSIZE, ys[1]-ys[0],
            MPFQ_DONE);

    /* Hmmm. This would deserve better thought. Surely we don't need 64
     * in the prime case. Anything which makes checks relevant will do.
     * For the binary case, we used to work with 64 as a constant, but
     * for the prime case we want to make this tunable (or maybe 1 ?)
     */
    mpfq_vbase Ac;
    mpfq_vbase_oo_field_init_byfeatures(Ac,
            MPFQ_PRIME_MPZ, bw->p,
            MPFQ_GROUPSIZE, nchecks,
            MPFQ_DONE);

    pi_datatype_ptr Ac_pi = pi_alloc_mpfq_datatype(pi, Ac);

    block_control_signals();

    matmul_top_init(mmt, A, pi, pl, bw->dir);
    mmt_vec_init(mmt,0,0, y,   bw->dir, /* shared ! */ 1, mmt->n[bw->dir]);
    mmt_vec_init(mmt,0,0, my, !bw->dir,                0, mmt->n[!bw->dir]);

    unsigned int unpadded = MAX(mmt->n0[0], mmt->n0[1]);

    serialize(pi->m);
    
    mmt_full_vec_set_zero(my);

    uint32_t * gxvecs = NULL;
    unsigned int nx = 0;
    if (!fake) {
        load_x(&gxvecs, bw->m, &nx, pi);
    } else {
        set_x_fake(&gxvecs, bw->m, &nx, pi);
    }

    indices_twist(mmt, gxvecs, nx * bw->m, bw->dir);

    /* let's be generous with interleaving protection. I don't want to be
     * bothered, really */
    pi_interleaving_flip(pi);
    pi_interleaving_flip(pi);
    matmul_top_comm_bench(mmt, bw->dir);
    pi_interleaving_flip(pi);
    pi_interleaving_flip(pi);

    char * v_name = NULL;
    int rc = asprintf(&v_name, V_FILE_BASE_PATTERN, ys[0], ys[1]);
    ASSERT_ALWAYS(rc >= 0);
    if (!fake) {
        if (tcan_print) { printf("Loading %s...", v_name); fflush(stdout); }
        mmt_vec_load(y, v_name, bw->start, unpadded);
        if (tcan_print) { printf("done\n"); }
    } else {
        gmp_randstate_t rstate;
        gmp_randinit_default(rstate);
#if 0
        /* This is for setting the source vector to something consistent
         * across mappings, so that given a fixed (fake, here) matrix, any
         * splitting will give the same source vector. This is just a
         * copy of the mechanism which exists in prep for doing exactly
         * this. Alas, in what we denote as a fake situation, there is
         * no chance of course that two different splittings lesad to
         * identical matrices ! Hence, we'd rather not bother with
         * generating something consistent.
         */
        if (pi->m->trank == 0 && !bw->seed) {
            bw->seed = time(NULL);
            MPI_Bcast(&bw->seed, 1, MPI_INT, 0, pi->m->pals);
        }
        serialize_threads(pi->m);
        gmp_randseed_ui(rstate, bw->seed);
        if (tcan_print) {
            printf("// Random generator seeded with %d\n", bw->seed);
        }
        if (tcan_print) { printf("Creating fake %s...", v_name); fflush(stdout); }
        mmt_vec_set_random_through_file(mmt, NULL, v_name, bw->dir, bw->start, unpadded, rstate);
        if (tcan_print) { printf("done\n"); }
#else
        unsigned long g = pi->m->jrank * pi->m->ncores + pi->m->trank;
        gmp_randseed_ui(rstate, bw->seed + g);
        mmt_vec_set_random_inconsistent(y, rstate);
        mmt_vec_allreduce(y);
#endif
        gmp_randclear(rstate);
    }

    mmt_vec check_vector;
    void * ahead = NULL;

    mpfq_vbase_tmpl AxAc;
    mpfq_vbase_oo_init_templates(AxAc, A, Ac);

    if (!bw->skip_online_checks) {
        /* We do the dot product by working on the local vector chunks.
         * Therefore, we must really understand the check vector as
         * playing a role in the very same direction of the y vector!
         */
        mmt_vec_init(mmt, Ac, Ac_pi,
                check_vector, bw->dir, THREAD_SHARED_VECTOR, mmt->n[bw->dir]);
        if (tcan_print) { printf("Loading check vector..."); fflush(stdout); }
        mmt_vec_load(check_vector, CHECK_FILE_BASE, bw->interval,  mmt->n0[bw->dir]);
        if (tcan_print) { printf("done\n"); }
    }

    if (!bw->skip_online_checks) {
        cheating_vec_init(A, &ahead, nchecks);
    }

    /* We'll store all xy matrices locally before doing reductions. Given
     * the small footprint of these matrices, it's rather innocuous.
     */
    void * xymats;

    if (tcan_print) {
        printf("Each thread allocates %zd kb for the A matrices\n",
                A->vec_elt_stride(A, bw->m*bw->interval) >> 10);
    }
    cheating_vec_init(A, &xymats, bw->m*bw->interval);
    
    if (bw->end == 0) {
        /* Decide on an automatic ending value */
        unsigned int length;
        /* The padded dimension is not the important one */
        length = MAX(mmt->n0[0], mmt->n0[1]);
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

#if 0
    /* FIXME -- that's temporary ! only for debugging */
    pi_log_init(pi->m);
    pi_log_init(pi->wr[0]);
    pi_log_init(pi->wr[1]);
#endif

    timing_init(timing, 4 * mmt->nmatrices, bw->start, bw->interval * iceildiv(bw->end, bw->interval));
    for(int i = 0 ; i < mmt->nmatrices; i++) {
        timing_set_timer_name(timing, 4*i, "CPU%d", i);
        timing_set_timer_items(timing, 4*i, mmt->matrices[i]->mm->ncoeffs);
        timing_set_timer_name(timing, 4*i+1, "cpu-wait%d", i);
        timing_set_timer_name(timing, 4*i+2, "COMM%d", i);
        timing_set_timer_name(timing, 4*i+3, "comm-wait%d", i);
    }

    pi_interleaving_flip(pi);
    pi_interleaving_flip(pi);

    for(int s = bw->start ; s < bw->end ; s += bw->interval ) {
        // Plan ahead. The check vector is here to predict the final A matrix.
        // Note that our share of the dot product is determined by the
        // intersections of the i0..i1 intervals on both sides.
        
        /* Note that the check vector is always stored untwisted in
         * memory */

        if (!bw->skip_online_checks) {
            A->vec_set_zero(A, ahead, nchecks);
            AxAc->dotprod(A->obj, Ac->obj, ahead,
                    mmt_my_own_subvec(check_vector),
                    mmt_my_own_subvec(y),
                    mmt_my_own_size_in_items(y));
        }

        /* Create an empty slot in program execution, so that we don't
         * impose strong constraints on twist/untwist_vector being free of
         * MPI calls.
         */
        pi_interleaving_flip(pi);
        pi_interleaving_flip(pi);
        mmt_vec_twist(mmt, y);

        A->vec_set_zero(A, xymats, bw->m*bw->interval);
        serialize(pi->m);
        for(int i = 0 ; i < bw->interval ; i++) {
            /* The first part of this loop must be guaranteed to be free
             * of any mpi calls */
            pi_interleaving_flip(pi);

            /* Compute the product by x */
            x_dotprod(A->vec_subvec(A, xymats, i * bw->m),
                    gxvecs, bw->m, nx, y, 1);

            /* copy code from matmul_top_mul() ; this should really be
             * factored out */

            /* Timers take, in sequence, the 4n values defined in the
             * definition of the timing structure. The first is timer 0.
             */
            mmt_vec_ptr x[2] = { y, my };
            if (y->d == 0) {
                for(int mx = 0, k = 0 ; mx < mmt->nmatrices ; mx++, k^=1) {
                    mmt_vec_ptr src = x[k];
                    mmt_vec_ptr dst = x[!k];

                    matmul_top_mul_cpu(mmt, mx, dst, src);

                    timing_next_timer(timing);
                    pi_interleaving_flip(pi);
                    serialize(pi->m);
                    timing_next_timer(timing);

                    /* Now we can resume MPI communications. */
                    if (mx == (mmt->nmatrices - 1) && (mmt->nmatrices & 1) == 1) {
                        ASSERT_ALWAYS(src == y);
                        matmul_top_mul_comm(src, dst);
                    } else {
                        mmt_vec_allreduce(dst);
                    }

                    timing_next_timer(timing);
                    serialize(pi->m);
                    timing_next_timer(timing);
                }
            } else {
                /* same, in reverse order */
                for(int mx = mmt->nmatrices, k = 0 ; mx-- > 0 ; k^=1) {
                    mmt_vec_ptr src = x[k];
                    mmt_vec_ptr dst = x[!k];
                    ASSERT_ALWAYS(src->consistency == 2);

                    matmul_top_mul_cpu(mmt, mx, dst, src);
                    ASSERT_ALWAYS(dst->consistency == 0);

                    timing_next_timer(timing);
                    pi_interleaving_flip(pi);
                    serialize(pi->m);
                    timing_next_timer(timing);

                    if (mx == 0 && (mmt->nmatrices & 1) == 1) {
                        ASSERT_ALWAYS(src == y);
                        matmul_top_mul_comm(src, dst);
                    } else {
                        mmt_vec_allreduce(dst);
                    }
                    timing_next_timer(timing);
                    serialize(pi->m);
                    timing_next_timer(timing);
                }
            }

            timing_check(pi, timing, s+i+1, tcan_print);
        }
        serialize(pi->m);

        /* See remark above. */
        pi_interleaving_flip(pi);
        pi_interleaving_flip(pi);

        if (!bw->skip_online_checks) {
            /* Last dot product. This must cancel ! */
            x_dotprod(ahead, gxvecs, nchecks, nx, y, -1);

            pi_allreduce(NULL, ahead, nchecks, mmt->pitype, BWC_PI_SUM, pi->m);
            if (!A->vec_is_zero(A, ahead, nchecks)) {
                printf("Failed check at iteration %d\n", s + bw->interval);
                exit(1);
            }
        }

        mmt_vec_untwist(mmt, y);

        /* Now (and only now) collect the xy matrices */
        pi_allreduce(NULL, xymats,
                bw->m * bw->interval,
                mmt->pitype, BWC_PI_SUM, pi->m);

        if (pi->m->trank == 0 && pi->m->jrank == 0) {
            char * tmp;
            int rc;
            rc = asprintf(&tmp, A_FILE_PATTERN, ys[0], ys[1], s, s+bw->interval);
            FILE * f = fopen(tmp, "wb");
            rc = fwrite(xymats, A->vec_elt_stride(A, 1), bw->m*bw->interval, f);
            if (rc != bw->m*bw->interval) {
                fprintf(stderr, "Ayee -- short write\n");
                // make sure our input data won't be deleted -- this
                // chunk will have to be redone later, maybe the disk
                // failure is temporary (?)
            }
            fclose(f);
            free(tmp);
        }

        mmt_vec_save(y, v_name, s + bw->interval, mmt->n0[bw->dir]);

        if (pi->m->trank == 0 && pi->m->jrank == 0)
            keep_rolling_checkpoints(v_name, s + bw->interval);

        serialize(pi->m);

        // reached s + bw->interval. Count our time on cpu, and compute the sum.
        timing_disp_collective_oneline(pi, timing, s + bw->interval, tcan_print, "krylov");
    }

    timing_final_tally(pi, timing, tcan_print, "krylov");

#if 0
    pi_log_clear(pi->m);
    pi_log_clear(pi->wr[0]);
    pi_log_clear(pi->wr[1]);
#endif

    if (tcan_print) {
        printf("Done krylov.\n");
    }
    serialize(pi->m);

    cheating_vec_clear(A, &xymats, bw->m*bw->interval);

    if (!bw->skip_online_checks) {
        mmt_vec_clear(mmt, check_vector);
        cheating_vec_clear(A, &ahead, nchecks);
    }

    free(gxvecs);
    free(v_name);

    mmt_vec_clear(mmt, y);
    mmt_vec_clear(mmt, my);

    matmul_top_clear(mmt);
    pi_free_mpfq_datatype(pi, Ac_pi);
    A->oo_field_clear(A);
    Ac->oo_field_clear(Ac);

    timing_clear(timing);

    return NULL;
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
    /* interpret our parameters */
    if (bw->ys[0] < 0) { fprintf(stderr, "no ys value set\n"); exit(1); }
    catch_control_signals();

    if (param_list_warn_unused(pl)) {
        param_list_print_usage(pl, bw->original_argv[0], stderr);
        exit(EXIT_FAILURE);
    }

    pi_go(krylov_prog, pl, 0);

    parallelizing_info_finish();
    param_list_clear(pl);
    bw_common_clear_new(bw);

    return 0;
}

