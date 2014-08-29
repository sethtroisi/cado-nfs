#include "cado.h"

#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <limits.h>
#include <dirent.h>
#include <errno.h>
#include "bwc_config.h"
#include "parallelizing_info.h"
#include "matmul_top.h"
#include "select_mpi.h"

#include "params.h"
#include "xvectors.h"
#include "portability.h"
#include "misc.h"
#include "bw-common-mpi.h"
#include "filenames.h"
#include "balancing.h"
#include "mpfq/mpfq.h"
#include "mpfq/mpfq_vbase.h"

struct sfile_info {
    unsigned int n0,n1;
    unsigned int iter;
};

struct sfiles_list {
    unsigned int s0,s1;
    int sfiles_alloc;
    int nsfiles;
    struct sfile_info (*sfiles)[1];
};

struct sols_list {
    int sols_alloc;
    int nsols;
    struct sfiles_list (*sols)[1];
};

int exitcode = 0;

static void prelude(parallelizing_info_ptr pi, struct sols_list * sl)
{
    sl->sols_alloc=0;
    sl->sols = NULL;
    sl->nsols=0;
    serialize_threads(pi->m);
    const char * spat = S_FILE_BASE_PATTERN ".%u" "%n";
    if (pi->m->jrank == 0 && pi->m->trank == 0) {
        /* It's our job to collect the directory data.
         */
        DIR * dir;
        dir = opendir(".");
        struct dirent * de;
        for( ; (de = readdir(dir)) != NULL ; ) {
            int k;
            int rc;
            unsigned int n0,n1;
            unsigned int s0,s1;
            unsigned int iter;
            rc = sscanf(de->d_name, spat, &s0, &s1, &n0, &n1, &iter, &k);
            if (rc < 5 || k != (int) strlen(de->d_name)) {
                continue;
            }
            int si = 0;
            for( ; si < sl->nsols ; si++) {
                if (sl->sols[si]->s0 == s0 && sl->sols[si]->s1 == s1)
                    break;
            }
            if (si == sl->nsols) {
                if (si >= sl->sols_alloc) {
                    sl->sols_alloc += 8 + sl->sols_alloc/2;
                    sl->sols = realloc(sl->sols, sl->sols_alloc * sizeof(struct sfiles_list));
                }
                memset(sl->sols[si], 0, sizeof(sl->sols[si]));
                sl->sols[si]->s0 = s0;
                sl->sols[si]->s1 = s1;
                si = sl->nsols++;
            }

            struct sfiles_list * s = sl->sols[si];

            /* append this file to the list of files for this solution set */
            if (s->nsfiles >= s->sfiles_alloc) {
                for( ; s->nsfiles >= s->sfiles_alloc ; s->sfiles_alloc += 8 + s->sfiles_alloc/2);
                s->sfiles = realloc(s->sfiles, s->sfiles_alloc * sizeof(struct sfile_info));
            }
            s->sfiles[s->nsfiles]->n0 = n0;
            s->sfiles[s->nsfiles]->n1 = n1;
            s->sfiles[s->nsfiles]->iter = iter;
            if (bw->interval && iter % bw->interval != 0) {
                fprintf(stderr,
                        "Warning: %s is not a checkpoint at a multiple of "
                        "the interval value %d -- this might indicate a "
                        "severe bug with leftover data, likely to corrupt "
                        "the final computation\n", de->d_name, bw->interval);
            }
            s->nsfiles++;
        }
        closedir(dir);
    }
    /* Note that it's not necessary to care about the file names -- in
     * practice, the file name is only relevant to the job/thread doing
     * actual I/O.
     */
    serialize(pi->m);
    global_broadcast(pi->m, sl, sizeof(struct sols_list), 0, 0);
    if (pi->m->jrank || pi->m->trank) {
        sl->sols = malloc(sl->sols_alloc * sizeof(struct sfiles_list));
    }
    for(int i = 0 ; i < sl->nsols ; i++) {
        struct sfiles_list * s = sl->sols[i];
        global_broadcast(pi->m, s, sizeof(struct sfiles_list), 0, 0);
        if (pi->m->jrank || pi->m->trank) {
            s->sfiles = malloc(s->sfiles_alloc * sizeof(struct sfile_info));
        }
        global_broadcast(pi->m, s->sfiles, s->nsfiles * sizeof(struct sfile_info), 0, 0);
    }
    serialize(pi->m);
}

int agree_on_flag(pi_wiring_ptr w, int v)
{
    static pthread_mutex_t mutex[1] = {PTHREAD_MUTEX_INITIALIZER};
    int * ptr = &v;
    thread_broadcast(w, (void**) &ptr, sizeof(void*), 0);
    for(unsigned int i = 0 ; i < w->ncores ; i++) {
        int ptrc = pthread_mutex_lock(mutex);
        ASSERT_ALWAYS(ptrc == 0);
        * ptr &= v;
        ptrc = pthread_mutex_unlock(mutex);
        ASSERT_ALWAYS(ptrc == 0);
    }
    serialize_threads(w);

    if (w->trank == 0) {
        int err = MPI_Allreduce(MPI_IN_PLACE, &v, 1, MPI_INT, MPI_LAND, w->pals);
        ASSERT_ALWAYS(!err);
    }
    serialize_threads(w);
    if (w->trank != 0) {
        v = *ptr;
    }
    /* serialize again because the leader must not release its stack
     * before all pals have read it ! */
    serialize_threads(w);
    return v;
}

void * gather_prog(parallelizing_info_ptr pi, param_list pl, void * arg MAYBE_UNUSED)
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

    mpfq_vbase A;
    mpfq_vbase_oo_field_init_byfeatures(A, 
            MPFQ_PRIME_MPZ, bw->p,
            MPFQ_GROUPSIZE, bw->nsolvecs,
            MPFQ_DONE);


    matmul_top_init(mmt, A, pi, flags, pl, bw->dir);
    unsigned int unpadded = MAX(mmt->n0[0], mmt->n0[1]);

    mmt_wiring_ptr mcol = mmt->wr[bw->dir];
    mmt_wiring_ptr mrow = mmt->wr[!bw->dir];

    struct sols_list sl[1];
    prelude(pi, sl);
    if (sl->nsols == 0) {
        if (tcan_print) {
            fprintf(stderr, "Found zero S files. Problem with command line ?\n");
            exitcode = 1;
        }
        serialize_threads(pi->m);
        return NULL;
    }

    pi_wiring_ptr picol = mmt->pi->wr[bw->dir];

    unsigned int ii0, ii1;
    unsigned int di = mcol->i1 - mcol->i0;

    ii0 = mcol->i0 + di * (picol->jrank * picol->ncores + picol->trank) /
        (picol->njobs * picol->ncores);
    ii1 = mcol->i0 + di * (picol->jrank * picol->ncores + picol->trank + 1) /
        (picol->njobs * picol->ncores);

    mmt_vec svec;
    mmt_vec tvec;
    vec_init_generic(mmt->pi->m, A, svec, 0, ii1-ii0);
    vec_init_generic(mmt->pi->m, A, tvec, 0, ii1-ii0);

    const char * rhs_name = param_list_lookup_string(pl, "rhs");
    FILE * rhs = NULL;
    mmt_vec rhsvec;
    if (rhs_name) {
        if (!pi->m->trank && !pi->m->jrank)
            rhs = fopen(rhs_name, "rb");
        matmul_top_vec_init_generic(mmt, A, rhsvec, bw->dir, 0);
    }

    for(int s = 0 ; s < sl->nsols ; s++) {
        struct sfiles_list * sf = sl->sols[s];

        char * kprefix;

        int rc = asprintf(&kprefix, K_FILE_BASE_PATTERN, sf->s0, sf->s1);

        ASSERT_ALWAYS(rc >= 0);

        if (tcan_print) {
            printf("Trying to build solutions %u..%u\n", sf->s0, sf->s1);
        }

        A->vec_set_zero(A, mrow->v->v, mrow->i1 - mrow->i0);

        /* This array receives the bw->nrhs coefficients affecting the
         * rhs * columns in the identities obtained */
        void * rhscoeffs;

        if (rhs_name) {
            A->vec_init(A, &rhscoeffs, bw->nrhs);

            if (tcan_print) {
                printf("Reading rhs coefficients for solution %u\n", sf->s0);
            }
            if (rhs) {      /* main thread */
                for(int j = 0 ; j < bw->nrhs ; j++) {
                    /* We're implicitly supposing that elements are stored
                     * alone, not grouped à la SIMD (which does happen for
                     * GF(2), of course). */
                    ASSERT_ALWAYS(sf->s1 == sf->s0 + 1);
                    rc = fseek(rhs, A->vec_elt_stride(A, j * bw->n + sf->s0), SEEK_SET);
                    ASSERT_ALWAYS(rc >= 0);
                    rc = fread(A->vec_coeff_ptr(A, rhscoeffs, j), A->vec_elt_stride(A,1), 1, rhs);
                    if (A->is_zero(A, A->vec_coeff_ptr_const(A, rhscoeffs, j))) {
                        printf("Notice: coefficient for vector " V_FILE_BASE_PATTERN " in solution %d is zero\n", j, j+1, sf->s0);
                    }
                    ASSERT_ALWAYS(rc == 1);
                }
            }
            global_broadcast(pi->m, rhscoeffs, A->vec_elt_stride(A,bw->nrhs), 0, 0);

            A->vec_set_zero(A, tvec->v, ii1 - ii0);
            for(int j = 0 ; j < bw->nrhs ; j++) {
                /* Add c_j times v_j to tvec */
                ASSERT_ALWAYS(sf->s1 == sf->s0 + 1);
                char * tmp;
                int rc = asprintf(&tmp, V_FILE_BASE_PATTERN, j, j + 1);
                ASSERT_ALWAYS(rc >= 0);

                pi_load_file_2d(pi, bw->dir, tmp, 0, svec->v, A->vec_elt_stride(A, ii1 - ii0), A->vec_elt_stride(A, unpadded));
                A->vec_scal_mul(A, svec->v, svec->v, A->vec_coeff_ptr(A, rhscoeffs, j), ii1 - ii0);
                A->vec_add(A, tvec->v, tvec->v, svec->v, ii1 - ii0);
                free(tmp);
            }

            /* Now save the sum to a temp file, which we'll read later on
             */
            {
                char * tmp;
                int rc = asprintf(&tmp, R_FILE_BASE_PATTERN, sf->s0, sf->s1);
                ASSERT_ALWAYS(rc >= 0);
                pi_save_file_2d(pi, bw->dir, tmp, 0, tvec->v, A->vec_elt_stride(A, ii1 - ii0), A->vec_elt_stride(A, unpadded));
                free(tmp);
            }
        }

        /* Let now tvec be the sum of the LHS contributions */
        A->vec_set_zero(A, tvec->v, ii1 - ii0);
        for(int i = 0 ; i < sf->nsfiles ; i++) {
            char * tmp;
            int rc = asprintf(&tmp, S_FILE_BASE_PATTERN,
                    sf->s0, sf->s1,
                    sf->sfiles[i]->n0, sf->sfiles[i]->n1);
            ASSERT_ALWAYS(rc >= 0);

            if (tcan_print) {
                printf("loading %s.%u\n", tmp, sf->sfiles[i]->iter);
            }
            pi_load_file_2d(pi, bw->dir, tmp, sf->sfiles[i]->iter, svec->v, A->vec_elt_stride(A, ii1 - ii0), A->vec_elt_stride(A, unpadded));
            free(tmp);
            A->vec_add(A, tvec->v, tvec->v, svec->v, ii1 - ii0);
        }

        /* As before, we're simply using file I/O as a means of
         * broadcast. It's much easier, since we'll be doing I/O
         * anyway...  */
        pi_save_file_2d(pi, bw->dir, kprefix, 0, tvec->v, A->vec_elt_stride(A, ii1 - ii0),  A->vec_elt_stride(A, unpadded));
        


        int is_zero = 0;

        serialize(pi->m);
        matmul_top_load_vector(mmt, kprefix, bw->dir, 0, unpadded);
        matmul_top_twist_vector(mmt, bw->dir);

        unsigned int how_many;
        unsigned int offset_c;
        unsigned int offset_v;
        how_many = intersect_two_intervals(&offset_c, &offset_v,
                    mrow->i0, mrow->i1,
                    mcol->i0, mcol->i1);

        void * check_area = SUBVEC(mcol->v, v, offset_v);
        is_zero = A->vec_is_zero(A, check_area, how_many);
        is_zero = agree_on_flag(pi->m, is_zero);

        if (is_zero) {
            fprintf(stderr, "Found zero vector. Most certainly a bug. "
                    "No solution found.\n");
            continue;
        }
        serialize(pi->m);

        /* Note that for the inhomogeneous case, we'll do the loop only
         * once, since we end with a break. */
        for(int i = 1 ; i < 10 ; i++) {
            serialize(pi->m);

            matmul_top_mul(mmt, bw->dir);

            matmul_top_untwist_vector(mmt, bw->dir);

            if (rhs_name) {
                char * tmp;
                int rc = asprintf(&tmp, R_FILE_BASE_PATTERN, sf->s0, sf->s1);
                ASSERT_ALWAYS(rc >= 0);
                matmul_top_load_vector_generic(mmt, rhsvec, tmp, bw->dir, 0, unpadded);
                free(tmp);
                if (rhs) {
                    /* Now get rid of the R file, we don't need it */
                    rc = asprintf(&tmp, R_FILE_BASE_PATTERN ".%u", sf->s0, sf->s1, 0);
                    ASSERT_ALWAYS(rc >= 0);
                    rc = unlink(tmp);
                    ASSERT_ALWAYS(rc == 0);
                    free(tmp);
                }

                A->vec_add(A, check_area, check_area,
                        SUBVEC(rhsvec, v, offset_v), how_many);
            }

            is_zero = A->vec_is_zero(A, check_area, how_many);
            is_zero = agree_on_flag(pi->m, is_zero);
            serialize(pi->m);
            if (agree_on_flag(pi->m, is_zero)) {
                if (tcan_print) {
                    if (rhs_name) {
                        printf("M^%u * V + R is zero\n", i);
                        // printf("M^%u * V + R is zero [" K_FILE_BASE_PATTERN ".%u contains M^%u * V, R.sols%u-%u contains R]!\n", i, sf->s0, sf->s1, i-1, i-1, sf->s0, sf->s1);
                    } else {
                        printf("M^%u * V is zero [" K_FILE_BASE_PATTERN ".%u contains M^%u * V]!\n", i, sf->s0, sf->s1, i-1, i-1);
                    }
                }
                break;
            }
            if (rhs_name) break;

            matmul_top_save_vector(mmt, kprefix, bw->dir, i, unpadded);
            matmul_top_twist_vector(mmt, bw->dir);
        }
        if (!is_zero) {
            if (tcan_print) {
                printf("Solution range %u..%u: no solution found, most probably a bug\n", sf->s0, sf->s1);
            }
            exitcode=1;
            return NULL;
        }
        if (rhs) {
            /* We now append the rhs coefficients to the vector we
             * have just saved. Of course only the I/O thread does this. */
            char * tmp;
            int rc = asprintf(&tmp, "%s.%u", kprefix, 0);
            ASSERT_ALWAYS(rc >= 0);
            printf("Expanding %s so as to include the coefficients for the %d RHS columns\n", tmp, bw->nrhs);
            FILE * f = fopen(tmp, "ab");
            rc = fseek(f, 0, SEEK_END);
            ASSERT_ALWAYS(rc >= 0);
            rc = fwrite(rhscoeffs, A->vec_elt_stride(A, 1), bw->nrhs, f);
            ASSERT_ALWAYS(rc == bw->nrhs);
            fclose(f);
            printf("%s is now a right nullspace vector for (M|RHS).\n", tmp);
            free(tmp);
        }

        serialize(pi->m);
        free(kprefix);
        if (rhs_name) A->vec_clear(A, &rhscoeffs, bw->nrhs);
    }

    if (rhs) fclose(rhs);
    if (rhs_name) matmul_top_vec_clear_generic(mmt, rhsvec, bw->dir);

    vec_clear_generic(mmt->pi->m, svec, ii1-ii0);
    vec_clear_generic(mmt->pi->m, tvec, ii1-ii0);

    matmul_top_clear(mmt);
    A->oo_field_clear(A);

    for(int i = 0 ; i < sl->nsols ; i++) {
        free(sl->sols[i]->sfiles);
    }
    free(sl->sols);

    return NULL;
}


void usage()
{
    fprintf(stderr, "Usage: ./gather <options>\n");
    fprintf(stderr, "%s", bw_common_usage_string());
    fprintf(stderr, "Relevant options here: wdir cfg m n mpi thr matrix interval\n");
    fprintf(stderr, "Note: data files must be found in wdir !\n");
    fprintf(stderr, "All S* files are taken in wdir. Spurious leftover data may corrupt result !\n");
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

    pi_go(gather_prog, pl, 0);

    param_list_clear(pl);
    bw_common_clear_mpi(bw);
    return exitcode;
}

