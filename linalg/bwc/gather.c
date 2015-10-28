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
#include "bw-common.h"
#include "filenames.h"
#include "balancing.h"
#include "mpfq/mpfq.h"
#include "mpfq/mpfq_vbase.h"
#include "cheating_vec_init.h"

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
    pi_bcast(sl, sizeof(struct sols_list), BWC_PI_BYTE, 0, 0, pi->m);
    if (pi->m->jrank || pi->m->trank) {
        sl->sols = malloc(sl->sols_alloc * sizeof(struct sfiles_list));
    }
    for(int i = 0 ; i < sl->nsols ; i++) {
        struct sfiles_list * s = sl->sols[i];
        pi_bcast(s, sizeof(struct sfiles_list), BWC_PI_BYTE, 0, 0, pi->m);
        if (pi->m->jrank || pi->m->trank) {
            s->sfiles = malloc(s->sfiles_alloc * sizeof(struct sfile_info));
        }
        pi_bcast(s->sfiles,
                s->nsfiles * sizeof(struct sfile_info), BWC_PI_BYTE,
                0, 0,
                pi->m);
    }
    serialize(pi->m);
}

void * gather_prog(parallelizing_info_ptr pi, param_list pl, void * arg MAYBE_UNUSED)
{
    /* Interleaving does not make sense for this program. So the second
     * block of threads just leave immediately */
    if (pi->interleaved && pi->interleaved->idx)
        return NULL;

    int tcan_print = bw->can_print && pi->m->trank == 0;
    matmul_top_data mmt;
    mmt_vec y, my;

    int nsolvecs = bw->nsolvecs;
    unsigned int multi = 1;
    unsigned int nsolvecs_pervec = nsolvecs;

    if (mpz_cmp_ui(bw->p,2) != 0 && nsolvecs > 1) {
        if (tcan_print) {
            fprintf(stderr,
"Note: the present code is a quick hack.\n"
"Some SIMD operations like\n"
"    w_i += c_i * v (with v a constant vector, c_i scalar)\n"
"are performed with a simple-minded loop on i, while there is an acknowledged\n"
"potential for improvement by using more SIMD-aware code.\n");
        }
        multi = nsolvecs;
        nsolvecs_pervec = 1;
    }

    ASSERT_ALWAYS(multi * nsolvecs_pervec == (unsigned int) nsolvecs);

    mpfq_vbase A;
    mpfq_vbase_oo_field_init_byfeatures(A, 
            MPFQ_PRIME_MPZ, bw->p,
            MPFQ_GROUPSIZE, nsolvecs_pervec,
            MPFQ_DONE);

    matmul_top_init(mmt, A, pi, pl, bw->dir);
    mmt_vec_init(mmt,0,0, y,   bw->dir, /* shared ! */ 1, mmt->n[bw->dir]);
    mmt_vec_init(mmt,0,0, my, !bw->dir,                0, mmt->n[!bw->dir]);


    /* this is really a misnomer, because in the typical case, M is
     * rectangular, and then the square matrix does induce some padding.
     * This is however not the padding we are interested in. The padding
     * we're referring to in the naming of this variable is the one which
     * is related to the number of jobs and threads: internal dimensions
     * are arranged to be multiples */
    unsigned int unpadded = MAX(mmt->n0[0], mmt->n0[1]);
    unsigned int nrhs = 0;

    struct sols_list sl[1];
    prelude(pi, sl);
    if (sl->nsols == 0) {
        if (tcan_print) {
            fprintf(stderr, "Found zero S files. Problem with command line ?\n");
            pthread_mutex_lock(pi->m->th->m);
            exitcode=1;
            pthread_mutex_unlock(pi->m->th->m);
        }
        serialize_threads(pi->m);
        return NULL;
    }

    size_t eblock = mmt_my_own_size_in_items(y);
    void * svec;
    void * tvec;
    cheating_vec_init(A, &svec, eblock);
    cheating_vec_init(A, &tvec, eblock);

    const char * rhs_name = param_list_lookup_string(pl, "rhs");
    if (rhs_name != NULL) {
        if (pi->m->jrank == 0 && pi->m->trank == 0)
            get_rhs_file_header(rhs_name, NULL, &nrhs, NULL);
        pi_bcast(&nrhs, 1, BWC_PI_UNSIGNED, 0, 0, pi->m);
    }

    if (tcan_print && nrhs) {
        if (nrhs) {
            printf("** Informational note about GF(p) inhomogeneous system:\n");
            printf("   Original matrix dimensions: %" PRIu32" %" PRIu32"\n", mmt->n0[0], mmt->n0[1]);
            printf("   We expect to obtain a vector of size %" PRIu32"\n",
                    mmt->n0[!bw->dir] + nrhs);
            printf("   which we hope will be a kernel vector for (M_square||RHS).\n"
                   "   We will discard the coefficients which correspond to padding columns\n"
                   "   This entails keeping coordinates in the intervals\n"
                   "   [0..%" PRIu32"[ and [%" PRIu32"..%" PRIu32"[\n"
                   "   in the result.\n"
                   "** end note.\n",
            mmt->n0[bw->dir], mmt->n0[!bw->dir], mmt->n0[!bw->dir] + nrhs);
        }
    }

    const char * rhscoeffs_name = param_list_lookup_string(pl, "rhscoeffs");
    FILE * rhscoeffs_file = NULL;
    mmt_vec rhscoeffs_vec;
    if (rhscoeffs_name) {
        if (!pi->m->trank && !pi->m->jrank)
            rhscoeffs_file = fopen(rhscoeffs_name, "rb");
        mmt_vec_init(mmt, A, mmt->pitype, rhscoeffs_vec, bw->dir, 0, mmt->n[bw->dir]);
    }

    for(int s = 0 ; s < sl->nsols ; s++) {
        struct sfiles_list * sf = sl->sols[s];

        char * kprefix;

        int rc = asprintf(&kprefix, K_FILE_BASE_PATTERN, sf->s0, sf->s1);

        ASSERT_ALWAYS(rc >= 0);

        if (tcan_print) {
            printf("Trying to build solutions %u..%u\n", sf->s0, sf->s1);
        }

        mmt_full_vec_set_zero(my);

        /* This array receives the nrhs coefficients affecting the
         * rhs * columns in the identities obtained */
        void * rhscoeffs = NULL;

        if (rhscoeffs_name) {
            cheating_vec_init(A, &rhscoeffs, nrhs);

            if (tcan_print) {
                printf("Reading rhs coefficients for solution %u\n", sf->s0);
            }
            if (rhscoeffs_file) {      /* main thread */
                for(unsigned int j = 0 ; j < nrhs ; j++) {
                    /* We're implicitly supposing that elements are stored
                     * alone, not grouped à la SIMD (which does happen for
                     * GF(2), of course). */
                    ASSERT_ALWAYS(sf->s1 == sf->s0 + 1);
                    rc = fseek(rhscoeffs_file, A->vec_elt_stride(A, j * bw->n + sf->s0), SEEK_SET);
                    ASSERT_ALWAYS(rc >= 0);
                    rc = fread(A->vec_coeff_ptr(A, rhscoeffs, j), A->vec_elt_stride(A,1), 1, rhscoeffs_file);
                    if (A->is_zero(A, A->vec_coeff_ptr_const(A, rhscoeffs, j))) {
                        printf("Notice: coefficient for vector " V_FILE_BASE_PATTERN " in solution %d is zero\n", j, j+1, sf->s0);
                    }
                    ASSERT_ALWAYS(rc == 1);
                }
            }
            pi_bcast(rhscoeffs, nrhs, mmt->pitype, 0, 0, pi->m);

            A->vec_set_zero(A, tvec, eblock);
            for(unsigned int j = 0 ; j < nrhs ; j++) {
                /* Add c_j times v_j to tvec */
                ASSERT_ALWAYS(sf->s1 == sf->s0 + 1);
                char * tmp;
                int rc = asprintf(&tmp, V_FILE_BASE_PATTERN ".0", j, j + 1);
                ASSERT_ALWAYS(rc >= 0);

                pi_file_handle f;
                rc = pi_file_open(f, pi, bw->dir, tmp, "rb");
                if (tcan_print && !rc) fprintf(stderr, "%s: not found\n", tmp);
                ASSERT_ALWAYS(rc);
                ssize_t s = pi_file_read(f, svec, A->vec_elt_stride(A, eblock), A->vec_elt_stride(A, unpadded));
                ASSERT_ALWAYS(s >= 0 && s == A->vec_elt_stride(A, unpadded));
                pi_file_close(f);
                free(tmp);

                A->vec_scal_mul(A, svec, svec, A->vec_coeff_ptr(A, rhscoeffs, j), eblock);
                A->vec_add(A, tvec, tvec, svec, eblock);
            }

            /* Now save the sum to a temp file, which we'll read later on
             */
            {
                char * tmp;
                int rc = asprintf(&tmp, R_FILE_BASE_PATTERN ".0", sf->s0, sf->s1);
                ASSERT_ALWAYS(rc >= 0);
                pi_file_handle f;
                rc = pi_file_open(f, pi, bw->dir, tmp, "wb");
                if (tcan_print && !rc) fprintf(stderr, "%s: not found\n", tmp);
                ASSERT_ALWAYS(rc);
                ssize_t s = pi_file_write(f, tvec, A->vec_elt_stride(A, eblock), A->vec_elt_stride(A, unpadded));
                ASSERT_ALWAYS(s >= 0 && s == A->vec_elt_stride(A, unpadded));
                pi_file_close(f);
                free(tmp);
            }
        }

        /* Collect now the sum of the LHS contributions */
        mmt_full_vec_set_zero(y);
        void * sv = mmt_my_own_subvec(y);

        for(int i = 0 ; i < sf->nsfiles ; i++) {
            char * tmp;
            int rc = asprintf(&tmp, S_FILE_BASE_PATTERN ".%u",
                    sf->s0, sf->s1,
                    sf->sfiles[i]->n0, sf->sfiles[i]->n1,
                    sf->sfiles[i]->iter);
            ASSERT_ALWAYS(rc >= 0);

            if (tcan_print && verbose_enabled(CADO_VERBOSE_PRINT_BWC_LOADING_MKSOL_FILES)) {
                printf("loading %s\n", tmp);
            }
            pi_file_handle f;
            rc = pi_file_open(f, pi, bw->dir, tmp, "rb");
            if (tcan_print && !rc) fprintf(stderr, "%s: not found\n", tmp);
            ASSERT_ALWAYS(rc);
            ssize_t s = pi_file_read(f, svec, A->vec_elt_stride(A, eblock), A->vec_elt_stride(A, unpadded));
            ASSERT_ALWAYS(s >= 0 && s == A->vec_elt_stride(A, unpadded));
            pi_file_close(f);
            free(tmp);

            A->vec_add(A, sv, sv, svec, eblock);
        }
        y->consistency = 1;

        mmt_vec_broadcast(y);

        mmt_vec_unapply_T(mmt, y);
        mmt_vec_save(y, kprefix, 0, unpadded);
        mmt_vec_apply_T(mmt, y);
                
        serialize(pi->m);
        mmt_vec_twist(mmt, y);

        int is_zero = A->vec_is_zero(A,
                mmt_my_own_subvec(y), mmt_my_own_size_in_items(y));
        pi_allreduce(NULL, &is_zero, 1, BWC_PI_INT, BWC_PI_MIN, pi->m);

        if (is_zero) {
            fprintf(stderr, "Found zero vector. Most certainly a bug. "
                    "No solution found.\n");
            exitcode=1;
            continue;
        }
        serialize(pi->m);

        /* Note that for the inhomogeneous case, we'll do the loop only
         * once, since we end with a break. */
        for(int i = 1 ; i < 10 ; i++) {
            serialize(pi->m);

            matmul_top_mul(mmt, my, y);

            mmt_vec_untwist(mmt, y);

            if (rhscoeffs_name) {
                char * tmp;
                int rc = asprintf(&tmp, R_FILE_BASE_PATTERN, sf->s0, sf->s1);
                ASSERT_ALWAYS(rc >= 0);
                mmt_vec_load(rhscoeffs_vec, tmp, 0, unpadded);
                free(tmp);
                if (rhscoeffs_file) {
                    /* Now get rid of the R file, we don't need it */
                    rc = asprintf(&tmp, R_FILE_BASE_PATTERN ".%u", sf->s0, sf->s1, 0);
                    ASSERT_ALWAYS(rc >= 0);
                    rc = unlink(tmp);
                    ASSERT_ALWAYS(rc == 0);
                    free(tmp);
                }

                A->vec_add(A,
                        mmt_my_own_subvec(y),
                        mmt_my_own_subvec(y),
                        mmt_my_own_subvec(rhscoeffs_vec),
                        mmt_my_own_size_in_items(y));
            }

            is_zero = A->vec_is_zero(A,
                    mmt_my_own_subvec(y), mmt_my_own_size_in_items(y));
            pi_allreduce(NULL, &is_zero, 1, BWC_PI_INT, BWC_PI_MIN, pi->m);

            if (is_zero) {
                if (tcan_print) {
                    if (rhscoeffs_name) {
                        printf("M^%u * V + R is zero\n", i);
                        // printf("M^%u * V + R is zero [" K_FILE_BASE_PATTERN ".%u contains M^%u * V, R.sols%u-%u contains R]!\n", i, sf->s0, sf->s1, i-1, i-1, sf->s0, sf->s1);
                    } else {
                        printf("M^%u * V is zero [" K_FILE_BASE_PATTERN ".%u contains M^%u * V]!\n", i, sf->s0, sf->s1, i-1, i-1);
                    }
                }
                break;
            }
            if (rhscoeffs_name) break;

            mmt_vec_unapply_T(mmt, y);
            mmt_vec_save(y, kprefix, i, unpadded);
            mmt_vec_apply_T(mmt, y);
            mmt_vec_twist(mmt, y);
        }
        if (!is_zero) {
            if (tcan_print) {
                printf("Solution range %u..%u: no solution found, most probably a bug\n", sf->s0, sf->s1);
            }
            pthread_mutex_lock(pi->m->th->m);
            exitcode=1;
            pthread_mutex_unlock(pi->m->th->m);
            return NULL;
        }
        if (rhscoeffs_file) {   /* leader node only ! */
            /* We now append the rhs coefficients to the vector we
             * have just saved. Of course only the I/O thread does this. */
            char * tmp;
            int rc = asprintf(&tmp, "%s.%u", kprefix, 0);
            ASSERT_ALWAYS(rc >= 0);
            printf("Expanding %s so as to include the coefficients for the %d RHS columns\n", tmp, nrhs);

            FILE * f = fopen(tmp, "ab");
            ASSERT_ALWAYS(f);
            rc = fseek(f, 0, SEEK_END);
            ASSERT_ALWAYS(rc >= 0);
            rc = fwrite(rhscoeffs, A->vec_elt_stride(A, 1), nrhs, f);
            ASSERT_ALWAYS(rc == (int) nrhs);
            fclose(f);
            printf("%s is now a right nullspace vector for (M|RHS) (for a square M of dimension %" PRIu32"x%" PRIu32").\n", tmp, unpadded, unpadded);

            char * tmp2;
            rc = asprintf(&tmp2, "%s.%u.truncated.txt", kprefix, 0);
            ASSERT_ALWAYS(rc >= 0);

            f = fopen(tmp, "rb");

            FILE *f2 = fopen(tmp2, "w");
            ASSERT_ALWAYS(f2);
            void * data = malloc( A->vec_elt_stride(A, 1));
            for(uint32_t i = 0 ; i < mmt->n0[bw->dir] ; i++) {
                size_t rc = fread(data, A->vec_elt_stride(A, 1), 1, f);
                ASSERT_ALWAYS(rc == 1);
                A->fprint(A, f2, data);
                fprintf(f2, "\n");
            }
            free(data);
            for(uint32_t i = 0 ; i < nrhs ; i++) {
                A->fprint(A, f2, A->vec_coeff_ptr(A, rhscoeffs, i));
                fprintf(f2, "\n");
            }
            fclose(f2);
            printf("%s (in ascii) is now a right nullspace vector for (M|RHS) (for the original M, of dimension %" PRIu32"x%" PRIu32").\n", tmp2, mmt->n0[0], mmt->n0[1]);
            free(tmp2);
            free(tmp);
        }

        serialize(pi->m);
        free(kprefix);
        if (rhscoeffs_name) cheating_vec_clear(A, &rhscoeffs, nrhs);
    }

    if (rhscoeffs_file) fclose(rhscoeffs_file);
    if (rhscoeffs_name) mmt_vec_clear(mmt, rhscoeffs_vec);

    cheating_vec_clear(A, &svec, eblock);
    cheating_vec_clear(A, &tvec, eblock);

    mmt_vec_clear(mmt, y);
    mmt_vec_clear(mmt, my);

    matmul_top_clear(mmt);

    A->oo_field_clear(A);

    for(int i = 0 ; i < sl->nsols ; i++) {
        free(sl->sols[i]->sfiles);
    }
    free(sl->sols);

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
    param_list_decl_usage(pl, "rhs",
            "file with the right-hand side vectors for inhomogeneous systems mod p (only the header is read by this program)");
    param_list_decl_usage(pl, "rhscoeffs",
            "for the solution vector(s), this corresponds to the contribution(s) on the columns concerned by the rhs");

    bw_common_parse_cmdline(bw, pl, &argc, &argv);

    bw_common_interpret_parameters(bw, pl);
    parallelizing_info_lookup_parameters(pl);
    matmul_top_lookup_parameters(pl);
    /* interpret our parameters: none here (so far). */
    param_list_lookup_string(pl, "rhs");
    param_list_lookup_string(pl, "rhscoeffs");

    if (param_list_warn_unused(pl)) {
        param_list_print_usage(pl, argv[0], stderr);
        exit(EXIT_FAILURE);
    }

    pi_go(gather_prog, pl, 0);

    parallelizing_info_finish();
    param_list_clear(pl);
    bw_common_clear_new(bw);

    return exitcode;
}

