#include "cado.h"

#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <limits.h>
#include <dirent.h>
#include <errno.h>
#include <vector>
#include <algorithm>
#include "bwc_config.h"
#include "parallelizing_info.h"
#include "matmul_top.h"
#include "select_mpi.h"

#include "params.h"
#include "xvectors.h"
#include "portability.h"
#include "misc.h"
#include "bw-common.h"
#include "balancing.h"
#include "mpfq/mpfq.h"
#include "mpfq/mpfq_vbase.h"
#include "cheating_vec_init.h"

using namespace std;

struct sfile_info {
    unsigned int iter0;
    unsigned int iter1;
    char name[NAME_MAX];
    static bool match(sfile_info& v, const char * name, unsigned int sol0, unsigned int sol1) {
        int k;
        int rc;
        unsigned int s0,s1;
        unsigned int iter0, iter1;
        rc = sscanf(name, "S.sols%u-%u.%u-%u%n", &s0, &s1, &iter0, &iter1, &k);
        if (rc < 4 || k != (int) strlen(name))
            return false;
        if (s0 != sol0 || s1 != sol1)
            return false;
        v.iter0 = iter0;
        v.iter1 = iter1;
        size_t res = strlcpy(v.name, name, NAME_MAX);
        ASSERT_ALWAYS(res < NAME_MAX);
        return true;
    }
};

static bool operator<(sfile_info const & a, sfile_info const & b)
{
    return a.iter1 < b.iter1;
}

int exitcode = 0;

vector<sfile_info> prelude(parallelizing_info_ptr pi)
{
    vector<sfile_info> res;
    serialize_threads(pi->m);
    if (pi->m->jrank == 0 && pi->m->trank == 0) {
        /* It's our job to collect the directory data.  */
        DIR * dir = opendir(".");
        struct dirent * de;
        for( ; (de = readdir(dir)) != NULL ; ) {
            sfile_info sf;
            if (!sfile_info::match(sf, de->d_name, bw->solutions[0], bw->solutions[1]))
                continue;
            res.push_back(sf);
            if (bw->interval && sf.iter1 % bw->interval != 0) {
                fprintf(stderr,
                        "Warning: %s is not a checkpoint at a multiple of "
                        "the interval value %d -- this might indicate a "
                        "severe bug with leftover data, likely to corrupt "
                        "the final computation\n", sf.name, bw->interval);
            }
        }
        closedir(dir);

        sort(res.begin(), res.end());

        unsigned int prev_iter = 0;
        for(size_t i = 0 ; i < res.size() ; i++) {
            sfile_info const & cur(res[i]);
            if (cur.iter0 != prev_iter) {
                fprintf(stderr, "Within the set of S files, file "
                        "%s seems to be the first "
                        "to come after iteration %u, therefore there is "
                        "a gap for the range %u..%u\n",
                        cur.name, prev_iter, prev_iter, cur.iter0);
                exit(EXIT_FAILURE);
            }
            prev_iter = cur.iter1;
        }
    }

    serialize(pi->m);
    unsigned long s = res.size();
    pi_bcast(&s, 1, BWC_PI_UNSIGNED_LONG, 0, 0, pi->m);
    if (pi->m->jrank || pi->m->trank) {
        res.assign(s, sfile_info());
    }
    pi_bcast(&(res[0]), s * sizeof(sfile_info), BWC_PI_BYTE, 0, 0, pi->m);
    serialize(pi->m);
    return res;
}

void * gather_prog(parallelizing_info_ptr pi, param_list pl, void * arg MAYBE_UNUSED)
{
    ASSERT_ALWAYS(!pi->interleaved);

    int tcan_print = bw->can_print && pi->m->trank == 0;
    matmul_top_data mmt;

    unsigned int solutions[2] = { bw->solutions[0], bw->solutions[1], };
    int char2 = mpz_cmp_ui(bw->p, 2) == 0;
    int splitwidth = char2 ? 64 : 1;

    /* Define and initialize our arithmetic back-ends. More or less the
     * same deal as for mksol (for the "solutions" part). See comments
     * there. */

    /* {{{ main arithmetic backend: for the solutions we compute. */
    // unsigned int A_multiplex MAYBE_UNUSED = 1;
    unsigned int A_width = solutions[1]-solutions[0];
    if ((char2 && (A_width != 64 && A_width != 128 && A_width != 256))
            || (!char2 && A_width > 1))
    {
        fprintf(stderr,
                "We cannot support computing %u solutions at a time "
                "with one single Spmv operation, given the currently "
                "implemented code\n",
                A_width);
        exit(EXIT_FAILURE);
    }
    mpfq_vbase A;
    mpfq_vbase_oo_field_init_byfeatures(A,
            MPFQ_PRIME_MPZ, bw->p,
            MPFQ_GROUPSIZE, A_width,
            MPFQ_DONE);

    /* }}} */

    /* {{{ This arithmetic back-end will be used only for reading the rhs
     * vectors from the V files, and do an addmul later on.
     */
    mpfq_vbase Av;
    unsigned int Av_width = splitwidth;
    mpfq_vbase_oo_field_init_byfeatures(Av,
            MPFQ_PRIME_MPZ, bw->p,
            MPFQ_GROUPSIZE, Av_width,
            MPFQ_DONE);
    pi_datatype_ptr Av_pi = pi_alloc_mpfq_datatype(pi, Av);
    /* }}} */

    /* {{{ ... and the combined operations */
    mpfq_vbase_tmpl AvxA;
    mpfq_vbase_oo_init_templates(AvxA, Av, A);
    /* }}} */

    matmul_top_init(mmt, A, pi, pl, bw->dir);

    mmt_vec ymy[2];
    mmt_vec_ptr y = ymy[0];
    mmt_vec_ptr my = ymy[1];

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

    vector<sfile_info> sl = prelude(pi);
    if (sl.empty()) {
        if (tcan_print) {
            fprintf(stderr, "Found zero S files for solution range %u..%u. "
                    "Problem with command line ?\n",
                    solutions[0], solutions[1]);
            pthread_mutex_lock(pi->m->th->m);
            exitcode=1;
            pthread_mutex_unlock(pi->m->th->m);
        }
        serialize_threads(pi->m);
        return NULL;
    }

    size_t eblock = mmt_my_own_size_in_items(y);

    const char * rhs_name = param_list_lookup_string(pl, "rhs");
    void * rhscoeffs = NULL;

    if (rhs_name) {
        ASSERT_ALWAYS(!char2);
        ASSERT_ALWAYS(A_width == 1);
        if (pi->m->jrank == 0 && pi->m->trank == 0)
            get_rhs_file_header(rhs_name, NULL, &nrhs, NULL);
        pi_bcast(&nrhs, 1, BWC_PI_UNSIGNED, 0, 0, pi->m);

        if (tcan_print) {
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
        cheating_vec_init(A, &rhscoeffs, nrhs);
    }

    char * kprefix;

    if (tcan_print)
        printf("Trying to build solutions %u..%u\n", solutions[0], solutions[1]);

    int rc = asprintf(&kprefix, "K.sols%u-%u", solutions[0], solutions[1]);
    ASSERT_ALWAYS(rc >= 0);

    mmt_full_vec_set_zero(my);

    mmt_full_vec_set_zero(y);

    { /* {{{ Collect now the sum of the LHS contributions */
        mmt_vec svec;
        mmt_vec_init(mmt,0,0, svec,bw->dir, /* shared ! */ 1, mmt->n[bw->dir]);
        for(size_t i = 0 ; i < sl.size() ; i++) {
            if (tcan_print && verbose_enabled(CADO_VERBOSE_PRINT_BWC_LOADING_MKSOL_FILES)) {
                printf("loading %s\n", sl[i].name);
            }
            mmt_vec_load(svec, sl[i].name, unpadded);

            A->vec_add(A,
                    mmt_my_own_subvec(y), 
                    mmt_my_own_subvec(y),
                    mmt_my_own_subvec(svec),
                    eblock);
        }
        y->consistency = 1;
        mmt_vec_broadcast(y);
        mmt_vec_clear(mmt, svec);
    } /* }}} */

    /* Note that for the inhomogeneous case, we'll do the loop only
     * once, since we end with a break. */
    int winning_iter = 0;

    /* This is used as an out-of-loop check for an exceptional condition.
     * This is slighly clumsy, I should do it better */
    int is_zero;

    for(int i = 1 ; nrhs ? (i<=1) : (i < 10) ; i++) {

        mmt_vec_unapply_T(mmt, y);
        /* Here, we want to make sure that we have something non-zero in
         * the **input coordinate space**. It matters little to us if we
         * found a solution of (M||zero-pad) * x = 0 with x having
         * non-zero coordinates only in the zero-pad part.
         *
         * Therefore, the following check is unfortunately not good
         * enough:
         
        int is_zero = A->vec_is_zero(A,
                mmt_my_own_subvec(y), mmt_my_own_size_in_items(y));

         * instead, we want to check only up to index mmt->n0[bw->dir]
         * (and bw->dir is y->d). (this is valid because at this point, y
         * is untwisted and has T unapplied).
         */
        size_t my_input_coordinates;
        size_t my_pad_coordinates;
        ASSERT_ALWAYS(y->d == bw->dir);
        size_t my_i0 = y->i0 + mmt_my_own_offset_in_items(y);
        size_t my_i1 = my_i0 + mmt_my_own_size_in_items(y);

        if (my_i0 >= mmt->n0[bw->dir]) {
            my_input_coordinates = 0;
            my_pad_coordinates =  mmt_my_own_size_in_items(y);
        } else if (my_i1 >= mmt->n0[bw->dir]) {
            my_input_coordinates = mmt->n0[bw->dir] - my_i0;
            my_pad_coordinates =  my_i1 - mmt->n0[bw->dir];
        } else {
            my_input_coordinates = mmt_my_own_size_in_items(y);
            my_pad_coordinates = 0;
        }
        int input_is_zero = A->vec_is_zero(A,
                mmt_my_own_subvec(y),
                my_input_coordinates);
        int pad_is_zero = A->vec_is_zero(A,
                A->vec_subvec(A, mmt_my_own_subvec(y), my_input_coordinates),
                my_pad_coordinates);

        pi_allreduce(NULL, &input_is_zero, 1, BWC_PI_INT, BWC_PI_MIN, pi->m);
        pi_allreduce(NULL, &pad_is_zero, 1, BWC_PI_INT, BWC_PI_MIN, pi->m);

        if (input_is_zero) {
            if (tcan_print)
                fprintf(stderr, "Found zero vector. (coordinates on the padding part are %s zero). Most certainly a bug. "
                        "No solution found.\n",
                        pad_is_zero ? "also" : "NOT");
            serialize(pi->m);
            pthread_mutex_lock(pi->m->th->m);
            exitcode=1;
            pthread_mutex_unlock(pi->m->th->m);
            return NULL;
        }
        {
            char * tmp;
            int rc = asprintf(&tmp, "%s.%d", kprefix, i-1);
            ASSERT_ALWAYS(rc >= 0);
            mmt_vec_save(y, tmp, unpadded);
            free(tmp);
        }
        mmt_vec_apply_T(mmt, y);


        serialize(pi->m);
        mmt_vec_twist(mmt, y);
        matmul_top_mul(mmt, ymy, NULL);
        mmt_vec_untwist(mmt, y);

        /* Add the contributions from the right-hand side vectors, to see
         * whether that makes the sum equal to zero */

        if (nrhs) {
            mmt_vec vi;
            
            mmt_vec_init(mmt,Av,Av_pi, vi,bw->dir, /* shared ! */ 1, mmt->n[bw->dir]);
            for(unsigned int j = 0 ; j < nrhs ; j++) {
                char * tmp;
                int rc = asprintf(&tmp, "V%u-%u.0", j, j + 1);
                ASSERT_ALWAYS(rc >= 0);
                if (tcan_print && verbose_enabled(CADO_VERBOSE_PRINT_BWC_LOADING_MKSOL_FILES)) {
                    printf("loading %s\n", tmp);
                }
                mmt_vec_load(vi, tmp, unpadded);
                free(tmp);

                if (pi->m->trank == 0 && pi->m->jrank == 0) {
                    rc = asprintf(&tmp, "F.sols%u-%u.%u-%u.rhs", solutions[0], solutions[1], j, j+Av_width);
                    if (verbose_enabled(CADO_VERBOSE_PRINT_BWC_LOADING_MKSOL_FILES)) {
                        printf("loading %s\n", tmp);
                    }
                    ASSERT_ALWAYS(rc >= 0);
                    FILE * f = fopen(tmp, "rb");
                    rc = fread(
                            A->vec_subvec(A, rhscoeffs, j),
                            A->vec_elt_stride(A,1),
                            1, f);
                    ASSERT_ALWAYS(rc == 1);
                    if (A->is_zero(A, A->vec_coeff_ptr_const(A, rhscoeffs, j))) {
                        printf("Notice: coefficient for vector V%u-%u in file %s is zero\n", j, j+1, tmp);
                    }
                    fclose(f);
                    free(tmp);
                }
                pi_bcast(A->vec_subvec(A, rhscoeffs, j), 1, mmt->pitype, 0, 0, pi->m);

                AvxA->addmul_tiny(Av, A, 
                                mmt_my_own_subvec(y),
                                mmt_my_own_subvec(vi),
                                A->vec_subvec(A, rhscoeffs, j),
                                eblock);
            }
            y->consistency = 1;
            mmt_vec_broadcast(y);
            mmt_vec_reduce_mod_p(y);
            mmt_vec_clear(mmt, vi);
        }

        /* This "is zero" check is also valid on the padded matrix of
         * course, so we don't heave the same headache as above */
        is_zero = A->vec_is_zero(A,
                mmt_my_own_subvec(y), mmt_my_own_size_in_items(y));
        pi_allreduce(NULL, &is_zero, 1, BWC_PI_INT, BWC_PI_MIN, pi->m);

        if (tcan_print) {
            if (nrhs) {
                const char * strings[] = {
                    "V * M^%u + R is %s\n", /* unsupported anyway */
                    "M^%u * V + R is %s\n",};
                printf(strings[bw->dir], i, is_zero ? "zero" : "NOT zero");
            } else {
                const char * strings[] = {
                    "V * M^%u is %s [K.sols%u-%u.%u contains V * M^%u]!\n",
                    "M^%u * V is %s [K.sols%u-%u.%u contains M^%u * V]!\n",
                };
                printf(strings[bw->dir], i,
                        is_zero ? "zero" : "NOT zero",
                        solutions[0], solutions[1], i-1, i-1);
            }
        }
        if (is_zero) {
            winning_iter = i-1;
            break;
        }
    }

    /* we exit with an untwisted vector. */

    if (!is_zero) {
        int nb_nonzero_coeffs=0;
        for(unsigned int i = 0 ; i < mmt_my_own_size_in_items(y) ; i++) {
            nb_nonzero_coeffs += !
                A->vec_is_zero(A,
                        A->vec_subvec(A, mmt_my_own_subvec(y), i),
                        1);
        }
        pi_allreduce(NULL, &nb_nonzero_coeffs, 1, BWC_PI_INT, BWC_PI_SUM, pi->m);
        if (tcan_print) {
            printf("Solution range %u..%u: no solution found [%d non zero coefficients in result], most probably a bug\n", solutions[0], solutions[1], nb_nonzero_coeffs);
            if (nb_nonzero_coeffs < bw->n) {
                printf("There is some likelihood that by combining %d different solutions (each entailing a separate mksol run), a full solution can be recovered. Ask for support.\n", nb_nonzero_coeffs + 1);
            }
        }
        pthread_mutex_lock(pi->m->th->m);
        exitcode=1;
        pthread_mutex_unlock(pi->m->th->m);
        return NULL;
    }

    if (pi->m->jrank == 0 && pi->m->trank == 0) {   /* leader node only ! */
        /* About the ASSERT below: The 0-characteristic subspace
         * nilpotent part is not really a topic of interest in the case
         * of inhomogenous linear systems.  Exiting the loop with i>1 in
         * the inhomogenous case is an error, it seems, and I believe
         * that the solution provided was never really a solution. On the
         * contrary, the i>1 case in homogenous situations is perfectly
         * valid, and we're handling it now instead of using
         * winning_iter==0 always.
         */
        ASSERT_ALWAYS(winning_iter == 0 || rhs_name == NULL);

        /* We now append the rhs coefficients to the vector we
         * have just saved. Of course only the I/O thread does this. */
        char * tmp;
        int rc = asprintf(&tmp, "%s.%u", kprefix, winning_iter);
        ASSERT_ALWAYS(rc >= 0);
        if (nrhs)
        printf("Expanding %s so as to include the coefficients for the %d RHS columns\n", tmp, nrhs);

        FILE * f = fopen(tmp, "ab");
        ASSERT_ALWAYS(f);
        rc = fseek(f, 0, SEEK_END);
        ASSERT_ALWAYS(rc >= 0);
        rc = fwrite(rhscoeffs, A->vec_elt_stride(A, 1), nrhs, f);
        ASSERT_ALWAYS(rc == (int) nrhs);
        fclose(f);
        printf("%s is now a %s nullspace vector for %s (for a square M of dimension %" PRIu32"x%" PRIu32").\n", tmp, bw_dirtext[bw->dir], nrhs ? "(M|RHS)" : "M", unpadded, unpadded);


        /* write an ascii version while we're at it */
        char * tmp2;
        rc = asprintf(&tmp2, "%s.%u.txt", kprefix, winning_iter);
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
        printf("%s (in ascii) is now a %s nullspace vector for %s (for the original M, of dimension %" PRIu32"x%" PRIu32").\n", tmp2, bw_dirtext[bw->dir], nrhs ? "(M|RHS)" : "M", mmt->n0[0], mmt->n0[1]);
        free(tmp2);
        free(tmp);
    }

    serialize(pi->m);
    free(kprefix);


    if (rhs_name)
        cheating_vec_clear(A, &rhscoeffs, nrhs);
    mmt_vec_clear(mmt, y);
    mmt_vec_clear(mmt, my);

    matmul_top_clear(mmt);

    pi_free_mpfq_datatype(pi, Av_pi);
    A->oo_field_clear(A);
    Av->oo_field_clear(Av);

    return NULL;
}


int main(int argc, char * argv[])
{
    param_list pl;

    bw_common_init(bw, &argc, &argv);
    param_list_init(pl);
    parallelizing_info_init();

    bw_common_decl_usage(pl);
    parallelizing_info_decl_usage(pl);
    matmul_top_decl_usage(pl);
    /* declare local parameters and switches */
    param_list_decl_usage(pl, "rhs",
            "file with the right-hand side vectors for inhomogeneous systems mod p (only the header is read by this program, while the actual contents are recovered from the V*.0 files)");
    param_list_decl_usage(pl, "rhscoeffs",
            "for the solution vector(s), this corresponds to the contribution(s) on the columns concerned by the rhs");

    bw_common_parse_cmdline(bw, pl, &argc, &argv);

    param_list_remove_key(pl, "interleaving");

    bw_common_interpret_parameters(bw, pl);
    parallelizing_info_lookup_parameters(pl);
    matmul_top_lookup_parameters(pl);
    /* interpret our parameters */

    ASSERT_ALWAYS(!param_list_lookup_string(pl, "ys"));
    ASSERT_ALWAYS(param_list_lookup_string(pl, "solutions"));

    param_list_lookup_string(pl, "rhs");
    param_list_lookup_string(pl, "rhscoeffs");

    if (param_list_warn_unused(pl)) {
        param_list_print_usage(pl, argv[0], stderr);
        exit(EXIT_FAILURE);
    }

    pi_go(gather_prog, pl, 0);

    parallelizing_info_finish();
    param_list_clear(pl);
    bw_common_clear(bw);

    return exitcode;
}

