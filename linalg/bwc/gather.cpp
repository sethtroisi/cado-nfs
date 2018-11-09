#include "cado.h"

#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <limits.h>
#include <dirent.h>
#include <errno.h>
#include <vector>
#include <set>
#include <map>
#include <tuple>
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

struct sfile_info {/*{{{*/
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
};/*}}}*/

static bool operator<(sfile_info const & a, sfile_info const & b)/*{{{*/
{
    return a.iter1 < b.iter1;
}/*}}}*/

int exitcode = 0;

vector<sfile_info> prelude(parallelizing_info_ptr pi)/*{{{*/
{
    int leader = pi->m->jrank == 0 && pi->m->trank == 0;
    vector<sfile_info> res;
    serialize_threads(pi->m);
    if (leader) {
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
    if (!leader)
        res.assign(s, sfile_info());
    pi_bcast(&(res[0]), s * sizeof(sfile_info), BWC_PI_BYTE, 0, 0, pi->m);
    serialize(pi->m);
    return res;
}/*}}}*/

void fprint_signed(FILE * f, mpfq_vbase_ptr A, void * x)
{
    char * ss;
    void * minus;
    A->vec_init(A, &minus, 1);
    A->neg(A, minus, x);
    if (A->cmp(A, minus, x) < 0) {
        A->asprint(A, &ss, minus);
        fprintf(f, "-%s", ss);
    } else {
        A->asprint(A, &ss, x);
        fprintf(f, "%s", ss);
    }
    A->vec_clear(A, &minus, 1);
    free(ss);
}

void allgather(std::vector<unsigned int>& v, pi_comm_ptr wr)/*{{{*/
{
    /* want to collectively merge all vectors "v". boost mpi
     * would be great for that, really */

    /* first merge to one leader per node */

    std::vector<unsigned int> *mainv;
    mainv = &v;
    pi_thread_bcast(&mainv, sizeof(mainv), BWC_PI_BYTE, 0, wr);

    for(unsigned int j = 1 ; j < wr->ncores ; ++j) {
        if (wr->trank == j)
            mainv->insert(mainv->end(), v.begin(), v.end());
        serialize_threads(wr);
    }

    if (wr->trank == 0) {
        std::vector<unsigned int> allv;
        std::vector<int> sizes(wr->njobs, 0);
        std::vector<int> displs(wr->njobs, 0);
        sizes[wr->jrank] = v.size();
        int total = v.size();
        MPI_Allreduce(MPI_IN_PLACE, &total, 1, MPI_INT, MPI_SUM,
                wr->pals);
        MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
                &sizes[0], 1, MPI_INT, wr->pals);
        for(unsigned int i = 1 ; i < wr->njobs ; i++)
            displs[i] = displs[i-1] + sizes[i-1];
        allv.assign(total, 0);
        std::copy(v.begin(), v.end(), allv.begin() + displs[wr->jrank]);
        MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
                &allv[0], &sizes[0], &displs[0], MPI_UNSIGNED,
                wr->pals);
        std::swap(v, allv);
    }

    serialize_threads(wr);
    v = *mainv;
    serialize_threads(wr);
}/*}}}*/

void broadcast(std::vector<unsigned int>& v, pi_comm_ptr wr)/*{{{*/
{
    int total = v.size();
    pi_bcast(&total, 1, BWC_PI_INT, 0, 0, wr);
    if (wr->jrank || wr->trank) v.assign(total, 0);
    pi_bcast(&v[0], total, BWC_PI_UNSIGNED, 0, 0, wr);
    serialize(wr);
}/*}}}*/

std::vector<unsigned int> indices_of_zero_or_nonzero_values(mmt_vec_ptr y, unsigned int maxidx, int want_nonzero)/*{{{*/
{
    mpfq_vbase_ptr A = y->abase;
    parallelizing_info_ptr pi = y->pi;

    std::vector<unsigned int> myz;

    if (pi->wr[y->d]->trank == 0 && pi->wr[y->d]->jrank == 0) {
        for(unsigned int i = 0 ; i < maxidx ; i++) {
            if (y->i0 <= i && i < y->i1) {
                if (!!want_nonzero == !A->is_zero(A, A->vec_coeff_ptr_const(A, y->v, i - y->i0))) {
                    myz.push_back(i);
                }
            }
        }

        /* in fact, a single gather at node 0 thread 0 would do */
        allgather(myz, pi->wr[!y->d]);
    }

    broadcast(myz, pi->m);    /* And broadcast that to everyone as well. */

    return myz;
}/*}}}*/

std::vector<unsigned int> indices_of_zero_values(mmt_vec_ptr y, unsigned int maxidx)/*{{{*/
{
    return indices_of_zero_or_nonzero_values(y, maxidx, 0);
}/*}}}*/

std::vector<unsigned int> indices_of_nonzero_values(mmt_vec_ptr y, unsigned int maxidx)/*{{{*/
{
    return indices_of_zero_or_nonzero_values(y, maxidx, 1);
}/*}}}*/

std::vector<unsigned int> get_possibly_wrong_columns(matmul_top_data_ptr mmt)/*{{{*/
{
    parallelizing_info_ptr pi = mmt->pi;
    int tcan_print = bw->can_print && pi->m->trank == 0;

    std::vector<unsigned int> allz;

    if (mmt->n0[bw->dir] >= mmt->n0[!bw->dir]) return allz;

    /*
       if (tcan_print) {
       const char * name[2] = { "row", "column" };
       printf("// The original matrix has more %ss than %ss.\n// With nullspace=%s,this carries a slight danger that some parasite elements of a degree 2 nilpotent nullspace mask the real solution.\n// Trying to detect this situation.\n", name[!bw->dir], name[bw->dir], bw_dirtext[bw->dir]);
       }
        */

    /* Comments below assume that we're in the typical case
     * bw->dir==1, nrows > ncols */

    gmp_randstate_t rstate;
    gmp_randinit_default(rstate);

    if (pi->m->trank == 0 && !bw->seed) {
        /* note that bw is shared between threads, thus only thread 0 should
         * test and update it here.
         * at pi->m->jrank > 0, we don't care about the seed anyway
         */
        bw->seed = time(NULL);
        MPI_Bcast(&bw->seed, 1, MPI_INT, 0, pi->m->pals);
    }
    serialize_threads(pi->m);

    gmp_randseed_ui(rstate, bw->seed);
    if (tcan_print)
        printf("// Random generator seeded with %d\n", bw->seed);

    /* Do that in the opposite direction compared to ymy */
    mmt_vec zmz[2];
    mmt_vec_ptr z = zmz[0];
    mmt_vec_ptr mz = zmz[1];

    mmt_vec_init(mmt,0,0, z,  !bw->dir, /* shared ! */ 1, mmt->n[!bw->dir]);
    mmt_vec_init(mmt,0,0, mz,  bw->dir,                0, mmt->n[bw->dir]);


    mmt_vec_set_random_inconsistent(z, rstate);
    mmt_vec_truncate_above_index(mmt, z, mmt->n0[bw->dir]);
    mmt_vec_apply_T(mmt, z);
    mmt_vec_twist(mmt, z);
    matmul_top_mul(mmt, zmz, NULL);
    mmt_vec_untwist(mmt, z);
    mmt_apply_identity(mz, z);
    mmt_vec_allreduce(mz);
    mmt_vec_unapply_T(mmt, mz);
    serialize(pi->m);

    /* Check the indices of columns in the principal part. We're
     * interested in column indices which are still zero. So it's
     * really a loop until mmt->n0[bw->dir] */

    ASSERT_ALWAYS(mz->d == bw->dir);
    allz = indices_of_zero_values(mz, mmt->n0[bw->dir]);


    /* Do a second check with the *FULL* vector. Coordinates that are
     * still zero are no reason to be worried */

    mmt_vec_set_random_inconsistent(z, rstate);
    mmt_vec_truncate_above_index(mmt, z, mmt->n0[!bw->dir]);
    mmt_vec_apply_T(mmt, z);
    mmt_vec_twist(mmt, z);
    matmul_top_mul(mmt, zmz, NULL);
    mmt_vec_untwist(mmt, z);
    mmt_apply_identity(mz, z);
    mmt_vec_allreduce(mz);
    mmt_vec_unapply_T(mmt, mz);
    serialize(pi->m);

    ASSERT_ALWAYS(mz->d == bw->dir);
    std::set<unsigned int> allz_set(allz.begin(), allz.end());
    for(auto j : indices_of_zero_values(mz, mmt->n0[bw->dir])) {
        allz_set.erase(j);
    }
    allz.assign(allz_set.begin(), allz_set.end());

    mmt_vec_clear(mmt, z);
    mmt_vec_clear(mmt, mz);

    gmp_randclear(rstate);

    return allz;
}/*}}}*/

/* Take a vector, a priori sparse, and pick its coefficients at indices
 * marked by [rows], and put them, in order, in the j-th column of the
 * matrix pointed to by matrix, which has cblocks blocks of
 * A->simd_groupsize(A) entries. The column number j is thus made of elements
 * of the blocks whose index is congruent to
 * (j/groupsize)-th block mod cblocks ; A is my->abase.
 */
void compress_vector_to_sparse(void * matrix, unsigned int j, unsigned int cblocks, mmt_vec_ptr my, std::vector<unsigned int> & rows)
{
    mpfq_vbase_ptr A = my->abase;

    unsigned int own_i0 = my->i0 + mmt_my_own_offset_in_items(my);
    unsigned int own_i1 = own_i0 + mmt_my_own_size_in_items(my);
    cxx_mpz v;
    unsigned int jq = j / A->simd_groupsize(A);
    // unsigned int jr = j % A->simd_groupsize(A);

    int char2 = mpz_cmp_ui(bw->p, 2) == 0;
    ASSERT_ALWAYS(char2 || A->simd_groupsize(A) == 1);

    for(unsigned int ii = 0 ; ii < rows.size() ; ii++) {
        unsigned int i = rows[ii];
        if (own_i0 <= i && i < own_i1) {
            const void * src = A->vec_coeff_ptr_const(A, my->v, i - my->i0);
            void * dst = A->vec_coeff_ptr(A, matrix, ii * cblocks + jq);
            A->set(A, dst, src);
        }
    }
}

struct abase_proxy {

    parallelizing_info_ptr pi;
    mpfq_vbase A;
    pi_datatype_ptr A_pi;

    abase_proxy(parallelizing_info_ptr pi, int width) : pi(pi) {
        mpfq_vbase_oo_field_init_byfeatures(A,
                MPFQ_PRIME_MPZ, bw->p,
                MPFQ_SIMD_GROUPSIZE, width,
                MPFQ_DONE);
        A_pi = pi_alloc_mpfq_datatype(pi, A);
    }
    static abase_proxy most_natural(parallelizing_info_ptr pi) {
        return abase_proxy(pi, mpz_cmp_ui(bw->p, 2) == 0 ? 64 : 1);
    }
    std::map<mpfq_vbase_ptr, mpfq_vbase_tmpl> tdict;
    mpfq_vbase_tmpl_ptr templates(mpfq_vbase_ptr A1) {
        auto it = tdict.find(A1);
        if (it == tdict.end())
            mpfq_vbase_oo_init_templates(tdict[A1], A, A1);
        return tdict[A1];
    }
    ~abase_proxy()
    {
        pi_free_mpfq_datatype(pi, A_pi);
        A->oo_field_clear(A);
    }
};

class parasite_fixer {/*{{{*/
    matmul_top_data_ptr mmt;
    mpfq_vbase_ptr A;
    parallelizing_info_ptr pi;

    typedef std::map<std::pair<unsigned int, unsigned int>, cxx_mpz> pre_matrix_t;

    public:
    bool attempt_to_fix;

    std::vector<unsigned int> cols;
    std::vector<unsigned int> rows;
    void * matrix;

    parasite_fixer(matmul_top_data_ptr mmt) : mmt(mmt), A(mmt->abase), pi(mmt->pi) {/*{{{*/
        matrix = NULL;

        int tcan_print = bw->can_print && pi->m->trank == 0;

        cols = get_possibly_wrong_columns(mmt);

        attempt_to_fix = !cols.empty() && cols.size() < 64;

        if (!cols.empty() && tcan_print) {
            printf("# %zu possibly wrong coordinates detected"
                    " in solution vector because the matrix"
                    " has a non-trivial nilpotent space.\n",
                    cols.size());
            if (attempt_to_fix) {
                printf("# Will try to fix\n");
            } else {
                printf("# Deliberately avoiding attempt to fix because"
                        " this number of vectors is large."
                        " It is rather an indication that something is wrong."
                        " Proceeding anyway\n");
            }
        }
        compute_matrix();
    }/*}}}*/

    parasite_fixer(parasite_fixer const&) = delete;

    /*{{{ row_coordinates_of_nonzero_cols(matmul_top_data_ptr mmt) */
    std::vector<unsigned int> row_coordinates_of_nonzero_cols(matmul_top_data_ptr mmt, std::vector<unsigned int> const& cols)
    {
        int tcan_print = bw->can_print && pi->m->trank == 0;

        mpfq_vbase_ptr A = mmt->abase;

        mmt_vec ymy[2];
        mmt_vec_ptr y = ymy[0];
        mmt_vec_ptr my = ymy[1];

        mmt_vec_init(mmt,0,0, y,   bw->dir, /* shared ! */ 1, mmt->n[bw->dir]);
        mmt_vec_init(mmt,0,0, my, !bw->dir,                0, mmt->n[!bw->dir]);

        /* Now try to see which indices are potentially affected */
        unsigned int B = A->simd_groupsize(A);
        for(unsigned int jjq = 0 ; jjq < cols.size() ; jjq+=B) {
            mmt_full_vec_set_zero(y);
            for(unsigned int jjr = 0 ; jjr < B && (jjq + jjr < cols.size()) ; jjr++) {
                unsigned int j = cols[jjq + jjr];
                mmt_vec_add_basis_vector_at(y, jjr, j);
            }
            mmt_vec_apply_T(mmt, y);
            mmt_vec_twist(mmt, y);
            matmul_top_mul(mmt, ymy, NULL);
            mmt_vec_untwist(mmt, y);
            /* Not entirely clear to me if I should unapply_T here or not
            */
            mmt_apply_identity(my, y);
            mmt_vec_allreduce(my);
            mmt_vec_unapply_T(mmt, my);
            serialize(pi->m);

            /* This reads the global list */
            std::vector<unsigned int> kk = indices_of_nonzero_values(my, mmt->n0[!bw->dir]);

            if (tcan_print) {
                printf("# List of the %zu non-zero coordinates attached to input coordinate(s)", kk.size());
                for(unsigned int jjr = 0 ; jjr < B && (jjq + jjr < cols.size()) ; jjr++) {
                    unsigned int j = cols[jjq + jjr];
                    printf(" %u", j);
                }
                printf("\n");
                for(auto k : kk)
                    printf("#\t%u\n", k);
            }
            rows.insert(rows.end(), kk.begin(), kk.end());
            serialize(pi->m);
        }
        std::sort(rows.begin(), rows.end());
        std::unique(rows.begin(), rows.end());
        mmt_vec_clear(mmt, y);
        mmt_vec_clear(mmt, my);

        return rows;
    }/*}}}*/

    /* this one could live outside the class */
    void pre_matrix_to_matrix(pre_matrix_t /* const */ & pre_matrix)/*{{{*/
    {
        /* we should have a const arg, but an error in the mpfq api
         * prevents this.
         */
        if (pi->wr[!bw->dir]->jrank == 0 && pi->wr[!bw->dir]->trank == 0) {
            unsigned int kk = 0;
            for(unsigned int ii = 0 ; ii < rows.size() ; ii++) {
                unsigned int i = rows[ii];
                for(unsigned int jj = 0 ; jj < cols.size() ; jj++) {
                    unsigned int j = cols[jj];
                    auto it = pre_matrix.find(std::make_pair(i,j));
                    if (it != pre_matrix.end()) {
                        /* FIXME: Api error in mpfq. */
                        A->set_mpz(A, A->vec_coeff_ptr(A, matrix, kk), (mpz_ptr) it->second);
                    }
                    kk++;
                }
            }
        }

        pi_allreduce(NULL, matrix,
                rows.size() * cols.size(),
                mmt->pitype, BWC_PI_SUM, pi->m);
    }/*}}}*/

    void debug_print_local_matrix(void * nz = NULL)/*{{{*/
    {
        printf("# Dump of the full matrix%s as seen by J%uT%u\n",
                nz ? " (with coefficients of the vector encountered)" : "",
                pi->m->jrank, pi->m->trank);
        unsigned int kk = 0;
        unsigned int cblocks = iceildiv(cols.size(), A->simd_groupsize(A));
        for(unsigned int ii = 0 ; ii < rows.size() ; ii++) {
            printf("#\t");
            for(unsigned int jj = 0 ; jj < cblocks ; jj++) {
                printf(" ");
                fprint_signed(stdout, A, A->vec_coeff_ptr(A, matrix, kk));
                kk++;
            }
            if (nz) {
                printf(" ");
                fprint_signed(stdout, A, A->vec_coeff_ptr(A, nz, ii));
            }
            printf("\n");
        }
    }/*}}}*/

    void debug_print_all_local_matrices(void * nz = NULL)/*{{{*/
    {
        for(unsigned int jr = 0 ; jr < pi->m->njobs ; jr++) {
            for(unsigned int tr = 0 ; tr < pi->m->ncores ; tr++) {
                if (pi->m->jrank == jr && pi->m->trank == tr)
                    debug_print_local_matrix(nz);
                serialize(pi->m);
            }
        }
    }/*}}}*/

    void compute_matrix() {/*{{{*/
        if (!attempt_to_fix) return;

        int tcan_print = bw->can_print && pi->m->trank == 0;
        int leader = pi->m->jrank == 0 && pi->m->trank == 0;

        rows = row_coordinates_of_nonzero_cols(mmt, cols);

        if (tcan_print) {
            printf("# Candidate solutions will be checked for accordance with a matrix of dimension %zu*%zu\n", rows.size(), cols.size());
        }

        mpfq_vbase_ptr A = mmt->abase;

        unsigned int cblocks = iceildiv(cols.size(), A->simd_groupsize(A));

        cheating_vec_init(A, &matrix, cblocks * rows.size());
        A->vec_set_zero(A, matrix, cblocks * rows.size());

        /* code is similar to row_coordinates_of_nonzero_cols() */
        mmt_vec ymy[2];
        mmt_vec_ptr y = ymy[0];
        mmt_vec_ptr my = ymy[1];

        mmt_vec_init(mmt,0,0, y,   bw->dir, /* shared ! */ 1, mmt->n[bw->dir]);
        mmt_vec_init(mmt,0,0, my, !bw->dir,                0, mmt->n[!bw->dir]);

        /* Now try to see which indices are potentially affected */
        unsigned int B = A->simd_groupsize(A);
        for(unsigned int jjq = 0 ; jjq < cols.size() ; jjq+=A->simd_groupsize(A)) {
            mmt_full_vec_set_zero(y);
            for(unsigned int jjr = 0 ; jjr < B && (jjq + jjr < cols.size()) ; jjr++) {
                unsigned int j = cols[jjq + jjr];
                mmt_vec_add_basis_vector_at(y, jjr, j);
            }
            mmt_vec_apply_T(mmt, y);
            mmt_vec_twist(mmt, y);
            matmul_top_mul(mmt, ymy, NULL);
            mmt_vec_untwist(mmt, y);
            /* Not entirely clear to me if I should unapply_T here or not
            */
            mmt_apply_identity(my, y);
            mmt_vec_allreduce(my);
            mmt_vec_unapply_T(mmt, my);
            serialize(pi->m);

            compress_vector_to_sparse(matrix, jjq, cblocks, my, rows);

            serialize(pi->m);
        }
        mmt_vec_clear(mmt, y);
        mmt_vec_clear(mmt, my);

        pi_allreduce(NULL, matrix,
                rows.size() * cblocks,
                mmt->pitype, BWC_PI_SUM, pi->m);

        if (leader) debug_print_local_matrix();

        serialize(pi->m);
    }/*}}}*/

    ~parasite_fixer() {/*{{{*/
        if (matrix) cheating_vec_clear(A, &matrix, rows.size() * cols.size());
    }/*}}}*/

    bool attempt(mmt_vec_ptr my, mmt_vec_ptr y)/*{{{*/
    {
        int tcan_print = bw->can_print && pi->m->trank == 0;
        int leader = pi->m->jrank == 0 && pi->m->trank == 0;

        std::vector<unsigned int> nz_pos;

        /* Need to get the indices with respect to !bw->dir...  */
        mmt_apply_identity(my, y);
        mmt_vec_allreduce(my);
        mmt_vec_unapply_T(mmt, my);
        serialize(pi->m);

        nz_pos = indices_of_nonzero_values(my, mmt->n0[!bw->dir]);

        if (!std::includes(rows.begin(), rows.end(), nz_pos.begin(), nz_pos.end())) return false;

        if (tcan_print)
            printf("# Note: all the non-zero coordinates are included in the output of the \"possibly wrong\" columns, which is a good sign\n");

        void * nz;

        cheating_vec_init(A, &nz, rows.size());

        ASSERT_ALWAYS(my->abase == mmt->abase);

        compress_vector_to_sparse(nz, 0, 1, my, rows);

        // debug_print_all_local_matrices(nz);

        if (leader) debug_print_local_matrix(nz);

        cheating_vec_clear(A, &nz, rows.size());

        /* TODO: change that to a "return true" once we've done the fix
         * for good.
         */
        return false;
    }/*}}}*/
};/*}}}*/

std::tuple<int, int> check_zero_and_padding(mmt_vec_ptr y, unsigned int maxidx)/*{{{*/
{

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

    if (my_i0 >= maxidx) {
        my_input_coordinates = 0;
        my_pad_coordinates =  mmt_my_own_size_in_items(y);
    } else if (my_i1 >= maxidx) {
        my_input_coordinates = maxidx - my_i0;
        my_pad_coordinates =  my_i1 - maxidx;
    } else {
        my_input_coordinates = mmt_my_own_size_in_items(y);
        my_pad_coordinates = 0;
    }
    int input_is_zero = y->abase->vec_is_zero(y->abase,
            mmt_my_own_subvec(y),
            my_input_coordinates);
    int pad_is_zero = y->abase->vec_is_zero(y->abase,
            y->abase->vec_subvec(y->abase,
                mmt_my_own_subvec(y), my_input_coordinates),
            my_pad_coordinates);

    pi_allreduce(NULL, &input_is_zero, 1, BWC_PI_INT, BWC_PI_MIN, y->pi->m);
    pi_allreduce(NULL, &pad_is_zero, 1, BWC_PI_INT, BWC_PI_MIN, y->pi->m);

    return std::make_tuple(input_is_zero, pad_is_zero);
}/*}}}*/

struct rhs /*{{{*/ {
    matmul_top_data_ptr mmt;
    mpfq_vbase_ptr A;
    unsigned int nrhs;
    void * rhscoeffs;

    rhs(rhs const&) = delete;

    rhs(matmul_top_data_ptr mmt, const char * rhs_name, unsigned int solutions[2]) : mmt(mmt), A(mmt->abase) /* {{{ */
    {
        nrhs = 0;
        rhscoeffs = NULL;

        if (!rhs_name) return;

        parallelizing_info_ptr pi = mmt->pi;
        int tcan_print = bw->can_print && pi->m->trank == 0;
        int leader = pi->m->jrank == 0 && pi->m->trank == 0;

        int char2 = mpz_cmp_ui(bw->p, 2) == 0;
        ASSERT_ALWAYS(!char2);
        ASSERT_ALWAYS(A->simd_groupsize(A) == 1);

        if (leader)
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

        /* everyone allocates something */
        cheating_vec_init(A, &rhscoeffs, nrhs);

        if (leader) {
            int splitwidth = char2 ? 64 : 1;
            unsigned int Av_width = splitwidth;
            for(unsigned int j = 0 ; j < nrhs ; j++) {
                char * tmp;
                int rc = asprintf(&tmp, "F.sols%u-%u.%u-%u.rhs", solutions[0], solutions[1], j, j+Av_width);
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
        }
        pi_bcast(rhscoeffs, nrhs, mmt->pitype, 0, 0, pi->m);
    }/*}}}*/

    ~rhs() {/*{{{*/
        if (rhscoeffs)
            cheating_vec_clear(A, &rhscoeffs, nrhs);
    }/*}}}*/
    void fwrite_rhs_coeffs(FILE * f) /* {{{ */
    {
        int rc = fwrite(rhscoeffs, A->vec_elt_stride(A, 1), nrhs, f);
        ASSERT_ALWAYS(rc == (int) nrhs);
    }/*}}}*/
    void fprint_rhs_coeffs(FILE * f2) /* {{{ */
    {
        for(uint32_t i = 0 ; i < nrhs ; i++) {
            A->fprint(A, f2, A->vec_coeff_ptr(A, rhscoeffs, i));
            fprintf(f2, "\n");
        }
    }/*}}}*/
    operator bool() const { return nrhs; }

    void add_contribution(mmt_vec_ptr y)/*{{{*/
    {
        if (!nrhs) return;

        parallelizing_info_ptr pi = mmt->pi;
        mpfq_vbase_ptr A = mmt->abase;
        ASSERT_ALWAYS(y->abase == A);
        int tcan_print = bw->can_print && pi->m->trank == 0;
        unsigned int unpadded = MAX(mmt->n0[0], mmt->n0[1]);
        size_t eblock = mmt_my_own_size_in_items(y);

        abase_proxy natural = abase_proxy::most_natural(pi);
        mpfq_vbase_ptr Av = natural.A;
        pi_datatype_ptr Av_pi = natural.A_pi;

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

            natural.templates(A)->addmul_tiny(Av, A, 
                    mmt_my_own_subvec(y),
                    mmt_my_own_subvec(vi),
                    A->vec_subvec(A, rhscoeffs, j),
                    eblock);
        }

        mmt_vec_clear(mmt, vi);

        y->consistency = 1;
        mmt_vec_broadcast(y);
        mmt_vec_reduce_mod_p(y);
    }/*}}}*/
};
/*}}}*/

std::tuple<int, int, int> test_one_vector(matmul_top_data_ptr mmt, mmt_vec * ymy, rhs & R)
{
    mpfq_vbase_ptr A = mmt->abase;
    parallelizing_info_ptr pi = mmt->pi;

    mmt_vec_ptr y = ymy[0];

    int input_is_zero;
    int pad_is_zero;
    int hamming_out = -1;
    std::tie(input_is_zero, pad_is_zero) = check_zero_and_padding(y, mmt->n0[bw->dir]);
    if (!input_is_zero) {
        mmt_vec_apply_T(mmt, y);
        serialize(pi->m);
        mmt_vec_twist(mmt, y);
        matmul_top_mul(mmt, ymy, NULL);
        mmt_vec_untwist(mmt, y);
        mmt_vec_unapply_T(mmt, y);
        serialize(pi->m);
        /* Add the contributions from the right-hand side vectors, to see
         * whether that makes the sum equal to zero */
        R.add_contribution(y);

        /* This "is zero" check is also valid on the padded matrix of
         * course, so we don't heave the same headache as above */
        int is_zero = A->vec_is_zero(A,
                mmt_my_own_subvec(y), mmt_my_own_size_in_items(y));
        pi_allreduce(NULL, &is_zero, 1, BWC_PI_INT, BWC_PI_MIN, pi->m);

        hamming_out = is_zero ? 0 : mmt_vec_hamming_weight(y);
    }
    return std::make_tuple(input_is_zero, pad_is_zero, hamming_out);
}


void * gather_prog(parallelizing_info_ptr pi, param_list pl, void * arg MAYBE_UNUSED)
{
    ASSERT_ALWAYS(!pi->interleaved);

    int tcan_print = bw->can_print && pi->m->trank == 0;
    int leader = pi->m->jrank == 0 && pi->m->trank == 0;

    matmul_top_data mmt;

    unsigned int solutions[2] = { bw->solutions[0], bw->solutions[1], };
    int char2 = mpz_cmp_ui(bw->p, 2) == 0;

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

    abase_proxy abase_solutions(pi, A_width);
    mpfq_vbase_ptr A = abase_solutions.A;

    /* }}} */

    matmul_top_init(mmt, A, pi, pl, bw->dir);

    parasite_fixer pfixer(mmt);

    mmt_vec ymy[2];
    mmt_vec_ptr y = ymy[0];
    mmt_vec_ptr my = ymy[1];

    mmt_vec_init(mmt,0,0, y,   bw->dir, /* shared ! */ 1, mmt->n[bw->dir]);
    mmt_vec_init(mmt,0,0, my, !bw->dir,                0, mmt->n[!bw->dir]);
    mmt_full_vec_set_zero(y);
    mmt_full_vec_set_zero(my);


    /* this is really a misnomer, because in the typical case, M is
     * rectangular, and then the square matrix does induce some padding.
     * This is however not the padding we are interested in. The padding
     * we're referring to in the naming of this variable is the one which
     * is related to the number of jobs and threads: internal dimensions
     * are arranged to be multiples */
    unsigned int unpadded = MAX(mmt->n0[0], mmt->n0[1]);

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

    rhs R(mmt, param_list_lookup_string(pl, "rhs"), solutions);

    char * kprefix;

    if (tcan_print)
        printf("Trying to build solutions %u..%u\n", solutions[0], solutions[1]);

    int rc = asprintf(&kprefix, "K.sols%u-%u", solutions[0], solutions[1]);
    ASSERT_ALWAYS(rc >= 0);


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
                    mmt_my_own_size_in_items(y));
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

    /* The mmt_vec_unapply_T below is disturbing. Does that mean that the
     * vectors we feed are with T applied ? All of this is slightly
     * hairy...
     * 
     * I _think_ that this may be related to "apply_T" and "unapply_T"
     * being named in a way that reflects a direction which is not what
     * we have in mind all the time. "apply" for column vectors means
     * v<-T^-1*v !!
     *
     * We have, in order
     *
     * a vector read from disk
     * unapply_T
     *
     * apply_T
     * twist
     * matmul_top_mul
     * untwist
     * unapply_T
     *
     * which, for bw->dir==1, handles the following data (for a column
     * vector V on the left, for the same vector written as a row vector
     * on the right)
     *
     *  load :  T^-1 * V
     *  unapply : V
     *
     *  apply   : T^-1 * V
     *  twist   : Sc * T^-1 * V
     *  mul     : Sc*Mpad*T*Sc^-1 * T^-1 * V = Sc*Mpad*V
     *  untwist : Mpad*V
     *  unapply : T*Mpad*V
     *
     *  apply   : Mpad*V
     *  twist   : Sc*Mpad*V
     *  mul     : Sc*(Mpad*T)^2*T^-1*V
     *  untwist : (Mpad*T)^2*T^-1*V
     *  unapply : T*(Mpad*T)^2*T^-1*V
     *
     *  apply   : (Mpad*T)^2*T^-1*V
     *  (and so on)
     */
    mmt_vec_unapply_T(mmt, y);

    for(int i = 1 ; R ? (i<=1) : (i < 10) ; i++) {

        int input_is_zero;
        int pad_is_zero;
        int hamming_out;
        std::tie(input_is_zero, pad_is_zero, hamming_out) = test_one_vector(mmt, ymy, R);
        is_zero = hamming_out == 0;

        if (!is_zero && pfixer.attempt_to_fix) {
            if (hamming_out >= sqrt(mmt->n0[!bw->dir])) {
                printf("# not attempting to fix possibly wrong coordinates, too many non-zero coordinates in output\n");
            } else {
                pfixer.attempt(my, y);
                std::tie(input_is_zero, pad_is_zero, hamming_out) = test_one_vector(mmt, ymy, R);
                is_zero = hamming_out == 0;
            }
        }

        /* {{{ save file. If input is zero, bail out */
        {
            char * tmp;
            int rc = asprintf(&tmp, "%s%s.%d", input_is_zero ? "zero" : "",
                    kprefix, i-1);
            ASSERT_ALWAYS(rc >= 0);
            mmt_vec_save(y, tmp, unpadded);
            if (input_is_zero) {
                if (tcan_print)
                    fprintf(stderr,
                            "Using %sM^%u%s as input: Found zero vector."
                            " (coordinates on the padding part are %s zero)."
                            "\n"
                            "No solution found."
                            " Most certainly a bug."
                            "\n",
                            bw->dir ? "" : "V * ",
                            i-1,
                            bw->dir ? " * V" : "",
                            pad_is_zero ? "also" : "NOT");
                if (tcan_print && !pad_is_zero)
                    fprintf(stderr,
                            "For reference, this useless vector (non-zero out, zero in) is stored in %s."
                            "\n",
                            tmp);
                serialize(pi->m);
                pthread_mutex_lock(pi->m->th->m);
                exitcode=1;
                pthread_mutex_unlock(pi->m->th->m);
                return NULL;
            }
            free(tmp);
        }
        /* }}} */

        if (tcan_print) {
            if (R) {
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
            if (hamming_out && tcan_print)
                printf("Hamming weight is %d\n", hamming_out);
        }

        if (is_zero) {
            winning_iter = i-1;
            break;
        }

    }

    /* we exit with an untwisted vector. */

    if (!is_zero) {
        int nb_nonzero_coeffs=mmt_vec_hamming_weight(y);
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

    if (leader) {
        /* About the ASSERT below: The 0-characteristic subspace
         * nilpotent part is not really a topic of interest in the case
         * of inhomogenous linear systems.  Exiting the loop with i>1 in
         * the inhomogenous case is an error, it seems, and I believe
         * that the solution provided was never really a solution. On the
         * contrary, the i>1 case in homogenous situations is perfectly
         * valid, and we're handling it now instead of using
         * winning_iter==0 always.
         */
        ASSERT_ALWAYS(winning_iter == 0 || !R);

        char * tmp;
        int rc = asprintf(&tmp, "%s.%u", kprefix, winning_iter);
        ASSERT_ALWAYS(rc >= 0);
        /* {{{ append the RHS coefficients if relevant */
        if (R) {
            printf("Expanding %s so as to include the coefficients for the %d RHS columns\n", tmp, R.nrhs);
            FILE * f = fopen(tmp, "ab");
            ASSERT_ALWAYS(f);
            rc = fseek(f, 0, SEEK_END);
            ASSERT_ALWAYS(rc >= 0);
            R.fwrite_rhs_coeffs(f);
            fclose(f);
        }

        printf("%s is now a %s nullspace vector for %s"
                " (for a square M of dimension %" PRIu32"x%" PRIu32").\n",
                tmp,
                bw_dirtext[bw->dir],
                R ? "(M|RHS)" : "M",
                unpadded, unpadded);
        /* }}} */

        char * tmp2;
        rc = asprintf(&tmp2, "%s.%u.txt", kprefix, winning_iter);
        ASSERT_ALWAYS(rc >= 0);
        /* {{{ write an ascii version while we're at it */
        {
            FILE * f  = fopen(tmp, "rb");
            FILE * f2 = fopen(tmp2, "w");
            ASSERT_ALWAYS(f);
            ASSERT_ALWAYS(f2);
            void * data = malloc(A->vec_elt_stride(A, 1));
            for(uint32_t i = 0 ; i < mmt->n0[bw->dir] ; i++) {
                size_t rc = fread(data, A->vec_elt_stride(A, 1), 1, f);
                ASSERT_ALWAYS(rc == 1);
                A->fprint(A, f2, data);
                fprintf(f2, "\n");
            }
            free(data);
            R.fprint_rhs_coeffs(f2);
            fclose(f2);
        }

        printf("%s (in ascii) is now a %s nullspace vector for %s"
                " (for the original M, of dimension %" PRIu32"x%" PRIu32").\n",
                tmp2,
                bw_dirtext[bw->dir],
                R ? "(M|RHS)" : "M",
                mmt->n0[0], mmt->n0[1]);
        /* }}} */
        free(tmp2);
        free(tmp);
    }

    serialize(pi->m);
    free(kprefix);

    mmt_vec_clear(mmt, y);
    mmt_vec_clear(mmt, my);

    matmul_top_clear(mmt);

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

