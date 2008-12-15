#define _POSIX_C_SOURCE 200112L
#define _GNU_SOURCE

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>

// we're doing open close mmap truncate...
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>


#include "manu.h"

#include "matmul.h"
#include "matmul_top.h"
#include "select_mpi.h"
#include "intersections.h"
#include "info_file.h"

#if 0
abt * extra_svec_init(matmul_top_data_srcptr mmt) {
    return abinit(mmt->abase, mmt->i1 - mmt->i0);
}
void extra_svec_clear(matmul_top_data_srcptr mmt, abt * x) {
    abclear(mmt->abase, x, mmt->i1 - mmt->i0);
}
void extra_svec_update(matmul_top_data_srcptr mmt, abt * x) {
    abcopy(mmt->abase, x, mmt->v, mmt->i1 - mmt->i0);
}
void extra_svec_recall(matmul_top_data_srcptr mmt, const abt * x) {
    abcopy(mmt->abase, mmt->v, x, mmt->i1 - mmt->i0);
}
#endif

#ifndef CONJUGATED_PERMUTATIONS
#error "Do you really, really want to use arbitrary left and right sigmas ?"
#endif

void matmul_top_init(matmul_top_data_ptr mmt,
        abobj_ptr abase,
        matmul_ptr mm,
        parallelizing_info_ptr pi,
        const char * filename)
{
    memset(mmt, 0, sizeof(*mmt));

    mmt->abase = abase;
    mmt->pi = pi;
    mmt->mm = mm;

    // n[]
    // ncoeffs_total, ncoeffs,
    // locfile,
    // fences[]
    // wr[]->*
    // are filled in later within read_info_file

    mmt->wr[0]->all_v = malloc(mmt->pi->wr[0]->ncores * sizeof(abt *));
    mmt->wr[1]->all_v = malloc(mmt->pi->wr[1]->ncores * sizeof(abt *));

    read_info_file(mmt, filename);
}

void matmul_top_read_matrix(matmul_top_data_ptr mmt)
{
    mmt->mm = matmul_reload_cache(mmt->abase, mmt->locfile);
    if (mmt->mm)
        return;

    fprintf(stderr, "Could not find cache file %s\n", mmt->locfile);
    mmt->mm = matmul_build(mmt->abase, mmt->locfile);
    matmul_save_cache(mmt->mm, mmt->locfile);
}


void matmul_top_clear(matmul_top_data_ptr mmt, abobj_ptr abase MAYBE_UNUSED)
{
    matmul_clear(mmt->mm);
    for(int d = 0 ; d < 2 ; d++) {
        abclear(mmt->abase, mmt->wr[d]->v, mmt->wr[d]->i1 - mmt->wr[d]->i0);
    }
    free(mmt->locfile);
}

/* XXX
 * ``across'' (horizontally) and ``down'' (vertically) are just here for
 * exposition. The binding of this operation to job/thread arrangement
 * is flexible, through argument d.
 * XXX
 */
static void
broadcast_down(matmul_top_data_ptr mmt, int d)
{
    // broadcasting down columns is for common agreement on a right
    // source vector, once a left output has been computed.
    //
    // broadcasting down a column is when d == 1.
#ifndef  CONJUGATED_PERMUTATIONS
    choke me;
#endif

    pi_wiring_ptr picol = mmt->pi->wr[d];
    pi_wiring_ptr pirow = mmt->pi->wr[!d];

    mmt_wiring_ptr mcol = mmt->wr[d];
    // mmt_wiring_ptr mrow = mmt->wr[!d];

    // unsigned int ncols = mmt->n[d];
    unsigned int nrows = mmt->n[!d];

    /* The intersection of our column-wise input range [i0..i1[ with
     * all the fences existing for the output ranges is described by
     * mcol->x ; note how this intersection is common to a
     * complete column.
     */
    /* If our vector data is not shared amongst threads, it has to be
     * broadcasted. We define the ``extra'' flag in this case.
     */
    unsigned int extra = ((mmt->flags[d] & THREAD_SHARED_VECTOR) != 0);
    for(unsigned int t = 0 ; t < pirow->ncores + extra ; t++) {
        serialize_threads(pirow);
        if (t == pirow->trank) {
            for(unsigned int u = 0 ; u < picol->ncores ; u++) {
                serialize_threads(picol);
                if (u != picol->trank) {
                    continue;
                }
                // our turn.
                for(unsigned int i = 0 ; i < mcol->xlen ; i++) {
                    struct isect_info * xx = &(mcol->x[i]);
                    if (xx->k % picol->ncores != picol->trank) {
                        // then we're not going to have any work to do.
                        continue;
                    }
                    unsigned int off = aboffset(mmt->abase, xx->offset_me);
                    void * ptr = mcol->v + off;
                    size_t siz = xx->count * sizeof(abt);
                    unsigned int root = xx->k / picol->ncores;
                    MPI_Bcast(ptr, siz, MPI_BYTE, root, picol->pals);
                }
                /* Note that at this point, it is not guaranteed that all
                 * our input range [i0,i1[ has been covered. In the case
                 * of rectangular matrices, we might have to cope with
                 * this case. This is done later on. Note though that
                 * rectangular matrices (= block lanczos case) yield an
                 * allreduce call instead of broadcast, so the point is
                 * moot.
                 */
            }
        } else if (t == pirow->trank + 1 && extra) {
            /* Then, while the next column is communicating, we'll do the
             * inter-thread bcast if needed. It does not matter when the
             * right vectors are shared, of course, since no matter which
             * thread received the data, it ends up in the right place. */
            for(unsigned int i = 0 ; i < mcol->xlen ; i++) {
                struct isect_info * xx = &(mcol->x[i]);
                unsigned int src = xx->k % picol->ncores;
                unsigned int off = aboffset(mmt->abase, xx->offset_me);
                void * ptr = mcol->v + off;
                void * sptr = mmt->wr[d]->all_v[src] + off;
                size_t siz = xx->count * sizeof(abt);
                if (picol->trank != src) {
                    memcpy(ptr, sptr, siz);
                }
            }
        }
    }

    if (mcol->i1 > nrows) {
        /* untested */
        BUG();
        if (extra || picol->trank == 0) {
            abzero(mmt->abase,
                    mcol->v + aboffset(mmt->abase, nrows-mcol->i0),
                    mcol->i1-nrows);
        }
        if (!extra) {
            /* of course all threads within the same column have common
             * i0 and i1, so no bug here apparently */
            serialize_threads(picol);
        }
    }
}

static void
reduce_across(matmul_top_data_ptr mmt, int d)
{
#ifndef  CONJUGATED_PERMUTATIONS
    choke me;
#endif

    /* reducing across a row is when d == 0 */

    // pi_wiring_ptr pirow = mmt->pi->wr[d];
    // pi_wiring_ptr picol = mmt->pi->wr[!d];

}



void matmul_top_fill_random_source(matmul_top_data_ptr mmt, int d)
{
    // In conjugation mode, it is possible to fill exactly the data chunk
    // that will eventually be relevant. However, it's easy enough to
    // fill our output vector with garbage, and do broadcast_down
    // afterwards...
    abrandom(mmt->abase, mmt->wr[d]->v, mmt->wr[d]->i1 - mmt->wr[d]->i0);

    // reconcile all cells which correspond to the same vertical block.
    broadcast_down(mmt, d);
}


/* Vector I/O is done by only one job, one thread. It incurs a
 * significant amount of memory allocation, but this is done relatively
 * rarely.
 *
 * Vectors, as handled within the core routines, are permuted. He have
 * coordinates (v_{\sigma(0)}...v_{\sigma(n-1)}), where \sigma is the
 * permutation given by the *_perm files. Notice that we assume that we
 * have conjugated permutations on the two sides. This means that no
 * data manipulation is required within the critical loops (this would
 * imply a fuzzy communication pattern, boiling down to essentially
 * using allreduce instead of reduce).
 */

/* As always, comments are written with one given point of view in mind.
 * The vector wich gets saved is always the source vector. Here, our
 * arbitrary description choice is that we iterate M^i times vector,
 * hence the source vector is on the right: d == 1.
 */

// comments, variable names and so on are names which for simplicity
// reflect the situation d == 1 (save right vector). The stable state we
// start with is after a matrix-times-vector product (or just before
// one): The left vector has been computed, and the data for indices
// [i0..i1[ is on the nodes seeing those indices vertically as well. Data
// is therefore found in the right vector area, as after the reduce step.
//
// Doing a broadcast_down columns will ensure that each row contains
// the complete data set for our vector.

/* The first function is the back-end, run only on the top row.
 */
static void save_vector_toprow(matmul_top_data_ptr mmt, int d, unsigned int index, unsigned int iter)
{
    // pi_wiring_ptr picol = mmt->pi->wr[d];
    pi_wiring_ptr pirow = mmt->pi->wr[!d];

    // The vector is saved by _one_ node, _one_ thread.
    void * recvbuf = NULL;
    int fd = -1;
    size_t siz = abbytes(mmt->abase, mmt->n[!d]);
    // the page size is always a multiple of two, so rounding to the next
    // multiple is easy.
    size_t wsiz = ((siz - 1) | (sysconf(_SC_PAGESIZE)-1)) + 1;

    if (pirow->jrank == 0) {
        if (pirow->trank == 0) {
            char * filename;
            asprintf(&filename, "V%u.%u.twisted", index, iter);
            fd = open(filename, O_RDWR|O_CREAT, 0666);
            if (fd < 0) {
                fprintf(stderr, "fopen(%s): %s\n",
                        filename, strerror(errno));
                fprintf(stderr, "WARNING: checkpointing failed\n");
                free(filename);
                // XXX what do we do here ?
                return;
            }

            ftruncate(fd, wsiz);
            recvbuf = mmap(NULL, wsiz, PROT_WRITE, MAP_SHARED, fd, 0);
            if (recvbuf == MAP_FAILED) {
                fprintf(stderr, "mmap(%s): %s\n",
                        filename, strerror(errno));
                fprintf(stderr, "WARNING: checkpointing failed\n");
                free(filename);
                close(fd);
                return;
            }
            free(filename);
        }
        /* Now all threads from job zero see the area mmaped by their
         * leader thread.
         */
        thread_agreement(pirow, &recvbuf, 0);
    }

    /* Prepare gatherv arguments */
    void * sendbuf = (void *) mmt->wr[d]->v;
    int sendcount = abbytes(mmt->abase, mmt->wr[d]->i1 - mmt->wr[d]->i0);

    int * displs = (int *) malloc(pirow->njobs * sizeof(int));
    int * recvcounts = (int *) malloc(pirow->njobs * sizeof(int));

    /* recall that ``fences'' must be understood as vertical when
     * visually they form a vertical line. Since we are interested in
     * column indices here, this means vertical fences */
    for(unsigned int k = 0 ; k < pirow->njobs ; k++) {
        int x = k * pirow->ncores + pirow->trank;
        unsigned int dx0 = mmt->fences[d][x];
        unsigned int dx1 = mmt->fences[d][x+1];

        displs[k] = abbytes(mmt->abase, dx0);
        recvcounts[k] = abbytes(mmt->abase, dx1 - dx0);
    }

    /* At this point, if ever we could have as many communications via
     * the MPI channel as we have wires, that would be cool: we would
     * simply drop the pthread barrier.  Lacking this possibility, we
     * have to serialize here (for now).
     */
    for(unsigned int t = 0 ; t < pirow->ncores ; t++) {
        serialize_threads(pirow);
        if (t == pirow->trank) { // our turn.
            ASSERT(sendcount == recvcounts[pirow->jrank]);
            MPI_Gatherv(sendbuf, sendcount, MPI_BYTE,
                    recvbuf, recvcounts, displs, MPI_BYTE,
                    0, pirow->pals);
        }
    }
    // what a pain in the Xss. Because we don't exit with the mutex
    // locked, we can't do much... The main thread might end up calling
    // munmap before we complete the call to MPI_Gatherv.

    serialize_threads(pirow);
    
    free(recvcounts);
    free(displs);

    if (pirow->jrank == 0 && pirow->trank == 0) {
        munmap(recvbuf, wsiz);
        ftruncate(fd, siz);
        close(fd);
    }
}

void matmul_top_save_vector(matmul_top_data_ptr mmt, int d, unsigned int index, unsigned int iter)
{
    // just to make sure...
    broadcast_down(mmt, d);

    // we'll do some funky stuff, so serialize in order to avoid jokes.
    // Penalty is insignificant anyway since I/O are rare.
    serialize(mmt->pi->m);
    serialize_threads(mmt->pi->m);

    pi_wiring_ptr picol = mmt->pi->wr[d];
    // pi_wiring_ptr pirow = mmt->pi->wr[!d];

    // Now every job row has the complete vector.  We'll do I/O from only
    // one row, that's easier. Pick the topmost one.
    if (picol->jrank == 0 && picol->trank == 0) {
        save_vector_toprow(mmt, d, index, iter);
    }

    serialize(mmt->pi->m);
    serialize_threads(mmt->pi->m);

}

// now the exact opposite.

/* The backend is exactly dual to save_vector_toprow.
 */
static void load_vector_toprow(matmul_top_data_ptr mmt, int d, unsigned int index, unsigned int iter)
{
    // pi_wiring_ptr picol = mmt->pi->wr[d];
    pi_wiring_ptr pirow = mmt->pi->wr[!d];
    
    void * sendbuf = NULL;
    int fd = -1;
    size_t siz = abbytes(mmt->abase, mmt->n[!d]);
    size_t wsiz = ((siz - 1) | (sysconf(_SC_PAGESIZE)-1)) + 1;
    if (pirow->jrank == 0) {
        if (pirow->trank == 0) {
            char * filename;
            asprintf(&filename, "V%u.%u.twisted", index, iter);
            fd = open(filename, O_RDONLY, 0666);
            DIE_ERRNO_DIAG(fd < 0, "fopen", filename);
            sendbuf = mmap(NULL, wsiz, PROT_READ, MAP_SHARED, fd, 0);
            DIE_ERRNO_DIAG(sendbuf == MAP_FAILED, "mmap", filename);
            free(filename);
        }
        /* Now all threads from job zero see the area mmaped by their
         * leader thread.
         */
        thread_agreement(pirow, &sendbuf, 0);
    }

    void * recvbuf = (void *) mmt->wr[d]->v;
    int recvcount = abbytes(mmt->abase, mmt->wr[d]->i1 - mmt->wr[d]->i0);

    int * displs = (int *) malloc(pirow->njobs * sizeof(int));
    int * sendcounts = (int *) malloc(pirow->njobs * sizeof(int));

    for(unsigned int k = 0 ; k < pirow->njobs ; k++) {
        int x = k * pirow->ncores + pirow->trank;
        unsigned int dx0 = mmt->fences[d][x];
        unsigned int dx1 = mmt->fences[d][x+1];

        displs[k] = abbytes(mmt->abase, dx0);
        sendcounts[k] = abbytes(mmt->abase, dx1 - dx0);
    }

    for(unsigned int t = 0 ; t < pirow->ncores ; t++) {
        serialize_threads(pirow);
        if (t == pirow->trank) { // our turn.
            ASSERT(recvcount == sendcounts[pirow->jrank]);
            MPI_Scatterv(sendbuf, sendcounts, displs, MPI_BYTE,
                    recvbuf, recvcount, MPI_BYTE,
                    0, pirow->pals);
        }
    }

    serialize_threads(pirow);

    free(sendcounts);
    free(displs);

    if (pirow->jrank == 0 && pirow->trank == 0) {
        munmap(sendbuf, wsiz);
        close(fd);
    }
}

void matmul_top_load_vector(matmul_top_data_ptr mmt, int d, unsigned int index, unsigned int iter)
{
    serialize(mmt->pi->m);
    serialize_threads(mmt->pi->m);

    pi_wiring_ptr picol = mmt->pi->wr[d];
    pi_wiring_ptr pirow = mmt->pi->wr[!d];

    if (picol->jrank == 0 && picol->trank == 0) {
        load_vector_toprow(mmt, d, index, iter);
    }

    serialize(mmt->pi->m);
    serialize_threads(mmt->pi->m);

    // after loading, we need to broadcast the data. As in other
    // occcasions, this has to be one one thread at a time.

    size_t siz = abbytes(mmt->abase, mmt->wr[d]->i1 - mmt->wr[d]->i0);
    if (picol->trank == 0) {
        // first, down on columns.
        for(unsigned int t = 0 ; t < pirow->ncores ; t++) {
            serialize_threads(pirow);
            if (t == pirow->trank) {
                void * ptr = mmt->wr[d]->v;
                MPI_Bcast(ptr, siz, MPI_BYTE, 0, picol->pals);
            }
        }
    }

    serialize(mmt->pi->m);

    // then maybe amongst threads if they're not sharing data.
    if ((mmt->flags[d] & THREAD_SHARED_VECTOR) == 0) {
        // we have picol->ncores threads, which would be delighted to get
        // access to some common data. The data is known to sit in thread
        // zero, so that's relatively easy.
        
        if (picol->trank != 0) {
            memcpy(mmt->wr[d]->v, mmt->wr[d]->all_v[0], siz);
        }
    }
}
