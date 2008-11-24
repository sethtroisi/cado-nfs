#define _POSIX_C_SOURCE 200112L

#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cerrno>

// we're doing open close mmap truncate...
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>


#include <string>


#include "matmul.h"
#include "matmul_top.h"
#include "select_mpi.h"
#include "intersections.hpp"
#include "info_file.h"

using namespace std;

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

void matmul_top_init(matmul_top_data_ptr mmt,
        abobj_ptr abase,
        parallelizing_info_ptr pi,
        const char * filename)
{
    memset(mmt, 0, sizeof(*mmt));

    mmt->abase = abase;
    for(int d = 0 ; d < 2 ; d++) {
        mmt->wr[d]->v = abinit(mmt->abase, mmt->wr[d]->i1 - mmt->wr[d]->i0);
    }

    mmt->pi = pi;

    read_info_file(mmt, filename);

}

void matmul_top_read_matrix(matmul_top_data_ptr mmt)
{
    mmt->mm = matmul_reload_cache(mmt->abase, mmt->filename);
    if (mmt->mm)
        return;

    fprintf(stderr, "Could not find cache file %s\n", mmt->filename);
    mmt->mm = matmul_build(mmt->abase, mmt->filename);
    matmul_save_cache(mmt->mm, mmt->filename);
}


void matmul_top_clear(matmul_top_data_ptr mmt, abobj_ptr abase MAYBE_UNUSED)
{
    matmul_clear(mmt->mm);
    for(int d = 0 ; d < 2 ; d++) {
        abclear(mmt->abase, mmt->wr[d]->v, mmt->wr[d]->i1 - mmt->wr[d]->i0);
    }
    free(mmt->filename);
}

static void
broadcast_across(matmul_top_data_ptr mmt, int d)
{
    // broadcasting across columns is for common agreement on a column
    // vector.
#ifdef  CONJUGATED_PERMUTATIONS
    for(unsigned int i = 0 ; i < mmt->wr[d]->xlen ; i++) {
        void * ptr = mmt->wr[d]->v + aboffset(mmt->abase, mmt->wr[d]->x[i].offset_me);
        MPI_Bcast(ptr, mmt->wr[d]->x[i].count * sizeof(abt),
                MPI_BYTE, mmt->wr[d]->x[i].k, mmt->pi->wr[d]->pals);
    }
    if (mmt->wr[d]->i1 > mmt->wr[!d]->n) {
        abzero(mmt->abase,
                mmt->wr[d]->v + aboffset(mmt->abase, mmt->wr[!d]->n-mmt->wr[d]->i0),
                mmt->wr[d]->i1-mmt->wr[!d]->n);
    }
#else
    choke me;
#endif
}

void matmul_top_fill_random(matmul_top_data_ptr mmt, int d)
{
    // In conjugation mode, it is possible to fill exactly the data chunk
    // that will eventually be relevant. However, it's easy enough to
    // fill our output vector with garbage, and do broadcast_across_column
    // afterwards...

    abrandom(mmt->abase, mmt->wr[d]->v, mmt->wr[d]->i1 - mmt->wr[d]->i0);

    // reconcile all cells which correspond to the same vertical block.
    broadcast_across(mmt, d);
}

// comments are written assuming that d == 1 (save right vector). The
// stable state we start with is after a matrix-times-vector product (or
// just before one): The left vector has been computed, and the data for
// indices [i0..i1[ is on the nodes seeing those indices vertically
// as well. Data is therefore found in the right vector area, as after
// the reduce step.
//
// Doing a broadcast_across columns will ensure that each row contains
// the complete data set for our vector.
void save_vector(matmul_top_data_ptr mmt, int d, unsigned int index, unsigned int iter)
{
    // just to make sure...
    broadcast_across(mmt, !d);

    // Now in our row, everybody has the vector. We'll do I/O only from
    // one row, that's easier. Pick the one whose rank is zero in the
    // column communicator.
    if (mmt->pi->wr[d]->jrank != 0)
        return;

    // The vector is saved by _one_ node.
    int fd;

    if (mmt->pi->wr[!d]->jrank == 0) {
        string filename = "V";
        filename += index;
        filename += ".";
        filename += iter;
        filename += ".shuffled";
        fd = open(filename.c_str(), O_WRONLY|O_CREAT, 0666);
        if (fd < 0) {
            fprintf(stderr, "fopen(%s): %s\n",
                    filename.c_str(), strerror(errno));
            // XXX what do we do here ?
            return;
        }
        ftruncate(fd, abbytes(mmt->abase, mmt->wr[d]->n));
    }
        void * sendbuf = (void *) mmt->wr[d]->v;
        int sendcount = abbytes(mmt->abase, mmt->wr[d]->i1 - mmt->wr[d]->i0);
        void * recvbuf = mmap(NULL, abbytes(mmt->abase, mmt->wr[d]->n),
                PROT_WRITE, 0, fd, 0);
        int len = mmt->pi->wr[d]->ncores * mmt->pi->wr[d]->njobs;


    my_pthread_barrier_t bar[1];
    my_pthread_barrier_init(bar, NULL, mmt->pi->wr[d]->ncores);

    /* At this point, if ever we could have as many communications via
     * the MPI channel as we have wires, that would be cool (IOW, we
     * would simply drop the pthread barrier, + pay some attention to our
     * helper variables).
     */
    int * recvcounts = (int *) malloc(mmt->pi->wr[d]->njobs * sizeof(int));
    int * displs = (int *) malloc(mmt->pi->wr[d]->njobs * sizeof(int));

    /* Because of the mispleasure mentioned above, we have to serialize
     * here (for now).
     */
    for(unsigned int t = 0 ; t < mmt->pi->wr[d]->ncores ; t++) {
        my_pthread_barrier_wait(bar);
        if (t != mmt->pi->wr[d]->trank)
            continue;

        int dis = 0;
        for(int k = 0 ; k < len ; k++) {
            int v = abbytes(mmt->abase, mmt->wr[d]->fences[k+1] - mmt->wr[d]->fences[k]);
            if (k % mmt->pi->wr[d]->ncores == t) {
                displs[k/mmt->pi->wr[d]->ncores] = dis;
                recvcounts[k/mmt->pi->wr[d]->ncores] = v;
            }
            dis += v;
        }

        // Do gatherv directly into the mmaped area.
        MPI_Gatherv(sendbuf, sendcount, MPI_BYTE, recvbuf, recvcounts, displs, MPI_BYTE, 0, mmt->pi->wr[d]->pals);

        munmap(recvbuf, abbytes(mmt->abase, mmt->wr[d]->n));
    }

    free(recvcounts);
    free(displs);

    my_pthread_barrier_destroy(bar);

    if (mmt->pi->wr[!d]->jrank == 0) {
        close(fd);
    }
}

// now the exact opposite.
void load_vector(matmul_top_data_ptr mmt, int d, unsigned int index, unsigned int iter)
{
    // just to make sure...
    broadcast_across(mmt, !d);

    // Now in our row, everybody has the vector. We'll do I/O only from
    // one row, that's easier. Pick the one whose rank is zero in the
    // column communicator.
    if (mmt->pi->wr[d]->jrank != 0)
        return;

    // The vector is saved by _one_ node.
    int fd;

    if (mmt->pi->wr[!d]->jrank == 0) {
        string filename = "V";
        filename += index;
        filename += ".";
        filename += iter;
        filename += ".shuffled";
        fd = open(filename.c_str(), O_RDONLY);
        if (fd < 0) {
            fprintf(stderr, "fopen(%s): %s\n",
                    filename.c_str(), strerror(errno));
            // XXX what do we do here ?
            return;
        }
    }
        void * recvbuf = (void *) mmt->wr[d]->v;
        int recvcount = abbytes(mmt->abase, mmt->wr[d]->i1 - mmt->wr[d]->i0);
        void * sendbuf = mmap(NULL, abbytes(mmt->abase, mmt->wr[d]->n),
                PROT_READ, 0, fd, 0);
        int len = mmt->pi->wr[d]->ncores * mmt->pi->wr[d]->njobs;


    my_pthread_barrier_t bar[1];
    my_pthread_barrier_init(bar, NULL, mmt->pi->wr[d]->ncores);

    /* At this point, if ever we could have as many communications via
     * the MPI channel as we have wires, that would be cool (IOW, we
     * would simply drop the pthread barrier, + pay some attention to our
     * helper variables).
     */
    int * sendcounts = (int *) malloc(mmt->pi->wr[d]->njobs * sizeof(int));
    int * displs = (int *) malloc(mmt->pi->wr[d]->njobs * sizeof(int));

    /* Because of the mispleasure mentioned above, we have to serialize
     * here (for now).
     */
    for(unsigned int t = 0 ; t < mmt->pi->wr[d]->ncores ; t++) {
        my_pthread_barrier_wait(bar);
        if (t != mmt->pi->wr[d]->trank)
            continue;

        int dis = 0;
        for(int k = 0 ; k < len ; k++) {
            int v = abbytes(mmt->abase, mmt->wr[d]->fences[k+1] - mmt->wr[d]->fences[k]);
            if (k % mmt->pi->wr[d]->ncores == t) {
                displs[k/mmt->pi->wr[d]->ncores] = dis;
                sendcounts[k/mmt->pi->wr[d]->ncores] = v;
            }
            dis += v;
        }

        // Do scatterv directly from the mmaped area.
        MPI_Scatterv(sendbuf, sendcounts, displs, MPI_BYTE, recvbuf, recvcount, MPI_BYTE, 0, mmt->pi->wr[d]->pals);

        munmap(sendbuf, abbytes(mmt->abase, mmt->wr[d]->n));
    }

    free(sendcounts);
    free(displs);

    my_pthread_barrier_destroy(bar);

    if (mmt->pi->wr[!d]->jrank == 0) {
        close(fd);
    }
}

