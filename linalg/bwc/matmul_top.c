#include "cado.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>             // truncate()
#include <sys/types.h>          // truncate()
#include <errno.h>
#include <sys/stat.h>
#include "bwc_config.h"
#include "matmul.h"
#include "matmul_top.h"
#include "select_mpi.h"
#include "intersections.h"
#include "filenames.h"
#include "balancing_workhorse.h"
#include "portability.h"
#include "misc.h"
#include "random_matrix.h"

/* Our innermost communication routines are essentially all-gather and
 * reduce-scatter, following the MPI terminology. We provide several
 * choices for doing this. The compile-time definitions here allow to
 * change which is used. The "best" should in theory always be the one
 * provided by the MPI implementation. Unfortunately, it is not always
 * that clear.
 *
 * Note that because we have several threads wanting to do all these
 * operations one after the other, we also have some interest in using
 * non-blockng collectives, or even also some sort of
 * communication-computation overlap.
 */

/* choices for reduce-scatter */
#define RS_CHOICE_STOCK_RS           1
#define RS_CHOICE_STOCK_RSBLOCK      (MPI_VERSION_ATLEAST(2,2) ? 2 : -2)
#define RS_CHOICE_STOCK_IRSBLOCK     (MPI_VERSION_ATLEAST(3,0) ? 3 : -3)
#define RS_CHOICE_MINE               4
#define RS_CHOICE_MINE_DROP_IN       (MPI_VERSION_ATLEAST(3,0) ? 5 : -5)
#define RS_CHOICE_MINE_PARALLEL      6
#define RS_CHOICE_MINE_OVERLAPPING   7  /* TODO */
/* choices for all-gather */
#define AG_CHOICE_STOCK_AG           1
#define AG_CHOICE_STOCK_IAG          (MPI_VERSION_ATLEAST(3,0) ? 2 : -2)

/* The actual performance is affected by the communicator size, the chunk
 * size, and the operation. When we're doing bxor, we have the following
 * measurements for RS.
 *
 * n=2, chunk=174 MB, openmpi-1.8.3, IB FDR.
 * RS_CHOICE_MINE_PARALLEL .25
 * RS_CHOICE_MINE .45
 * RS_CHOICE_MINE_DROP_IN  .48
 * RS_CHOICE_STOCK_RS      .77
 * RS_CHOICE_STOCK_RSBLOCK 1.63    
 * RS_CHOICE_STOCK_IRSBLOCK        => seg fault
 * RS_CHOICE_MINE_OVERLAPPING      => to be implemented
 */

/* _this_ part can be configured */
#define RS_CHOICE RS_CHOICE_MINE_PARALLEL
#define AG_CHOICE (MPI_VERSION_ATLEAST(3,0) ? AG_CHOICE_STOCK_IAG : AG_CHOICE_STOCK_AG)

/* some sanity checking */
#if RS_CHOICE < 0
#error "Choice for reduce-scatter strategy is invalid or not supported"
#endif
#if AG_CHOICE < 0
#error "Choice for reduce-scatter strategy is invalid or not supported"
#endif

/* we no longer use this.
#ifdef  HAVE_MPI3_API
#define MPI_LIBRARY_NONBLOCKING_COLLECTIVES
#endif
*/

///////////////////////////////////////////////////////////////////
/* Start with stuff that does not depend on abase at all -- this
 * provides a half-baked interface */

/* At some point we had this. Not sure it's still useful. */
#define ABASE_UNIVERSAL_READAHEAD_ITEMS 8

static permutation_data_ptr permutation_data_alloc();
static void permutation_data_free(permutation_data_ptr a);
static void permutation_data_push(permutation_data_ptr a, uint32_t u[2]);

void matmul_top_decl_usage(param_list_ptr pl)
{
    param_list_decl_usage(pl, "matrix",
            "the matrix file (binary)");
    param_list_decl_usage(pl, "balancing",
            "the matrix balancing file, as computed by mf_bal");
    param_list_decl_usage(pl, "random_matrix",
            "characteristics of a random matrix to be used for staged runs.");

    param_list_decl_usage(pl, "rebuild_cache",
            "force rebuilding matrix caches");
    param_list_decl_usage(pl, "export_cachelist",
            "print the per-node needed cache files and exit");
    param_list_decl_usage(pl, "save_submatrices",
            "after dispatching, save a copy of the local uncompressed matrix before creating the cache file");
    param_list_decl_usage(pl, "sequential_cache_build",
            "build the cache files sequentially on each node");
    param_list_decl_usage(pl, "sequential_cache_read",
            "read the cache files sequentially on each node");
    balancing_decl_usage(pl);
    matmul_decl_usage(pl);
}

void matmul_top_lookup_parameters(param_list_ptr pl)
{
    param_list_lookup_string(pl, "matrix");
    param_list_lookup_string(pl, "balancing");
    param_list_lookup_string(pl, "random_matrix");
    param_list_lookup_string(pl, "rebuild_cache");
    param_list_lookup_string(pl, "export_cachelist");
    param_list_lookup_string(pl, "save_submatrices");
    param_list_lookup_string(pl, "sequential_cache_build");
    param_list_lookup_string(pl, "sequential_cache_read");
    balancing_lookup_parameters(pl);
    matmul_lookup_parameters(pl);
}

/* {{{ vector init/clear */
void mmt_vec_init(matmul_top_data_ptr mmt, mpfq_vbase_ptr abase, pi_datatype_ptr pitype, mmt_vec_ptr v, int d, int flags)
{
    if (v == NULL) v = mmt->wr[d]->v;
    if (abase == NULL) abase = mmt->abase;
    if (pitype == NULL) pitype = mmt->pitype;
    unsigned int n = mmt->wr[d]->i1 - mmt->wr[d]->i0;
    matmul_aux(mmt->mm, MATMUL_AUX_GET_READAHEAD, &n);
    pi_comm_ptr picol = mmt->pi->wr[d];
    v->abase = abase;
    v->pitype = pitype;
    v->flags = flags;

    v->all_v = shared_malloc(picol, picol->ncores * sizeof(void *));

    /* Because we sometimes provide extra readahead space, we need to
     * zero out the allocated area. Strictly speaking, it would suffice
     * to zero out the readahead zone only, but it's easier this way.
     */
    if (flags & THREAD_SHARED_VECTOR) {
        if (picol->trank == 0) {
            abase->vec_init(abase, &v->v, n + ABASE_UNIVERSAL_READAHEAD_ITEMS);
            abase->vec_set_zero(abase, v->v, n + ABASE_UNIVERSAL_READAHEAD_ITEMS);
            for(unsigned int t = 0 ; t < picol->ncores ; t++)
                v->all_v[t] = v->v;
        }
        serialize_threads(picol);
        v->v = v->all_v[picol->trank];
    } else {
        abase->vec_init(abase, &v->v, n + ABASE_UNIVERSAL_READAHEAD_ITEMS);
        abase->vec_set_zero(abase, v->v, n + ABASE_UNIVERSAL_READAHEAD_ITEMS);
        v->all_v[picol->trank] = v->v;
        serialize_threads(picol);
    }
}

void mmt_vec_clear(matmul_top_data_ptr mmt, mmt_vec_ptr v, int d)
{
    if (v == NULL) v = mmt->wr[d]->v;
    unsigned int n = mmt->wr[d]->i1 - mmt->wr[d]->i0;
    matmul_aux(mmt->mm, MATMUL_AUX_GET_READAHEAD, &n);
    pi_comm_ptr picol = mmt->pi->wr[d];
    if (v->flags & THREAD_SHARED_VECTOR) {
        if (picol->trank == 0)
            v->abase->vec_clear(v->abase, &v->v, n + ABASE_UNIVERSAL_READAHEAD_ITEMS);
    } else {
        v->abase->vec_clear(v->abase, &v->v, n + ABASE_UNIVERSAL_READAHEAD_ITEMS);
    }
    shared_free(picol, v->all_v);
}
/* }}} */

/* Each job and thread is really the "owner" of one particular data area
 * (which is narrower than the data it's interested in for doing the
 * matrix times vector product, because that one is actually spread
 * across more threads and jobs).
 */
size_t mmt_my_own_offset_in_items(matmul_top_data_ptr mmt, int d)
{
    mmt_comm_ptr mcol = mmt->wr[d];
    pi_comm_ptr picol = mmt->pi->wr[d];
    /* set up pointers for our local share of the vectors to be
     * considered. This is with respect to vectors on size mcol.  */
    ASSERT_ALWAYS((mcol->i1 - mcol->i0) % picol->totalsize == 0);
    size_t eblock = (mcol->i1 - mcol->i0) /  picol->totalsize;
    /* the "column vector" (the one which is good for M*v) has size which
     * is the number of columns, so the size of rows... */
    int pos = picol->jrank * picol->ncores + picol->trank;
    return pos * eblock;
}

size_t mmt_my_own_offset_in_bytes(matmul_top_data_ptr mmt, mmt_vec_ptr v, int d)
{
    mmt_comm_ptr mcol = mmt->wr[d];
    if (v == NULL) v = mcol->v;
    return v->abase->vec_elt_stride(v->abase, mmt_my_own_offset_in_items(mmt, d));
}


void * mmt_my_own_subvec(matmul_top_data_ptr mmt, mmt_vec_ptr v, int d)
{
    mmt_comm_ptr mcol = mmt->wr[d];
    if (v == NULL) v = mcol->v;
    return SUBVEC(v, v, mmt_my_own_offset_in_items(mmt, d));
}


size_t mmt_my_own_size_in_items(matmul_top_data_ptr mmt, int d)
{
    mmt_comm_ptr mcol = mmt->wr[d];
    pi_comm_ptr picol = mmt->pi->wr[d];
    ASSERT_ALWAYS((mcol->i1 - mcol->i0) % picol->totalsize == 0);
    size_t eblock = (mcol->i1 - mcol->i0) /  picol->totalsize;
    return eblock;
}

size_t mmt_my_own_size_in_bytes(matmul_top_data_ptr mmt, mmt_vec_ptr v, int d)
{
    mmt_comm_ptr mcol = mmt->wr[d];
    if (v == NULL) v = mcol->v;
    return v->abase->vec_elt_stride(v->abase, mmt_my_own_size_in_items(mmt, d));
}

/* This copies **ONLY** our locally owned data from v to w. No provision
 * is made to handle the case where we would want a broader copy. For
 * that, a priori we want all threads to contribute.
 */
void mmt_vec_set(matmul_top_data_ptr mmt, mmt_vec_ptr w, mmt_vec_ptr v, int d)
{
    mmt_comm_ptr mcol = mmt->wr[d];
    if (w == NULL) w = mcol->v;
    if (v == NULL) v = mcol->v;
    ASSERT_ALWAYS(v->abase == w->abase);
    v->abase->vec_set(v->abase,
            mmt_my_own_subvec(mmt, w, d),
            mmt_my_own_subvec(mmt, v, d),
            mmt_my_own_size_in_items(mmt, d));
}

void mmt_vec_set_zero(matmul_top_data_ptr mmt, mmt_vec_ptr v, int d)
{
    pi_comm_ptr picol = mmt->pi->wr[d];
    mmt_comm_ptr mcol = mmt->wr[d];
    if (v == NULL) v = mcol->v;
    if (v->flags & THREAD_SHARED_VECTOR) {
        if (picol->trank == 0)
            v->abase->vec_set_zero(v->abase, v->v, mcol->i1 - mcol->i0);
        serialize_threads(picol);
    } else {
        v->abase->vec_set_zero(v->abase, v->v, mcol->i1 - mcol->i0);
    }
}


/* {{{ broadcast_down (generic interface) */
/* broadcast_down reads data in mmt->wr[d]->v, and broadcasts it across the
 * communicator mmt->pi->wr[d] ; eventually everybody on the communicator
 * mmt->pi->wr[d] has the data.
 *
 * Note that the combination of reduce_across + broadcast_down is not the
 * identity (because of the shuffled_product).
 */

/* XXX
 * ``across'' (horizontally) and ``down'' (vertically) are just here for
 * exposition. The binding of this operation to job/thread arrangement
 * is flexible, through argument d.
 * XXX
 */
void
broadcast_down(matmul_top_data_ptr mmt, mmt_vec_ptr v, int d)
{
    if (v == NULL) v = mmt->wr[d]->v;

    // broadcasting down columns is for common agreement on a right
    // source vector, once a left output has been computed.
    //
    // broadcasting down a column is when d == 1.
    int err;

    mmt_comm_ptr mcol = mmt->wr[d];
    // mmt_comm_ptr mrow = mmt->wr[!d];
    
    pi_comm_ptr picol = mmt->pi->wr[d];
#ifndef MPI_LIBRARY_MT_CAPABLE
    pi_comm_ptr pirow = mmt->pi->wr[!d];
#endif

    /* The intersection of our column-wise input range [i0..i1[ with
     * all the fences existing for the output ranges is described by
     * mcol->x ; note how this intersection is common to a
     * complete column.
     */

    pi_log_op(mmt->pi->m, "[%s] enter first loop", __func__);
    /* Make sure that no thread on the column is wandering in other
     * places -- when we're leaving reduce, this is important. */
    serialize_threads(picol);

    /* This loop suffers from two-dimensional serializing, so the
     * simplistic macros SEVERAL_THREADS_PLAY_MPI_BEGIN and END do not
     * help here.
     */

#ifndef MPI_LIBRARY_MT_CAPABLE
    if ((v->flags & THREAD_SHARED_VECTOR) == 0) {
        size_t eblock = (mcol->i1 - mcol->i0) / picol->totalsize;
        if (picol->trank) {
            int pos = picol->jrank * picol->ncores + picol->trank;
            int offset = pos * eblock;
            v->abase->vec_set(v->abase,
                    SUBVEC(v, all_v[0], offset),
                    SUBVEC(v, v, offset),
                    eblock);
        }
        serialize_threads(picol);
    }
    if (picol->trank == 0) {
#if AG_CHOICE == AG_CHOICE_STOCK_IAG
        MPI_Request * req = shared_malloc(pirow, pirow->ncores * sizeof(MPI_Request));
#endif  /* AG_CHOICE == AG_CHOICE_STOCK_IAG */

        for(unsigned int t = 0 ; t < pirow->ncores ; t++) {
            pi_log_op(pirow, "[%s] serialize_threads", __func__);
            serialize_threads(pirow);
            pi_log_op(pirow, "[%s] serialize_threads done", __func__);
            if (t != pirow->trank)
                continue;   // not our turn.
            // although the openmpi man page looks funny, I'm assuming that
            // MPI_Allgather wants MPI_IN_PLACE as a sendbuf argument.
            pi_log_op(picol, "[%s] MPI_Allgather", __func__);
            size_t eblock = mmt_my_own_size_in_items(mmt, d);
#if AG_CHOICE == AG_CHOICE_STOCK_IAG
            err = MPI_Iallgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, v->v, v->abase->vec_elt_stride(v->abase, 1) * eblock * picol->ncores, MPI_BYTE, picol->pals, &req[t]);
#elif AG_CHOICE == AG_CHOICE_STOCK_AG
            err = MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, v->v, v->abase->vec_elt_stride(v->abase, 1) * eblock * picol->ncores, MPI_BYTE, picol->pals);
#else   /* AG_CHOICE */
#error "Bad AG_CHOICE setting"
#endif  /* AG_CHOICE */
            pi_log_op(picol, "[%s] MPI_Allgather done", __func__);
            ASSERT_ALWAYS(!err);
        }

#if AG_CHOICE == AG_CHOICE_STOCK_IAG
        serialize_threads(pirow);
        if (pirow->trank == 0) {
            for(unsigned int t = 0 ; t < pirow->ncores ; t++) {
                MPI_Wait(&req[t], MPI_STATUS_IGNORE);
            }
        }
        shared_free(pirow, req);
#endif  /* AG_CHOICE == AG_CHOICE_STOCK_IAG */
    }
    if ((v->flags & THREAD_SHARED_VECTOR) == 0) {
        serialize_threads(picol);
        if (picol->trank) {
            v->abase->vec_set(v->abase,
                    v->v, v->all_v[0], mcol->i1 - mcol->i0);
        }
    }
#else   /* MPI_LIBRARY_MT_CAPABLE */
    /* Code deleted 20110119, as I've never been able to have enough
     * trust in an MPI implementation to check this */
    ASSERT_ALWAYS(0);
#endif  /* MPI_LIBRARY_MT_CAPABLE */

    pi_log_op(mmt->pi->m, "[%s] trailer", __func__);
}
/* }}} */

/* {{{ generic interfaces for load/save */
/* {{{ load (generic) */
void mmt_vec_load(matmul_top_data_ptr mmt, mmt_vec_ptr v, const char * name, int d, unsigned int iter, unsigned int itemsondisk)
{
    if (v == NULL) v = mmt->wr[d]->v;

    serialize(mmt->pi->m);

    size_t sizeondisk = v->abase->vec_elt_stride(v->abase, itemsondisk);
    void * mychunk = mmt_my_own_subvec(mmt, v, d);
    size_t mysize = mmt_my_own_size_in_bytes(mmt, v, d);

    char * filename;
    int rc = asprintf(&filename, "%s.%u", name, iter);
    ASSERT_ALWAYS(rc >= 0);

    pi_file_handle f;
    if (!pi_file_open(f, mmt->pi, d, filename, "rb"))
        goto mmt_vec_load_error;

    ssize_t s = pi_file_read(f, mychunk, mysize, sizeondisk);
    pi_file_close(f);

    if (s < 0 || (size_t) s < sizeondisk) {
        goto mmt_vec_load_error;
    }

    free(filename);

    /* not clear it's useful, but well. */
    broadcast_down(mmt, v, d);
    serialize_threads(mmt->pi->m);
    return;
mmt_vec_load_error:
    if (mmt->pi->m->trank == 0 && mmt->pi->m->jrank == 0) {
        fprintf(stderr, "ERROR: failed to load %s\n", filename);
    }
    if (mmt->pi->m->trank == 0)
        MPI_Abort(mmt->pi->m->pals, EXIT_FAILURE);
    serialize_threads(mmt->pi->m);
}
/* }}} */
/* {{{ save (generic) */
void mmt_vec_save(matmul_top_data_ptr mmt, mmt_vec_ptr v, const char * name, int d, unsigned int iter, unsigned int itemsondisk)
{
    if (v == NULL) v = mmt->wr[d]->v;

    /* XXX in shuffled mode, this should be useless. in non-shuffled
     * mode, it's not completely clear. */
    // we want row 0 to have everything.
    broadcast_down(mmt, v, d);
    serialize_threads(mmt->pi->m);

    size_t sizeondisk = v->abase->vec_elt_stride(v->abase, itemsondisk);
    void * mychunk = mmt_my_own_subvec(mmt, v, d);
    size_t mysize = mmt_my_own_size_in_bytes(mmt, v, d);

    char * filename;
    int rc = asprintf(&filename, "%s.%u", name, iter);
    ASSERT_ALWAYS(rc >= 0);

    pi_file_handle f;
    pi_file_open(f, mmt->pi, d, filename, "wb");
    ssize_t s = pi_file_write(f, mychunk, mysize, sizeondisk);
    pi_file_close(f);

    if (s < 0 || (size_t) s < sizeondisk) {
        if (mmt->pi->m->trank == 0 && mmt->pi->m->jrank == 0) {
            fprintf(stderr, "WARNING: failed to save %s\n", filename);
            unlink(filename);
        }
    }

    free(filename);
}
/* }}} */
/* }}} */

//////////////////////////////////////////////////////////////////////////

#if 0/*{{{*/
/* no longer used -- was only used by prep.
 * It's not buggy, but making this work in a context where we have
 * multiple threads is tricky.
 */
void matmul_top_fill_random_source_generic(matmul_top_data_ptr mmt, size_t stride, mmt_vec_ptr v, int d)
{
    if (v == NULL) v = mmt->wr[d]->v;

    // In conjugation mode, it is possible to fill exactly the data chunk
    // that will eventually be relevant. However, it's easy enough to
    // fill our output vector with garbage, and do broadcast_down
    // afterwards...
    if ((v->flags & THREAD_SHARED_VECTOR) == 0 || mmt->pi->wr[d]->trank == 0)
        mpfq_generic_random(stride, v->v, mmt->wr[d]->i1 - mmt->wr[d]->i0);

    // reconcile all cells which correspond to the same vertical block.
    broadcast_down(mmt, v, d);
}
#endif/*}}}*/

/* {{{ reduce_across */
/* {{{ various reduce_scatter implementations */
/* {{{ alternative_reduce_scatter */
#if RS_CHOICE == RS_CHOICE_MINE
void alternative_reduce_scatter(matmul_top_data_ptr mmt, mmt_vec_ptr v, int eitems, int d)
{
    pi_comm_ptr wr = mmt->pi->wr[d];
    int njobs = wr->njobs;
    int rank = wr->jrank;
    MPI_Datatype t = v->pitype->datatype;

    /* If the rsbuf[] buffers have not yet been allocated, it is time to
     * do so now. We also take the opportunity to possibly re-allocate
     * them if because of a larger abase, the corresponding storage has
     * to be expanded.
     */
    size_t needed = v->abase->vec_elt_stride(v->abase, 
        mmt->n[d] / mmt->pi->m->totalsize * mmt->pi->wr[d]->ncores);
    if (mmt->wr[d]->rsbuf_size < needed) {
        mmt->wr[d]->rsbuf[0] = malloc(needed);
        mmt->wr[d]->rsbuf[1] = malloc(needed);
        mmt->wr[d]->rsbuf_size = needed;
    }

    void *b[2];
    b[0] = mmt->wr[d]->rsbuf[0];
    b[1] = mmt->wr[d]->rsbuf[1];
    v->abase->vec_set_zero(v->abase, b[0], eitems);

    int l = (rank + 1) % njobs;
    int srank = (rank + 1) % njobs;
    int drank = (rank + njobs - 1) % njobs;

    for (int i = 0; i < njobs; i++) {
        int j0, j1;
        j0 = l * eitems; 
        j1 = j0 +  eitems;
        v->abase->vec_add(v->abase, b[0], b[0], SUBVEC(v, all_v[0], j0), j1-j0);
        if (i == njobs - 1)  
            break;
        MPI_Sendrecv(b[0], eitems, t, drank, (i<<16) + rank,
                     b[1], eitems, t, srank, (i<<16) + srank,
                     wr->pals, MPI_STATUS_IGNORE);
        void * tb = b[0];   b[0] = b[1];  b[1] = tb; 
        l = (l + 1) % njobs;
    }
    v->abase->vec_set(v->abase, v->all_v[0], b[0], eitems);
}
#endif  /* RS_CHOICE == RS_CHOICE_MINE */
/* }}} */
/* {{{ alternative_reduce_scatter_parallel */
#if RS_CHOICE == RS_CHOICE_MINE_PARALLEL
/* Example data for a factoring matrix (rsa100) of size 135820*135692,
 * split over 2x3 mpi jobs, and 7x5 threads.
 *
 * 2 == mmt->pi->wr[1]->njobs (number of jobs encountered on a vertical axis).
 * 7 == mmt->pi->wr[1]->ncores (number of jobs per core on a vertical axis).
 * 3 == mmt->pi->wr[0]->njobs (number of jobs encountered on an horiz. axis).
 * 5 == mmt->pi->wr[0]->ncores (number of jobs per core on an horiz. axis).
 *
 * matrix is padded to a multiple of 210 = 2*3*5*7, which is * N=135870=210*647
 *
 * we'll call 647 the "small chunk" size.
 *
 * for all jobs/threads, the following relations hold:
 *      mmt->wr[0]->i1 - mmt->wr[0]->i0 == N/14 == 9705 == 15 * 647
 *      mmt->wr[1]->i1 - mmt->wr[1]->i0 == N/15 == 9058 == 14 * 647
 *
 * a reduce_across operation, in the context of factoring, is with d==1
 * below. Hence, in fact, we're doing a reduction down a column.
 *
 * the value eitems fed to this function is mmt->pi->wr[d]->ncores (here,
 * 7) times the small chunk size. Here 7*647 == 4529.
 */
/* all threads in mmt->wr[!d], one after another a priori, are going to
 * do alternative_reduce_scatter on their vector v[i]
 */
void alternative_reduce_scatter_parallel(matmul_top_data_ptr mmt, mmt_vec_ptr * vs, int eitems, int d)
{
    /* we write all data counts below with comments indicating the typical
     * size in our toy example above */

    mpfq_vbase_ptr ab = vs[0]->abase;
    pi_comm_ptr xr = mmt->pi->wr[!d]; /* 3 jobs, 5 cores */
    pi_comm_ptr wr = mmt->pi->wr[d];  /* 2 jobs, 7 cores */
    /* what we're going to do will happen completely in parallel over
     * xr->njobs==3 communicators. We're simulating here the split into 5
     * communicators, even though those collide into a unique
     * communicator as far as MPI is concerned.
     */
    int njobs = wr->njobs;      /* 2 */
    /* note that we no longer care at this point about what happened at
     * the thread level in our dimension. This is already done, period.
     */
    int rank = wr->jrank;
    MPI_Datatype t = vs[0]->pitype->datatype;

    /* If the rsbuf[] buffers have not yet been allocated, it is time to
     * do so now. We also take the opportunity to possibly re-allocate
     * them if because of a larger abase, the corresponding storage has
     * to be expanded.
     */
    size_t needed_per_rsthread = ab->vec_elt_stride(ab, 
        mmt->n[d] / mmt->pi->m->totalsize * wr->ncores);
    if (mmt->wr[d]->rsbuf_size < needed_per_rsthread) {
        mmt->wr[d]->rsbuf[0] = malloc(needed_per_rsthread);
        mmt->wr[d]->rsbuf[1] = malloc(needed_per_rsthread);
        mmt->wr[d]->rsbuf_size = needed_per_rsthread;
    }

    /* are we sure ? */
    /* caller sets:
     * eitems == (mmt->wr[1]->i1-mmt->wr[1]->i0) / mmt->pi->wr[1]->njobs
     *        == mmt->n[1] / mmt->pi->wr[0]->totalsize / mmt->pi->wr[1]->njobs
     *        == mmt->n[1] / mmt->pi->m->totalsize * mmt->pi->wr[1]->ncores
     *
     * as long as the identity
     *  mmt->wr[1]->(i1-i0) == mmt->n[1] / mmt->pi->wr[0]->totalsize
     * is maintained throughout the program execution, we are safe on the
     * assert below.
     */
    ASSERT_ALWAYS(needed_per_rsthread == (size_t) ab->vec_elt_stride(ab, eitems));       /* 4529 */


    void *(*foo)[2] = shared_malloc(xr, xr->ncores /* 5 */ * sizeof(*foo));
    void ** mine = foo[xr->trank];
    ab->vec_set_zero(ab, mmt->wr[d]->rsbuf[0], eitems);
    mine[0] = mmt->wr[d]->rsbuf[0];
    mine[1] = mmt->wr[d]->rsbuf[1];
    serialize_threads(xr);

    int srank = (rank + 1) % njobs;
    int drank = (rank + njobs - 1) % njobs;

    /* We describe the algorithm for one of the xr->ncores==5 threads.
     * local vector areas [9058] split into [wr->njobs==2] sub-areas of
     * size [eitems (==4529)].
     *
     * on each job, each thread has 9058 == 2*4529 items
     *          more generally: njobs * eitems items.
     *
     * each has a place where it wants to receive data, let's call that
     * "his" zone, even though everyone has data at all places.
     *
     */
    /* scheme for 4 jobs: (input, output. All chunks have size eitmes. AA
     * is the sum a3+a2+a1+a0)
     *
     * j0: a0 b0 c0 d0      j0: AA .. .. ..
     * j1: a1 b1 c1 d1      j1: .. BB .. ..
     * j2: a2 b2 c2 d2      j2: .. .. CC ..
     * j3: a3 b3 c3 d3      j3: .. .. .. DD
     *
     * loop 0: j0: sends       b0 to j3, receives       c1 from j1
     * loop 1: j0: sends    c1+c0 to j3, receives    d2+d1 from j1
     * loop 2: j0: sends d2+d1+d0 to j3, receives a3+a2+a1 from j1
     * loop 3: we only have to adjust.
     */

    for (int i = 0, s = 0; i < njobs; i++, s=!s) {
        int l = (rank + 1 + i) % njobs;
        int j0 = l * eitems; 
        int j1 = j0 +  eitems;
        ab->vec_add(ab, mine[s], mine[s],
                ab->vec_subvec(ab, vs[xr->trank]->all_v[0], j0),
                j1-j0);
        serialize_threads(xr);

        if (i == njobs - 1)  
            break;

        if (xr->trank == 0) {
            MPI_Request * r = NULL;
            r = malloc(2 * xr->ncores * sizeof(MPI_Request));
            for(unsigned int w = 0 ; w < xr->ncores ; w++) {
                MPI_Request * rs = r + 2*w;
                MPI_Request * rr = r + 2*w + 1;
                MPI_Isend(foo[w][s],  eitems, t, drank, 0xb00+w, wr->pals, rs);
                MPI_Irecv(foo[w][!s], eitems, t, srank, 0xb00+w, wr->pals, rr);
                /*
                MPI_Sendrecv(foo[w][s], eitems, t, drank, 0xbeef,
                        foo[w][!s], eitems, t, srank, 0xbeef,
                        wr->pals, MPI_STATUS_IGNORE);
                        */
                // MPI_Waitall(2, r + 2*w, MPI_STATUSES_IGNORE);
            }
            MPI_Waitall(2 * xr->ncores, r, MPI_STATUSES_IGNORE);
            free(r);
        }
        serialize_threads(xr);
    }

    ab->vec_set(ab, vs[xr->trank]->all_v[0], mine[(njobs-1)&1], eitems);

    shared_free(xr, foo);
    pi_log_op(wr, "[%s] MPI_Reduce_scatter done", __func__);
}
#endif  /* RS_CHOICE == RS_CHOICE_MINE_PARALLEL */
/* }}} */
/* {{{ my_MPI_Reduce_scatter_block */
#if RS_CHOICE == RS_CHOICE_MINE_DROP_IN
int my_MPI_Reduce_scatter_block(void *sendbuf, void *recvbuf, int recvcount,
                MPI_Datatype datatype, MPI_Op op, MPI_Comm wr)
{
    int njobs;
    int rank;
    MPI_Comm_size(wr, &njobs);
    MPI_Comm_rank(wr, &rank);

    int tsize;
    MPI_Type_size(datatype, &tsize);

    assert(sendbuf == MPI_IN_PLACE);
    void * v = recvbuf;
    
    /* This is a deliberate leak. Note that we expect to be serailized
     * here, so there is no concurrency issue with the static data. */
    static size_t rsbuf_size = 0;
    static unsigned long * rsbuf[2];

    size_t needed = recvcount * tsize;

    if (rsbuf_size < needed) {
        rsbuf[0] = malloc(needed);
        rsbuf[1] = malloc(needed);
        rsbuf_size = needed;
    }

    memset(rsbuf[0], 0, recvcount * tsize);

    int srank = (rank + 1) % njobs;
    int drank = (rank + njobs - 1) % njobs;

    for (int i = 0, w = 0; i < njobs; i++, w^=1) {
        int j0 = ((rank + i + 1) % njobs) * recvcount;
        void * share = pointer_arith_const(v, j0 * tsize);
#if MPI_VERSION_ATLEAST(2,2)
        MPI_Reduce_local(share, rsbuf[w], recvcount, datatype, op);
#else
        {
            ASSERT_ALWAYS(datatype == MPI_UNSIGNED_LONG);
            ASSERT_ALWAYS(op == MPI_BXOR);
            unsigned long * a = share;
            unsigned long * b = rsbuf[w];
            for(int k = 0 ; k < recvcount ; k++) {
                b[k] ^= a[k];
            }
        }
#endif
        if (i == njobs - 1) {
            memcpy(v, rsbuf[w], recvcount * tsize);
            break;
        }
        MPI_Sendrecv(rsbuf[w],  recvcount, datatype, drank, (i<<16) + rank,
                     rsbuf[!w], recvcount, datatype, srank, (i<<16) + srank,
                     wr, MPI_STATUS_IGNORE);
    }
    return 0;
}
#endif
/* }}} */
/* }}} */

/* reduce_across_inner reads data in v (vector for side d), sums it up
 * across the communicator mmt->pi->wr[d], and collects the results in
 * vector v again, except that it's put in thread0's buffer (counting
 * threads in direction d, of course), *AT THE BEGINNING* of the data
 * area (which is surprising).
 *
 * reduce_across completes the work by saving the resulting data in vector w
 * (vector in dimension !d).
 *
 * Note that the combination of reduce_across + broadcast_down is not the
 * identity (because of the shuffled_product).
 */
void
reduce_across_inner(matmul_top_data_ptr mmt, mmt_vec_ptr v, int d)
{
    if (v == NULL) v = mmt->wr[d]->v;

    /* reducing across a row is when d == 0 */
    pi_comm_ptr pirow = mmt->pi->wr[d];
#ifndef MPI_LIBRARY_MT_CAPABLE
    pi_comm_ptr picol = mmt->pi->wr[!d];
#endif

    mmt_comm_ptr mrow = mmt->wr[d];

    unsigned int z = mmt->n[d] / mmt->pi->m->totalsize;
    unsigned int z1 = z * pirow->njobs;

    pi_log_op(mmt->pi->m, "[%s] enter first loop", __func__);

    // I don't think that the case of shared vectors has been tested
    // correctly for reduction. Well, to start with, I doubt it really
    // makes a lot of sense anyhow.
    // ASSERT_ALWAYS((v->flags & THREAD_SHARED_VECTOR) == 0);

    if (pirow->ncores > 1 && !(v->flags & THREAD_SHARED_VECTOR)) {
        /* row threads have to sum up their data. Of course it's
         * irrelevant when there is only one such thread...
         *
         * Concerning locking, we have to make sure that everybody on the
         * row has finished its computation task, but besides that,
         * there's no locking until we start mpi stuff.
         */
        pi_log_op(pirow, "[%s] serialize_threads", __func__);
        serialize_threads(pirow);
        pi_log_op(pirow, "[%s] serialize_threads done", __func__);

        /* Our [i0,i1[ range is split into pirow->ncores parts. This range
         * represent coordinates which are common to all threads on
         * pirow. Corresponding data has to be summed. Amongst the
         * pirow->ncores threads, thread k takes the responsibility of
         * summing data in the k-th block (that is, indices
         * i0+k*(i1-i0)/pirow->ncores and further). As a convention,
         * the thread which eventually owns the final data is thread 0.
         *
         * Note that one should consider that data in threads other than
         * the destination thread may be clobbered by the operation
         * (although in the present implementation it is not).
         */
        unsigned int ii0 = mrow->i0 + pirow->trank * z1;
        unsigned int ii1 = ii0 + z1;
        void * dptr = SUBVEC(v, all_v[0], ii0 - mrow->i0);
        for(uint32_t w = 1 ; w < pirow->ncores ; w++) {
            const void * sptr = SUBVEC(v, all_v[w], ii0 - mrow->i0);
            v->abase->vec_add(v->abase, dptr, dptr, sptr, ii1-ii0);
        }
        pi_log_op(pirow, "[%s] thread reduction done", __func__);
    }

    /* Good. Now on each node, thread 0 has the reduced data for the
     * range [i0..i1[ (well, actually, this corresponds to indices
     * distributed in a funny manner in the original matrix, but as far
     * as we care here, it's really the range mrow->i0..mrow->i1
     *
     * These data areas must now be reduced or reduce-scattered across
     * nodes. Note that in the case where the MPI library does not
     * support MT operation, we are performing the reduction in place,
     * and copy to the buffer in the other direction is done in a
     * secondary step -- which will be among the reasons which imply a
     * thread serialization before anything interesting can be done after
     * this function.
     */

    /* XXX We have to serialize threads here. At least row-wise, so that
     * the computations above are finished. Easily seen, that's just a
     * couple of lines above.
     *
     * Unfortunately, that's not the whole story. We must also serialize
     * column-wise. Because we're writing on the right vector buffers,
     * which are used as input for the matrix multiplication code --
     * which may be still running at this point because there has been no
     * column serialization in the current procedure yet.
     */

    pi_log_op(mmt->pi->m, "[%s] secondary loop", __func__);

    pi_log_op(mmt->pi->m, "[%s] serialize_threads", __func__);
    serialize_threads(pirow);
    pi_log_op(mmt->pi->m, "[%s] serialize_threads done", __func__);

#ifndef MPI_LIBRARY_MT_CAPABLE
    if (pirow->trank == 0) {
        /* openmpi-1.8.2 does not seem to have a working non-blocking
         * reduce_scatter, at least not a very efficient one. All
         * I've been able to do is to run MPI_Ireduce_scatter with
         * MPI_UNSIGNED_LONG and MPI_BXOR. With custom types it seems
         * to crash. And anyway, it's very inefficient.
         */

        ASSERT((mrow->i1 - mrow->i0) % pirow->totalsize == 0);
#if RS_CHOICE == RS_CHOICE_STOCK_RS 
        int z = (mrow->i1 - mrow->i0) / pirow->njobs;
        void * dptr = v->all_v[0];
        SEVERAL_THREADS_PLAY_MPI_BEGIN(picol) {
            // all recvcounts are equal
            int * rc = malloc(pirow->njobs * sizeof(int));
            for(unsigned int k = 0 ; k < pirow->njobs ; rc[k++]=z) ;
            err = MPI_Reduce_scatter(dptr, dptr, rc,
                    v->pitype->datatype,
                    BWC_PI_SUM->custom,
                    pirow->pals);
            free(rc);
            ASSERT_ALWAYS(!err);
        }
        SEVERAL_THREADS_PLAY_MPI_END();
#elif RS_CHOICE == RS_CHOICE_STOCK_RSBLOCK
        int z = (mrow->i1 - mrow->i0) / pirow->njobs;
        void * dptr = v->all_v[0];
        SEVERAL_THREADS_PLAY_MPI_BEGIN(picol) {
            err = MPI_Reduce_scatter_block(dptr, dptr, z,
                    v->pitype->datatype,
                    BWC_PI_SUM->custom,
                    pirow->pals);
            ASSERT_ALWAYS(!err);
        }
        SEVERAL_THREADS_PLAY_MPI_END();
#elif RS_CHOICE == RS_CHOICE_MINE_DROP_IN
        int z = (mrow->i1 - mrow->i0) / pirow->njobs;
        void * dptr = v->all_v[0];
        SEVERAL_THREADS_PLAY_MPI_BEGIN(picol) {
            err = my_MPI_Reduce_scatter_block(dptr, dptr, z,
                    v->pitype->datatype,
                    BWC_PI_SUM->custom,
                    pirow->pals);
            ASSERT_ALWAYS(!err);
        }
        SEVERAL_THREADS_PLAY_MPI_END();
#elif RS_CHOICE == RS_CHOICE_STOCK_IRSBLOCK
        int z = (mrow->i1 - mrow->i0) / pirow->njobs;
        void * dptr = v->all_v[0];
        MPI_Request * req = shared_malloc(picol, picol->ncores * sizeof(MPI_Request));
        SEVERAL_THREADS_PLAY_MPI_BEGIN(picol) {
            err = MPI_Ireduce_scatter(dptr, dptr, z,
                    v->pitype->datatype,
                    BWC_PI_SUM->custom,
                    pirow->pals, &req[t__]);
            ASSERT_ALWAYS(!err);
            pi_log_op(pirow, "[%s] MPI_Reduce_scatter done", __func__);
        }
        SEVERAL_THREADS_PLAY_MPI_END();
        serialize_threads(picol);
        if (picol->trank == 0) {
            for(unsigned int t = 0 ; t < picol->ncores ; t++) {
                MPI_Wait(&req[t], MPI_STATUS_IGNORE);
            }
        }
        shared_free(picol, req);
#elif RS_CHOICE == RS_CHOICE_MINE
        int z = (mrow->i1 - mrow->i0) / pirow->njobs;
        /* This strategy exposes code which is really similar to
         * RS_CHOICE_MINE_DROP_IN, with the only exception that we
         * have a slightly different interface. There's no reason for
         * both to stay in the long run.
         */
        SEVERAL_THREADS_PLAY_MPI_BEGIN(picol) {
            alternative_reduce_scatter(mmt, v, z, d);
        }
        SEVERAL_THREADS_PLAY_MPI_END();
#elif RS_CHOICE == RS_CHOICE_MINE_PARALLEL
        int z = (mrow->i1 - mrow->i0) / pirow->njobs;
        mmt_vec_ptr * vs = shared_malloc(picol, picol->ncores * sizeof(mmt_vec_ptr));
        vs[picol->trank] = v;
        serialize_threads(picol);
        alternative_reduce_scatter_parallel(mmt, vs, z, d);
        shared_free(picol, vs);
#elif RS_CHOICE == RS_CHOICE_MINE_OVERLAPPING
#error "not implemented, but planned"
#endif
    }
    serialize_threads(pirow);
}
void
reduce_across(matmul_top_data_ptr mmt, mmt_vec_ptr w, mmt_vec_ptr v, int d)
{
    reduce_across_inner(mmt, v, d);
    if (v == NULL) v = mmt->wr[d]->v;
    if (w == NULL) w = mmt->wr[!d]->v;
    // row threads pick what they're interested in in thread0's reduced
    // buffer. Writes are non-overlapping in the mcol buffer here.
    // Different row threads always have different mcol buffers, and
    // sibling col threads write to different locations in their
    // (generally shared) mcol buffer, depending on which row they
    // intersect.

    // Notice that results are being written inside v->all_v[0],
    // just packed at the beginning.

    // row job rj, row thread rt has a col buffer containing
    // picol->totalsize blocks of size eblock.  col job cj, col
    // thread ct has in its row buffer pirow->totalsize blocks, but
    // only virtually. Because of the reduce_scatter operation, only
    // pirow->ncores blocks are here.
    //
    // Thus the k-th block in the row buffer is rather understood as
    // the one of indek rj * pirow->ncores + k in the data the row
    // threads have collectively computed.
    //
    // This block is of interest to the row thread of index k of
    // course, thus we restrict to rt==k
    //
    // Now among the picol->totalsize blocks of the col buffer, this
    // will go to position cj * picol->ncores + ct
 
    size_t eblock = mmt_my_own_size_in_items(mmt, d);
    ASSERT_ALWAYS(mmt_my_own_size_in_items(mmt, !d) == eblock);

    v->abase->vec_set(v->abase,
            mmt_my_own_subvec(mmt, w, !d),
            /* Note: reduce-scatter packs everything at the beginning in
             * all_v[0], which is why we don't have our usual offset
             * here. */
            SUBVEC(v, all_v[0], mmt->pi->wr[d]->trank * eblock),
            eblock);
#else   /* MPI_LIBRARY_MT_CAPABLE */
    /* Code deleted 20110119, as I've never been able to have enough
     * trust in an MPI implementation to check this */
    ASSERT_ALWAYS(0);
#endif

    // as usual, we do not serialize on exit. Up to the next routine to
    // do so if needed.
    // 
    // what _is_ guaranteed is that in each column, the _leader_ thread
    // has all the necessary data to begin the column broadcast.
    //
    // In most cases, threads must be prevented from starting computation
    // before the leader thread has finished importing data. This means
    // that a column thread serialization is probably needed in most
    // circumstances after this step.
}

/* This small variant is used only for the transposition routines. We do
 * not want, in that case, to rely on the data moving back and forth
 * between the left and right vectors, because in full generality, we
 * have no guarantee that they are the same size.
 *
 * Therefore we're merely building upon reduce_across_inner, only with
 * the weirdo "pack-at-the-beginning" behaviour removed.
 */
void
reduce_across_sameside(matmul_top_data_ptr mmt, mmt_vec_ptr v, int d)
{
    reduce_across_inner(mmt, v, d);
    if (v == NULL) v = mmt->wr[d]->v;
    size_t eblock = mmt_my_own_size_in_items(mmt, d);
    v->abase->vec_set(v->abase,
            mmt_my_own_subvec(mmt, v, d),
            SUBVEC(v, all_v[0], mmt->pi->wr[d]->trank * eblock),
            eblock);
    if (mmt->pi->wr[d]->jrank) {
        v->abase->vec_set_zero(v->abase,
                SUBVEC(v, all_v[0], mmt->pi->wr[d]->trank * eblock),
                eblock);
    }
    serialize_threads(mmt->pi->m);
}
/* }}} */

/* {{{ allreduce_across */
/* This is only for convenience now. Eventually this will be relevant for
 * block Lanczos.  Note that allreduce is conceptually much simpler.
 * There is no funny permutation to be considered.
 */
void
allreduce_across(matmul_top_data_ptr mmt, mmt_vec_ptr v, int d)
{
    /* reducing across a row is when d == 0 */
    pi_comm_ptr pirow = mmt->pi->wr[d];

    mmt_comm_ptr mrow = mmt->wr[d];
    if (v == NULL) v = mrow->v;

    unsigned int z = mmt->n[d] / mmt->pi->m->totalsize;
    unsigned int z1 = z * pirow->njobs;

    pi_log_op(mmt->pi->m, "[%s] enter first loop", __func__);

    serialize_threads(mmt->pi->m);
    /* sum up row threads, so that only one thread on each row is used
     * for communication */
    unsigned int ii0 = mrow->i0 + pirow->trank * z1;
    unsigned int ii1 = ii0 + z1;
    if (!(v->flags & THREAD_SHARED_VECTOR)) {
        for(unsigned int k = 1 ; k < pirow->ncores ; k++) {
            void * dv = SUBVEC(v, v, ii0 - mrow->i0);
            void * sv = SUBVEC(v, all_v[(pirow->trank+k) % pirow->ncores], ii0 - mrow->i0);
            v->abase->vec_add(v->abase, dv, dv, sv, ii1 - ii0);
        }
    }
    SEVERAL_THREADS_PLAY_MPI_BEGIN(mmt->pi->m) {
        MPI_Allreduce(MPI_IN_PLACE,
                SUBVEC(v, v, ii0 - mrow->i0),
                z1,
                v->pitype->datatype,
                BWC_PI_SUM->custom,
                pirow->pals);
    }
    SEVERAL_THREADS_PLAY_MPI_END();
    if (!(v->flags & THREAD_SHARED_VECTOR)) {
        for(unsigned int k = 1 ; k < pirow->ncores ; k++) {
            void * sv = SUBVEC(v, v, ii0 - mrow->i0);
            void * dv = SUBVEC(v, all_v[(pirow->trank+k)%pirow->ncores], ii0 - mrow->i0);
            v->abase->vec_set(v->abase, dv, sv, ii1 - ii0);
        }
    }
}
/* }}} */

/**********************************************************************/
/* bench code */

void matmul_top_comm_bench_helper(int * pk, double * pt,
                                  void (*f) (matmul_top_data_ptr, mmt_vec_ptr, int),
				  matmul_top_data_ptr mmt, int d)
{
    int k;
    double t0, t1;
    int cont;
    t0 = wct_seconds();
    for (k = 0;; k++) {
	t1 = wct_seconds();
	cont = t1 < t0 + 0.25;
	cont = cont && (t1 < t0 + 1 || k < 100);
        pi_allreduce(&cont, &cont, 1, BWC_PI_INT, BWC_PI_MIN, mmt->pi->m);
	if (!cont)
	    break;
	(*f) (mmt, NULL, d);
    }
    int target = 10 * k / (t1 - t0);
    ASSERT_ALWAYS(target >= 0);
    if (target > 100)
	target = 100;
    if (target == 0)
        target = 1;
    pi_bcast(&target, 1, BWC_PI_INT, 0, 0, mmt->pi->m);
    t0 = wct_seconds();
    for (k = 0; k < target; k++) {
        pi_log_op(mmt->pi->m, "[%s] iter%d/%d", __func__, k, target);
        (*f) (mmt, NULL, d);
    }
    serialize(mmt->pi->m);
    t1 = wct_seconds();
    *pk = k;
    *pt = t1 - t0;
}


void matmul_top_comm_bench(matmul_top_data_ptr mmt, int d)
{
    /* like matmul_top_mul_comm, we'll call reduce_across with !d, and
     * broadcast_down with d */
    int k;
    double dt;

    void (*funcs[2])(matmul_top_data_ptr, mmt_vec_ptr, int) = {
        broadcast_down,
        reduce_across_inner,
    };
    const char * text[2] = { "bd", "ra" };

    mpfq_vbase_ptr abase = mmt->abase;

    mmt_comm_ptr mrow = mmt->wr[!d];
    mmt_comm_ptr mcol = mmt->wr[d];
    pi_comm_ptr pirow = mmt->pi->wr[!d];
    pi_comm_ptr picol = mmt->pi->wr[d];

    /* within each row, all jobs are concerned with the same range
     * mrow->i0 to mrow->i1. This is split into m=pirow->njobs chunks,
     * and reduce_scatter has m-1 communication rounds where all of the m
     * nodes output and receive one such chunk. All this happens for all
     * column threads, so we multiply by picol->ncores, too.
     */
    size_t data_out_ra = abase->vec_elt_stride(abase, picol->ncores * (mrow->i1 - mrow->i0) / pirow->njobs * (pirow->njobs - 1));

    /* one way to do all-gather is to mimick this, except that each node
     * will output the same chunk at each round. Beyond that, the
     * calculation is similar, and we'll use it as a guide. Note of
     * course that if hardware-level multicast is used, our throughput
     * estimation is way off.
     */
    size_t data_out_ag = abase->vec_elt_stride(abase, pirow->ncores * (mcol->i1 - mcol->i0) / picol->njobs * (picol->njobs - 1));

    size_t datasize[2] = { data_out_ag, data_out_ra };

    for(int s = 0 ; s < 2 ; s++) {
        /* we have our axis, and the other axis */
        pi_comm_ptr wr = mmt->pi->wr[d ^ s];          /* our axis */
        pi_comm_ptr xr = mmt->pi->wr[d ^ s ^ 1];      /* other axis */
        /* our operation has operated on the axis wr ; hence, we must
         * display data relative to the different indices within the
         * communicator xr.
         */
        matmul_top_comm_bench_helper(&k, &dt, funcs[s], mmt, d^s);
        if (verbose_enabled(CADO_VERBOSE_PRINT_BWC_TIMING_GRIDS)) {
            for(unsigned int z = 0 ; z < xr->njobs ; z++) {
                if (xr->jrank == z && wr->jrank == 0) {
                    for(unsigned int w = 0 ; w < xr->ncores ; w++) {
                        if (xr->trank == w && wr->trank == 0) {
                            char buf[16];
                            printf("%s %2d/%d, %s: %2d in %.1fs ; one: %.2fs (xput: %s/s)\n",
                                    wr->th->desc,
                                    xr->jrank * xr->ncores + xr->trank,
                                    xr->totalsize,
                                    text[s],
                                    k, dt, dt/k,
                                    size_disp(datasize[d^s] * k/dt, buf));
                        }
                    }
                }
                serialize(mmt->pi->m);
            }
        }
        serialize(mmt->pi->m);
    }
}


/**********************************************************************/
/* Utility stuff for applying standard permutations to vectors.
 *
 * Extra documentation about the permutations S and P can be found in
 * linalg/bwc/README.
 *
 * The permutations Sr and Sc are the row and columns permutations
 * computed by mf_bal (which are forced to be equal in the block
 * Wiedemann case).
 *
 * The matmul_top_mul_comm routine, which is part of matmul_top_mul, does
 * not only a reduction+broadcast, but induces multiplication by a
 * permutation matrix P (operating on the vector area which is both the
 * source and destination of matmul_top_mul). This permutation is known
 * statically, but depends on the splitting.
 *
 * The ``twisted'' vectors are as follows:
 *
 * Square matrices, with Sr==Sc (block Wiedemann case).
 *
 * [TODO: say how I define the matrix associated with a permutation]
 *
 *   The blocks dispatched across the different jobs and threads are not
 *   the blocks of M, nor the blocks of Sr*M*Sc^-1, but those of the
 *   Mt==P*Sr*M*Sc^-1. This way, we are working with iterates of the
 *   matrix P*Sc*M*(P*Sr)^-1 for nullspace=left, and Sr*M*Sc^-1 for
 *   nullspace=right. For Sc==Sr, these are conjugates of M in both
 *   cases.
 *    
 *     Case of nullspace=left, vector times matrix, the "left" vector
 *     mmt->wr[0]->v is represented as a row vector.
 *
 *       [TODO: hmmm. not clear it's P^-1, not P]
 *       The vector times matrix operation also multiplies by P^-1
 *
 *       We have:
 *            Ytwisted = Y * (P*Sc)^-1
 *            Ytwisted * Mt * P^-1 = Y * (P*Sc)^-1 * P*Sr*M*Sc^-1 * P^-1
 *                                 = Y * M * (P*Sc)^-1
 *                                 = (Y*M)twisted
 *
 *     Case of nullspace=right, matrix times vector, the "right" vector is
 *     represented as a column vector:
 *       The matrix times vector operation also multiplies by P.
 *       We have:
 *          Ytwisted = Sc * Y
 *          P * Mt * Ytwisted = Sr*M*Y = (M*Y)_twisted
 *    
 *   Therefore it is possible to twist and untwist files consistently.
 *
 * For rectangular matrices, we have no permutation P, because
 * reduce_across is not called. Therefore the situation is somewhat
 * simpler, except that of course we must be careful to distinguish left
 * and right vectors properly.

 * The functions below apply the transformations to the input vector.
 *
 * None of these functions is optimized for speed. Their purpose is only
 * normalization w.r.t permutations being used.
 */

/* {{{ utility: matmul_top_zero_vec_area */
void matmul_top_zero_vec_area(matmul_top_data_ptr mmt, int d)
{
    mmt_comm_ptr mdst = mmt->wr[d];
    pi_comm_ptr pidst = mmt->pi->wr[d];
    if (mdst->v->flags & THREAD_SHARED_VECTOR) {
        serialize_threads(pidst);
        if (pidst->trank == 0)
            mdst->v->abase->vec_set_zero(mdst->v->abase, mdst->v->v, mdst->i1 - mdst->i0);
        serialize_threads(pidst);
    } else {
        mdst->v->abase->vec_set_zero(mdst->v->abase, mdst->v->v, mdst->i1 - mdst->i0);
    }
}
/* }}} */
/* {{{ utility: matmul_top_restrict_vec_area */
/* On a shared vector which is assumed to be commonly known by all
 * nodes/threads, select only one portion to remain, and zero out the
 * rest (for later use by e.g.  matmul_top_mul_comm or other integrated
 * function).
 */
static void matmul_top_restrict_vec_area(matmul_top_data_ptr mmt, int d)
{
    mmt_comm_ptr mdst = mmt->wr[d];
    pi_comm_ptr pidst = mmt->pi->wr[d];
    if (mdst->v->flags & THREAD_SHARED_VECTOR)
        serialize_threads(pidst);
    if (pidst->trank == 0 || !(mdst->v->flags & THREAD_SHARED_VECTOR)) {
        if (pidst->jrank == 0 && pidst->trank == 0) {
            /* ok, we keep the data */
        } else {
            mdst->v->abase->vec_set_zero(mdst->v->abase, mdst->v->v, mdst->i1 - mdst->i0);
        }
    }
}
/* }}} */

/* {{{ apply_balancing_permutation -- this applies the permutation S */
/* This is sort of an equivalent to matmul_top_cpu, although since it is
 * meant to work on permutation matrices, it works equally well even if
 * output buffers are shared.  */
void apply_balancing_permutation(matmul_top_data_ptr mmt, int d)
{
    mmt_comm_ptr mrow = mmt->wr[0];
    mmt_comm_ptr mcol = mmt->wr[1];

    matmul_top_zero_vec_area(mmt, !d);
    serialize_threads(mmt->pi->m);

    for(uint32_t i = 0 ; i < mmt->mm->ntwists ; i++) {
        uint32_t u = mmt->mm->twist[i][0];
        uint32_t v = mmt->mm->twist[i][1];
        ASSERT_ALWAYS(u >= mrow->i0);
        ASSERT_ALWAYS(u < mrow->i1);
        ASSERT_ALWAYS(v >= mcol->i0);
        ASSERT_ALWAYS(v < mcol->i1);
        if (d==0) {
            mcol->v->abase->vec_set(mcol->v->abase,
                    SUBVEC(mcol->v, v, v - mcol->i0),
                    SUBVEC(mrow->v, v, u - mrow->i0),
                    1);
        } else {
            mrow->v->abase->vec_set(mrow->v->abase,
                    SUBVEC(mrow->v, v, u - mrow->i0),
                    SUBVEC(mcol->v, v, v - mcol->i0),
                    1);
        }
    }
}
/* }}} */

/* {{{ apply_decorrelating_permutation -- this applies the fixed
 * permutation which we use unconditionally in bwc to avoid correlation
 * of row and column weights. Because of this permutation T, we never
 * work with the matrix M the user thinks, but rather with the matrix MT.
 */
/* This is sort of an equivalent to matmul_top_cpu, although since it is
 * meant to work on permutation matrices, it works equally well even if
 * output buffers are shared.  */
void apply_decorrelating_permutation(matmul_top_data_ptr mmt, int d)
{
    /* See matmul_top_apply_T and matmul_top_unapply_T below. Because we
     * use the same mechanism as for matmul_top_apply_S and friends, we
     * must pay attention that the real operation we're concerned with is
     * only with vector for direction d==1 (column vector).
     *
     * At the point where we get here, the situation where the action of
     * T is asked for the vector in direction 0 has already been evicted
     * -- there's an early return in matmul_top_unapply_T. Therefore we
     * must use d (which, by convention, indicates the source argument)
     * to indicate whether we're at the beginning or the end of the
     * chain. If at the beginning, by convention, we're unapplying.
     */

    mmt_comm_ptr msrc = mmt->wr[d];
    mmt_comm_ptr mdst = mmt->wr[!d];

    matmul_top_zero_vec_area(mmt, !d);
    serialize_threads(mmt->pi->m);

    mpfq_vbase_ptr A = mmt->abase;

    for(uint32_t u = msrc->i0 ; u < msrc->i1 ; u++) {
        uint32_t v;
        if (d == 1) {
            /* source is the col vector, hence we're at the beginning of the
             * chain */
            v = balancing_pre_unshuffle(mmt->bal, u);
        } else {
            v = balancing_pre_shuffle(mmt->bal, u);
        }
        if (v < mdst->i0) continue;
        if (v >= mdst->i1) continue;
        A->vec_set(A,
                A->vec_coeff_ptr(A, mdst->v->v, v - mdst->i0),
                A->vec_coeff_ptr_const(A, msrc->v->v, u - msrc->i0),
                1);
    }
}
/* }}} */

/* {{{ This applies the identity matrix (thus copies the vector from [d]
 * to [!d]) */
void apply_identity(matmul_top_data_ptr mmt, int d)
{
    mmt_comm_ptr msrc = mmt->wr[d];
    mmt_comm_ptr mdst = mmt->wr[!d];

    matmul_top_zero_vec_area(mmt, !d);
    serialize_threads(mmt->pi->m);

    for(unsigned int i = msrc->i0 ; i < msrc->i1 ; i++) {
        if (i >= mdst->i0 && i < mdst->i1) {
            mdst->v->abase->vec_set(mdst->v->abase,
                    SUBVEC(mdst->v, v, i - mdst->i0),
                    SUBVEC(msrc->v, v, i - msrc->i0),
                    1);
        }
    }
}
/* }}} */

/* {{{ matmul_top_unapply_S_unapply_P */
// for nullspace=left, bw->dir == 0 : vectors are represented as
// horizontal blocks n*N
// d == 0       v <- v * Sr * Pr
// d == 1       v <- v * Sr^-1 * Pr^-1
//
// for nullspace=right, bw->dir == 1 : vectors are represented as
// vertical blocks N*n
// d == 1       v <- Pr * Sr * v;
// d == 0       v <- Pr^-1 * Sr^-1 * v;
void matmul_top_unapply_S_unapply_P(matmul_top_data_ptr mmt, int d)
{
    // if (mmt->pi->m->jrank == 0 && mmt->pi->m->trank == 0) printf("[%s] d == %d\n", __func__, d);
    apply_balancing_permutation(mmt, d);
    matmul_top_mul_comm(mmt, d);
}
/* }}} */

/* {{{ matmul_top_apply_P_apply_S */
// for nullspace=left, bw->dir == 0 : vectors are represented as
// horizontal blocks n*N
// d == 0       v <- v * Pr * Sr
// d == 1       v <- v * Pr^-1 * Sr^-1
//
// for nullspace=right, bw->dir == 1 : vectors are represented as
// vertical blocks N*n
// d == 1       v <- Sr * Pr * v;
// d == 0       v <- Sr^-1 * Pr^-1 * v;
void matmul_top_apply_P_apply_S(matmul_top_data_ptr mmt, int d)
{
    // if (mmt->pi->m->jrank == 0 && mmt->pi->m->trank == 0) printf("[%s] d == %d\n", __func__, d);
    pi_comm_ptr picol = mmt->pi->wr[d];
    matmul_top_restrict_vec_area(mmt, d);
    serialize_threads(picol);
    matmul_top_mul_comm(mmt, !d);
    apply_balancing_permutation(mmt, !d);
    allreduce_across(mmt, NULL, d);
}
/* }}} */

/* {{{ matmul_top_apply_P */
void matmul_top_apply_P(matmul_top_data_ptr mmt, int d)
{
    // if (mmt->pi->m->jrank == 0 && mmt->pi->m->trank == 0) printf("[%s] d == %d\n", __func__, d);
    matmul_top_restrict_vec_area(mmt, d);
    matmul_top_mul_comm(mmt, !d);
    apply_identity(mmt, !d);
    allreduce_across(mmt, NULL, d);
}
/* }}} */

/* {{{ matmul_top_unapply_P */
void matmul_top_unapply_P(matmul_top_data_ptr mmt, int d)
{
    apply_identity(mmt, d);
    matmul_top_mul_comm(mmt, d);
}
/* }}} */

/* {{{ matmul_top_unapply_S */
void matmul_top_unapply_S(matmul_top_data_ptr mmt, int d)
{
    apply_balancing_permutation(mmt, d);
    allreduce_across(mmt, NULL, !d);
    apply_identity(mmt, !d);
    allreduce_across(mmt, NULL, d);
}
/* }}} */

/* {{{ matmul_top_apply_S */
void matmul_top_apply_S(matmul_top_data_ptr mmt, int d)
{
    apply_identity(mmt, d);
    allreduce_across(mmt, NULL, !d);
    apply_balancing_permutation(mmt, !d);
    allreduce_across(mmt, NULL, d);
}
/* }}} */

/* {{{ matmul_top_unapply_T */
void matmul_top_unapply_T(matmul_top_data_ptr mmt, int d)
{
    if (d == 0) return;
    apply_decorrelating_permutation(mmt, d);
    allreduce_across(mmt, NULL, !d);
    apply_identity(mmt, !d);
    allreduce_across(mmt, NULL, d);
}
/* }}} */

/* {{{ matmul_top_apply_T */
void matmul_top_apply_T(matmul_top_data_ptr mmt, int d)
{
    if (d == 0) return;
    apply_identity(mmt, d);
    allreduce_across(mmt, NULL, !d);
    apply_decorrelating_permutation(mmt, !d);
    allreduce_across(mmt, NULL, d);
}
/* }}} */

/* {{{ matmul_top_twist_vector */
void matmul_top_twist_vector(matmul_top_data_ptr mmt, int d)
{
    serialize(mmt->pi->m);
    if (d == 1)
        matmul_top_apply_S(mmt, d);
    else
        matmul_top_unapply_S_unapply_P(mmt, d);
    // I know this looks odd, but I have witnessed inconsistencies when
    // matmul_top_mul is called directly after twist_vector, for the case
    // of shared vectors. twist_vector() is rare enough to allow an
    // extra serializing call.
    serialize(mmt->pi->m);
}
/* }}} */

/* {{{ matmul_top_untwist_vector */
void matmul_top_untwist_vector(matmul_top_data_ptr mmt, int d)
{
    serialize(mmt->pi->m);
    if (d == 1)
        matmul_top_unapply_S(mmt, d);
    else
        matmul_top_apply_P_apply_S(mmt, d);
    // I know this looks odd, but I have witnessed inconsistencies when
    // matmul_top_mul is called directly after twist_vector, for the case
    // of on-shared vectors. twist_vector() is rare enough to allow an
    // extra serializing call.
    serialize(mmt->pi->m);
}
/* }}} */

/* {{{ Application of permutations to indices */
/* FIXME: this is stupid and does not belong here. We should instead do
 * this as a parallelizing_info collective operation */
static void share_u32_table(pi_comm_ptr wr, uint32_t * g, unsigned int n)
{
    uint32_t ** allg = shared_malloc(wr, wr->ncores * sizeof(uint32_t *));
    allg[wr->trank] = g;
    serialize_threads(wr);
    if (wr->trank == 0) {
        for(unsigned int t = 1 ; t < wr->ncores ; t++) {
            for(unsigned int k = 0 ; k < n ; k++) {
                g[k] += allg[t][k];
            }
        }
        /* We expect one given table cell to contain a non-zero value
         * only at one location. Thus it doesn't really matter whether we
         * do +, XOR, OR */
        int err = MPI_Allreduce(MPI_IN_PLACE, g, n*sizeof(uint32_t), MPI_BYTE, MPI_BXOR, wr->pals);
        ASSERT_ALWAYS(!err);
    }
    serialize_threads(wr);
    if (wr->trank)
        memcpy(g, allg[0], n * sizeof(uint32_t));
    serialize_threads(wr);
    shared_free(wr, allg);
}

#if 0
void indices_apply_S(matmul_top_data_ptr mmt, uint32_t * xs, unsigned int n, int d)
{
    mmt_comm_ptr mrow = mmt->wr[!d];
    mmt_comm_ptr mcol = mmt->wr[d];
    unsigned int ncols = mmt->n[d];

    /* rebuild the complete permutation S. Slightly expensive, but if we
     * try to do otherwise, we'll go nuts */
    uint32_t * perm = malloc(ncols * sizeof(uint32_t));
    memset(perm, 0, ncols * sizeof(uint32_t));
    for(uint32_t i = 0 ; i < mmt->mm->ntwists ; i++) {
        uint32_t u = mmt->mm->twist[i][0];
        uint32_t v = mmt->mm->twist[i][1];
        ASSERT_ALWAYS(v >= mrow->i0);
        ASSERT_ALWAYS(v < mrow->i1);
        ASSERT_ALWAYS(u >= mcol->i0);
        ASSERT_ALWAYS(u < mcol->i1);
        perm[u] = v;
    }

    share_u32_table(mmt->pi->m, perm, ncols);

    for(unsigned int j = 0 ; j < n ; j++) {
        xs[j] = perm[xs[j]];
    }
    free(perm);
}
#endif

void indices_apply_S(matmul_top_data_ptr mmt, uint32_t * xs, unsigned int n, int d MAYBE_UNUSED)
{
    mmt_comm_ptr mrow = mmt->wr[0];
    mmt_comm_ptr mcol = mmt->wr[1];

    uint32_t * r = NULL;
    r = malloc((mrow->i1 - mrow->i0) * sizeof(uint32_t));
    memset(r, 0, (mrow->i1 - mrow->i0) * sizeof(uint32_t));

    for(uint32_t i = 0 ; i < mmt->mm->ntwists ; i++) {
        uint32_t u = mmt->mm->twist[i][1];
        uint32_t v = mmt->mm->twist[i][0];
        ASSERT_ALWAYS(v >= mrow->i0);
        ASSERT_ALWAYS(v < mrow->i1);
        ASSERT_ALWAYS(u >= mcol->i0);
        ASSERT_ALWAYS(u < mcol->i1);
        r[v - mrow->i0] = u;
    }
    for(unsigned int j = 0 ; j < n ; j++) {
        uint32_t i = xs[j];
        if (i >= mrow->i0 && i < mrow->i1)
            xs[j] = r[i - mrow->i0];
        else
            xs[j] = 0;
    }
    free(r);
    share_u32_table(mmt->pi->m, xs, n);
}

void indices_unapply_S(matmul_top_data_ptr mmt, uint32_t * xs, unsigned int n, int d MAYBE_UNUSED)
{
    mmt_comm_ptr mrow = mmt->wr[0];
    mmt_comm_ptr mcol = mmt->wr[1];

    uint32_t * c = NULL;
    c = malloc((mcol->i1 - mcol->i0) * sizeof(uint32_t));
    memset(c, 0, (mcol->i1 - mcol->i0) * sizeof(uint32_t));

    for(uint32_t i = 0 ; i < mmt->mm->ntwists ; i++) {
        uint32_t u = mmt->mm->twist[i][1];
        uint32_t v = mmt->mm->twist[i][0];
        ASSERT_ALWAYS(v >= mrow->i0);
        ASSERT_ALWAYS(v < mrow->i1);
        ASSERT_ALWAYS(u >= mcol->i0);
        ASSERT_ALWAYS(u < mcol->i1);
        c[u - mcol->i0] = v;
    }
    for(unsigned int j = 0 ; j < n ; j++) {
        uint32_t i = xs[j];
        if (i >= mcol->i0 && i < mcol->i1)
            xs[j] = c[i - mcol->i0];
        else
            xs[j] = 0;
    }
    free(c);
    share_u32_table(mmt->pi->m, xs, n);
}

void indices_apply_P(matmul_top_data_ptr mmt, uint32_t * xs, unsigned int n, int d)
{
    pi_comm_ptr picol = mmt->pi->wr[d];
    pi_comm_ptr pirow = mmt->pi->wr[!d];
    unsigned int ncols = mmt->n[d];
    unsigned int nv = pirow->totalsize;
    unsigned int nh = picol->totalsize;

    for(unsigned int j = 0 ; j < n ; j++) {
        uint32_t v = xs[j];
        /* Now apply the permutation P to v */
        uint32_t z = ncols / (nh*nv);
        uint32_t q = v / z;
        uint32_t r = v % z;
        uint32_t q1 = q / nv;
        uint32_t q0 = q % nv;
        xs[j] = (q0 * nh + q1) * z + r;
    }
}
/* }}} */

/**********************************************************************/
#if 0
/* no longer used -- was only used by prep.
 * It's not buggy, but making this work in a context where we have
 * multiple threads is tricky.
 */
void matmul_top_fill_random_source(matmul_top_data_ptr mmt, int d)
{
    matmul_top_fill_random_source_generic(mmt, mmt->vr->stride, NULL, d);
}
#endif

/* Takes data in mmt->wr[d]->v, and compute the corresponding partial result in
 * mmt->wr[!d]->v.
 */
void matmul_top_mul_cpu(matmul_top_data_ptr mmt, int d)
{
#ifndef NDEBUG
    unsigned int di_in = mmt->wr[d]->i1 - mmt->wr[d]->i0;
    unsigned int di_out = mmt->wr[!d]->i1 - mmt->wr[!d]->i0;
    ASSERT(mmt->mm->dim[1] == (d ? di_in : di_out));
    ASSERT(mmt->mm->dim[0] == (d ? di_out : di_in));
#endif
    mmt_comm_ptr mrow = mmt->wr[!d];
    mmt_comm_ptr mcol = mmt->wr[d];

    ASSERT_ALWAYS((mrow->v->flags & THREAD_SHARED_VECTOR) == 0);

    pi_log_op(mmt->pi->m, "[%s] enter matmul_mul", __func__);

    /* Note that matmul_init copies the calling abase argument to the
     * lower-level mm structure. It can quite probably be qualified as a
     * flaw.
     */
    matmul_mul(mmt->mm, mrow->v->v, mcol->v->v, d);
}

/* This takes partial results in the areas mmt->wr[!d]->v, and puts the
 * collected and re-broadcasted results in the areas mmt->wd[d]->v
 *
 * Note that for the shuffled product, this is not equivalent to a trivial
 * operation.
 */
void matmul_top_mul_comm(matmul_top_data_ptr mmt, int d)
{
    pi_comm_ptr picol = mmt->pi->wr[d];
    mmt_comm_ptr mcol = mmt->wr[d];

    pi_log_op(mmt->pi->m, "[%s] enter reduce_across", __func__);
    reduce_across(mmt, NULL, NULL, !d);
    pi_log_op(mmt->pi->m, "[%s] enter broadcast_down", __func__);
    broadcast_down(mmt, NULL, d);

    /* If we have shared input data for the column threads, then we'd
     * better make sure it has arrived completely, because while all
     * threads will need the data, only one is actually importing it.
     */
    if (mcol->v->flags & THREAD_SHARED_VECTOR) {
        pi_log_op(picol, "[%s] serialize threads", __func__);
        serialize_threads(picol);
    }
}


/* Vector I/O is done by only one job, one thread. It incurs a
 * significant amount of memory allocation, but this is done relatively
 * rarely.
 *
 * Vectors, as handled within the core routines, are permuted. He have
 * coordinates (v_{\sigma(0)}...v_{\sigma(n-1)}), where \sigma is the
 * permutation given by the *_perm files. For block Wiedemann, we assume
 * that we have conjugated permutations on the two sides. This means that
 * no data manipulation is required within the critical loops (this would
 * imply a fuzzy communication pattern, boiling down to essentially using
 * allreduce instead of reduce).
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

void mmt_vec_set_random_through_file(matmul_top_data_ptr mmt, mmt_vec_ptr v, const char * name, int d, unsigned int iter, unsigned int itemsondisk, gmp_randstate_t rstate)
{
    /* FIXME: this generates the complete vector on rank 0, saves it, and
     * loads it again. But I'm a bit puzzled by the choice of saving a
     * number of items which is n0[d]. Seems to me that this is in fact
     * incorrect, we want n0[!d] here.
     */
    mpfq_vbase_ptr A = mmt->abase;
    pi_datatype_ptr A_pi = mmt->pitype;
    parallelizing_info_ptr pi = mmt->pi;

    if (pi->m->trank == 0) {
        char * filename;
        int rc;
        rc = asprintf(&filename, "%s.%d", name, iter);
        ASSERT_ALWAYS(rc >= 0);

        void * y;
        A->vec_init(A, &y, mmt->n[d]);
        A->vec_set_zero(A, y, mmt->n[d]);
        /* Again, important. Generate zero coordinates for padding !
         * This provides reproducibility of random choices.
         */
        if (pi->m->jrank == 0)
            A->vec_random(A, y, mmt->n0[d], rstate);
        int err = MPI_Bcast(y, mmt->n[d], A_pi->datatype, 0, pi->m->pals);
        ASSERT_ALWAYS(!err);
        FILE * f = fopen(filename, "wb");
        ASSERT_ALWAYS(f);
        rc = fwrite(y, A->vec_elt_stride(A,1), itemsondisk, f);
        ASSERT_ALWAYS(rc == (int) itemsondisk);
        fclose(f);
        A->vec_clear(A, &y, mmt->n[d]);
        free(filename);
    }
    mmt_vec_load(mmt, v, name, d, iter, itemsondisk);
}

void mmt_vec_set_random_inconsistent(matmul_top_data_ptr mmt, mmt_vec_ptr v, int d, gmp_randstate_t rstate)
{
    if (v == NULL) v = mmt->wr[d]->v;
    mmt_comm_ptr wr = mmt->wr[d];
    mmt->abase->vec_random(mmt->abase, v->v, wr->i1 - wr->i0, rstate);
}



/**********************************************************************/
static void mmt_finish_init(matmul_top_data_ptr mmt, int const *, param_list_ptr pl, int optimized_direction);
static void matmul_top_read_submatrix(matmul_top_data_ptr mmt, param_list_ptr pl, int optimized_direction);

static void mmt_fill_fields_from_balancing(matmul_top_data_ptr mmt, param_list_ptr pl)
{
    mmt->n[0] = mmt->bal->trows;
    mmt->n[1] = mmt->bal->tcols;
    mmt->n0[0] = mmt->bal->h->nrows;
    mmt->n0[1] = mmt->bal->h->ncols;
   
    // nslices[0] is the number of slices in a column -- this is the number of
    // ``horizontal slices'', in balance.c parlance.
    //
    // nslices[0] is also the number of big pieces into which a left
    // vector is split (eventually, it is split in nh * nv pieces of
    // course, but the topmost split level is nh).
    unsigned int nslices[2];
    nslices[0] = mmt->bal->h->nh;
    nslices[1] = mmt->bal->h->nv;

    int ok = 1;

    for(int d = 0 ; d < 2 ; d++) {
        /* the column communicator must be split in exactly as many
         * pieces as we use to split the left vector */
        ok = ok && mmt->pi->wr[d]->totalsize == nslices[!d];
    }

    if (!ok) {
        fprintf(stderr, "Configured split %ux%ux%ux%u does not match "
                "split %ux%u found in balancing file.\n",
                mmt->pi->wr[1]->njobs,
                mmt->pi->wr[0]->njobs,
                mmt->pi->wr[1]->ncores,
                mmt->pi->wr[0]->ncores,
                nslices[0],
                nslices[1]);
        serialize(mmt->pi->m);
        exit(1);
    }

    // mmt->ncoeffs_total = 0;
    unsigned int ix[2];
    // ix[1] is something between 0 and nslices[0]-1, i.e. between 0 and nh-1
    // ix[0] is something between 0 and nslices[1]-1, i.e. between 0 and nv-1

    // fences[0] has 1+nslices[0] stop, which is nh+1
    // fences[1] has 1+nslices[1] stop, which is nv+1
    //
    // ix[1] and ix[1]+1 are thus indices within fences[0]
    for(int d = 0 ; d < 2 ; d++)  {
        unsigned int * f = malloc((1 + nslices[d])*sizeof(unsigned int));
        // we use UINT_MAX markers for sanity checking.
        memset(f, -1, (1+nslices[d]) * sizeof(unsigned int));
        f[nslices[d]] = mmt->n[d];
        uint32_t quo = mmt->n[d] / nslices[d];
        unsigned int rem = mmt->n[d] % nslices[d];
        for(unsigned int k = 0 ; k < nslices[d] ; k++) {
            f[k] = k * quo + MIN(rem, k);
        }
        mmt->fences[d] = f;
    }

    for(int d = 0 ; d < 2 ; d++)  {
        ix[!d] = mmt->pi->wr[!d]->jrank * mmt->pi->wr[!d]->ncores;
        ix[!d]+= mmt->pi->wr[!d]->trank;
        mmt->wr[d]->i0 = mmt->fences[d][ix[!d]];
        mmt->wr[d]->i1 = mmt->fences[d][ix[!d] + 1];
        ASSERT_ALWAYS(mmt->wr[d]->i0 < mmt->wr[d]->i1);
        ASSERT_ALWAYS(mmt->wr[d]->i1 <= mmt->n[d]);
        ASSERT_ALWAYS((mmt->wr[d]->i1 - mmt->wr[d]->i0) % mmt->pi->wr[d]->totalsize == 0);

        // the mm layer has an ncoeffs field as well, which is used,
        // while to the contrary the one in the mmt structure is totally
        // unused (and is now gone)
        // mmt->ncoeffs = ncoeffs;
        // mmt->ncoeffs_total += ncoeffs;
    }

    char * base;
    {
        const char * mname;
        // get basename of matrix file. Remove potential suffix, .txt or .bin
        mname = param_list_lookup_string(pl, "matrix");
        ASSERT_ALWAYS(mname);
        char * last_slash = strrchr(mname, '/');
        if (last_slash) {
            mname = last_slash + 1;
        }
        const char * suffixes[] = { ".txt", ".bin" };
        unsigned int l = strlen(mname);
        for(unsigned int k = 0 ; k < sizeof(suffixes) / sizeof(suffixes[0]) ; k++) {
            unsigned int ls = strlen(suffixes[k]);
            if (ls > l) continue;
            if (strcmp(mname+l-ls,suffixes[k]) == 0)
                l -= ls;
        }
        base = strndup(mname, l);
    }

    /* The name choice is only indicative. Depending on the value of
     * the nullspace argument, the matrices in h?.v? may in fact be
     * transposed because this is what the cache building code prefers to
     * receive.
     */
    int rc = asprintf(&mmt->locfile, "%s.%dx%d.%08" PRIx32 ".h%d.v%d", base, nslices[0], nslices[1], mmt->bal->h->checksum, ix[1], ix[0]);
    free(base);
    ASSERT_ALWAYS(rc >= 0);
}

void matmul_top_init(matmul_top_data_ptr mmt,
        mpfq_vbase_ptr abase,
        pi_datatype_ptr pitype,
        /* matmul_ptr mm, */
        parallelizing_info_ptr pi,
        int const * flags,
        param_list_ptr pl,
        int optimized_direction)
{
    memset(mmt, 0, sizeof(*mmt));

    mmt->abase = abase;
    mmt->pitype = pitype;
    mmt->pi = pi;
    mmt->mm = NULL;

    serialize_threads(mmt->pi->m);

    // n[]
    // ncoeffs_total, ncoeffs,
    // locfile,
    // fences[]
    // wr[]->*
    // are filled in later within read_info_file
    //

    const char * tmp, * rtmp;
    tmp = param_list_lookup_string(pl, "balancing");
    rtmp = param_list_lookup_string(pl, "random_matrix");
    if (tmp && rtmp) {
        fprintf(stderr, "balancing= and random_matrix= are incompatible\n");
        exit(1);
    }
    if (!tmp && !rtmp) {
        fprintf(stderr, "Missing parameter balancing=\n");
        exit(1);
    }
    if (pi->m->jrank == 0 && pi->m->trank == 0) {
        if (tmp) {
            balancing_read_header(mmt->bal, tmp);
        } else {
            random_matrix_fill_fake_balancing_header(mmt->bal, pi, rtmp);
        }
    }

    pi_bcast(mmt->bal, sizeof(mmt->bal), BWC_PI_BYTE, 0, 0, mmt->pi->m);

    // after that, we need to check for every node if the cache file can
    // be found. If none is found, rebuild the cache files
    mmt_finish_init(mmt, flags, pl, optimized_direction);
}

/* Here, we get a copy of the rowperm and colperm.
 *
 * For each (job,thread), two pairs of intervals are defined.
 *
 * for row indices: [i0[0], i1[0][ = [mrow->i0, mrow->i1[ and [Xi0[0], Xi1[0][
 * for col indices: [i0[1], i1[1][ = [mcol->i0, mcol->i1[ and [Xi0[1], Xi1[1][
 *
 * For each job (common among threads), two other pairs are defined
 * (those are only of interest to the internals of this function):
 *
 * for row indices: [i0_job[0], i1_job[0][ and [Xi0_job[0], Xi1_job[0][
 * for col indices: [i0_job[1], i1_job[1][ and [Xi0_job[1], Xi1_job[1][
 *
 *
 * The parts we get are:
 *      [i, rowperm[i]] for    i in [mrow->i0, mrow->i1[
 *                         and rowperm[i] in [X0, X1[
 *      [j, colperm[j]] for    j in [mcol->i0, mcol->i1[
 *                         and colperm[i] in [Y0, Y1[
 * where X0 is what would be mcol->i0 if the matrix as many columns as
 * rows, and ditto for X1,Y0,Y1.
 */
static void mmt_get_local_permutation_data(matmul_top_data_ptr mmt, param_list_ptr pl)
{
    mmt_comm_ptr mwr[2];
    pi_comm_ptr piwr[2];
    size_t eblock[2];
    int pos[2];
    unsigned int i0[2], i1[2];
    unsigned int Xi0[2], Xi1[2];
    // unsigned int i0_job[2], i1_job[2];
    // unsigned int Xi0_job[2], Xi1_job[2];

    for(int d = 0 ; d < 2 ; d++)  {
        mwr[d] = mmt->wr[d];
        piwr[d] = mmt->pi->wr[d];
        eblock[d] = (mwr[d]->i1 - mwr[d]->i0) /  piwr[d]->totalsize;
        pos[d] = piwr[d]->jrank * piwr[d]->ncores + piwr[d]->trank;
    }
    /* Now that we've understood this quite well, decide whiche are our
     * intervals for both the left and right vectors. This is used for
     * the Xi0 and Xi1 things.
     */
    for(int d = 0 ; d < 2 ; d++)  {
        size_t e = eblock[d];
        // int pos_job = piwr[!d]->jrank * piwr[!d]->ncores;
        // int Xpos_job = piwr[d]->jrank * piwr[d]->ncores;
        i0[d] = e * pos[!d] * piwr[d]->totalsize;
        i1[d] = e * (pos[!d]+1) * piwr[d]->totalsize;
        Xi0[d] = e * pos[d] * piwr[!d]->totalsize;
        Xi1[d] = e * (pos[d]+1) * piwr[!d]->totalsize;
        // i0_job[d] = e *  pos_job    * piwr[d]->totalsize;
        // i1_job[d] = e * (pos_job+piwr[!d]->ncores) * piwr[d]->totalsize;
        // Xi0_job[d] = e *  Xpos_job * piwr[!d]->totalsize;
        // Xi1_job[d] = e * (Xpos_job+piwr[d]->ncores) * piwr[!d]->totalsize;
        ASSERT_ALWAYS(mwr[d]->i0 == i0[d]);
        ASSERT_ALWAYS(mwr[d]->i1 == i1[d]);
    }
    const char * text[2] = { "left", "right", };
    /*
    for(int d = 0 ; d < 2 ; d++)  {
        SEVERAL_THREADS_PLAY_MPI_BEGIN(mmt->pi->m) {
            fprintf(stderr, "[%s] J%uT%u: on %s vectors: read [%u,%u[ virtual-write [%u,%u[\n",
                    mmt->pi->nodenumber_s,
                    mmt->pi->m->jrank, mmt->pi->m->trank,
                    text[d],
                    i0[d], i1[d],
                    Xi0[d], Xi1[d]);
            if (mmt->pi->m->trank == 0) {
                fprintf(stderr, "[%s] J%u*: on %s vectors: read [%u,%u[ virtual-write [%u,%u[\n",
                        mmt->pi->nodenumber_s,
                        mmt->pi->m->jrank,
                        text[d],
                        i0_job[d], i1_job[d],
                        Xi0_job[d], Xi1_job[d]);
            }
        }
        SEVERAL_THREADS_PLAY_MPI_END();
    }
    */

    /* we're going to read the balancing file completely, on the leader
     * node. Make sure there's nothing there yet.
     */
    ASSERT_ALWAYS(mmt->bal->rowperm == NULL);
    ASSERT_ALWAYS(mmt->bal->colperm == NULL);
    unsigned int rowperm_items=0;
    unsigned int colperm_items=0;

    /* Define a complete structure for the balancing which is shared
     * among threads, so that we'll be able to access it from all threads
     * simultaneously. We will put things in bal_tmp->rowperm and
     * bal_tmp->colperm, but beyond that, the header part will be wrong
     * at non-root nodes.
     */
    balancing_ptr bal_tmp = shared_malloc_set_zero(mmt->pi->m, sizeof(balancing));

    if (mmt->pi->m->jrank == 0 && mmt->pi->m->trank == 0) {
        const char * tmp = param_list_lookup_string(pl, "balancing");
        if (tmp)
            balancing_read(bal_tmp, tmp);
        /* It's fine if we have nothing. This just means that we'll have
         * no balancing to deal with. */
        rowperm_items = bal_tmp->rowperm != NULL ? bal_tmp->trows : 0;
        colperm_items = bal_tmp->colperm != NULL ? bal_tmp->tcols : 0;
    }
    pi_bcast(&rowperm_items, 1, BWC_PI_UNSIGNED, 0, 0, mmt->pi->m);
    pi_bcast(&colperm_items, 1, BWC_PI_UNSIGNED, 0, 0, mmt->pi->m);

    if (mmt->pi->m->trank == 0) {
        if (rowperm_items) {
            if (mmt->pi->m->jrank != 0)
                bal_tmp->rowperm = malloc(mmt->bal->trows * sizeof(uint32_t));
            MPI_Bcast(bal_tmp->rowperm, mmt->bal->trows * sizeof(uint32_t), MPI_BYTE, 0, mmt->pi->m->pals);
        }
        if (colperm_items) {
            if (mmt->pi->m->jrank != 0)
                bal_tmp->colperm = malloc(mmt->bal->tcols * sizeof(uint32_t));
            MPI_Bcast(bal_tmp->colperm, mmt->bal->tcols * sizeof(uint32_t), MPI_BYTE, 0, mmt->pi->m->pals);
        }

#if 0
        /* now restrict the number of rowperm_items to the number we care
         * about locally */
        if (i0_block[0] != 0) {
            memcpy(bal_tmp->rowperm, bal_tmp->rowperm + i0_block[0], (i1_block[0] - i0_block[0]) * sizeof(uint32_t));
        }
        rowperm_items = i1_block[0] - i0_block[0];
        bal_tmp->rowperm = realloc(bal_tmp->rowperm, (i1_block[0] - i0_block[0]) * sizeof(uint32_t));
#endif
    }
    serialize_threads(mmt->pi->m);      /* important ! */

    uint32_t * balperm[2] = { bal_tmp->rowperm, bal_tmp->colperm };

    for(int d = 0 ; d < 2 ; d++)  {
        if (!balperm[d]) continue;

        mmt->perm[d] = permutation_data_alloc();

        /* now create the really local permutation */
        for(unsigned int i = i0[d] ; i < i1[d] ; i++) {
            uint32_t ij[2] = { i, balperm[d][i] };
            if (ij[1] >= Xi0[d] && ij[1] < Xi1[d]) {
                permutation_data_push(mmt->perm[d], ij);
            }
        }
        printf("[%s] J%uT%u does %u/%u permutation pairs for %s vectors\n",
                mmt->pi->nodenumber_s,
                mmt->pi->m->jrank, mmt->pi->m->trank,
                mmt->perm[d]->n, d ? mmt->bal->tcols : mmt->bal->trows,
                text[d]);
    }

    serialize_threads(mmt->pi->m);      /* important ! */

    if (mmt->pi->m->trank == 0) {
        if (bal_tmp->colperm) free(bal_tmp->colperm);
        if (bal_tmp->rowperm) free(bal_tmp->rowperm);
    }

    shared_free(mmt->pi->m, bal_tmp);
}




/* Some work has to be done in order to fill the remaining fields in the
 * matmul_top structure.
 */
static void mmt_finish_init(matmul_top_data_ptr mmt, int const * flags, param_list_ptr pl, int optimized_direction)
{
    mmt_fill_fields_from_balancing(mmt, pl);

    mmt_get_local_permutation_data(mmt, pl);

    // TODO: reuse the ../bw-matmul/matmul/matrix_base.cpp things for
    // displaying communication info.

    if (!mmt->pi->interleaved) {
        matmul_top_read_submatrix(mmt, pl, optimized_direction);
    } else {
        /* Interleaved threads will share their matrix data. The first
         * thread to arrive will do the initialization, and the second
         * one will just grab the pointer. The trick is to be able to
         * pick the pointer in the right location ! We add a generic
         * pointer dictionary feature in the parallelizing_info
         * interface for this purpose.
         */

#define MMT_MM_MAGIC_KEY        0xaa000000UL

        if (mmt->pi->interleaved->idx == 0) {
            matmul_top_read_submatrix(mmt, pl, optimized_direction);
            pi_store_generic(mmt->pi, MMT_MM_MAGIC_KEY, mmt->pi->m->trank, mmt->mm);
        } else {
            mmt->mm = pi_load_generic(mmt->pi, MMT_MM_MAGIC_KEY, mmt->pi->m->trank);
        }
    }

    /* NOTE: flags default to zero */
    mmt_vec_init(mmt, NULL, NULL, NULL, 0, flags ? flags[0] : 0);
    mmt_vec_init(mmt, NULL, NULL, NULL, 1, flags ? flags[1] : 0);
}

static int export_cache_list_if_requested(matmul_top_data_ptr mmt, param_list_ptr pl)
{
    const char * cachelist = param_list_lookup_string(pl, "export_cachelist");
    if (!cachelist) return 0;

    char * myline = NULL;
    int rc;
    rc = asprintf(&myline, "%s %s", mmt->pi->nodename, mmt->mm->cachefile_name);
    ASSERT_ALWAYS(rc >= 0);
    ASSERT_ALWAYS(myline != NULL);
    char ** tlines = NULL;
    tlines = shared_malloc_set_zero(mmt->pi->m, mmt->pi->m->ncores * sizeof(const char *));
    tlines[mmt->pi->m->trank] = myline;
    serialize_threads(mmt->pi->m);

    /* Also, just out of curiosity, try to see what we have currently */
    struct stat st[1];
    int * has_cache = shared_malloc_set_zero(mmt->pi->m, mmt->pi->m->totalsize * sizeof(int));
    rc = stat(mmt->mm->cachefile_name, st);
    unsigned int mynode = mmt->pi->m->ncores * mmt->pi->m->jrank;
    has_cache[mynode + mmt->pi->m->trank] = rc == 0;
    serialize_threads(mmt->pi->m);

    int len = 0;
    if (mmt->pi->m->trank == 0) {
    MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
                has_cache, mmt->pi->m->ncores, MPI_INT, mmt->pi->m->pals);

        for(unsigned int j = 0 ; j < mmt->pi->m->ncores ; j++) {
            int s = strlen(tlines[j]);
            if (s >= len) len = s;
        }
        MPI_Allreduce(MPI_IN_PLACE, &len, 1, MPI_INT, MPI_MAX, mmt->pi->m->pals);
        char * info = malloc(mmt->pi->m->totalsize * (len + 1));
        char * mybuf = malloc(mmt->pi->m->ncores * (len+1));
        memset(mybuf, 0, mmt->pi->m->ncores * (len+1));
        for(unsigned int j = 0 ; j < mmt->pi->m->ncores ; j++) {
            memcpy(mybuf + j * (len+1), tlines[j], strlen(tlines[j])+1);
        }
        MPI_Allgather(mybuf, mmt->pi->m->ncores * (len+1), MPI_BYTE,
                info, mmt->pi->m->ncores * (len+1), MPI_BYTE,
                mmt->pi->m->pals);
        if (mmt->pi->m->jrank == 0) {
            FILE * f = fopen(cachelist, "wb");
            DIE_ERRNO_DIAG(f == NULL, "fopen", cachelist);
            for(unsigned int j = 0 ; j < mmt->pi->m->njobs ; j++) {
                unsigned int j0 = j * mmt->pi->m->ncores;
                fprintf(f, "get-cache ");
                for(unsigned int k = 0 ; k < mmt->pi->m->ncores ; k++) {
                    char * t = info + (j0 + k) * (len + 1);
                    char * q = strchr(t, ' ');
                    ASSERT_ALWAYS(q);
                    if (!k) {
                        *q='\0';
                        fprintf(f, "%s", t);
                        *q=' ';
                    }
                    fprintf(f, "%s", q);
                }
                fprintf(f, "\n");
            }
            for(unsigned int j = 0 ; j < mmt->pi->m->njobs ; j++) {
                unsigned int j0 = j * mmt->pi->m->ncores;
                fprintf(f, "has-cache ");
                for(unsigned int k = 0 ; k < mmt->pi->m->ncores ; k++) {
                    char * t = info + (j0 + k) * (len + 1);
                    char * q = strchr(t, ' ');
                    ASSERT_ALWAYS(q);
                    if (!k) {
                        *q='\0';
                        fprintf(f, "%s", t);
                        *q=' ';
                    }
                    if (!has_cache[mmt->pi->m->ncores * j + k]) continue;
                    fprintf(f, "%s", q);
                }
                fprintf(f, "\n");
            }
            fclose(f);
        }
        free(info);
        free(mybuf);
    }
    serialize_threads(mmt->pi->m);
    shared_free(mmt->pi->m, tlines);
    shared_free(mmt->pi->m, has_cache);
    serialize_threads(mmt->pi->m);
    free(myline);
    serialize(mmt->pi->m);

    return 1;
}


static void matmul_top_read_submatrix(matmul_top_data_ptr mmt, param_list_ptr pl, int optimized_direction)
{
    int rebuild = 0;
    param_list_parse_int(pl, "rebuild_cache", &rebuild);

    mmt->mm = matmul_init(mmt->abase,
            mmt->wr[0]->i1 - mmt->wr[0]->i0,
            mmt->wr[1]->i1 - mmt->wr[1]->i0,
            mmt->locfile, NULL  /* means: choose mm_impl from pl */,
            pl, optimized_direction); 

    mmt->mm->nslices[1] = mmt->bal->h->nh;
    mmt->mm->nslices[0] = mmt->bal->h->nv;

    // *IF* we need to do a collective read of the matrix, we need to
    // provide the pointer *now*.
    unsigned int sqread = 0;
    param_list_parse_uint(pl, "sequential_cache_read", &sqread);

    int cache_loaded = 0;

    if (export_cache_list_if_requested(mmt, pl)) {
        /* This has the effect of terminating the program if we are being
         * called from dispatch */
        return;
    }

    if (!rebuild) {
        if (mmt->pi->m->jrank == 0 && mmt->pi->m->trank == 0) {
            printf("Now trying to load matrix cache files\n");
        }
        if (sqread) {
            for(unsigned int j = 0 ; j < mmt->pi->m->ncores ; j++) {
                serialize_threads(mmt->pi->m);
                if (j == mmt->pi->m->trank) {
                    cache_loaded = matmul_reload_cache(mmt->mm);
                }
            }
        } else {
            cache_loaded = matmul_reload_cache(mmt->mm);
        }
    }

    if (!pi_data_eq(&cache_loaded, 1, BWC_PI_INT, mmt->pi->m)) {
        if (mmt->pi->m->trank == 0) {
        if (mmt->pi->m->jrank == 0) {
            fprintf(stderr, "Fatal error: cache files not present at expected locations\n");
        }
        }
        SEVERAL_THREADS_PLAY_MPI_BEGIN(mmt->pi->m) {
            fprintf(stderr, "[%s] J%uT%u: cache %s: %s\n",
                    mmt->pi->nodenumber_s,
                    mmt->pi->m->jrank,
                    mmt->pi->m->trank,
                    mmt->mm->cachefile_name,
                    cache_loaded ? "ok" : "not ok");
        }
        SEVERAL_THREADS_PLAY_MPI_END();
        serialize(mmt->pi->m);
        abort();
    }

    unsigned int sqb = 0;
    param_list_parse_uint(pl, "sequential_cache_build", &sqb);

    matrix_u32 m;
    memset(m, 0, sizeof(matrix_u32));
    /* see remark in raw_matrix_u32.h about data ownership for type
     * matrix_u32 */

    if (!cache_loaded) {
        m->mfile = param_list_lookup_string(pl, "matrix");
        m->bfile = param_list_lookup_string(pl, "balancing");
        // the mm layer is informed of the higher priority computations
        // that will take place. Depending on the implementation, this
        // may cause the direct or transposed ordering to be preferred.
        // Thus we have to read this back from the mm structure.
        m->transpose = mmt->mm->store_transposed;
        {
            mpz_t p;
            mpz_init(p);
            mmt->abase->field_characteristic(mmt->abase, p);
            m->withcoeffs = mpz_cmp_ui(p, 2) != 0;
            mpz_clear(p);
        }
        if (param_list_lookup_string(pl, "random_matrix")) {
            if (mmt->pi->m->jrank == 0 && mmt->pi->m->trank == 0) {
                printf("Begin creation of fake matrix data in parallel\n");
            }
            random_matrix_get_u32(mmt->pi, pl, m);
            /* Fill in ntwists. It's important for BL, since otherwise
             * with fake matrices the untwisting will in effect multiply
             * by zero...
             */
            mmt_comm_ptr mrow = mmt->wr[0];
            mmt_comm_ptr mcol = mmt->wr[1];
            unsigned int nx;
            unsigned int row_offset;
            unsigned int col_offset;
            nx = intersect_two_intervals(&row_offset, &col_offset, mrow->i0, mrow->i1, mcol->i0, mcol->i1);
            ASSERT_ALWAYS(m->ntwists == 0);
            ASSERT_ALWAYS(m->twist == NULL);
            m->ntwists = nx;
            m->twist = malloc(nx * sizeof(uint32_t[2]));
            for(unsigned int i = 0 ; i < nx ; i++) {
                m->twist[i][0] = mrow->i0 + row_offset + i;
                m->twist[i][1] = mrow->i0 + row_offset + i;
            }
        } else {
            if (mmt->pi->m->jrank == 0 && mmt->pi->m->trank == 0) {
                printf("Matrix dispatching starts\n");
            }
            balancing_get_matrix_u32(mmt->pi, pl, m);
        }

        int ssm = 0;
        param_list_parse_int(pl, "save_submatrices", &ssm);
        if (ssm) {
            fprintf(stderr, "DEBUG: creating %s\n", mmt->locfile);
            FILE * f = fopen(mmt->locfile, "wb");
            fwrite(m->p, sizeof(uint32_t), m->size, f);
            fclose(f);
        }
    }


    if (!sqb) {
        if (!cache_loaded) {
            // everybody does it in parallel
            if (verbose_enabled(CADO_VERBOSE_PRINT_BWC_CACHE_MAJOR_INFO))
                printf("[%s] J%uT%u building cache for %s\n",
                        mmt->pi->nodenumber_s,
                        mmt->pi->m->jrank,
                        mmt->pi->m->trank,
                        mmt->locfile);
            matmul_build_cache(mmt->mm, m);
            matmul_save_cache(mmt->mm);
        }
    } else {
        for(unsigned int j = 0 ; j < mmt->pi->m->ncores + 1 ; j++) {
            serialize_threads(mmt->pi->m);
            if (cache_loaded) continue;
            if (j == mmt->pi->m->trank) {
                if (verbose_enabled(CADO_VERBOSE_PRINT_BWC_CACHE_MAJOR_INFO))
                    printf("[%s] J%uT%u building cache for %s\n",
                            mmt->pi->nodenumber_s,
                            mmt->pi->m->jrank,
                            mmt->pi->m->trank,
                            mmt->locfile);
                matmul_build_cache(mmt->mm, m);
            } else if (j == mmt->pi->m->trank + 1) {
                matmul_save_cache(mmt->mm);
            }
        }
    }
    /* see remark in raw_matrix_u32.h about data ownership for type
     * matrix_u32 */

    if (verbose_enabled(CADO_VERBOSE_PRINT_BWC_CACHE_MAJOR_INFO)) {
        my_pthread_mutex_lock(mmt->pi->m->th->m);
        printf("[%s] J%uT%u uses cache file %s\n",
                mmt->pi->nodenumber_s,
                mmt->pi->m->jrank, mmt->pi->m->trank,
                /* cache for mmt->locfile, */
                mmt->mm->cachefile_name);
        my_pthread_mutex_unlock(mmt->pi->m->th->m);
    }
}


void matmul_top_clear(matmul_top_data_ptr mmt)
{
    serialize_threads(mmt->pi->m);
    serialize(mmt->pi->m);
    serialize_threads(mmt->pi->m);

    mmt_vec_clear(mmt,NULL,0);
    mmt_vec_clear(mmt,NULL,1);
#if RS_CHOICE == RS_CHOICE_MINE || RS_CHOICE == RS_CHOICE_MINE_PARALLEL
    for(int d = 0 ; d < 2 ; d++)  {
        permutation_data_free(mmt->perm[d]);
        free(mmt->wr[d]->rsbuf[0]); mmt->wr[d]->rsbuf[0] = NULL;
        free(mmt->wr[d]->rsbuf[1]); mmt->wr[d]->rsbuf[1] = NULL;
        mmt->wr[d]->rsbuf_size = 0;
        // both are expected to hold storage for:
        // (mmt->n[d] / mmt->pi->m->totalsize * mmt->pi->wr[d]->ncores))
        // elements, corresponding to the largest abase encountered.
    }
#endif  /* RS_CHOICE == RS_CHOICE_MINE || RS_CHOICE == RS_CHOICE_MINE_PARALLEL */
    serialize(mmt->pi->m);
    if (!mmt->pi->interleaved) {
        matmul_clear(mmt->mm);
    } else if (mmt->pi->interleaved->idx == 1) {
        matmul_clear(mmt->mm);
        /* group 0 is the first to leave, thus it doesn't to freeing.
         */
    }
    free(mmt->fences[0]);
    free(mmt->fences[1]);
    free(mmt->locfile);
}

static permutation_data_ptr permutation_data_alloc()
{
    permutation_data_ptr a = malloc(sizeof(permutation_data));
    a->n = 0;
    a->alloc = 0;
    a->x = NULL;
    return a;
}

static void permutation_data_free(permutation_data_ptr a)
{
    if (!a) return;
    a->n = 0;
    a->alloc = 0;
    if (a->x) free(a->x);
    free(a);
}

static void permutation_data_push(permutation_data_ptr a, uint32_t u[2])
{
    if (a->n >= a->alloc) {
        a->alloc = a->n + 16 + a->alloc / 4;
        a->x = realloc(a->x, a->alloc * sizeof(uint32_t[2]));
    }
    a->x[a->n][0] = u[0];
    a->x[a->n][1] = u[1];
    a->n++;
}
