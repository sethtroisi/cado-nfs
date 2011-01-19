#include "cado.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>
#include "bwc_config.h"
#include "matmul.h"
#include "matmul_top.h"
#include "select_mpi.h"
#include "intersections.h"
#include "debug.h"
#include "filenames.h"
#include "balancing_workhorse.h"

#ifndef CONJUGATED_PERMUTATIONS
#error "Do you really, really want to use arbitrary left and right sigmas ?"
#endif

///////////////////////////////////////////////////////////////////
/* Start with stuff that does not depend on abase at all -- this
 * provides a half-baked interface */

void vec_init_generic(pi_wiring_ptr picol, size_t stride, mmt_generic_vec_ptr v, int flags, unsigned int n)
{
    v->flags = flags;

    // all_v pointers are stored twice, in a rotating buffer manner, so
    // that it's possible for one thread to address the n other threads
    // without having to compute a modulo
    size_t allocsize = 2 * picol->ncores * sizeof(void *);
    v->all_v = alignable_malloc(allocsize);

    /* Because we sometimes provide extra readahead space, we need to
     * zero out the allocated area. Strictly speaking, it would suffice
     * to zero out the readahead zone only, but it's easier this way.
     */
    if (flags & THREAD_SHARED_VECTOR) {
        void * r;
        if (picol->trank == 0) {
            r = abase_generic_init(stride, n);
            abase_generic_zero(stride, r, n);
        }
        thread_broadcast(picol, &r, 0);
        v->v = r;
        for(unsigned int t = 0 ; t < picol->ncores ; t++) {
            v->all_v[t] = v->v;
        }
    } else {
        v->v = abase_generic_init(stride, n);
        abase_generic_zero(stride, v->v, n);
        v->all_v[picol->trank] = v->v;
        for(unsigned int t = 0 ; t < picol->ncores ; t++) {
            serialize_threads(picol);
            void * r;
            r = v->v;
            thread_broadcast(picol, &r, t);
            v->all_v[t] = r;
        }
    }
    for(unsigned int t = 0 ; t < picol->ncores ; t++) {
        v->all_v[picol->ncores + t] = v->all_v[t];
    }

}

void vec_clear_generic(pi_wiring_ptr picol, size_t stride MAYBE_UNUSED, mmt_generic_vec_ptr v, unsigned int n MAYBE_UNUSED)
{
    if (v->flags & THREAD_SHARED_VECTOR) {
        if (picol->trank == 0)
            abase_generic_clear(stride, v->v, n);
    } else {
        abase_generic_clear(stride, v->v, n);
    }
    size_t allocsize MAYBE_UNUSED;
    allocsize = 2 * picol->ncores * sizeof(void *);
    alignable_free(v->all_v, allocsize);
}

void matmul_top_vec_init_generic(matmul_top_data_ptr mmt, size_t stride, mmt_generic_vec_ptr v, int d, int flags)
{
    if (v == NULL) v = (mmt_generic_vec_ptr) mmt->wr[d]->v;
    unsigned int n1 = mmt->wr[d]->i1 - mmt->wr[d]->i0;
    matmul_aux(mmt->mm, MATMUL_AUX_GET_READAHEAD, &n1);
    vec_init_generic(mmt->pi->wr[d], stride, v, flags, n1);
}

void matmul_top_vec_clear_generic(matmul_top_data_ptr mmt, size_t stride, mmt_generic_vec_ptr v, int d)
{
    if (v == NULL) v = (mmt_generic_vec_ptr) mmt->wr[d]->v;
    unsigned int n1 = mmt->wr[d]->i1 - mmt->wr[d]->i0;
    matmul_aux(mmt->mm, MATMUL_AUX_GET_READAHEAD, &n1);
    vec_clear_generic(mmt->pi->wr[d], stride, v, n1);
}

/* broadcast_down reads data in mmt->wr[d]->v, and broadcasts it across the
 * communicator mmt->pi->wr[d] ; eventually everybody on the communicator
 * mmt->pi->wr[d] has the data.
 *
 * Note that for shuffled product, the combination of reduce_across +
 * broadcast_down is not the identity.
 */
static void
broadcast_down_generic(matmul_top_data_ptr mmt, size_t stride, mmt_generic_vec_ptr v, int d)
{
    // broadcasting down columns is for common agreement on a right
    // source vector, once a left output has been computed.
    //
    // broadcasting down a column is when d == 1.
#ifndef  CONJUGATED_PERMUTATIONS
    choke me;
#endif
    int err;
    if (v == NULL) v = (mmt_generic_vec_ptr) mmt->wr[d]->v;

    pi_wiring_ptr picol = mmt->pi->wr[d];
#ifndef MPI_LIBRARY_MT_CAPABLE
    pi_wiring_ptr pirow = mmt->pi->wr[!d];
#endif

    mmt_wiring_ptr mcol = mmt->wr[d];
    // mmt_wiring_ptr mrow = mmt->wr[!d];

    unsigned int ncols = mmt->n[d];
    unsigned int nrows = mmt->n[!d];

    /* The intersection of our column-wise input range [i0..i1[ with
     * all the fences existing for the output ranges is described by
     * mcol->x ; note how this intersection is common to a
     * complete column.
     */
    /* We used to support (prior to 2011/01/19) the source vector *not*
     * being shared among threads. This has never been actively checked,
     * so despite the fact that the code probably used to work at some
     * point, it was most probably buggy and has been removed.
     * Furthermore, the corresponding code has not been ported to support
     * shuffled mult. The corresponding code has thus been removed.
     */
    ASSERT_ALWAYS(v->flags & THREAD_SHARED_VECTOR);

    pi_log_op(mmt->pi->m, "[%s] enter first loop", __func__);
    /* Make sure that no thread on the column is wandering in other
     * places -- when we're leaving reduce, this is important. */
    serialize_threads(picol);

    /* This loop suffers from two-dimensional serializing, so the
     * simplistic macros SEVERAL_THREADS_PLAY_MPI_BEGIN and END do not
     * help here.
     */
#ifndef MPI_LIBRARY_MT_CAPABLE
    if (mmt->bal->h->flags & FLAG_SHUFFLED_MUL) {
        if (picol->trank == 0) {
            for(unsigned int t = 0 ; t < pirow->ncores ; t++) {
                pi_log_op(pirow, "[%s] serialize_threads", __func__);
                serialize_threads(pirow);
                pi_log_op(pirow, "[%s] serialize_threads done", __func__);
                if (t != pirow->trank)
                    continue;   // not our turn.
                // although the openmpi man page looks funny, I'm assuming that
                // MPI_Allgather wants MPI_IN_PLACE as a sendbuf argument.
                pi_log_op(picol, "[%s] MPI_Allgatherv", __func__);
                ASSERT((mcol->i1 - mcol->i0) % picol->totalsize == 0);
                size_t eblock = (mcol->i1 - mcol->i0) /  picol->totalsize;
                err = MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, mcol->v->v, stride * eblock * picol->ncores, MPI_BYTE, picol->pals);
                pi_log_op(picol, "[%s] MPI_Allgatherv done", __func__);
                ASSERT_ALWAYS(!err);
            }
        }
    } else {
        for(unsigned int t = 0 ; t < pirow->ncores ; t++) {
            pi_log_op(pirow, "[%s] serialize_threads", __func__);
            serialize_threads(pirow);
            pi_log_op(pirow, "[%s] serialize_threads done", __func__);
            if (t != pirow->trank) 
                continue;
            for(unsigned int u = 0 ; u < picol->ncores ; u++) {
                serialize_threads(picol);
                if (u != picol->trank)
                    continue;
                // our turn.
                ASSERT(mcol->i0 % picol->totalsize == 0);
                ASSERT(mcol->i1 % picol->totalsize == 0);
                unsigned int z = ncols / mmt->pi->m->totalsize;
                ASSERT(nrows == ncols);
                ASSERT(mcol->i0 % z == 0);
                ASSERT(mcol->i1 % z == 0);
                unsigned int nv = pirow->totalsize;
                unsigned int nh = picol->totalsize;
                // unsigned int cz = ncols / nv; ASSERT(cz == z * nh);
                unsigned int rz = nrows / nh; ASSERT(rz == z * nv);
                for(unsigned int ii = mcol->i0 ; ii < mcol->i1 ; ) {
                    unsigned int znum = ii / z;
                    unsigned int k = znum / nv;
                    unsigned int offset_there = (znum % nv) * z;
                    unsigned int count = MIN(mcol->i1 - ii, rz - offset_there);
                    unsigned int offset_me = ii - mcol->i0;
                    /*
                       ASSERT(l < mcol->xlen);
                       ASSERT(mcol->x[l].k == (int) k);
                       ASSERT(mcol->x[l].offset_me == offset_me);
                       ASSERT(mcol->x[l].offset_there == offset_there);
                       ASSERT(mcol->x[l].count == count);
                       */
                    if (k % picol->ncores == picol->trank) {
                        size_t off = offset_me * stride;
                        abase_generic_ptr ptr = abase_generic_ptr_add(v->v, off);
                        size_t siz = count * stride;
                        unsigned int root = k / picol->ncores;
                        pi_log_op(picol, "[%s] MPI_Bcast", __func__);
                        err = MPI_Bcast(ptr, siz, MPI_BYTE, root, picol->pals);
                        pi_log_op(picol, "[%s] MPI_Bcast done", __func__);
                        ASSERT_ALWAYS(!err);
                    }
                    ii += count;
                }
            }
        }
        }
#else   /* MPI_LIBRARY_MT_CAPABLE */
    /* Code deleted 20110119, as I've never been able to have enough
     * trust in an MPI implementation to check this */
    ASSERT_ALWAYS(0);
#endif  /* MPI_LIBRARY_MT_CAPABLE */

    pi_log_op(mmt->pi->m, "[%s] trailer", __func__);

    ASSERT_ALWAYS(mcol->i1 <= nrows);
}

#if 0
/* no longer used -- was only used by prep.
 * It's not buggy, but making this work in a context where we have
 * multiple threads is tricky.
 */
void matmul_top_fill_random_source_generic(matmul_top_data_ptr mmt, size_t stride, mmt_generic_vec_ptr v, int d)
{
    if (v == NULL) v = (mmt_generic_vec_ptr) mmt->wr[d]->v;

    // In conjugation mode, it is possible to fill exactly the data chunk
    // that will eventually be relevant. However, it's easy enough to
    // fill our output vector with garbage, and do broadcast_down
    // afterwards...
    if ((v->flags & THREAD_SHARED_VECTOR) == 0 || mmt->pi->wr[d]->trank == 0)
        abase_generic_random(stride, v->v, mmt->wr[d]->i1 - mmt->wr[d]->i0);

    // reconcile all cells which correspond to the same vertical block.
    broadcast_down_generic(mmt, stride, v, d);
}
#endif

/* It really something relevant to the pirow communicator. Turn it so */
static void save_vector_toprow_generic(matmul_top_data_ptr mmt, size_t stride, mmt_generic_vec_ptr v, const char * name, int d, unsigned int iter)
{
    if (v == NULL) v = (mmt_generic_vec_ptr) mmt->wr[d]->v;

    pi_wiring_ptr pirow = mmt->pi->wr[!d];
    mmt_wiring_ptr mcol = mmt->wr[d];

    char * filename;
    int rc;
    rc = asprintf(&filename, COMMON_VECTOR_ITERATE_PATTERN, name, iter, mmt->bal->h->checksum);
    FATAL_ERROR_CHECK(rc < 0, "out of memory");

    size_t mysize = stride * (mcol->i1 - mcol->i0);

    rc = pi_save_file(pirow, filename, (void*) v->v, mysize);

    if (rc == 0) {
        fprintf(stderr, "WARNING: checkpointing failed\n");
    }

    free(filename);
}
// now the exact opposite.

/* The backend is exactly dual to save_vector_toprow_generic.
*/
void load_vector_toprow_generic(matmul_top_data_ptr mmt, size_t stride, mmt_generic_vec_ptr v, const char * name, int d, unsigned int iter)
{
    if (v == NULL) v = (mmt_generic_vec_ptr) mmt->wr[d]->v;

    pi_wiring_ptr pirow = mmt->pi->wr[!d];
    mmt_wiring_ptr mcol = mmt->wr[d];

    char * filename;
    int rc;
    rc = asprintf(&filename, COMMON_VECTOR_ITERATE_PATTERN, name, iter, mmt->bal->h->checksum);
    FATAL_ERROR_CHECK(rc < 0, "out of memory");

    size_t mysize = stride * (mcol->i1 - mcol->i0);

    rc = pi_load_file(pirow, filename, (void*) v->v, mysize);

    if (rc == 0) {
        fprintf(stderr, "WARNING: checkpointing failed\n");
    }

    free(filename);
}

void matmul_top_load_vector_generic(matmul_top_data_ptr mmt, size_t stride, mmt_generic_vec_ptr v, const char * name, int d, unsigned int iter)
{
    int err;
    if (v == NULL) v = (mmt_generic_vec_ptr) mmt->wr[d]->v;

    serialize(mmt->pi->m);
    serialize_threads(mmt->pi->m);

    pi_wiring_ptr picol = mmt->pi->wr[d];
#ifndef MPI_LIBRARY_MT_CAPABLE
    pi_wiring_ptr pirow = mmt->pi->wr[!d];
#endif
    mmt_wiring_ptr mcol = mmt->wr[d];
    // mmt_wiring_ptr mrow = mmt->wr[!d];

    if (picol->jrank == 0 && picol->trank == 0) {
        load_vector_toprow_generic(mmt, stride, v, name, d, iter);
    }

    if (picol->jrank == 0) {
        serialize_threads(mmt->pi->m);
    }
    if (picol->trank == 0) {
        SEVERAL_THREADS_PLAY_MPI_BEGIN(pirow) {
            int err = MPI_Barrier(picol->pals);
            ASSERT_ALWAYS(!err);
        }
        SEVERAL_THREADS_PLAY_MPI_END;
    }
    serialize(mmt->pi->m);
    serialize_threads(mmt->pi->m);

    // after loading, we need to broadcast the data. As in other
    // occcasions, this has to be one one thread at a time unless the MPI
    // library can handle concurrent calls along parallel communicators.

    size_t siz = stride * (mcol->i1 - mcol->i0);
    if (picol->trank == 0) {
        // first, down on columns.
        SEVERAL_THREADS_PLAY_MPI_BEGIN(pirow) {
            void * ptr = v->v;
            err = MPI_Bcast(ptr, siz, MPI_BYTE, 0, picol->pals);
            ASSERT_ALWAYS(!err);
        }
        SEVERAL_THREADS_PLAY_MPI_END;
    }

    serialize(mmt->pi->m);

    // then maybe amongst threads if they're not sharing data.
    if ((v->flags & THREAD_SHARED_VECTOR) == 0) {
        // we have picol->ncores threads, which would be delighted to get
        // access to some common data. The data is known to sit in thread
        // zero, so that's relatively easy.

        if (picol->trank != 0) {
            memcpy(v->v, v->all_v[0], siz);
        }
    }
}

void matmul_top_save_vector_generic(matmul_top_data_ptr mmt, size_t stride, mmt_generic_vec_ptr v, const char * name, int d, unsigned int iter)
{
    if (v == NULL) v = (mmt_generic_vec_ptr) mmt->wr[d]->v;

    // we want row 0 to have everything.
    broadcast_down_generic(mmt, stride, v, d);

    // we'll do some funky stuff, so serialize in order to avoid jokes.
    // Penalty is insignificant anyway since I/O are rare.
    serialize(mmt->pi->m);
    serialize_threads(mmt->pi->m);

    pi_wiring_ptr picol = mmt->pi->wr[d];
    // pi_wiring_ptr pirow = mmt->pi->wr[!d];

    // Now every job row has the complete vector.  We'll do I/O from only
    // one row, that's easier. Pick the topmost one.
    if (picol->jrank == 0 && picol->trank == 0) {
        save_vector_toprow_generic(mmt, stride, v, name, d, iter);
    }

    serialize(mmt->pi->m);
    serialize_threads(mmt->pi->m);
}

//////////////////////////////////////////////////////////////////////////

void matmul_top_vec_init(matmul_top_data_ptr mmt, int d, int flags)
{
    matmul_top_vec_init_generic(mmt, abbytes(mmt->abase, 1), NULL, d, flags);
}

static void mmt_finish_init(matmul_top_data_ptr mmt, int const *, param_list pl, int optimized_direction);
static void matmul_top_read_submatrix(matmul_top_data_ptr mmt, param_list pl, int optimized_direction);

static void mmt_fill_fields_from_balancing(matmul_top_data_ptr mmt, param_list pl)
{
    mmt->n[0] = mmt->bal->trows;
    mmt->n[1] = mmt->bal->tcols;
   
    // nslices[0] is the number of slices in a row -- this is the number of
    // ``vertical slices'', in balance.c parlance.
    unsigned int nslices[2];
    nslices[1] = mmt->bal->h->nh;
    nslices[0] = mmt->bal->h->nv;

    int ok = 1;

    for(int d = 0 ; d < 2 ; d++)
        ok = ok && mmt->pi->wr[d]->totalsize == nslices[d];

    if (!ok) {
        fprintf(stderr, "Configured split %ux%ux%ux%u does not match "
                "on-disk data %ux%u\n",
                mmt->pi->wr[1]->njobs,
                mmt->pi->wr[0]->njobs,
                mmt->pi->wr[1]->ncores,
                mmt->pi->wr[0]->ncores,
                nslices[1],
                nslices[0]);
        serialize(mmt->pi->m);
        exit(1);
    }

    // mmt->ncoeffs_total = 0;
    unsigned int ix[2];

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
        mmt->fences[d^1] = f;
    }

    for(int d = 0 ; d < 2 ; d++)  {
        ix[d] = mmt->pi->wr[1^d]->jrank * mmt->pi->wr[1^d]->ncores;
        ix[d]+= mmt->pi->wr[1^d]->trank;
        mmt->wr[d]->i0 = mmt->fences[d][ix[d]];
        mmt->wr[d]->i1 = mmt->fences[d][ix[d] + 1];
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
        base = cado_strndup(mname, l);
    }

    int rc = asprintf(&mmt->locfile, "%s.%dx%d.%08"PRIx32".h%d.v%d", base, nslices[1], nslices[0], mmt->bal->h->checksum, ix[1], ix[0]);
    free(base);
    ASSERT_ALWAYS(rc >= 0);
}

void matmul_top_init(matmul_top_data_ptr mmt,
        abobj_ptr abase,
        /* matmul_ptr mm, */
        parallelizing_info_ptr pi,
        int const * flags,
        param_list pl,
        int optimized_direction)
{
    memset(mmt, 0, sizeof(*mmt));

    mmt->abase = abase;
    mmt->pi = pi;
    mmt->mm = NULL;

    // n[]
    // ncoeffs_total, ncoeffs,
    // locfile,
    // fences[]
    // wr[]->*
    // are filled in later within read_info_file
    //

    const char * tmp = param_list_lookup_string(pl, "balancing");
    if (!tmp) {
        fprintf(stderr, "Missing parameter balancing=\n");
        exit(1);
    }

    balancing_read_header(mmt->bal, tmp);
    serialize(mmt->pi->m);

    // after that, we need to check for every node if the cache file can
    // be found. If none is found, rebuild the cache files with an
    // equivalent of mf_dobal.
    mmt_finish_init(mmt, flags, pl, optimized_direction);
}

/* Some work has to be done in order to fill the remaining fields in the
 * matmul_top structure.
 */
static void mmt_finish_init(matmul_top_data_ptr mmt, int const * flags, param_list pl, int optimized_direction)
{
    mmt_fill_fields_from_balancing(mmt, pl);
#ifndef  CONJUGATED_PERMUTATIONS
    choke me;
#endif
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
    matmul_top_vec_init(mmt, 0, flags ? flags[0] : 0);
    matmul_top_vec_init(mmt, 1, flags ? flags[1] : 0);

}

static void matmul_top_read_submatrix(matmul_top_data_ptr mmt, param_list pl, int optimized_direction)
{
    const char * impl = param_list_lookup_string(pl, "mm_impl");
    int rebuild = 0;
    unsigned int cache_nbys = 0;
    param_list_parse_int(pl, "rebuild_cache", &rebuild);
    param_list_parse_uint(pl, "cache_nbys", &cache_nbys);

    unsigned int normal_nbys = abnbits(mmt->abase);

    if (cache_nbys) {
        abobj_set_nbys(mmt->abase, cache_nbys);
    }

    mmt->mm = matmul_init(mmt->abase,
            mmt->wr[0]->i1 - mmt->wr[0]->i0,
            mmt->wr[1]->i1 - mmt->wr[1]->i0,
            mmt->locfile, impl, pl, optimized_direction); 

    // *IF* we need to do a collective read of the matrix, we need to
    // provide the pointer *now*.
    unsigned int sqread = 0;
    param_list_parse_uint(pl, "sequential_cache_read", &sqread);

    int cache_loaded = 0;

    if (!rebuild) {
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
    if (!global_data_eq(mmt->pi, &cache_loaded, sizeof(int))) {
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
        SEVERAL_THREADS_PLAY_MPI_END;
        serialize(mmt->pi->m);
        abort();
    }

    unsigned int sqb = 0;
    param_list_parse_uint(pl, "sequential_cache_build", &sqb);

    uint32_t * matrix_ptr = NULL;

    if (!cache_loaded) {
        matrix_u32 m;
        if (mmt->pi->m->jrank == 0 && mmt->pi->m->trank == 0) {
            fprintf(stderr, "Matrix dispatching starts\n");
        }
        m->mfile = param_list_lookup_string(pl, "matrix");
        m->bfile = param_list_lookup_string(pl, "balancing");
        // the mm layer is informed of the higher priority computations
        // that will take place. Depending on the implementation, this
        // may cause the direct or transposed ordering to be preferred.
        // Thus we have to read this back from the mm structure.
        m->transpose = mmt->mm->store_transposed;
        balancing_get_matrix_u32(mmt->pi, pl, m);
        matrix_ptr = m->p;
        // note that this pointer will not necessarily be freed
        // immediately. Ownership is transferred to the mm layer, which
        // may decide to keep the pointer for the computations. This is
        // what the basic layer does. In contrast, the bucket layer
        // discards this data, as the post-processed matrix form is more
        // efficient eventually.

        int ssm = 0;
        param_list_parse_int(pl, "save_submatrices", &ssm);
        if (ssm) {
            fprintf(stderr, "DEBUG: creating %s\n", mmt->locfile);
            FILE * f = fopen(mmt->locfile, "w");
            fwrite(m->p, sizeof(uint32_t), m->size, f);
            fclose(f);
        }
    }


    if (!sqb) {
        if (!cache_loaded) {
            // everybody does it in parallel
            fprintf(stderr,"[%s] J%uT%u building cache for %s\n",
                    mmt->pi->nodenumber_s,
                    mmt->pi->m->jrank,
                    mmt->pi->m->trank,
                    mmt->locfile);
            matmul_build_cache(mmt->mm, matrix_ptr);
            matmul_save_cache(mmt->mm);
        }
    } else {
        for(unsigned int j = 0 ; j < mmt->pi->m->ncores + 1 ; j++) {
            serialize_threads(mmt->pi->m);
            if (cache_loaded) continue;
            if (j == mmt->pi->m->trank) {
                fprintf(stderr,"[%s] J%uT%u building cache for %s\n",
                        mmt->pi->nodenumber_s,
                        mmt->pi->m->jrank,
                        mmt->pi->m->trank,
                        mmt->locfile);
                matmul_build_cache(mmt->mm, matrix_ptr);
            } else if (j == mmt->pi->m->trank + 1) {
                matmul_save_cache(mmt->mm);
            }
        }
    }

    my_pthread_mutex_lock(mmt->pi->m->th->m);
    fprintf(stderr, "[%s] J%uT%u uses cache file %s\n",
            mmt->pi->nodenumber_s,
            mmt->pi->m->jrank, mmt->pi->m->trank,
            /* cache for mmt->locfile, */
            mmt->mm->cachefile_name);
    my_pthread_mutex_unlock(mmt->pi->m->th->m);

    abobj_set_nbys(mmt->abase, normal_nbys);
}

void matmul_top_vec_clear(matmul_top_data_ptr mmt, int d)
{
    matmul_top_vec_clear_generic(mmt, abbytes(mmt->abase, 1), NULL, d);
}


void matmul_top_clear(matmul_top_data_ptr mmt, abobj_ptr abase MAYBE_UNUSED)
{
    serialize(mmt->pi->m);
    matmul_top_vec_clear(mmt,0);
    matmul_top_vec_clear(mmt,1);
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

/* XXX
 * ``across'' (horizontally) and ``down'' (vertically) are just here for
 * exposition. The binding of this operation to job/thread arrangement
 * is flexible, through argument d.
 * XXX
 */
static void
broadcast_down(matmul_top_data_ptr mmt, int d)
{
    broadcast_down_generic(mmt, abbytes(mmt->abase, 1), NULL, d);
}

/* input: targ[i]=w_part[ydim*i+ycoord],
   output: rop[i]=(w[i] on this node) */
#if 0
void collect_row_r(vector rop, mpi_sparse_matrix M, vector src, vector targ)
{
    MPI_Status commstat;
    u32_t i, j, k, l, j0, j1;
    int sr_count;
    vector buf1, buf2, buf;

    sr_count = M.veclen_per_host * ULL2UL_MULTIPLIER;
    memset(rs_buffer1, 0, M.veclen_per_host * sizeof(*rs_buffer1));
    buf1 = rs_buffer1;
    buf2 = rs_buffer2;
    if (M.my_x)
	l = M.my_x - 1;
    else
	l = M.hsize - 1;
    for (i = 0; i < M.hsize; i++) {	/* data for node (l,my_y) */
	j0 = l * M.veclen_per_host;
	j1 = j0 + M.veclen_per_host;
	if (j1 > M.rows_per_host)
	    j1 = M.rows_per_host;
	for (j = j0, k = 0; j < j1; j++, k++)
	    buf1[k] ^= targ[j];
	if (i == M.hsize - 1)
	    break;
	MPI_Sendrecv(buf1, sr_count, MPI_U32_T, M.x_dest, TAG_CR0,
		      buf2, sr_count, MPI_U32_T, M.x_src, TAG_CR0,
		      M.hcomm, &commstat);
	buf = buf1;
	buf1 = buf2;
	buf2 = buf;
	if (l)
	    l--;
	else
	    l = M.hsize - 1;
    }
    memcpy(rop, buf1, M.veclen_per_host * sizeof(*rop));
}
#endif

/* Note that because reduce_across does additions, it is _NOT_ generic.
 *
 * reduce_across reads data in mmt->wr[d]->v, sums it up across the
 * communicator mmt->pi->wr[d], and collects the results in the areas pointed
 * to by mmt->wr[!d]->v
 *
 * Note that for shuffled product, the combination of reduce_across +
 * broadcast_down is not the identity.
 */
static void
reduce_across(matmul_top_data_ptr mmt, int d)
{
#ifndef  CONJUGATED_PERMUTATIONS
    choke me;
#endif

    int err;

    /* reducing across a row is when d == 0 */
    pi_wiring_ptr pirow = mmt->pi->wr[d];
#ifndef MPI_LIBRARY_MT_CAPABLE
    pi_wiring_ptr picol = mmt->pi->wr[!d];
#endif

    mmt_wiring_ptr mrow = mmt->wr[d];
    mmt_wiring_ptr mcol = mmt->wr[!d]; unsigned int nrows = mmt->n[d];
    unsigned int ncols = mmt->n[!d];

    pi_log_op(mmt->pi->m, "[%s] enter first loop", __func__);

    // I don't think that the case of shared vectors has been tested
    // correctly for reduction.
    assert((mmt->wr[d]->v->flags & THREAD_SHARED_VECTOR) == 0);

    if (pirow->ncores > 1) {
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
         * i0+k*(i1-i0)/pirow->ncores and further). As a convention, for
         * the shuffled product, the thread which eventually owns the
         * final data is thread 0. In the ``normal case'', the exact
         * thread which eventually owns the final data is not necessarily
         * the same for the whole range, nor does it carry any
         * relationship with the thread doing the computation !
         *
         * Note that one should consider that data in threads other than
         * the destination thread may be clobbered by the operation
         * (although in the present implementation it is not).
         */
        ASSERT(mrow->i0 % pirow->ncores == 0);
        ASSERT(mrow->i1 % pirow->ncores == 0);
        unsigned int z = nrows / mmt->pi->m->totalsize;
        unsigned int z1 = z * pirow->njobs;
        ASSERT(nrows == ncols);
        ASSERT(mrow->i0 % z1 == 0);
        ASSERT(mrow->i1 % z1 == 0);
        unsigned int nv = pirow->totalsize;
        unsigned int nh = picol->totalsize;
        unsigned int cz = ncols / nv; ASSERT(cz == z * nh);
        // unsigned int rz = nrows / nh; ASSERT(rz == z * nv);
        unsigned int ii0 = mrow->i0 + pirow->trank * z1;
        unsigned int ii1 = ii0 + z1;
        if (mmt->bal->h->flags & FLAG_SHUFFLED_MUL) {
            abt * dptr = mrow->v->all_v[0];
            for(uint32_t w = 1 ; w < pirow->ncores ; w++) {
                const abt * sptr = mrow->v->all_v[w];
                for(uint32_t ii = ii0 ; ii < ii1 ; ii++) {
                    abadd(mmt->abase, dptr + aboffset(mmt->abase, ii - mrow->i0),
                            sptr + aboffset(mmt->abase, ii - mrow->i0));
                    // TODO: for non-binary fields, here we would have an
                    // unreduced add (of probably already unreduced data
                    // as they come out of the matrix-vector
                    // multiplication), followed later on by a reduction
                    // in dptr.
                }
            }
        } else {
            for(unsigned int ii = ii0 ; ii < ii1 ; ) {
                unsigned int znum = ii / z;
                unsigned int k = znum / nh;
                unsigned int offset_there = (znum % nh) * z;
                unsigned int count = MIN(ii1 - ii, cz - offset_there);
                unsigned int offset_me = ii - mrow->i0;
                /*
                ASSERT(l < mrow->ylen);
                ASSERT(mrow->y[l].k == (int) k);
                ASSERT(mrow->y[l].offset_me == offset_me);
                ASSERT(mrow->y[l].offset_there == offset_there);
                ASSERT(mrow->y[l].count == count);
                */

                unsigned int dst = k % pirow->ncores;
                /* Note that ``offset_there'' does not count here, since
                 * at this point we're not yet flipping to the buffer in
                 * the other direction.
                 */
                unsigned int off = aboffset(mmt->abase, offset_me);
                abt * dptr = mrow->v->all_v[dst] + off;
                // size_t siz = abbytes(mmt->abase, count);
                ASSERT(offset_me + count <= mrow->i1 - mrow->i0);
                /* j indicates a thread number offset -- 0 means
                 * destination, and so on. all_v has been allocated with
                 * wraparound pointers precisely for accomodating the
                 * hack here. */
                for(unsigned int j = 1 ; j < pirow->ncores ; j++) {
                    /* A given data offset is summed by only one of the row
                     * threads. Otherwise, we're clearly creating rubbish
                     * since sptr and dptr lie within mrow->v
                     */
                    const abt * sptr = mrow->v->all_v[dst+j] + off;
                    for(unsigned int k = 0 ; k < count ; k++) {
                        abadd(mmt->abase, dptr + aboffset(mmt->abase, k),
                                sptr + aboffset(mmt->abase, k));
                        // TODO: for non-binary fields, here we would
                        // have an unreduced add (of probably already
                        // unreduced data as they come out of the
                        // matrix-vector multiplication), followed later
                        // on by a reduction in dptr.
                    }
                }
                ii += count;
            }
        }
        pi_log_op(pirow, "[%s] thread reduction done", __func__);
    }

    /* Good. Now on each node, there's one thread (thread 0 in the
     * shuffled case) whose buffer contains the reduced data for the
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
    if (mmt->bal->h->flags & FLAG_SHUFFLED_MUL) {
        serialize_threads(pirow);
    } else {
        serialize_threads(mmt->pi->m);
    }
    pi_log_op(mmt->pi->m, "[%s] serialize_threads done", __func__);

#ifndef MPI_LIBRARY_MT_CAPABLE
    if (mmt->bal->h->flags & FLAG_SHUFFLED_MUL) {
        if (pirow->trank == 0) {
            for(unsigned int t = 0 ; t < picol->ncores ; t++) {
                pi_log_op(picol, "[%s] serialize_threads", __func__);
                serialize_threads(picol);
                pi_log_op(picol, "[%s] serialize_threads done", __func__);
                if (t != picol->trank)
                    continue;   // not our turn.

                abase_generic_ptr dptr = mrow->v->all_v[0];
                pi_log_op(pirow, "[%s] MPI_Reduce_scatter", __func__);
                // all recvcounts are now equal
                int * rc = malloc(pirow->njobs * sizeof(int));
                ASSERT((mrow->i1 - mrow->i0) % pirow->totalsize == 0);
                for(unsigned int k = 0 ; k < pirow->njobs ; k++) {
                    rc[k] = (mrow->i1 - mrow->i0) / pirow->njobs;
                    rc[k] = abbytes(mmt->abase, rc[k]);
                }
                err = MPI_Reduce_scatter(MPI_IN_PLACE, dptr, rc, MPI_BYTE, MPI_BXOR, pirow->pals);
                free(rc);
                pi_log_op(pirow, "[%s] MPI_Reduce_scatter done", __func__);
                ASSERT_ALWAYS(!err);
            }
        }
        serialize_threads(pirow);
        // row threads pick what they're interested in in thread0's reduced
        // buffer. Writes are non-overlapping in the mcol buffer here.
        // Different row threads always have different mcol buffers, and
        // sibling col threads write to different locations in their
        // (generally shared) mcol buffer, depending on which row they
        // intersect.
        
        // Notice that results are being written inside mrow->v->all_v[0],
        // just packed at the beginning.

        size_t eblock = (mrow->i1 - mrow->i0) / pirow->totalsize;

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
        // Now amond the picol->totalsize blocks of the col buffer, this
        // will go to position cj * picol->ncores + ct
        size_t stride = abbytes(mmt->abase, 1);
        abase_generic_ptr sptr = abase_generic_ptr_add(mrow->v->all_v[0], stride * (pirow->trank * eblock));
        abase_generic_ptr dptr = abase_generic_ptr_add(mcol->v->v, stride * ((picol->jrank * picol->ncores + picol->trank) * eblock));
        abase_generic_copy(stride, dptr, sptr, eblock);
    } else {
        /* We need two levels of locking :-( */
        for(unsigned int t = 0 ; t < picol->ncores ; t++) {
            pi_log_op(picol, "[%s] serialize_threads", __func__);
            serialize_threads(picol);
            pi_log_op(picol, "[%s] serialize_threads done", __func__);
            if (t != picol->trank)
                continue;   // not our turn.

            // all row threads are likely to engage in a collective operation
            // at one moment or another.
            ASSERT(mrow->i0 % pirow->totalsize == 0);
            ASSERT(mrow->i1 % pirow->totalsize == 0);
            unsigned int z = nrows / mmt->pi->m->totalsize;
            ASSERT(nrows == ncols);
            ASSERT(mrow->i0 % z == 0);
            ASSERT(mrow->i1 % z == 0);
            unsigned int nv = pirow->totalsize;
            unsigned int nh = picol->totalsize;
            unsigned int cz = ncols / nv; ASSERT(cz == z * nh);
            // unsigned int rz = nrows / nh; ASSERT(rz == z * nv);
            for(unsigned int ii = mrow->i0 ; ii < mrow->i1 ; ) {
                unsigned int znum = ii / z;
                unsigned int k = znum / nh;
                unsigned int offset_there = (znum % nh) * z;
                unsigned int count = MIN(mrow->i1 - ii, cz - offset_there);
                unsigned int offset_me = ii - mrow->i0;
                /*
                   ASSERT(l < mrow->xlen);
                   ASSERT(mrow->x[l].k == (int) k);
                   ASSERT(mrow->x[l].offset_me == offset_me);
                   ASSERT(mrow->x[l].offset_there == offset_there);
                   ASSERT(mrow->x[l].count == count);
                   */
                unsigned int tdst = k % pirow->ncores;
                serialize_threads(pirow);
                if (pirow->trank == tdst) {
                    unsigned int jdst = k / pirow->ncores;
                    // MPI_Reduce does not know about const-ness !
                    abt * sptr = mrow->v->v + aboffset(mmt->abase, offset_me);
                    abt * dptr = mcol->v->v + aboffset(mmt->abase, offset_there);
                    size_t siz = abbytes(mmt->abase, count);
                    pi_log_op(pirow, "[%s] MPI_Reduce", __func__);
                    err = MPI_Reduce(sptr, dptr, siz, MPI_BYTE, MPI_BXOR, jdst, pirow->pals);
                    pi_log_op(pirow, "[%s] MPI_Reduce done", __func__);
                    ASSERT_ALWAYS(!err);
                }
                ii+=count;
            }
        }
    }
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

#if 0
/* no longer used -- was only used by prep.
 * It's not buggy, but making this work in a context where we have
 * multiple threads is tricky.
 */
void matmul_top_fill_random_source(matmul_top_data_ptr mmt, int d)
{
    matmul_top_fill_random_source_generic(mmt, abbytes(mmt->abase, 1), NULL, d);
}
#endif

#if 0
static void mmt_debug_writeout(matmul_top_data_ptr mmt, int d, const char * name)
{
    // serialize(mmt->pi->m);
    debug_write(mmt->wr[d]->v->v,
            abbytes(mmt->abase, mmt->wr[d]->i1 - mmt->wr[d]->i0),
            "%s.j%u.t%u", name, mmt->pi->m->jrank, mmt->pi->m->trank);
    // serialize(mmt->pi->m);
}
#endif

/* Takes data in mmt->wd[d]->v, and compute the corresponding partial result in
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
    mmt_wiring_ptr mrow = mmt->wr[!d];

    mmt_wiring_ptr mcol = mmt->wr[d];

    ASSERT_ALWAYS((mrow->v->flags & THREAD_SHARED_VECTOR) == 0);

    pi_log_op(mmt->pi->m, "[%s] enter matmul_mul", __func__);
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
    pi_wiring_ptr picol = mmt->pi->wr[d];
    mmt_wiring_ptr mcol = mmt->wr[d];

    pi_log_op(mmt->pi->m, "[%s] enter reduce_across", __func__);
    reduce_across(mmt, !d);
    pi_log_op(mmt->pi->m, "[%s] enter broadcast_down", __func__);
    broadcast_down(mmt, d);

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
void matmul_top_save_vector(matmul_top_data_ptr mmt, const char * name, int d, unsigned int iter)
{
    matmul_top_save_vector_generic(mmt, abbytes(mmt->abase, 1), NULL, name, d, iter);
}

void matmul_top_load_vector(matmul_top_data_ptr mmt, const char * name, int d, unsigned int iter)
{
    matmul_top_load_vector_generic(mmt, abbytes(mmt->abase, 1), NULL, name, d, iter);
}
