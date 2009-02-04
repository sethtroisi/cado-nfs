#define _POSIX_C_SOURCE 200112L

#include "xymats.h"
#include "manu.h"

/* This is used as a very simple reduction routine, mostly for the small
 * m*n matrices. It has nothing to do with the more involved reduce
 * operations relevant to vectors. Here, the data we're dealing with is
 * tiny.
 *
 * This operation serializes threads.
 */
void reduce_generic_threadlevel(abobj_t abase, mmt_vec_ptr v, pi_wiring_ptr wr, unsigned int n)
{
    if (v->flags & THREAD_SHARED_VECTOR)
        return; /* If the vector is shared, reduction is pointless */

    /* Do it star-like. Easy enough. */
    serialize_threads(wr);
    if (wr->trank == 0) {
        for(unsigned int t = 1 ; t < wr->ncores ; t++) {
            for(unsigned int k = 0 ; k < n ; k++) {
                abadd(abase, v->v + aboffset(abase,k),
                        v->all_v[t] + aboffset(abase,k));
            }
        }
    }
    serialize_threads(wr);
    if (wr->trank > 0) {
        abcopy(abase, v->v, v->all_v[0], n);
    }
}

void reduce_generic_mpilevel(abobj_t abase, mmt_vec_ptr v, pi_wiring_ptr wr, unsigned int n)
{
    if (wr->trank == 0) {
        size_t siz = abbytes(abase, n);
        int err = MPI_Allreduce(MPI_IN_PLACE, v->v, siz, MPI_BYTE, MPI_BXOR, wr->pals);
        BUG_ON(err);
    }
    serialize_threads(wr);

    if (wr->trank > 0 && (v->flags & THREAD_SHARED_VECTOR) == 0) {
        /* recover the reduced data from thread zero */
        abcopy(abase, v->v, v->all_v[0], n);
    }
}

/* This combines both, with fewer thread locks */
void reduce_generic(abobj_t abase, mmt_vec_ptr v, pi_wiring_ptr wr, unsigned int n)
{
    if (v->flags & THREAD_SHARED_VECTOR) {
        reduce_generic_mpilevel(abase, v, wr, n);
    }

    /* Do it star-like. Easy enough. */
    serialize_threads(wr);
    if (wr->trank == 0) {
        for(unsigned int t = 1 ; t < wr->ncores ; t++) {
            for(unsigned int k = 0 ; k < n ; k++) {
                abadd(abase, v->v + aboffset(abase,k),
                        v->all_v[t] + aboffset(abase,k));
            }
        }
        size_t siz = abbytes(abase, n);
        int err = MPI_Allreduce(MPI_IN_PLACE, v->v, siz, MPI_BYTE, MPI_BXOR, wr->pals);
        BUG_ON(err);
    }
    serialize_threads(wr);
    if (wr->trank > 0) {
        abcopy(abase, v->v, v->all_v[0], n);
    }
}

#if 0
#ifndef MPI_LIBRARY_MT_CAPABLE
    for(unsigned int t = 0 ; t < wr->ncores ; t++) {
        serialize_threads(wr);
        if (t != wr->trank)
            continue;   // not our turn.
        abt * dptr = xy_mats2;
        size_t siz = abbytes(abase, m * NBITER);
        int err = MPI_Allreduce(MPI_IN_PLACE, dptr, siz, MPI_BYTE, MPI_BXOR, wr->pals);
        BUG_ON(err);
    }
#else   /* MPI_LIBRARY_MT_CAPABLE */
    {
        abt * dptr = xy_mats2;
        size_t siz = abbytes(abase, m * NBITER);
        int err = MPI_Allreduce(MPI_IN_PLACE, dptr, siz, MPI_BYTE, MPI_BXOR, wr->pals);
        BUG_ON(err);
    }
#endif
#endif
