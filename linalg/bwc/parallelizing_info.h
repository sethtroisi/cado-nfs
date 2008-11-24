#ifndef PARALLELIZING_INFO_H_
#define PARALLELIZING_INFO_H_

#include "select_mpi.h"

struct pthread_things {
    my_pthread_barrier_t b[1];
    my_pthread_mutex_t m[1];
    char * desc;
    // int count;
};

struct pi_wiring {
    unsigned int nslices;
    unsigned int njobs;
    unsigned int ncores;
    unsigned int jcommon;       // was j
    unsigned int jrank; // was rank
    unsigned int tcommon;
    unsigned int trank; // was c
    MPI_Comm pals;

    struct pthread_things * th;
    int th_count;
};

struct parallelizing_info_s {
    // row-wise, column-wise.
    struct pi_wiring wr[2][1];
    // main.
    struct pi_wiring m[1];
};

typedef struct parallelizing_info_s parallelizing_info[1];
typedef struct parallelizing_info_s * parallelizing_info_ptr;
typedef const struct parallelizing_info_s * parallelizing_info_srcptr;

#ifdef __cplusplus
extern "C" {
#endif

#if 0
static inline unsigned int pi_jid(parallelizing_info_srcptr pi)
{ return pi->wr[1]->j + pi->wr[1]->njobs * pi->wr[0]->j; }
static inline unsigned int pi_tid(parallelizing_info_srcptr pi)
{ return pi->wr[1]->c + pi->wr[1]->ncores * pi->wr[0]->c; }
static inline unsigned int pi_id(parallelizing_info_srcptr pi)
{ return pi_tid(pi) + pi->wr[0]->ncores * pi->wr[1]->ncores * pi_jid(pi); }
static inline unsigned int pi_src_chunk(parallelizing_info_srcptr pi)
{ return pi->wr[1]->c + pi->wr[1]->ncores * pi->wr[1]->j; }
static inline unsigned int pi_dst_chunk(parallelizing_info_srcptr pi)
{ return pi->wr[0]->c + pi->wr[0]->ncores * pi->wr[0]->j; }
#endif

/* pi_go is the main function. It is responsible of creating all the
 * parallelizing_info data structures, set up the different
 * inter-job and inter-thread conciliation toys (communicators, pthread
 * barriers, and so on), and eventually run the provided function.
 *
 * nhj, nvj denotes the number of jobs in each column, and in each row,
 * respectively.
 *
 * nhc, nvc are the same for threads (cores).
 */
extern void pi_go(void *(*fcn)(parallelizing_info_ptr),
        unsigned int nhj, unsigned int nvj,
        unsigned int nhc, unsigned int nvc);

extern void hello(parallelizing_info_ptr pi);

#ifdef __cplusplus
}
#endif

#endif	/* PARALLELIZING_INFO_H_ */
