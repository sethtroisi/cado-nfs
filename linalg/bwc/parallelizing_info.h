#ifndef PARALLELIZING_INFO_H_
#define PARALLELIZING_INFO_H_

#include "select_mpi.h"

/*
 * The main guy here is the parallelizing_info data type. It is commonly
 * declared as pi.
 *
 * Threads and jobs are arranged in a grid-like manner. Several logical
 * communication group are defined. One for the global grid, as well as
 * one for each column / row.
 *
 * pi->m denotes the global group. There is only one such group. All
 * jobs, all threads contribute to it. At the mpi level, the communicator
 * which is relevant is of course MPI_COMM_WORLD.
 *
 * pi->wr[0] denotes horizontal, a.k.a. ROW groups.
 * pi->wr[1] denotes vertical, a.k.a. COLUMN groups.
 *
 * When a matrix is mapped to a grid process of this kind, say a matrix
 * having been split in nhs * nvs slices, then there are exactly nhs ROW
 * groups, and nvs COLUMN groups. ROW groups consist of nvs (job,thread)
 * items (as many as one finds COL groups), and conversely.
 */

/* utility structure. It is stored in thread-shared memory */
struct pthread_things {
    my_pthread_barrier_t b[1];
    my_pthread_mutex_t m[1];
    char * desc;
    void * utility_ptr;
    // int count;
};

struct pi_wiring_s {
    /* njobs : number of mpi jobs concerned by this logical group */
    /* ncores : number of threads concerned by this logical group */
    unsigned int njobs;
    unsigned int ncores;

    /* product njobs * ncores */
    unsigned int totalsize;
    unsigned int jcommon;       // was j
    unsigned int jrank; // was rank
    unsigned int tcommon;
    unsigned int trank; // was c
    MPI_Comm pals;

    struct pthread_things * th;
#ifdef  CONCURRENCY_DEBUG
    int th_count;
#endif
};

typedef struct pi_wiring_s pi_wiring[1];
typedef struct pi_wiring_s * pi_wiring_ptr;
typedef const struct pi_wiring_s * pi_wiring_srcptr;


struct parallelizing_info_s {
    // row-wise, column-wise.
    pi_wiring wr[2];
    // main.
    pi_wiring m;
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
extern void pi_go(void *(*fcn)(parallelizing_info_ptr, void * arg),
        unsigned int nhj, unsigned int nvj,
        unsigned int nhc, unsigned int nvc,
        void * arg);

extern void hello(parallelizing_info_ptr pi);

/* This must be viewed as a companion to MPI_Bcast. Threads must agree on
 * a shared area. The area is the *ptr value as presented on input by
 * the thread having index i ; i ==0 is a reasonable start in general. In
 * any case, it must be an integer between 0 and ncores-1.
 *
 * This function uses the wr->utility_ptr field.
 */
extern void thread_agreement(pi_wiring_ptr wr, void ** ptr, unsigned int i);

/* This one is the higher level thingy on top of MPI_Bcast and the
 * previous one. It broadcast data computed by job j, thread t (e.g. 0,0)
 * to all jobs and threads.
 */
extern void complete_broadcast(pi_wiring_ptr wr, void * ptr, size_t size, unsigned int j, unsigned int t);

#define serialize(w)   serialize__(w, __FILE__, __LINE__)
extern int serialize__(pi_wiring_ptr, const char *, unsigned int);
#define serialize_threads(w)   serialize_threads__(w, __FILE__, __LINE__)
extern int serialize_threads__(pi_wiring_ptr, const char *, unsigned int);

#ifdef __cplusplus
}
#endif

#endif	/* PARALLELIZING_INFO_H_ */
