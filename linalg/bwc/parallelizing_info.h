#ifndef PARALLELIZING_INFO_H_
#define PARALLELIZING_INFO_H_

#include "select_mpi.h"

#ifdef  HOMEMADE_BARRIERS
#include "barrier.h"
#endif

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

/* don't enable this. Clutters output a lot */
#define xxxCONCURRENCY_DEBUG

/*
 * MPI_LIBRARY_MT_CAPABLE: do mpi communications in a furiously
 * concurrent manner.
 *
 * This isn't ready yet. In presence of an MPI library which correctly
 * supports MPI_THREAD_MULTIPLE, chances are that this gets close to
 * working ok, and even maybe improve performance quite a bit.
 *
 * At the moment, with openmpi 1.2.8 compiled with
 * --enable-mpi-threads, it does not work (crash observed at
 * MPI_Init_thread -- fairly weird).
 *
 * Mpich2 1.1a2 works, apparently. And it even seems to end up being
 * noticeably faster for my toy examples (= latency-wise).
 */

#if defined(MPICH2) && MPICH2_NUMVERSION >= 10100002
#define MPI_LIBRARY_MT_CAPABLE
/*
 * at present I know of no version of openmpi with MPI_THREAD_MULTIPLE
 * working.
#elif defined(OPEN_MPI) && OMPI_MAJOR_VERSION >= 123456789
#define MPI_LIBRARY_MT_CAPABLE
 */
#else
/* Assume it does not work */
#endif


/* utility structure. It is stored in thread-shared memory */
struct pthread_things {
#ifdef  HOMEMADE_BARRIERS
    barrier_t b[1];
#else
    my_pthread_barrier_t b[1];
#endif

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
    unsigned int jcommon;
    unsigned int jrank;
    unsigned int tcommon;
    unsigned int trank;
    MPI_Comm pals;

#ifdef  HOMEMADE_BARRIERS
    unsigned long shared_data;
#endif

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
 * the thread having index i ; i == 0 is a reasonable start in general. In
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

/* prints the given string in a ascii-form matrix. */
extern void grid_print(parallelizing_info_ptr pi, char * buf, size_t siz, int print);

#define serialize(w)   serialize__(w, __FILE__, __LINE__)
extern int serialize__(pi_wiring_ptr, const char *, unsigned int);
#define serialize_threads(w)   serialize_threads__(w, __FILE__, __LINE__)
extern int serialize_threads__(pi_wiring_ptr, const char *, unsigned int);

extern int pi_save_file(pi_wiring_ptr, const char *, void *, size_t);
extern int pi_load_file(pi_wiring_ptr, const char *, void *, size_t);

extern int pi_load_file_2d(parallelizing_info_ptr, int, const char *, void *, size_t);
extern int pi_save_file_2d(parallelizing_info_ptr, int, const char *, void *, size_t);

#ifdef __cplusplus
}
#endif

/* This provides a fairly typical construct */
#ifndef MPI_LIBRARY_MT_CAPABLE
#define SEVERAL_THREADS_PLAY_MPI_BEGIN(comm)				\
    for(unsigned int t__ = 0 ; t__ < comm->ncores ; t__++) {		\
        serialize_threads(comm);					\
        if (t__ != comm->trank) continue; // not our turn.
#define SEVERAL_THREADS_PLAY_MPI_END    }
#else
#define SEVERAL_THREADS_PLAY_MPI_BEGIN(comm)    /**/
#define SEVERAL_THREADS_PLAY_MPI_END            /**/
#endif

#endif	/* PARALLELIZING_INFO_H_ */
