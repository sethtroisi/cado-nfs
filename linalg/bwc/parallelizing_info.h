#ifndef PARALLELIZING_INFO_H_
#define PARALLELIZING_INFO_H_

#include <sys/time.h>   /* for struct timeval */

#include "params.h"
#include "select_mpi.h"

/*
 * The main guy here is the parallelizing_info data type. It is commonly
 * declared as pi.
 *
 * Threads and jobs are arranged in a grid-like manner. Several logical
 * communication groups are defined. One for the global grid, as well as
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
 *
 * As for the ``fake mpi'' implementation we have, well, it's up to the
 * reader to decide. All collectives are nops, so it's trivially capable.
 * OTOH, not setting the flags ensures that the rest of the code compiles
 * ok.
 */

#if defined(MPICH2) && MPICH2_NUMVERSION >= 10100002
/* In fact, even in this case we might consider disabling it. */
#define xxxMPI_LIBRARY_MT_CAPABLE
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
    barrier_t bh[1];
    my_pthread_barrier_t b[1];

    my_pthread_mutex_t m[1];
    char * desc;
    void * utility_ptr;
    // int count;
};

/* need to forward-define this for log entries */
struct pi_wiring_s;
typedef struct pi_wiring_s * pi_wiring_ptr;
typedef const struct pi_wiring_s * pi_wiring_srcptr;

struct pi_log_entry {
    struct timeval tv[1];
    char what[80];
};

#define PI_LOG_BOOK_ENTRIES     32
struct pi_log_book {
    struct pi_log_entry t[PI_LOG_BOOK_ENTRIES];
    int hsize;  // history size -- only a count, once the things wraps.
    int next;   // next free pointer.
};

struct pi_wiring_s {
    /* njobs : number of mpi jobs concerned by this logical group */
    /* ncores : number of threads concerned by this logical group */
    unsigned int njobs;
    unsigned int ncores;

    /* product njobs * ncores */
    unsigned int totalsize;
    unsigned int jrank;
    unsigned int trank;
    MPI_Comm pals;

    struct pthread_things * th;
#ifdef  CONCURRENCY_DEBUG
    int th_count;
#endif
    struct pi_log_book * log_book;
};

typedef struct pi_wiring_s pi_wiring[1];

struct pi_interleaving_s {
    int idx;                    /* 0 or 1 */
    my_pthread_barrier_t * b;   /* not a 1-sized array on purpose --
                                   being next to index, it can't ! */
};
typedef struct pi_interleaving_s pi_interleaving[1];
typedef struct pi_interleaving_s * pi_interleaving_ptr;

struct pi_dictionary_entry_s;
typedef struct pi_dictionary_entry_s * pi_dictionary_entry_ptr;
struct pi_dictionary_entry_s {
    unsigned long key;
    unsigned long who;
    void * value;
    pi_dictionary_entry_ptr next;
};
typedef struct pi_dictionary_entry_s pi_dictionary_entry[1];

struct pi_dictionary_s {
    my_pthread_rwlock_t m[1];
    pi_dictionary_entry_ptr e;
};
typedef struct pi_dictionary_s pi_dictionary[1];
typedef struct pi_dictionary_s * pi_dictionary_ptr;

#define PI_NAMELEN      32
struct parallelizing_info_s {
    // row-wise, column-wise.
    pi_wiring wr[2];
    // main.
    pi_wiring m;
    pi_interleaving_ptr interleaved;
    pi_dictionary_ptr dict;
    char nodename[PI_NAMELEN];
    char nodeprefix[PI_NAMELEN];
    char nodenumber_s[PI_NAMELEN];
    /* This pointer is identical on all threads. It is non-null only in
     * case we happen to have sufficiently recent gcc, together with
     * sufficiently recent hwloc */
    void * cpubinding_info;
    int thr_orig[2];            /* when only_mpi is 1, this is what the
                                   thr parameter was set to originally.
                                   Otherwise we have {0,0} here. */
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
 * the param_list is checked for parameters mpi and thr, so as to define
 * the mpi and thr splttings.
 *
 * nhc, nvc are the same for threads (cores).
 */
extern void pi_go(void *(*fcn)(parallelizing_info_ptr, param_list pl, void * arg),
        param_list pl,
        void * arg);

extern void hello(parallelizing_info_ptr pi);

/* This must be viewed as a companion to MPI_Bcast. Threads must agree on
 * a shared area. The area is the *ptr value as presented on input by
 * the thread having index i ; i == 0 is a reasonable start in general. In
 * any case, it must be an integer between 0 and ncores-1.
 *
 * This function uses the wr->utility_ptr field.
 */
extern void thread_broadcast(pi_wiring_ptr wr, void * ptr, size_t size, unsigned int root);

typedef void (*thread_reducer_t)(void *, const void *, size_t);

extern void thread_reducer_int_min(void *, const void *, size_t size);
extern void thread_reducer_int_max(void *, const void *, size_t size);
extern void thread_reducer_int_sum(void *, const void *, size_t size);
/* dptr must be a shared area allocated by shared_malloc */
/* the area in ptr is clobbered, but its output results are undefined */
extern void thread_allreduce(pi_wiring_ptr wr, void * dptr, void * ptr, size_t size, thread_reducer_t f);

/* shared_malloc is like malloc, except that the pointer returned will be
 * equal on all threads (proper access will deserve proper locking of
 * course). shared_malloc_set_zero sets to zero too */
/* As a side-effect, all shared_* functions serialize threads */
extern void * shared_malloc(pi_wiring_ptr wr, size_t size);
extern void * shared_malloc_set_zero(pi_wiring_ptr wr, size_t size);
extern void shared_free(pi_wiring_ptr wr, void * ptr);


/* This one is the higher level thingy on top of MPI_Bcast and the
 * previous one. It broadcast data computed by job j, thread t (e.g. 0,0)
 * to all jobs and threads.
 */
extern void global_broadcast(pi_wiring_ptr wr, void * ptr, size_t size, unsigned int j, unsigned int t);

/* companions to the above, the two functions below compare a data area
 * between threads and/or mpi jobs, and collectively return the result.
 * Useful for deciding on a common way to go given a condition.
 * 
 * These assume different pointers on all threads. If equal (created with
 * shared_malloc), then another function may be called -- currently not
 * exposed because I couldn't come up with a satisfying name.
 */
extern int thread_data_eq(parallelizing_info_ptr pi, void *buffer, size_t sz);
extern int global_data_eq(parallelizing_info_ptr pi, void *buffer, size_t sz);

/* prints the given string in a ascii-form matrix. */
extern void grid_print(parallelizing_info_ptr pi, char * buf, size_t siz, int print);

#define serialize(w)   serialize__(w, __FILE__, __LINE__)
extern int serialize__(pi_wiring_ptr, const char *, unsigned int);
#define serialize_threads(w)   serialize_threads__(w, __FILE__, __LINE__)
extern int serialize_threads__(pi_wiring_ptr, const char *, unsigned int);

extern int pi_save_file(pi_wiring_ptr, const char *, unsigned int iter, void *, size_t, size_t);
extern int pi_load_file(pi_wiring_ptr, const char *, unsigned int iter, void *, size_t, size_t);

extern int pi_load_file_2d(parallelizing_info_ptr, int, const char *, unsigned int iter, void *, size_t, size_t);
extern int pi_save_file_2d(parallelizing_info_ptr, int, const char *, unsigned int iter, void *, size_t, size_t);

/* stuff related to log entry printing */
extern void pi_log_init(pi_wiring_ptr);
extern void pi_log_clear(pi_wiring_ptr);
extern void pi_log_op(pi_wiring_ptr, const char * fmt, ...);
extern void pi_log_print_all(parallelizing_info_ptr);
extern void pi_log_print(pi_wiring_ptr);

/* These are the calls for interleaving. The 2n threads are divided into
 * two grous. It is guaranteed that at a given point, the two groups of n
 * threads are separated on either size of the pi_interleaving_flip call.
 *
 * The called function must make sure that alternating blocks (delimited
 * by _flip) either do or don't contain mpi calls, IN TURN.
 *
 * _enter and _leave are called from pi_go, so although they are exposed,
 * one does not have to know about them.
 */
extern void pi_interleaving_flip(parallelizing_info_ptr);
extern void pi_interleaving_enter(parallelizing_info_ptr);
extern void pi_interleaving_leave(parallelizing_info_ptr);

extern void pi_store_generic(parallelizing_info_ptr, unsigned long, unsigned long, void *);
extern void * pi_load_generic(parallelizing_info_ptr, unsigned long, unsigned long);

extern void pi_dictionary_init(pi_dictionary_ptr);
extern void pi_dictionary_clear(pi_dictionary_ptr);

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
