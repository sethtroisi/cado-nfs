#ifndef BARRIER_H_
#define BARRIER_H_

#ifdef	__cplusplus
extern "C" {
#endif

/* This is the interface to barrier synchronization waits.
 *
 * A barrier object must first be initialized with barrier_init. The
 * second argument specifies the number of threads expected to join. The
 * last argument is better specified as NULL (it is a means to tweak the
 * cancel state, which looks unwise).
 *
 * barrier_wait implements the wait. The function only returns when all
 * the threads have reached it. The second and third arguments provide
 * additional capabilities:
 *
 * If the second argument is non-null, the waiting time is returned in
 * the pointed area.
 * 
 * If the third argument is non null, the function specified is called
 * for each thread with a mutex held. The int argument passed to this
 * function in an ordinal representing the order in which threads reach
 * the barrier. First thread to reach calls the function with this
 * argument set to 0, and last to reach with the argument set to (n-1).
 * The return value of this function *for the last arriving thread* is
 * used as return value from barrier_wait.
 *
 * If the third argument is null, the last thread to reach the barrier
 * returns with return code 1. Others are expected to return with code 0.
 * A negative value indicates failure.
 *
 * notes:
 *
 * It is safe to call barrier_wait several times in a row ; that is, in a
 * loop, you don't have to intermix barrier waits on two different
 * barriers. Even if a fast thread reaches the next barrier while the
 * slow one hardly escaped the previous one, the barrier data will not be
 * messed up.
 *
 * Negative return values indicate errors (opposite of error codes).
 * Therefore it is wise, if you use the third argument to barrier_wait,
 * to allow distinction with such errors.
 */

#ifdef ENABLE_PTHREADS
typedef struct barrier_tag {
    pthread_mutex_t     mutex;          /* Control access to barrier */
    pthread_cond_t      cv;             /* wait for barrier */
    int                 valid;          /* set when valid */
    int                 threshold;      /* number of threads required */
    int                 counter;        /* current number of threads */
    unsigned long       cycle;          /* count cycles */
    int			result;		/* Result of the companion function */
    int			cancelstate;	/* PTHREAD_CANCEL_DISABLE by default */
} barrier_t;

#define BARRIER_VALID   0xdbcafe
#define BARRIER_INITIALIZER(cnt) \
    {PTHREAD_MUTEX_INITIALIZER, PTHREAD_COND_INITIALIZER, \
    BARRIER_VALID, cnt, cnt, 0}
#else	/* ENABLE_PTHREADS */
typedef void * barrier_t;
#define BARRIER_INITIALIZER(n)	NULL
#endif	/* ENABLE_PTHREADS */

extern int barrier_init (barrier_t *, int, int *);
extern int barrier_destroy (barrier_t *);
extern int barrier_wait (barrier_t *, double *, int (*)(int));

#ifdef	__cplusplus
}
#endif

#endif	/* BARRIER_H_ */
