#ifndef WORKER_THREADS_H_
#define WORKER_THREADS_H_

#include <pthread.h>

#ifdef __cplusplus
extern "C" {
#endif

struct worker_threads_group;

typedef void (*worker_func_t)(struct worker_threads_group *, int, void *);

struct worker_thread_info {
    int i;
    struct worker_threads_group * tg;
};

struct worker_threads_group {
    /* The mutex can be freely used by callees */
    pthread_mutex_t mu;

    /* Total number of threads */
    unsigned int n;

    /* Number of worker threads at work (not idle) */
    unsigned int working;

    /* Number of threads which have completed their current task */
    unsigned int done;

    /* practically everything else is meant to be private and/or useless
     * anyway to callees */
    pthread_cond_t cond_in;
    pthread_cond_t cond_out;

    struct worker_thread_info * thread_table;
    pthread_t * threads;

    worker_func_t f;
    void * arg;
};


extern void worker_threads_do(struct worker_threads_group * tg, worker_func_t f, void * arg);
extern struct worker_threads_group * worker_threads_init(unsigned int n);
extern void worker_threads_clear(struct worker_threads_group * tg);

#ifdef __cplusplus
}
#endif

#endif	/* WORKER_THREADS_H_ */
