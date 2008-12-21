#define _POSIX_C_SOURCE 200112L

#define _BSD_SOURCE
/* for strdup */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "select_mpi.h"
#include "parallelizing_info.h"
#include "macros.h"
#include "manu.h"

#include <sys/time.h>   // gettimeofday

static inline void pi_wiring_init_pthread_things(pi_wiring_ptr w, const char * desc)
{
    struct pthread_things * res;

    res = malloc(sizeof(struct pthread_things));

    my_pthread_barrier_init(res->b, NULL, w->ncores);
    my_pthread_mutex_init(res->m, NULL);
    res->desc = strdup(desc);

    w->th = res;
}

static inline void pi_wiring_destroy_pthread_things(pi_wiring_ptr w)
{
    my_pthread_barrier_destroy(w->th->b);
    my_pthread_mutex_destroy(w->th->m);
    /* Beware ! Freeing mustn't happen more than once ! */
    free(w->th->desc);
    free(w->th);
    w->th = NULL;
}

static void print_several(unsigned int n1, unsigned int n2, char a, char b, unsigned int w)
{
    char * toto;
    toto = malloc(w+2);
    printf("%c",b);
    memset(toto, a, w+1);
    toto[w+1]='\0';
    for(unsigned int i = 0 ; i < n1 ; i++) {
        toto[w]=a;
        for(unsigned int j = 0 ; j < n2-1 ; j++) {
            printf(toto);
        }
        toto[w]=b;
        printf(toto);
    }
    printf("\n");
}

struct pi_go_helper_s {
    void *(*fcn)(parallelizing_info_ptr, void*);
    parallelizing_info_ptr p;
    void * arg;
};

void * pi_go_helper_func(struct pi_go_helper_s * s)
{
    return (s->fcn)(s->p, s->arg);
}

typedef void * (*pthread_callee_t)(void*);

/*
void pi_errhandler(MPI_Comm * comm, int * err, ...)
{
    int size;
    int rank;
    MPI_Comm_size(*comm, &size);
    MPI_Comm_rank(*comm, &rank);
    fprintf(stderr, "Fatal MPI error ;\n");
    fprintf(stderr, " job %d in a comm. of size %d ;\n", rank, size);
    if (err == NULL) {
        fprintf(stderr, " no error code.\n");
    } else {
        char buf[MPI_MAX_ERROR_STRING];
        int len = sizeof(buf);
        MPI_Error_string(*err, buf, &len);
        fprintf(stderr, " %s\n", buf);
    }
    abort();
}
*/

void pi_go(void *(*fcn)(parallelizing_info_ptr, void *),
        unsigned int nhj, unsigned int nvj,
        unsigned int nhc, unsigned int nvc,
        void * arg)
{
    /* used in several places for doing snprintf */
    char buf[20];
    int err;

    parallelizing_info pi;
    memset(pi, 0, sizeof(pi));

    pi->m->njobs = nhj * nvj;
    pi->m->ncores = nhc * nvc;
    pi->m->totalsize = pi->m->njobs * pi->m->ncores;

    MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN);

    // MPI_Errhandler my_eh;
    // MPI_Comm_create_errhandler(&pi_errhandler, &my_eh);
    // MPI_Comm_set_errhandler(MPI_COMM_WORLD, my_eh);

    MPI_Comm_dup(MPI_COMM_WORLD, & pi->m->pals);
    MPI_Comm_rank(pi->m->pals, (int*) & pi->m->jrank);

    int size;
    MPI_Comm_size(pi->m->pals, & size);
    if ((unsigned int) size != pi->m->njobs) {
        if (pi->m->jrank == 0) {
            fprintf(stderr, "Inconsistency -- exactly %u == %u * %u"
                    " MPI jobs are needed -- got %u\n",
                    pi->m->njobs, nhj, nvj, size);
        }
        exit(1);
    }

    pi->m->jcommon = 0;       // This can't be anything other than zero.
    pi->m->trank = 0;

    /* Init to the simple case */
    pi->wr[0]->njobs = nvj;
    pi->wr[0]->ncores = nvc;
    pi->wr[1]->njobs = nhj;
    pi->wr[1]->ncores = nhc;
    pi->wr[0]->totalsize = pi->wr[0]->njobs * pi->wr[0]->ncores;
    pi->wr[1]->totalsize = pi->wr[1]->njobs * pi->wr[1]->ncores;

    // Here we want the _common row number_ ; not to be confused with the
    // rank !
    pi->wr[0]->jcommon = pi->m->jrank / nvj;
    pi->wr[1]->jcommon = pi->m->jrank % nvj;

    pi->wr[0]->trank = 0;
    pi->wr[1]->trank = 0;

    for(int d = 0 ; d < 2 ; d++) {
        // A subgroup contains all process having the same jcommon value.
        // Therefore, we have as many horizontal MPI barriers set up as
        // one finds MPI job numbers across a column.
        err = MPI_Comm_split(pi->m->pals, pi->wr[d]->jcommon,
                pi->m->jrank, &pi->wr[d]->pals);
        BUG_ON(err);
        MPI_Comm_rank(pi->wr[d]->pals, (int*) &pi->wr[d]->jrank);
    }

    ASSERT_ALWAYS((unsigned int) pi->m->jrank == pi->wr[0]->jcommon * nvj +  pi->wr[0]->jrank);
    ASSERT_ALWAYS(pi->wr[0]->jcommon == (unsigned int) pi->wr[1]->jrank);

    ASSERT_ALWAYS((unsigned int) pi->m->jrank == pi->wr[1]->jrank * nvj +  pi->wr[1]->jcommon);
    ASSERT_ALWAYS(pi->wr[1]->jcommon == (unsigned int) pi->wr[0]->jrank);

    // OK. So far we have something reasonable except for per-thread
    // setup. We have to first duplicate our structure so as to
    // specialize the pi things, and then set up agreed barriers.

    // the global barrier is not too much work.
    pi_wiring_init_pthread_things(pi->m, "main");

    // column and row barriers are more tricky.
    parallelizing_info * grid;
    grid = (parallelizing_info *) malloc(nhc * nvc * sizeof(parallelizing_info));
    memset(grid, 0, nhc * nvc * sizeof(parallelizing_info));


    for(unsigned int k = 0 ; k < nhc * nvc ; k++) {
        parallelizing_info_ptr e = grid[k];
        memcpy(e, pi, sizeof(parallelizing_info));
        unsigned int k0 = k % nvc;
        unsigned int k1 = k / nvc;
        e->m->trank = k;
        e->wr[0]->trank = k0;
        e->wr[1]->trank = k1;
#ifdef  MPI_LIBRARY_MT_CAPABLE
        MPI_Comm_dup(pi->m->pals, &(e->m->pals));
        MPI_Comm_dup(pi->wr[0]->pals, &(e->wr[0]->pals));
        MPI_Comm_dup(pi->wr[1]->pals, &(e->wr[1]->pals));
#endif  /* MPI_LIBRARY_MT_CAPABLE */
    }

    // row barriers.
    // we've got pi->wr[1]->ncores cores working on separate rows (that's
    // the number of cores we encounter when walking down one column). So
    // each row leader is at index c * pi->wr[0]->ncores, because
    // pi->wr[0]->ncores is the number of cores working on separate
    // columns.
    //
    for(unsigned int c = 0 ; c < pi->wr[1]->ncores ; c++) {
        snprintf(buf, sizeof(buf), "r%u", 
                pi->wr[1]->jrank * pi->wr[1]->ncores + c);
        unsigned int Nc = c * pi->wr[0]->ncores;
        pi_wiring_ptr leader = grid[Nc]->wr[0];
        pi_wiring_init_pthread_things(leader, buf);
        // replicate.
        for(unsigned int k = 0 ; k < pi->wr[0]->ncores ; k++) {
            grid[k + Nc]->wr[0]->th = leader->th;
            grid[k + Nc]->wr[0]->tcommon = c;
        }
    }

    // column barriers.
    for(unsigned int c = 0 ; c < pi->wr[0]->ncores ; c++) {
        snprintf(buf, sizeof(buf), "c%u", 
                pi->wr[0]->jrank * pi->wr[0]->ncores + c);
        pi_wiring_ptr leader = grid[c]->wr[1];
        pi_wiring_init_pthread_things(leader, buf);
        // replicate.
        for(unsigned int k = 0 ; k < pi->wr[1]->ncores ; k++) {
            unsigned int Nk = k * pi->wr[0]->ncores;
            grid[c + Nk]->wr[1]->th = leader->th;
            grid[c + Nk]->wr[1]->tcommon = c;
        }
    }

    // Now do some printing, for debugging.
    for(unsigned int ij = 0 ; ij < pi->wr[1]->njobs ; ij++) {
        // not mt at the moment, so it's easy.
        err = MPI_Barrier(pi->m->pals);
        BUG_ON(err);

        if (pi->wr[1]->jrank != ij) continue;

        if (pi->m->jrank == 0)
            print_several(pi->wr[0]->njobs,pi->wr[0]->ncores, '-', '+', 12);

        for(unsigned int it = 0 ; it < pi->wr[1]->ncores ; it++) {
            if (pi->wr[0]->jrank == 0) printf("|");
            for(unsigned int jj = 0 ; jj < pi->wr[0]->njobs ; jj++) {
                err = MPI_Barrier(pi->wr[0]->pals);
                BUG_ON(err);
                for(unsigned int jt = 0 ; jt < pi->wr[0]->ncores ; jt++) {
                    parallelizing_info_srcptr tpi;
                    tpi = grid[it*pi->wr[0]->ncores+jt];
                    snprintf(buf, sizeof(buf), "J%uT%u%s%s",
                            tpi->m->jrank,
                            tpi->m->trank,
                            tpi->wr[0]->th->desc,
                            tpi->wr[1]->th->desc);
                    err = MPI_Bcast(buf, sizeof(buf),
                            MPI_BYTE, jj, pi->wr[0]->pals);
                    BUG_ON(err);
                    int fence = jt == pi->wr[0]->ncores - 1;
                    if (pi->wr[0]->jrank == 0)
                        printf("%-12s%c", buf, fence ? '|' : ' ');
                }
            }
            if (pi->wr[0]->jrank == 0) printf("\n");
        }
        if (pi->wr[1]->jrank == ij && pi->wr[0]->jrank == 0)
            print_several(pi->wr[0]->njobs,pi->wr[0]->ncores, '-', '+', 12);
    }

    err = MPI_Barrier(pi->m->pals);
    BUG_ON(err);
    if (pi->m->jrank == pi->m->njobs - 1)
        printf("going multithread now\n");
    err = MPI_Barrier(pi->m->pals);
    BUG_ON(err);

    // go mt.
    my_pthread_t * tgrid;
    tgrid = (my_pthread_t *) malloc(nhc * nvc * sizeof(my_pthread_t));

    struct pi_go_helper_s * helper_grid;
    helper_grid = malloc(nhc * nvc * sizeof(struct pi_go_helper_s));

    for(unsigned int k = 0 ; k < nhc * nvc ; k++) {
        helper_grid[k].fcn = fcn;
        helper_grid[k].arg = arg;
        helper_grid[k].p = (parallelizing_info_ptr) grid + k;
#if 0
        printf("Job %u Thread %u"
                " ROW (%u,%u) Job %u Thread %u"
                " COL (%u,%u) Job %u Thread %u\n",
                grid[k]->m->jrank, grid[k]->m->trank,
                grid[k]->wr[0]->jcommon,
                grid[k]->wr[0]->tcommon,
                grid[k]->wr[0]->jrank,
                grid[k]->wr[0]->trank,
                grid[k]->wr[1]->jcommon,
                grid[k]->wr[1]->tcommon,
                grid[k]->wr[1]->jrank,
                grid[k]->wr[1]->trank);
#endif
        my_pthread_create(tgrid+k, NULL,
                (pthread_callee_t) pi_go_helper_func, helper_grid+k);
#if 0
        printf("Thread (%u,%u)/(%u,%u) started\n",
                grid[k]->wr[0]->jrank, grid[k]->wr[1]->jrank,
                grid[k]->wr[0]->trank, grid[k]->wr[1]->trank);
#endif
    }

#if 0
    if (pi->m->jrank == 0) {
        printf("Done starting threads\n");
    }
#endif
    // nothing done here. We're the main job.
    
    for(unsigned int k = 0 ; k < nhc * nvc ; k++) {
        my_pthread_join(tgrid[k], NULL);
#if 0
        printf("Thread (%u,%u) finished\n", grid[k]->wr[0]->trank, grid[k]->wr[1]->trank);
#endif
    }
    free(helper_grid);
    free(tgrid);
    // destroy row barriers.
    for(unsigned int c = 0 ; c < pi->wr[1]->ncores ; c++) {
        unsigned int Nc = c * pi->wr[0]->ncores;
        pi_wiring_destroy_pthread_things(grid[Nc]->wr[0]);
    }
    // destroy column barriers
    for(unsigned int c = 0 ; c < pi->wr[0]->ncores ; c++) {
        pi_wiring_destroy_pthread_things(grid[c]->wr[1]);
    }
#ifdef  MPI_LIBRARY_MT_CAPABLE
    for(unsigned int k = 0 ; k < nhc * nvc ; k++) {
        MPI_Comm_free(&grid[k]->m->pals);
        MPI_Comm_free(&grid[k]->wr[0]->pals);
        MPI_Comm_free(&grid[k]->wr[1]->pals);
    }
#endif  /* MPI_LIBRARY_MT_CAPABLE */

    free(grid);

    pi_wiring_destroy_pthread_things(pi->m);

    for(int d = 0 ; d < 2 ; d++) {
        MPI_Comm_free(&pi->wr[d]->pals);
    }
    MPI_Comm_free(&pi->m->pals);

    // MPI_Errhandler_free(&my_eh);
}

void thread_agreement(pi_wiring_ptr wr, void ** ptr, unsigned int i)
{
    // FIXME: this design is flawed. Correct reentrancy is not achieved.
    // All the stuff around the barrier here should be done with the
    // mutex locked, otherwise disaster may occur.
    // A fix would incur either basically rewriting barrier code, so as
    // to place stuff at the right spot. Note that there _is_ some
    // barrier code in bw/barrier.[ch] which could easily be expanded for
    // this purpose.
    if (wr->trank == i) {
        wr->th->utility_ptr = *ptr;
    }
    my_pthread_barrier_wait(wr->th->b);
    *ptr = wr->th->utility_ptr;
}

void complete_broadcast(pi_wiring_ptr wr, void * ptr, size_t size, unsigned int j, unsigned int t)
{
    int err;

    ASSERT(j < wr->njobs);
    ASSERT(t < wr->ncores);
    if (wr->trank == t) {
        err = MPI_Bcast(ptr, size, MPI_BYTE, j, wr->pals);
        BUG_ON(err);
    }
    void * leader_ptr = ptr;
    thread_agreement(wr, &leader_ptr, 0);
    if (ptr != leader_ptr)
        memcpy(ptr, leader_ptr, size);
}


int serialize__(pi_wiring_ptr w, const char * s MAYBE_UNUSED, unsigned int l MAYBE_UNUSED)
{
    int err;
#ifdef  CONCURRENCY_DEBUG
    // note how w->th_count is normally a thread-private thing !
    w->th_count++;
    printf("(%s:%u) barrier #%d, %u/%u on %s [%p]\n", s, l, w->th_count,
            w->trank, w->ncores,
            w->th->desc, w->th->b);
#endif
    my_pthread_barrier_wait(w->th->b);
    if (w->trank == 0) {
        err = MPI_Barrier(w->pals);
        BUG_ON(err);
    }
    // struct timeval tv[1];
    // gettimeofday(tv, NULL);
    // printf("%.2f\n", tv->tv_sec + (double) tv->tv_usec / 1.0e6);
    // sleep(1);
    return w->jrank == 0 && w->trank == 0;
}

int serialize_threads__(pi_wiring_ptr w, const char * s MAYBE_UNUSED, unsigned int l MAYBE_UNUSED)
{
#ifdef  CONCURRENCY_DEBUG
    // note how w->th_count is normally a thread-private thing !
    w->th_count++;
    printf("(%s:%u) tbarrier #%d, %u/%u on %s [%p]\n", s, l, w->th_count,
            w->trank, w->ncores,
            w->th->desc, w->th->b);
#endif
    my_pthread_barrier_wait(w->th->b);
    // struct timeval tv[1];
    // gettimeofday(tv, NULL);
    // printf("%.2f\n", tv->tv_sec + (double) tv->tv_usec / 1.0e6);
    // sleep(1);
    return w->jrank == 0 && w->trank == 0;
}

// XXX Many serializing calls in the two routines below have no
// functional use, apart from the fact that they exert some pressure on
// the mpi+pthreads way of doing things.

void say_hello(pi_wiring_ptr w, parallelizing_info_ptr pi)
{
    serialize(w);
    for(unsigned int j = 0 ; j < w->njobs ; j++) {
        serialize(w);
        if ((unsigned int) w->jrank != j)
            continue;
        my_pthread_mutex_lock(w->th->m);
        printf("(%s) J%uT%u ; %s%s ; (%s:j%ut%u) (%s:j%ut%u)\n",
                w->th->desc,
                pi->m->jrank,
                pi->m->trank,
                pi->wr[0]->th->desc,
                pi->wr[1]->th->desc,
                pi->wr[0]->th->desc, pi->wr[0]->jrank, pi->wr[0]->trank,
                pi->wr[1]->th->desc, pi->wr[1]->jrank, pi->wr[1]->trank
              );
        my_pthread_mutex_unlock(w->th->m);
    }
}

void hello(parallelizing_info_ptr pi)
{
    if (serialize(pi->m)) {
        printf("Doing hello world loop\n");
    }
    say_hello(pi->m, pi);

    // before changing the grain of the barriers, we have to serialize at
    // the thread level (i.e. wait for the mpi serialization calls to
    // complete). Otherwise, inter-job communication occuring within
    // say_hello might clash with the inter-job calls done within
    // serialize().
    // Rule is: never let mpi-level synchronisations on a communicator
    // wander outside an execution period for which all threads sharing
    // this communicator are doing the same.
    for(unsigned int i = 0 ; i < pi->wr[0]->totalsize ; i++) {
        serialize(pi->m);
        serialize_threads(pi->m);
        if (i == pi->wr[0]->jrank * pi->wr[0]->ncores + pi->wr[0]->trank) {
            say_hello(pi->wr[1], pi);
        }
    }

    for(unsigned int i = 0 ; i < pi->wr[1]->totalsize ; i++) {
        serialize(pi->m);
        serialize_threads(pi->m);
        if (i == pi->wr[1]->jrank * pi->wr[1]->ncores + pi->wr[1]->trank) {
            say_hello(pi->wr[0], pi);
        }
    }

    if (serialize(pi->m)) {
        printf("OK: Finished hello world loop\n");
    }
}

