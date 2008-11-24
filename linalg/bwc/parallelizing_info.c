#define _POSIX_C_SOURCE 200112L

#define _BSD_SOURCE
/* for strdup */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "select_mpi.h"
#include "parallelizing_info.h"
#include "macros.h"

#include <sys/time.h>   // gettimeofday

#define xxxCONCURRENCY_DEBUG

static inline void pi_wiring_init_pthread_things(struct pi_wiring * w, const char * desc)
{
    struct pthread_things * res;

    res = malloc(sizeof(struct pthread_things));

    my_pthread_barrier_init(res->b, NULL, w->ncores);
    my_pthread_mutex_init(res->m, NULL);
    res->desc = strdup(desc);

    w->th = res;
}

static inline void pi_wiring_destroy_pthread_things(struct pi_wiring * w)
{
    my_pthread_barrier_destroy(w->th->b);
    my_pthread_mutex_destroy(w->th->m);
    /* Beware ! Freeing mustn't happen more than once ! */
    free(w->th->desc);
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

typedef void * (*pthread_callee_t)(void*);

void pi_go(void *(*fcn)(parallelizing_info_ptr),
        unsigned int nhj, unsigned int nvj,
        unsigned int nhc, unsigned int nvc)
{
    /* used in several places for doing snprintf */
    char buf[20];

    parallelizing_info pi;
    memset(pi, 0, sizeof(pi));

    pi->m->njobs = nhj * nvj;
    pi->m->ncores = nhc * nvc;
    MPI_Comm_dup(MPI_COMM_WORLD, & pi->m->pals);
    MPI_Comm_rank(pi->m->pals, (int*) & pi->m->jrank);

    int size;
    MPI_Comm_size(MPI_COMM_WORLD, & size);
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

    // Here we want the _common row number_ ; not to be confused with the
    // rank !
    pi->wr[0]->jcommon = pi->m->jrank / nvj;
    pi->wr[1]->jcommon = pi->m->jrank % nvj;

    pi->wr[0]->trank = 0;
    pi->wr[1]->trank = 0;

    for(int d = 0 ; d < 2 ; d++) {
        MPI_Comm_split(pi->m->pals, pi->wr[d]->jcommon,
                pi->m->jrank, &pi->wr[d]->pals);
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
        memcpy(grid+k, pi, sizeof(parallelizing_info));
        grid[k]->m->trank = k;
        grid[k]->wr[0]->trank = k % nvc;
        grid[k]->wr[1]->trank = k / nvc;
    }

    // row barriers.
    // we've got pi->wr[1]->ncores working on separate rows (that's the
    // number of cores we encounter when walking down one column). So
    // each row leader is at index c * pi->wr[0]->ncores, because
    // pi->wr[0]->ncores is the number of cores working on separate
    // columns.
    //
    for(unsigned int c = 0 ; c < pi->wr[1]->ncores ; c++) {
        snprintf(buf, sizeof(buf), "r%u", 
                pi->wr[1]->jrank * pi->wr[1]->ncores + c);
        unsigned int Nc = c * pi->wr[0]->ncores;
        struct pi_wiring * leader = grid[Nc]->wr[0];
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
        struct pi_wiring * leader = grid[c]->wr[1];
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
        MPI_Barrier(pi->m->pals);

        if (pi->wr[1]->jrank != ij) continue;

        if (pi->m->jrank == 0)
            print_several(pi->wr[0]->njobs,pi->wr[0]->ncores, '-', '+', 12);

        for(unsigned int it = 0 ; it < pi->wr[1]->ncores ; it++) {
            if (pi->wr[0]->jrank == 0) printf("|");
            for(unsigned int jj = 0 ; jj < pi->wr[0]->njobs ; jj++) {
                MPI_Barrier(pi->wr[0]->pals);
                for(unsigned int jt = 0 ; jt < pi->wr[0]->ncores ; jt++) {
                    parallelizing_info_srcptr tpi;
                    tpi = grid[it*pi->wr[0]->ncores+jt];
                    snprintf(buf, sizeof(buf), "J%uT%u%s%s",
                            tpi->m->jrank,
                            tpi->m->trank,
                            tpi->wr[0]->th->desc,
                            tpi->wr[1]->th->desc);
                    MPI_Bcast(buf, sizeof(buf), MPI_BYTE, jj, pi->wr[0]->pals);
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

    MPI_Barrier(pi->m->pals);
    if (pi->m->jrank == pi->m->njobs - 1)
        printf("going multithread now\n");
    MPI_Barrier(pi->m->pals);

    // go mt.
    my_pthread_t * tgrid;
    tgrid = (my_pthread_t *) malloc(nhc * nvc * sizeof(my_pthread_t));
    for(unsigned int k = 0 ; k < nhc * nvc ; k++) {
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
        my_pthread_create(tgrid+k, NULL, (pthread_callee_t) fcn, grid+k);
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

    free(grid);

    pi_wiring_destroy_pthread_things(pi->m);

    for(int d = 0 ; d < 2 ; d++) {
        MPI_Comm_free(&pi->wr[d]->pals);
    }
}

#define serialize(w)   serialize__(w, __LINE__)
int serialize__(struct pi_wiring * w, unsigned int l MAYBE_UNUSED)
{
#ifdef  CONCURRENCY_DEBUG
    // note how w->th_count is normally a thread-private thing !
    w->th_count++;
    printf("(line %u) barrier #%d, %u/%u on %s [%p]\n", l, w->th_count,
            w->trank, w->ncores,
            w->th->desc, w->th->b);
#endif
    if (w->trank == 0) {
        MPI_Barrier(w->pals);
    }
    my_pthread_barrier_wait(w->th->b);
    // struct timeval tv[1];
    // gettimeofday(tv, NULL);
    // printf("%.2f\n", tv->tv_sec + (double) tv->tv_usec / 1.0e6);
    // sleep(1);
    return w->jrank == 0 && w->trank == 0;
}

#define serialize_threads(w)   serialize_threads__(w, __LINE__)
int serialize_threads__(struct pi_wiring * w, unsigned int l MAYBE_UNUSED)
{
#ifdef  CONCURRENCY_DEBUG
    // note how w->th_count is normally a thread-private thing !
    w->th_count++;
    printf("(line %u) tbarrier #%d, %u/%u on %s [%p]\n", l, w->th_count,
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

void say_hello(struct pi_wiring * w, parallelizing_info_ptr pi)
{
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
    serialize(w);       //
}

void hello(parallelizing_info_ptr pi)
{
    if (serialize(pi->m)) {
        printf("Doing hello world loop\n");
    }
    say_hello(pi->m, pi);

    serialize(pi->m);   //
    for(unsigned int i = 0 ; i < pi->wr[0]->njobs * pi->wr[0]->ncores ; i++) {
        serialize(pi->m);       //
        if (i == pi->wr[0]->jrank * pi->wr[0]->ncores + pi->wr[0]->trank) {
            say_hello(pi->wr[1], pi);
        }
        // XXX we need to serialize __threads__ here, before attempting
        // anything at the job level ! Otherwise, me might have an
        // inter-job communication occuring within say_hello above, which
        // could clash with the main serialization below...
        serialize_threads(pi->m);
    }

    serialize(pi->m);   //
    for(unsigned int i = 0 ; i < pi->wr[1]->njobs * pi->wr[1]->ncores ; i++) {
        serialize(pi->m);       //
        if (i == pi->wr[1]->jrank * pi->wr[1]->ncores + pi->wr[1]->trank) {
            say_hello(pi->wr[0], pi);
        }
        serialize_threads(pi->m);
    }

    serialize(pi->m);   //
}

