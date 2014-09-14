#include "cado.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

// we're doing open close mmap truncate...
#ifdef HAVE_MMAN_H
#include <sys/mman.h>
#endif
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <errno.h>
#include <stdarg.h>
#include <stdint.h>
#include <inttypes.h>
#include <ctype.h>
#include "bwc_config.h"
#include "select_mpi.h"
#include "parallelizing_info.h"
#include "portability.h"
#include "macros.h"
#include "misc.h"
#include "memory.h"

#include <sys/time.h>   // gettimeofday
#ifdef  HAVE_UTSNAME_H
#include <sys/utsname.h>
#endif

#if defined(HAVE_HWLOC) && defined(HAVE_CXX11)
#include "cpubinding.h"
#endif  /* defined(HAVE_HWLOC) && defined(HAVE_CXX11) */

static inline void pi_wiring_init_pthread_things(pi_wiring_ptr w, const char * desc)
{
    struct pthread_things * res;

    res = malloc(sizeof(struct pthread_things));

    barrier_init(res->bh, w->ncores);
    my_pthread_barrier_init(res->b, NULL, w->ncores);
    my_pthread_mutex_init(res->m, NULL);
    res->desc = strdup(desc);

    w->th = res;
}

static inline void pi_wiring_destroy_pthread_things(pi_wiring_ptr w)
{
    barrier_destroy(w->th->bh);
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
            printf("%s", toto);
        }
        toto[w]=b;
        printf("%s", toto);
    }
    printf("\n");
    free(toto);
}

struct pi_go_helper_s {
    void *(*fcn)(parallelizing_info_ptr, param_list pl, void*);
    parallelizing_info_ptr p;
    param_list_ptr pl;
    void * arg;
};

void * pi_go_helper_func(struct pi_go_helper_s * s)
{
    pi_interleaving_enter(s->p);
#if defined(HAVE_HWLOC) && defined(HAVE_CXX11)
    cpubinding_do_pinning(s->p->cpubinding_info, s->p->wr[0]->trank, s->p->wr[1]->trank);
#endif /* defined(HAVE_HWLOC) && defined(HAVE_CXX11) */
    void * ret = (s->fcn)(s->p, s->pl, s->arg);
    pi_interleaving_flip(s->p);
    pi_interleaving_leave(s->p);
    return ret;
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

/* TODO: There would be a fairly easy generalization to that, and quite
 * handy. Modify the meaning of the ``siz'' argument to mean something
 * which doesn't have to be equal on all threads. Instead, pick the
 * maximum and use it for formatting.
 */

void grid_print(parallelizing_info_ptr pi, char * buf, size_t siz, int print)
{
    pi_wiring_ptr wr = pi->m;

    char * strings;
    strings = malloc(wr->njobs * wr->ncores * siz);

    int me = wr->jrank * wr->ncores + wr->trank;

    /* instead of doing memcpy, we align the stuff. */
    char * fmt;
    int rc = asprintf(&fmt, "%%-%zus", siz-1);
    ASSERT_ALWAYS(rc >= 0);
    snprintf(strings + me * siz, siz, fmt, buf);
    free(fmt);

    /* ceinture et bretelles */

    serialize(wr);
    for(unsigned int j = 0 ; j < wr->njobs ; j++) {
        for(unsigned int t = 0 ; t < wr->ncores ; t++) {
            global_broadcast(wr, strings + (j * wr->ncores + t) * siz,
                    siz, j, t);
        }
    }
    serialize(wr);

    char * ptr = strings;

    if (print) {
        /* There's also the null byte counted in siz. So that makes one
         * less.
         */
        unsigned int nj0 = pi->wr[0]->njobs;
        unsigned int nj1 = pi->wr[1]->njobs;
        unsigned int nt0 = pi->wr[0]->ncores;
        unsigned int nt1 = pi->wr[1]->ncores;
        print_several(nj0, nt0, '-', '+', siz-1);
        for(unsigned int j1 = 0 ; j1 < nj1 ; j1++) {
            for(unsigned int t1 = 0 ; t1 < nt1 ; t1++) {
                printf("|");
                for(unsigned int j0 = 0 ; j0 < nj0 ; j0++) {
                    for(unsigned int t0 = 0 ; t0 < nt0 ; t0++) {
                        int fence = t0 == nt0 - 1;
                        printf("%s%c", ptr, fence ? '|' : ' ');
                        ptr += siz;
                    }
                }
                printf("\n");
            }
            print_several(nj0, nt0, '-', '+', siz-1);
        }
    }
    serialize(wr);
    free(strings);
}

static void get_node_number_and_prefix(parallelizing_info_ptr pi)
{
    int len = strlen(pi->nodename);
    int minlen = len;
    int maxlen = len;
    MPI_Allreduce(MPI_IN_PLACE, &minlen, 1, MPI_INT, MPI_MIN, pi->m->pals);
    MPI_Allreduce(MPI_IN_PLACE, &maxlen, 1, MPI_INT, MPI_MAX, pi->m->pals);
    int j;
    for(j = 0 ; j < minlen ; j++) {
        int c = pi->nodename[j];
        if (isdigit(c)) {
            c = -1;
        }
        int cmin = c;
        int cmax = c;
        MPI_Allreduce(MPI_IN_PLACE, &cmin, 1, MPI_INT, MPI_MIN, pi->m->pals);
        MPI_Allreduce(MPI_IN_PLACE, &cmax, 1, MPI_INT, MPI_MAX, pi->m->pals);
        if (cmin != cmax || cmin < 0)
            break;
    }
    memcpy(pi->nodeprefix, pi->nodename, j);
    pi->nodeprefix[j]='\0';
    if (j && pi->nodeprefix[j-1] == '-')
        pi->nodeprefix[j-1] = '\0';
    memset(pi->nodenumber_s, ' ', maxlen - j);
    pi->nodenumber_s[maxlen-j]='\0';
    memcpy(pi->nodenumber_s + maxlen-len, pi->nodename + j, len-j);
    if (j == maxlen) {
        // all nodes have the same name. Most probably because we're only
        // on one node. Arrange for the display to look vaguely right.
        memcpy(pi->nodeprefix, pi->nodename, PI_NAMELEN);
        memcpy(pi->nodenumber_s, pi->nodename, PI_NAMELEN);
    }
}

static void display_process_grid(parallelizing_info_ptr pi)
{
    char * all_node_ids = malloc(PI_NAMELEN * pi->m->njobs);
    char * all_node_ids2 = malloc(PI_NAMELEN * pi->m->njobs);
    typedef int cpair[2];
    cpair my_coords = { pi->wr[0]->jrank, pi->wr[1]->jrank, };
    cpair * all_cpairs = malloc(sizeof(cpair) * pi->m->njobs);

    MPI_Allgather(pi->nodenumber_s, PI_NAMELEN, MPI_BYTE,
            all_node_ids, PI_NAMELEN, MPI_BYTE, pi->m->pals);
    MPI_Allgather(my_coords, 2, MPI_INT, all_cpairs, 2, MPI_INT, pi->m->pals);

    for(unsigned int i = 0 ; i < pi->m->njobs ; i++) {
        int x = all_cpairs[i][0];
        int y = all_cpairs[i][1];
        unsigned int j = y * pi->wr[0]->njobs + x;
        if (j >= pi->m->njobs) {
            fprintf(stderr, "uh ?\n");
            MPI_Abort(pi->m->pals, 42);
        }
        memcpy(all_node_ids2 + j * PI_NAMELEN,
                all_node_ids + i * PI_NAMELEN,
                PI_NAMELEN);
    }

    if (pi->m->jrank == 0) {
        if (strlen(pi->nodeprefix)) {
            printf("%d nodes on %s\n", pi->m->njobs, pi->nodeprefix);
        } else {
            printf("%d nodes, no common name prefix\n", pi->m->njobs);
        }
        char * pad = malloc(PI_NAMELEN);
        memset(pad, ' ', PI_NAMELEN);
        int node_id_len = strlen(pi->nodenumber_s);
        pad[node_id_len]='\0';
        if (node_id_len) pad[node_id_len/2]='.';
        for(unsigned int i = 0 ; i < pi->wr[1]->njobs ; i++) {       // i == y
        for(unsigned int it = 0 ; it < pi->wr[1]->ncores ; it++) {
            for(unsigned int j = 0 ; j < pi->wr[0]->njobs ; j++) {       // j == x
            for(unsigned int jt = 0 ; jt < pi->wr[0]->ncores ; jt++) {       // jt == x
                char * what = pad;
                if (!it && !jt) {
                    what = all_node_ids2 + PI_NAMELEN * (i * pi->wr[0]->njobs + j); 
                }
                printf(" %s", what);
            }
            }
            printf("\n");
        }
        }
        free(pad);
    }
    free(all_node_ids2);
    free(all_node_ids);
    free(all_cpairs);
    MPI_Barrier(pi->m->pals);
}

static void pi_init_mpilevel(parallelizing_info_ptr pi, param_list pl)
{
    int err;
    int mpi[2] = { 1, 1, };
    int thr[2] = { 1, 1, };

    param_list_parse_intxint(pl, "mpi", mpi);
    param_list_parse_intxint(pl, "thr", thr);

    unsigned int nhj = mpi[0];
    unsigned int nvj = mpi[1];
    unsigned int nhc = thr[0];
    unsigned int nvc = thr[1];

    memset(pi, 0, sizeof(parallelizing_info));

#ifdef HAVE_UTSNAME_H
    {
        struct utsname u[1];
        uname(u);
        memcpy(pi->nodename, u->nodename, MIN(PI_NAMELEN, sizeof(u->nodename)));
        pi->nodename[PI_NAMELEN-1]='\0';
        char * p = strchr(pi->nodename, '.');
        if (p) *p='\0';
    }
#else
    strncpy(pi->nodename, "unknown", sizeof(pi->nodename));
#endif
    pi->m->njobs = nhj * nvj;
    pi->m->ncores = nhc * nvc;
    pi->m->totalsize = pi->m->njobs * pi->m->ncores;

#ifndef NDEBUG
    /* Must make sure that we have proper ASSERTS after MPI calls then. */
    MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN);
#endif

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
    unsigned int jcommon[2];
    jcommon[0] = pi->m->jrank / nvj;
    jcommon[1] = pi->m->jrank % nvj;

    pi->wr[0]->trank = 0;
    pi->wr[1]->trank = 0;

    for(int d = 0 ; d < 2 ; d++) {
        // A subgroup contains all process having the same jcommon value.
        // Therefore, we have as many horizontal MPI barriers set up as
        // one finds MPI job numbers across a column.
        err = MPI_Comm_split(pi->m->pals, jcommon[d],
                pi->m->jrank, &pi->wr[d]->pals);
        ASSERT_ALWAYS(!err);
        MPI_Comm_rank(pi->wr[d]->pals, (int*) &pi->wr[d]->jrank);
    }

    get_node_number_and_prefix(pi);
    display_process_grid(pi);

#if defined(HAVE_HWLOC) && defined(HAVE_CXX11)
    /* prepare the cpu binding messages, and print the unique messages we
     * receive */
    char * cpubinding_messages;
    pi->cpubinding_info = cpubinding_get_info(&cpubinding_messages, pl, thr);
    int msgsize = 0;
    if (cpubinding_messages)
        msgsize = strlen(cpubinding_messages);
    MPI_Allreduce(MPI_IN_PLACE, &msgsize, 1, MPI_INT, MPI_MAX, pi->m->pals);
    if (msgsize == 0) {
        if (cpubinding_messages)
            free(cpubinding_messages);
    }
    msgsize++;
    int chunksize = PI_NAMELEN + msgsize;
    char * big_pool = malloc(pi->m->njobs * chunksize);
    memset(big_pool, 0, pi->m->njobs * chunksize);
    if (cpubinding_messages) {
        int rc;
        rc = strlcpy(big_pool + pi->m->jrank * chunksize, pi->nodename, PI_NAMELEN);
        ASSERT_ALWAYS(rc == (int) strlen(pi->nodename));
        rc = strlcpy(big_pool + pi->m->jrank * chunksize + PI_NAMELEN, cpubinding_messages, msgsize);
        ASSERT_ALWAYS(rc == (int) strlen(cpubinding_messages));
        free(cpubinding_messages);
    }

    MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
            big_pool, chunksize, MPI_BYTE, pi->m->pals);

    if (pi->m->jrank == 0) {
        const char * refnode = NULL;
        const char * ref = NULL;
        for(unsigned int i = 0 ; i < pi->m->njobs ; i++) {
            const char * node = big_pool + i * chunksize;
            const char * msg = node + PI_NAMELEN;
            if (!ref || strcmp(msg, ref) != 0) {
                printf("cpubinding messages on %s:\n%s", node, msg);
                ref = msg;
                refnode = node;
            } else {
                printf("cpubinding messages on %s: same as %s\n", node, refnode);
            }
        }
    }
    free(big_pool);
#endif /* defined(HAVE_HWLOC) && defined(HAVE_CXX11) */
}

    static parallelizing_info *
pi_grid_init(parallelizing_info_ptr pi)
{
    // unsigned int nvj = pi->wr[0]->njobs;
    // unsigned int nhj = pi->wr[1]->njobs;
    unsigned int nvc = pi->wr[0]->ncores;
    unsigned int nhc = pi->wr[1]->ncores;

    /* used in several places for doing snprintf */
    char buf[20];

    // OK. So far we have something reasonable except for per-thread
    // setup. We have to first duplicate our structure so as to
    // specialize the pi things, and then set up agreed barriers.

    // the global barrier is not too much work.
    pi_wiring_init_pthread_things(pi->m, "main");

    // column and row barriers are more tricky.
    parallelizing_info * grid;
    grid = (parallelizing_info *) malloc(nhc * nvc * sizeof(parallelizing_info));
    // setting to zero is important, since several fields are ensured as
    // null.
    memset(grid, 0, nhc * nvc * sizeof(parallelizing_info));

    char commname[64];
#ifndef MPI_LIBRARY_MT_CAPABLE
    MPI_Comm_set_name(pi->m->pals, "main");
    snprintf(commname, sizeof(commname), "rows%u-%u", 
            pi->wr[1]->jrank * pi->wr[1]->ncores,
            (pi->wr[1]->jrank + 1) * pi->wr[1]->ncores - 1);
    MPI_Comm_set_name(pi->wr[0]->pals, commname);
    snprintf(commname, sizeof(commname), "cols%u-%u", 
            pi->wr[0]->jrank * pi->wr[0]->ncores,
            (pi->wr[0]->jrank + 1) * pi->wr[0]->ncores - 1);
    MPI_Comm_set_name(pi->wr[1]->pals, commname);
#endif

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

        /* give a name */
        snprintf(commname, sizeof(commname), "main.%u", k);
        MPI_Comm_set_name(e->m->pals, commname);
        snprintf(commname, sizeof(commname), "row%u.%u", 
                e->wr[1]->jrank * e->wr[1]->ncores + e->wr[1]->trank,
                e->wr[0]->trank);
        MPI_Comm_set_name(e->wr[0]->pals, commname);
        snprintf(commname, sizeof(commname), "col%u.%u", 
                e->wr[0]->jrank * e->wr[0]->ncores + e->wr[0]->trank,
                e->wr[1]->trank);
        MPI_Comm_set_name(e->wr[1]->pals, commname);
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
        }
    }

    return grid;
}

#if 0
/* TODO: rewrite this using grid_print instead */
static void
pi_grid_print_sketch(parallelizing_info_ptr pi, parallelizing_info * grid)
{
    char buf[20];
    for(unsigned int ij = 0 ; ij < pi->wr[1]->njobs ; ij++) {
        // note that at this point, we're not yet MT.
        int err = MPI_Barrier(pi->m->pals);
        ASSERT_ALWAYS(!err);

        if (pi->wr[1]->jrank != ij) continue;

        if (pi->m->jrank == 0)
            print_several(pi->wr[0]->njobs,pi->wr[0]->ncores, '-', '+', 12);

        for(unsigned int it = 0 ; it < pi->wr[1]->ncores ; it++) {
            if (pi->wr[0]->jrank == 0) printf("|");
            for(unsigned int jj = 0 ; jj < pi->wr[0]->njobs ; jj++) {
                err = MPI_Barrier(pi->wr[0]->pals);
                ASSERT_ALWAYS(!err);
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
                    ASSERT_ALWAYS(!err);
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
}
#endif

/* given a grid of processes which has beed set up with pi_grid_init,
 * start the processes.
 *
 * The ngrids argument is here to accomodate the interleaved case.
 */
static void pi_go_mt_now(
        parallelizing_info ** grids,
        unsigned int ngrids,
        unsigned int n,
        void *(*fcn)(parallelizing_info_ptr, param_list pl, void *),
        param_list pl,
        void * arg)
{
    my_pthread_t * tgrid;
    struct pi_go_helper_s * helper_grid;

    tgrid = (my_pthread_t *) malloc(n * ngrids * sizeof(my_pthread_t));
    helper_grid = malloc(n * ngrids * sizeof(struct pi_go_helper_s));

    for(unsigned int k = 0 ; k < n * ngrids ; k++) {
        helper_grid[k].fcn = fcn;
        helper_grid[k].pl = pl;
        helper_grid[k].arg = arg;
        helper_grid[k].p = (parallelizing_info_ptr) grids[k/n] + (k%n);
        my_pthread_create(tgrid+k, NULL,
                (pthread_callee_t) pi_go_helper_func, helper_grid+k);
    }

    // nothing done here. We're the main job.
    for(unsigned int k = 0 ; k < n * ngrids ; k++) {
        my_pthread_join(tgrid[k], NULL);
    }
    free(helper_grid);
    free(tgrid);
}

static void pi_grid_clear(parallelizing_info_ptr pi, parallelizing_info * grid)
{
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
    for(unsigned int k = 0 ; k < pi->m->ncores ; k++) {
        MPI_Comm_free(&grid[k]->m->pals);
        MPI_Comm_free(&grid[k]->wr[0]->pals);
        MPI_Comm_free(&grid[k]->wr[1]->pals);
    }
#endif  /* MPI_LIBRARY_MT_CAPABLE */

    /* Don't do pi_log_clear, just like we haven't done pi_log_init. The
     * sub-threads may define it, in which case it's their responsibility
     * to clear the thing.
     */

    free(grid);
}

static void pi_clear_mpilevel(parallelizing_info_ptr pi)
{
#if defined(HAVE_HWLOC) && defined(HAVE_CXX11)
    int thr[2] = { pi->wr[0]->ncores, pi->wr[1]->ncores };
    cpubinding_free_info(pi->cpubinding_info, thr);
#endif /* defined(HAVE_HWLOC) && defined(HAVE_CXX11) */

    pi_wiring_destroy_pthread_things(pi->m);

    for(int d = 0 ; d < 2 ; d++) {
        MPI_Comm_free(&pi->wr[d]->pals);
    }
    MPI_Comm_free(&pi->m->pals);

    // MPI_Errhandler_free(&my_eh);
}

#if 0
static void shout_going_mt(pi_wiring_ptr m)
{
    int err = MPI_Barrier(m->pals);
    ASSERT_ALWAYS(!err);
    if (m->jrank == m->njobs - 1)
        printf("going multithread now\n");
    err = MPI_Barrier(m->pals);
    ASSERT_ALWAYS(!err);
}
#endif

static void pi_go_inner_not_interleaved(
        void *(*fcn)(parallelizing_info_ptr, param_list pl, void *),
        param_list pl,
        void * arg)
{
    parallelizing_info pi;
    parallelizing_info * grid;
    pi_init_mpilevel(pi, pl);
    grid = pi_grid_init(pi);
    // shout_going_mt(pi->m);
    // pi_grid_print_sketch(pi, grids[0]);
    pi_go_mt_now(&grid, 1, pi->m->ncores, fcn, pl, arg);
    pi_grid_clear(pi, grid);
    pi_clear_mpilevel(pi);
}

static void pi_go_inner_interleaved(
        void *(*fcn)(parallelizing_info_ptr, param_list pl, void *),
        param_list pl,
        void * arg)
{
    parallelizing_info pi[2];
    parallelizing_info * grids[2];

    pi_init_mpilevel(pi[0], pl);
    memcpy(pi[1], pi[0], sizeof(parallelizing_info));

    pi[0]->interleaved = malloc(sizeof(pi_interleaving));
    pi[0]->interleaved->idx = 0;

    pi[1]->interleaved = malloc(sizeof(pi_interleaving));
    pi[1]->interleaved->idx = 1;

    /* Now the whole point is that it's the _same_ barrier ! */
    my_pthread_barrier_t * b = malloc(sizeof(my_pthread_barrier_t));
    my_pthread_barrier_init(b, NULL, 2 * pi[0]->m->ncores);
    pi[0]->interleaved->b = b;
    pi[1]->interleaved->b = b;

    pi_dictionary d;
    pi_dictionary_init(d);
    pi[0]->dict = d;
    pi[1]->dict = d;

    grids[0] = pi_grid_init(pi[0]);
    grids[1] = pi_grid_init(pi[1]);

    // shout_going_mt(pi[0]->m);
    // pi_grid_print_sketch(pi, grids[0]);

    pi_go_mt_now(grids, 2, pi[0]->m->ncores, fcn, pl, arg);

    pi_grid_clear(pi[0], grids[0]);
    pi_grid_clear(pi[1], grids[1]);

    my_pthread_barrier_destroy(b);
    pi_dictionary_clear(d);

    free(pi[0]->interleaved);
    free(pi[1]->interleaved);

    pi_clear_mpilevel(pi[0]);
    /* pi[1] is a simple copy, it doesn't have to be cleared. */
}

void pi_dictionary_init(pi_dictionary_ptr d)
{
    my_pthread_rwlock_init(d->m, NULL);
    d->e = NULL;
}

void pi_dictionary_clear(pi_dictionary_ptr d)
{
    /* last grab of the lock -- altough it would be preferred if we
     * weren't mt of course ! */
    my_pthread_rwlock_wrlock(d->m);
    pi_dictionary_entry_ptr ne = d->e;
    for( ; ne ; ) {
        pi_dictionary_entry_ptr nne = ne->next;
        free(ne);
        ne = nne;
    }
    d->e = NULL;
    my_pthread_rwlock_unlock(d->m);
    my_pthread_rwlock_destroy(d->m);
}

void pi_store_generic(parallelizing_info_ptr pi, unsigned long key, unsigned long who, void * value)
{
    pi_dictionary_ptr d = pi->dict;
    ASSERT_ALWAYS(d != NULL);
    pi_dictionary_entry_ptr n = malloc(sizeof(pi_dictionary_entry));
    n->key = key;
    n->who = who;
    n->value = value;
    my_pthread_rwlock_wrlock(d->m);
    n->next = d->e;
    d->e = n;
    my_pthread_rwlock_unlock(d->m);
}

void * pi_load_generic(parallelizing_info_ptr pi, unsigned long key, unsigned long who)
{
    pi_dictionary_ptr d = pi->dict;
    ASSERT_ALWAYS(d != NULL);
    void * r = NULL;
    my_pthread_rwlock_rdlock(d->m);
    pi_dictionary_entry_ptr e = d->e;
    for( ; e ; e = e->next) {
        if (e->key == key && e->who == who) {
            r = e->value;
            break;
        }
    }
    my_pthread_rwlock_unlock(d->m);
    return r;
}

void pi_go(void *(*fcn)(parallelizing_info_ptr, param_list pl, void *),
        param_list pl,
        void * arg)
{
    int interleaving = 0;
#ifndef HAVE_MMAN_H
    fprintf(stderr, "Warning: copy-based mmaped I/O replacement has never been tested, and could very well be bogus\n");
#endif
    param_list_parse_int(pl, "interleaving", &interleaving);
    if (interleaving) {
        pi_go_inner_interleaved(fcn, pl, arg);
    } else {
        pi_go_inner_not_interleaved(fcn, pl, arg);
    }
}


void pi_log_init(pi_wiring_ptr wr)
{
    struct pi_log_book * lb = malloc(sizeof(struct pi_log_book));
    memset(lb, 0, sizeof(struct pi_log_book));
    wr->log_book = lb;

    char commname[MPI_MAX_OBJECT_NAME];
    int namelen = sizeof(commname);
    MPI_Comm_get_name(wr->pals, commname, &namelen);
    if (wr->jrank == 0 && wr->trank == 0) {
        printf("Enabled logging for %s\n", commname);
    }
}

void pi_log_clear(pi_wiring_ptr wr)
{
    if (wr->log_book)
        free(wr->log_book);
    wr->log_book = (void*) 0xdeadbeef;
}

void pi_log_op(pi_wiring_ptr wr, const char * fmt, ...)
{
    va_list ap;
    struct pi_log_book * lb = wr->log_book;
    if (!lb)
        return;

    va_start(ap, fmt);

    ASSERT_ALWAYS(lb->next < PI_LOG_BOOK_ENTRIES);

    struct pi_log_entry * e = &(lb->t[lb->next]);

    gettimeofday(e->tv, NULL);

    vsnprintf(e->what, sizeof(e->what), fmt, ap);
    va_end(ap);

    lb->next++;
    lb->next %= PI_LOG_BOOK_ENTRIES;
    lb->hsize++;
}

static void pi_log_print_backend(pi_wiring_ptr wr, const char * myname, char ** strings, int * n, int alloc)
{
    struct pi_log_book * lb = wr->log_book;
    if (!lb) return;

    ASSERT_ALWAYS(lb->next < PI_LOG_BOOK_ENTRIES);

    int h = lb->hsize;
    if (h >= PI_LOG_BOOK_ENTRIES) {
        h = PI_LOG_BOOK_ENTRIES;
    }
    int i0 = lb->next - h;
    if (i0 < 0)
        i0 += PI_LOG_BOOK_ENTRIES;
    char commname[MPI_MAX_OBJECT_NAME];
    int namelen = sizeof(commname);
    MPI_Comm_get_name(wr->pals, commname, &namelen);

    int rc;

    for(int i = i0 ; h-- ; ) {
        struct pi_log_entry * e = &(lb->t[i]);

        ASSERT_ALWAYS(*n < alloc);
        rc = asprintf(&(strings[*n]), "%" PRIu64 ".%06u (%s) %s %s", 
                (uint64_t) e->tv->tv_sec,
                (unsigned int) e->tv->tv_usec,
                myname,
                e->what,
                commname);
        ASSERT_ALWAYS(rc != -1);
        ++*n;

        i++;
        if (i == PI_LOG_BOOK_ENTRIES)
            i = 0;
    }
}

void pi_log_print(pi_wiring_ptr wr)
{
    /* FIXME stdout or stderr ? */
    int alloc = PI_LOG_BOOK_ENTRIES;
    char ** strings = malloc(alloc * sizeof(char*));
    int n = 0;
    pi_log_print_backend(wr, "", strings, &n, alloc);
    for(int i = 0 ; i < n ; i++) {
        puts(strings[i]);
        free(strings[i]);
    }
    fflush(stdout);
    free(strings);
}

typedef int (*sortfunc_t)(const void *, const void *);

int p_strcmp(char const * const * a, char const * const * b)
{
    /* Input is formatted so that sorting makes sense. */
    return strcmp(*a, *b);
}

void pi_log_print_all(parallelizing_info_ptr pi)
{
    int alloc = 3 * PI_LOG_BOOK_ENTRIES;
    char ** strings = malloc(alloc * sizeof(char*));
    int n = 0;
    char * myname;
    int rc;
    rc = asprintf(&myname, "%s%s", pi->wr[0]->th->desc, pi->wr[1]->th->desc);
    ASSERT_ALWAYS(rc != -1);

    pi_log_print_backend(pi->m, myname, strings, &n, alloc);
    pi_log_print_backend(pi->wr[0], myname, strings, &n, alloc);
    pi_log_print_backend(pi->wr[1], myname, strings, &n, alloc);
    qsort(strings, n, sizeof(char*), (sortfunc_t) &p_strcmp);

    for(int i = 0 ; i < n ; i++) {
        puts(strings[i]);
        free(strings[i]);
    }
    free(myname);
    free(strings);
}


struct thread_broadcast_arg {
    pi_wiring_ptr wr;
    /* this pointer is written from all threads. Its contents are
     * meaningful only for the root thread */
    void * ptr;
    unsigned int root;
    size_t size;
    thread_reducer_t f;         /* only for reduce */
    void * dptr;                /* only for reduce */
};

void thread_broadcast_in(int s MAYBE_UNUSED, struct thread_broadcast_arg * a)
{
    /* here we have a mutex locked */
    if (a->wr->trank == a->root) {
        a->wr->th->utility_ptr = a->ptr;
    }
}

void thread_broadcast_out(int s MAYBE_UNUSED, struct thread_broadcast_arg * a)
{
    /* here we know that exactly one thread has just filled the
     * a->wr->th->utility_ptr, so that a memcpy is possible, as long as
     * the owner of the data area has not touched it again */
    if (a->ptr != a->wr->th->utility_ptr)       /* placate valgrind */
        memcpy(a->ptr, a->wr->th->utility_ptr, a->size);
}

/* broadcast the data area pointed to by [ptr, ptr+size[ to all threads.
 * Pointers at calling threads must differ. */
void thread_broadcast(pi_wiring_ptr wr, void * ptr, size_t size, unsigned int i)
{
    struct thread_broadcast_arg a[1];
    a->wr = wr;
    a->ptr = ptr;
    a->root = i;
    a->size = size;
    barrier_wait(wr->th->bh,
            (void(*)(int,void*)) &thread_broadcast_in,
            (void(*)(int,void*)) &thread_broadcast_out,
            a);
    serialize_threads(wr);
}

void thread_reducer_int_min(void * mine, const void * other, size_t size MAYBE_UNUSED)
{
    int * imine = mine;
    const int * iother = other;
    if (*iother < *imine) *imine = *iother;
}

void thread_reducer_int_max(void * mine, const void * other, size_t size MAYBE_UNUSED)
{
    int * imine = mine;
    const int * iother = other;
    if (*iother > *imine) *imine = *iother;
}

void thread_reducer_int_sum(void * mine, const void * other, size_t size MAYBE_UNUSED)
{
    int * imine = mine;
    const int * iother = other;
    *imine += *iother;
}

/* this gets called with s = (count-1), (count-2), ..., 0 */
void thread_allreduce_in(int s, struct thread_broadcast_arg * a)
{
    if (s < (int) a->wr->ncores - 1) {
        (*a->f)(a->ptr, a->wr->th->utility_ptr, a->size);
    }
    a->wr->th->utility_ptr = a->ptr;
    if (!s)
        memcpy(a->dptr, a->wr->th->utility_ptr, a->size);
}

void thread_allreduce(pi_wiring_ptr wr, void * dptr, void * ptr, size_t size, thread_reducer_t f)
{
    struct thread_broadcast_arg a[1];
    a->wr = wr;
    a->ptr = ptr;
    a->root = 0;        /* meaningless */
    a->size = size;
    a->f = f;
    a->dptr = dptr;
    barrier_wait(wr->th->bh,
            (void(*)(int,void*)) &thread_allreduce_in,
            NULL,
            a);
    serialize_threads(wr);
}

void * shared_malloc(pi_wiring_ptr wr, size_t size)
{
    void * ptr = NULL;
    if (wr->trank == 0) ptr = malloc(size);
    thread_broadcast(wr, &ptr, sizeof(void*), 0);
    return ptr;
}

void * shared_malloc_set_zero(pi_wiring_ptr wr, size_t size)
{
    void * ptr = NULL;
    if (wr->trank == 0) {
        ptr = malloc(size);
        memset(ptr, 0, size);
    }
    thread_broadcast(wr, &ptr, sizeof(void*), 0);
    return ptr;
}

void shared_free(pi_wiring_ptr wr, void * ptr)
{
    serialize_threads(wr);
    if (wr->trank == 0) free(ptr);
}

void global_broadcast(pi_wiring_ptr wr, void * ptr, size_t size, unsigned int j, unsigned int t)
{
    int err;

    ASSERT(j < wr->njobs);
    ASSERT(t < wr->ncores);
    if (wr->trank == t) {
        err = MPI_Bcast(ptr, size, MPI_BYTE, j, wr->pals);
        ASSERT_ALWAYS(!err);
    }
    thread_broadcast(wr, ptr, size, t);
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
        ASSERT_ALWAYS(!err);
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
    return w->trank == 0;
}

// XXX Many serializing calls in the two routines below have no
// functional use, apart from the fact that they exert some pressure on
// the mpi+pthreads way of doing things.

void say_hello(pi_wiring_ptr w, parallelizing_info_ptr pi MAYBE_UNUSED)
{
    serialize(w);
    for(unsigned int j = 0 ; j < w->njobs ; j++) {
        serialize(w);
        if ((unsigned int) w->jrank != j)
            continue;
        my_pthread_mutex_lock(w->th->m);
#ifdef CONCURRENCY_DEBUG
        /* Make it less verbose -- if it ever hangs in there, then
         * we can re-enable it */
        printf("(%s) J%uT%u ; %s%s ; (%s:j%ut%u) (%s:j%ut%u)\n",
                w->th->desc,
                pi->m->jrank,
                pi->m->trank,
                pi->wr[0]->th->desc,
                pi->wr[1]->th->desc,
                pi->wr[0]->th->desc, pi->wr[0]->jrank, pi->wr[0]->trank,
                pi->wr[1]->th->desc, pi->wr[1]->jrank, pi->wr[1]->trank
              );
#endif
        my_pthread_mutex_unlock(w->th->m);
    }
}

void hello(parallelizing_info_ptr pi)
{
    if (serialize(pi->m)) {
#ifdef  CONCURRENCY_DEBUG
        printf("Doing hello world loop\n");
#endif
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
#ifdef CONCURRENCY_DEBUG
        printf("OK: Finished hello world loop\n");
#endif
    }
    serialize(pi->m);
}

static int get_counts_and_displacements(pi_wiring_ptr w, int my_size,
        int * displs, int * counts)
{
    counts[w->jrank] = my_size;
    SEVERAL_THREADS_PLAY_MPI_BEGIN(w) {
        int err = MPI_Allgather(MPI_IN_PLACE, 1, MPI_INT, counts, 1, MPI_INT, w->pals);
        ASSERT_ALWAYS(!err);
    }
    SEVERAL_THREADS_PLAY_MPI_END;

    /* Now all threads know the size argument of their friend threads. We
     * need to know the displacement values. These are obtained by
     * looking at the progressive sum over threads.
     *
     * We're serialization-happy here, since anyway I/O s occur rarely
     */
    /* The _void thingy makes gcc happy, and I'm not too comfortable
     * about deciding whether those pesky aliasing warnings are false
     * positives or not. Better shut them up for good. */
    void * allcounts_void;
    if (w->trank == 0) {
        allcounts_void = malloc(w->totalsize*sizeof(int));
    }
    thread_broadcast(w, &allcounts_void, sizeof(void*), 0);
    int * allcounts = (int *) allcounts_void;
    for(unsigned int k = 0 ; k < w->njobs ; k++) {
        allcounts[k * w->ncores + w->trank] = counts[k];
    }
    serialize_threads(w);
    int dis = 0;
    for(unsigned int k = 0 ; k < w->njobs ; k++) {
        unsigned int t;
        for(t = 0 ; t < w->trank ; t++) {
            dis += allcounts[k * w->ncores + t];
        }
        displs[k] = dis;
        for( ; t < w->ncores ; t++) {
            dis += allcounts[k * w->ncores + t];
        }
    }
    serialize_threads(w);
    if (w->trank == 0) {
        free(allcounts);
    }
    return dis;
}

/* d indicates whether we read the ranks row by row (d==0) or column by
 * column */
int get_counts_and_displacements_2d(parallelizing_info_ptr pi, int d,
        int my_size, int * displs, int * counts)
{
    pi_wiring_ptr w = pi->m;

    memset(displs, -1, w->njobs * sizeof(int));

    counts[w->jrank] = my_size;
    SEVERAL_THREADS_PLAY_MPI_BEGIN(w) {
        int err = MPI_Allgather(MPI_IN_PLACE, 1, MPI_INT, counts, 1, MPI_INT, w->pals);
        ASSERT_ALWAYS(!err);
    }
    SEVERAL_THREADS_PLAY_MPI_END;

    /* Now all threads know the size argument of their friend threads. We
     * need to know the displacement values. These are obtained by
     * looking at the progressive sum over threads.
     *
     * We're serialization-happy here, since anyway I/O s occur rarely
     */
    void * allcounts_void;
    if (w->trank == 0) {
        allcounts_void = malloc(w->njobs*w->ncores*sizeof(int));
    }
    thread_broadcast(w, &allcounts_void, sizeof(int*), 0);
    int * allcounts = (int *) allcounts_void;

    // Given an orientation value d, the numbering for the thread (it,jt)
    // on job (ij,jj) is:
    //
    // it+pi->wr[d]->ncores*(ij+pi->wr[d]->njobs*(jt+pi->wr[!d]->ncores*jj))
    // where ij ranges over pi->wr[d]->njobs and jj over pi->wr[!d]->njobs
    // (it and jt similar of course)
    //
    // Since we're filling the array in a multithreaded way, we're only
    // interested by the relevant value of it,jt ; the contribution
    // depending on ij,jj can be added later on. So we split the sum:
    //
    // it+pi->wr[d]->ncores*pi->wr[d]->njobs*jt
    // +
    // pi->wr[d]->ncores*(ij+pi->wr[d]->njobs*pi->wr[!d]->ncores*jj)
    // == 
    // it+pi->wr[d]->ncores*pi->wr[d]->njobs*jt
    // +
    // pi->wr[d]->ncores*ij + pi->m->ncores*pi->wr[d]->njobs*jj

    unsigned int v0;

    v0 = pi->wr[d]->trank+pi->wr[!d]->trank*pi->wr[d]->ncores*pi->wr[d]->njobs;
    for(unsigned int jj = 0 ; jj < pi->wr[!d]->njobs ; jj++) {
        unsigned int v = v0;
        /* d == 1 : we'll order offsets column by column */
        /* d == 0 : we'll order row column by row */
        for(unsigned int ij = 0 ; ij < pi->wr[d]->njobs ; ij++) {
            /* relative to pi->m, which job exactly corresponds to this
             * count ?
             * it's row_index * ncols + col_index. But we have to know
             * which is the row index.
             * d == 1 : ij is the row index
             * d == 0 : jj is the row index
             */
            unsigned int which = 0;
            which += ij * (d == 1 ? pi->wr[0]->njobs : 1);
            which += jj * (d == 0 ? pi->wr[0]->njobs : 1);
            ASSERT(which < pi->m->njobs);
            allcounts[v] = counts[which];
            v += pi->wr[d]->ncores;
        }
        v0 += pi->m->ncores*pi->wr[d]->njobs;
    }

    serialize_threads(w);

    int * alldisps = (int*) malloc(w->njobs*w->ncores*sizeof(int));
    int dis = 0;
    for(unsigned int k = 0 ; k < w->njobs * w->ncores ; k++) {
        alldisps[k] = dis;
        dis += allcounts[k];
    }

    v0 = pi->wr[d]->trank+pi->wr[!d]->trank*pi->wr[d]->ncores*pi->wr[d]->njobs;
    for(unsigned int jj = 0 ; jj < pi->wr[!d]->njobs ; jj++) {
        unsigned int v = v0;
        for(unsigned int ij = 0 ; ij < pi->wr[d]->njobs ; ij++) {
            unsigned int which = 0;
            which += ij * (d == 1 ? pi->wr[0]->njobs : 1);
            which += jj * (d == 0 ? pi->wr[0]->njobs : 1);
            ASSERT(which < pi->m->njobs);
            displs[which] = alldisps[v];
            v += pi->wr[d]->ncores;
        }
        v0 += pi->m->ncores*pi->wr[d]->njobs;
    }

    free(alldisps);

#if 0
    for(unsigned int i = 0 ; i < w->njobs ; i++) {
        my_pthread_mutex_lock(w->th->m);
        printf("J%uT%u %u @ %u\n", i, w->trank, counts[i], displs[i]);
        my_pthread_mutex_unlock(w->th->m);
    }
#endif

    serialize_threads(w);

#if 0
    ASSERT_ALWAYS(displs[w->jrank] >= 0);
    if (displs[w->jrank] == 0)
        ASSERT_ALWAYS(w->trank == 0);
#endif

    if (w->trank == 0) {
        free(allcounts);
    }

    return dis;
}

int truncate_counts_and_displacements(pi_wiring_ptr w, int * displs, int * counts, int sizeondisk)
{
    for(unsigned int k = 0 ; k < w->njobs ; k++) {
        if (displs[k] + counts[k] >= sizeondisk) {
            if (displs[k] >= sizeondisk) {
                displs[k] = sizeondisk;
                counts[k] = 0;
            } else {
                counts[k] = sizeondisk - displs[k];
            }
        }
    }
    return sizeondisk;
}

int area_is_zero(const void * src, ptrdiff_t offset0, ptrdiff_t offset1)
{
    const char * b0 = pointer_arith_const(src, offset0);
    const char * b1 = pointer_arith_const(src, offset1);
    for( ; b0 < b1 ; b0++) {
        if (*b0) return 0;
    }
    return 1;
}

/* That's a workalike to MPI_File_write_ordered, except that we work with
 * the finer-grained pi_wiring_ptr's.
 */
int pi_save_file(pi_wiring_ptr w, const char * name, unsigned int iter, void * buf, size_t mysize, size_t sizeondisk)
{
    int * displs = (int *) malloc(w->njobs * sizeof(int));
    int * recvcounts = (int *) malloc(w->njobs * sizeof(int));
    size_t siz;
    siz = get_counts_and_displacements(w, mysize, displs, recvcounts);

    // the page size is always a power of two, so rounding to the next
    // multiple is easy.
    size_t wsiz = ((siz - 1) | (pagesize()-1)) + 1;
    int leader = w->jrank == 0 && w->trank == 0;
    void * recvbuf = NULL;
    int fd = -1;        // only used by leader
    int rc;
    int err;

    if (leader) {
        char * filename;
        int rc;
        rc = asprintf(&filename, "%s.%u", name, iter);
        FATAL_ERROR_CHECK(rc < 0, "out of memory");
#ifdef HAVE_MINGW
        fd = open(filename, O_RDWR | O_CREAT | O_BINARY, 0666);
#else
        fd = open(filename, O_RDWR | O_CREAT, 0666);
#endif
        if (fd < 0) {
            fprintf(stderr, "fopen(%s): %s\n", filename, strerror(errno));
            goto pi_save_file_leader_init_done;
        }

        rc = ftruncate(fd, wsiz);
        if (rc < 0) {
            fprintf(stderr, "ftruncate(%s): %s\n",
                    filename, strerror(errno));
            close(fd);
            goto pi_save_file_leader_init_done;
        }
#ifdef HAVE_MMAN_H
        recvbuf = mmap(NULL, wsiz, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
        if (recvbuf == MAP_FAILED) {
            fprintf(stderr, "mmap(%s): %s\n",
                    filename, strerror(errno));
            recvbuf = NULL;
            close(fd);
            goto pi_save_file_leader_init_done;
        }
#else
        recvbuf = malloc(wsiz);
        FATAL_ERROR_CHECK(!recvbuf, "out of memory");
#endif
pi_save_file_leader_init_done:
        free(filename);
    }

    /* Now all threads from job zero see the area mmaped by their leader
     * thread. It does not make sense on other jobs of course, since for
     * these, recvbuf is unused. */
    if (w->jrank == 0) {
        thread_broadcast(w, &recvbuf, sizeof(void*), 0);
    }

    /* Rather unfortunate, but error checking requires some checking. As
     * mentioned earlier, we don't feel concerned a lot by this, since
     * in our context I/Os are rare enough
     */
    int ok = recvbuf != NULL;
    SEVERAL_THREADS_PLAY_MPI_BEGIN(w) {
        err = MPI_Bcast(&ok, 1, MPI_INT, 0, w->pals);
        ASSERT_ALWAYS(!err);
    }
    SEVERAL_THREADS_PLAY_MPI_END;
    if (!ok) return 0;

    SEVERAL_THREADS_PLAY_MPI_BEGIN(w) {
        err = MPI_Gatherv(buf, mysize, MPI_BYTE,
                recvbuf, recvcounts, displs, MPI_BYTE,
                0, w->pals);
        ASSERT_ALWAYS(!err);
    }
    SEVERAL_THREADS_PLAY_MPI_END;

    // what a pain in the Xss. Because we don't exit with the mutex
    // locked, we can't do much... The main thread might end up calling
    // munmap before we complete the call to MPI_Gatherv.

    serialize_threads(w);

    free(recvcounts);
    free(displs);

    if (leader) {
        ASSERT_ALWAYS(area_is_zero(recvbuf, sizeondisk, siz));
#ifdef HAVE_MMAN_H
        munmap(recvbuf, wsiz);
#else
        write(fd, recvbuf, wsiz);
        free(recvbuf);
#endif
        rc = ftruncate(fd, sizeondisk);
        if (rc < 0) {
            fprintf(stderr, "ftruncate(): %s\n", strerror(errno));
            /* If only ftruncate failed, don't return an error */
        }
        close(fd);
    }

    return 1;
}

/* d indicates whether we read the ranks row by row (d==0) or column by
 * column */
int pi_save_file_2d(parallelizing_info_ptr pi, int d, const char * name, unsigned int iter, void * buf, size_t mysize, size_t sizeondisk)
{
    pi_wiring_ptr w = pi->m;

    int * displs = (int *) malloc(w->njobs * sizeof(int));
    int * recvcounts = (int *) malloc(w->njobs * sizeof(int));
    size_t siz;
    siz = get_counts_and_displacements_2d(pi, d, mysize, displs, recvcounts);

    // the page size is always a power of two, so rounding to the next
    // multiple is easy.
    size_t wsiz = ((siz - 1) | (pagesize()-1)) + 1;
    int leader = w->jrank == 0 && w->trank == 0;
    void * recvbuf = NULL;
    int fd = -1;        // only used by leader
    int rc;
    int err;

    if (leader) {
        char * filename;
        int rc;
        rc = asprintf(&filename, "%s.%u", name, iter);
        FATAL_ERROR_CHECK(rc < 0, "out of memory");
#ifdef HAVE_MINGW
        fd = open(filename, O_RDWR | O_CREAT | O_BINARY, 0666);
#else
        fd = open(filename, O_RDWR | O_CREAT, 0666);
#endif
        if (fd < 0) {
            fprintf(stderr, "fopen(%s): %s\n", filename, strerror(errno));
            goto pi_save_file_2d_leader_init_done;
        }

        rc = ftruncate(fd, wsiz);
        if (rc < 0) {
            fprintf(stderr, "ftruncate(%s): %s\n",
                    filename, strerror(errno));
            close(fd);
            goto pi_save_file_2d_leader_init_done;
        }
#ifdef HAVE_MMAN_H
        recvbuf = mmap(NULL, wsiz, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
        if (recvbuf == MAP_FAILED) {
            fprintf(stderr, "mmap(%s): %s\n",
                    filename, strerror(errno));
            recvbuf = NULL;
            close(fd);
            goto pi_save_file_2d_leader_init_done;
        }
#else
        recvbuf = malloc(wsiz);
        FATAL_ERROR_CHECK(!recvbuf, "out of memory");
#endif
pi_save_file_2d_leader_init_done:
        free(filename);
    }

    /* Now all threads from job zero see the area mmaped by their leader
     * thread. It does not make sense on other jobs of course, since for
     * these, recvbuf is unused. */
    if (w->jrank == 0) {
        thread_broadcast(w, &recvbuf, sizeof(void*), 0);
    }

    /* Rather unfortunate, but error checking requires some checking. As
     * mentioned earlier, we don't feel concerned a lot by this, since
     * inour context I/Os are rare enough
     */
    int ok = recvbuf != NULL;
    SEVERAL_THREADS_PLAY_MPI_BEGIN(w) {
        err = MPI_Bcast(&ok, 1, MPI_INT, 0, w->pals);
        ASSERT_ALWAYS(!err);
    }
    SEVERAL_THREADS_PLAY_MPI_END;
    if (!ok) return 0;

    SEVERAL_THREADS_PLAY_MPI_BEGIN(w) {
        err = MPI_Gatherv(buf, mysize, MPI_BYTE,
                recvbuf, recvcounts, displs, MPI_BYTE,
                0, w->pals);
        ASSERT_ALWAYS(!err);
    }
    SEVERAL_THREADS_PLAY_MPI_END;

    // what a pain in the Xss. Because we don't exit with the mutex
    // locked, we can't do much... The main thread might end up calling
    // munmap before we complete the call to MPI_Gatherv.

    serialize_threads(w);

    free(recvcounts);
    free(displs);

    if (leader) {
        ASSERT_ALWAYS(area_is_zero(recvbuf, sizeondisk, siz));
#ifdef HAVE_MMAN_H
        munmap(recvbuf, wsiz);
#else
        write(fd, recvbuf, wsiz);
        free(recvbuf);
#endif
        rc = ftruncate(fd, sizeondisk);
        if (rc < 0) {
            fprintf(stderr, "ftruncate(): %s\n", strerror(errno));
            /* If only ftruncate failed, don't return an error */
        }
        close(fd);
    }

    return 1;
}

int pi_load_file(pi_wiring_ptr w, const char * name, unsigned int iter, void * buf, size_t mysize, size_t sizeondisk)
{
    int * displs = (int *) malloc(w->njobs * sizeof(int));
    int * sendcounts = (int *) malloc(w->njobs * sizeof(int));
    size_t siz;

    siz = get_counts_and_displacements(w, mysize, displs, sendcounts);
    truncate_counts_and_displacements(w, displs, sendcounts, sizeondisk);

    // the page size is always a power of two, so rounding to the next
    // multiple is easy.
    size_t wsiz = ((siz - 1) | (pagesize()-1)) + 1;
    int leader = w->jrank == 0 && w->trank == 0;
    void * sendbuf = NULL;
    int fd = -1;        // only used by leader
    int err;

    if (leader) {
        char * filename;
        int rc;
        rc = asprintf(&filename, "%s.%u", name, iter);
        FATAL_ERROR_CHECK(rc < 0, "out of memory");
#ifdef HAVE_MINGW
        fd = open(filename, O_RDONLY | O_BINARY, 0666);
#else
        fd = open(filename, O_RDONLY, 0666);
#endif
        DIE_ERRNO_DIAG(fd < 0, "fopen", filename);
#ifdef HAVE_MMAN_H
        sendbuf = mmap(NULL, wsiz, PROT_READ, MAP_SHARED, fd, 0);
        DIE_ERRNO_DIAG(sendbuf == MAP_FAILED, "mmap", filename);
#else
        sendbuf = malloc(wsiz);
        FATAL_ERROR_CHECK(!sendbuf, "out of memory");
        read(fd, sendbuf, siz);
#endif
        free(filename);
    }

    if (w->jrank == 0) {
        thread_broadcast(w, &sendbuf, sizeof(void*), 0);
    }

    /* For loading, in contrast to saving, we simply abort if reading
     * can't work. So no need to agree on an ok value. */

    SEVERAL_THREADS_PLAY_MPI_BEGIN(w) {
        // the send count might actually be a tiny bit smaller.
        // ASSERT_ALWAYS(mysize == (size_t) sendcounts[w->jrank]);
#if 1
        /* There seems to be a concurrency problem with the openmpi
         * implementation of MPI_Scatterv. We resort to an emulation
         * which uses MPI_Scatter. It's tailored-fit to our purpose,
         * since we know that the send counts will almost all be equal in
         * most cases. We also use somewhat jumbo frames.
         */
        int chunksize = 1 << 20;
        void * tsendbuf = NULL;
        void * trecvbuf = malloc(chunksize);
        if (w->jrank == 0) {
            tsendbuf = malloc(w->njobs * chunksize);
            memset(tsendbuf, 0, w->njobs * chunksize);
        }
        memset(trecvbuf, 0, chunksize);
        int cs[w->njobs];
        int cd[w->njobs];
        for(unsigned int j = 0 ; j < w->njobs ; j++) {
            cs[j] = sendcounts[j];
            cd[j] = displs[j];
        }
        err = 0;
        for(int spin = 0;;spin++) {
            int ts[w->njobs];
            int ts0=1;
            for(unsigned int j = 0 ; j < w->njobs ; j++) {
                ts[j] = MIN(chunksize, cs[j]);
                if (ts[j]) {
                    if (w->jrank == 0)
                        memcpy(tsendbuf + j * chunksize, sendbuf + cd[j], ts[j]);
                    ts0=0;
                }
            }
            if (ts0) break;
            // printf("J%uT%u start mpi_scatter\n", w->jrank, w->trank);
            err = MPI_Scatter(tsendbuf, chunksize, MPI_BYTE, trecvbuf,
                    chunksize, MPI_BYTE, 0, w->pals);
            // printf("J%uT%u end mpi_scatter\n", w->jrank, w->trank);
            memcpy(buf + cd[w->jrank] - displs[w->jrank], trecvbuf, ts[w->jrank]);
            for(unsigned int j = 0 ; j < w->njobs ; j++) {
                cs[j] -= ts[j];
                cd[j] += ts[j];
            }
        }
        free(trecvbuf);
        if (w->jrank == 0) {
            free(tsendbuf);
        }
#else
        err = MPI_Scatterv(sendbuf, sendcounts, displs, MPI_BYTE,
                buf, sendcounts[w->jrank], MPI_BYTE,
                0, w->pals);
#endif
        ASSERT_ALWAYS(!err);
    }
    SEVERAL_THREADS_PLAY_MPI_END;

    /* clear the tail if needed ! */
    if (mysize >= (size_t) sendcounts[w->jrank]) {
        memset(pointer_arith(buf, sendcounts[w->jrank]),
                0,
                mysize - sendcounts[w->jrank]);
    }

    serialize_threads(w);

    free(sendcounts);
    free(displs);

    if (leader) {
#ifdef HAVE_MMAN_H
        munmap(sendbuf, wsiz);
#else
        free(sendbuf);
#endif
        close(fd);
    }
    return 1;
}

int pi_load_file_2d(parallelizing_info_ptr pi, int d, const char * name, unsigned int iter, void * buf, size_t mysize, size_t sizeondisk)
{
    pi_wiring_ptr w = pi->m;

    int * displs = (int *) malloc(w->njobs * sizeof(int));
    int * sendcounts = (int *) malloc(w->njobs * sizeof(int));
    size_t siz;
    siz = get_counts_and_displacements_2d(pi, d, mysize, displs, sendcounts);
    truncate_counts_and_displacements(w, displs, sendcounts, sizeondisk);

    // the page size is always a power of two, so rounding to the next
    // multiple is easy.
    size_t wsiz = ((siz - 1) | (pagesize()-1)) + 1;
    int leader = w->jrank == 0 && w->trank == 0;
    void * sendbuf = NULL;
    int fd = -1;        // only used by leader
    int err;

    if (leader) {
        char * filename;
        int rc;
        rc = asprintf(&filename, "%s.%u", name, iter);
        FATAL_ERROR_CHECK(rc < 0, "out of memory");
#ifdef HAVE_MINGW
        fd = open(filename, O_RDONLY | O_BINARY, 0666);
#else
        fd = open(filename, O_RDONLY, 0666);
#endif
        DIE_ERRNO_DIAG(fd < 0, "fopen", filename);
#ifdef HAVE_MMAN_H
        sendbuf = mmap(NULL, wsiz, PROT_READ, MAP_SHARED, fd, 0);
        DIE_ERRNO_DIAG(sendbuf == MAP_FAILED, "mmap", filename);
#else
        sendbuf = malloc(wsiz);
        FATAL_ERROR_CHECK(!sendbuf, "out of memory");
        read(fd, sendbuf, siz);
#endif
        free(filename);
    }

    if (w->jrank == 0) {
        thread_broadcast(w, &sendbuf, sizeof(void*), 0);
    }

    /* For loading, in contrast to saving, we simply abort if reading
     * can't work. So no need to agree on an ok value. */

    SEVERAL_THREADS_PLAY_MPI_BEGIN(w) {
        // the send count might actually be a tiny bit smaller.
        // ASSERT_ALWAYS(mysize == (size_t) sendcounts[w->jrank]);
        err = MPI_Scatterv(sendbuf, sendcounts, displs, MPI_BYTE,
                buf, sendcounts[w->jrank], MPI_BYTE,
                0, w->pals);
        ASSERT_ALWAYS(!err);
    }
    SEVERAL_THREADS_PLAY_MPI_END;

    /* clear the tail if needed ! */
    if (mysize >= (size_t) sendcounts[w->jrank]) {
        memset(pointer_arith(buf, sendcounts[w->jrank]),
                0,
                mysize - sendcounts[w->jrank]);
    }

    serialize_threads(w);

    free(sendcounts);
    free(displs);

    if (leader) {
#ifdef HAVE_MMAN_H
        munmap(sendbuf, wsiz);
#else
        free(sendbuf);
#endif
        close(fd);
    }
    return 1;
}

void pi_interleaving_flip(parallelizing_info_ptr pi)
{
    if (!pi->interleaved)
        return;
    my_pthread_barrier_wait(pi->interleaved->b);
}

void pi_interleaving_enter(parallelizing_info_ptr pi)
{
    if (!pi->interleaved || pi->interleaved->idx == 0)
        return;

    my_pthread_barrier_wait(pi->interleaved->b);
}

void pi_interleaving_leave(parallelizing_info_ptr pi)
{
    if (!pi->interleaved || pi->interleaved->idx == 1)
        return;

    my_pthread_barrier_wait(pi->interleaved->b);
}

#if 0
/* XXX the name is not very informative, notably about the important
 * point that pointer must be identical.
 */
/* Considering a data area with a pointer identical on all calling
 * threads, tell whether the contents agree on all MPI jobs+threads
 * (returns 1 if this is the case)
 *
 * More exactly, this function only cares about the pointer on thread 0.
 */
int mpi_data_eq(parallelizing_info_ptr pi, void *buffer, size_t sz)
{
    /* TODO: think about allocating less */
    int cmp = 0;
    if (pi->m->trank == 0) {
        void * b_and = malloc(sz);
        void * b_or = malloc(sz);
        MPI_Allreduce(buffer, b_and, sz, MPI_BYTE, MPI_BAND, pi->m->pals);
        MPI_Allreduce(buffer, b_or, sz, MPI_BYTE, MPI_BOR, pi->m->pals);
        cmp = memcmp(b_and, b_or, sz);
        free(b_and);
        free(b_or);
    }
    thread_broadcast(pi->m, &cmp, sizeof(int), 0);
    return cmp == 0;
}
#endif

/* Here we assume pointers which are different on several calling
 * threads, and we compare the contents. Results may disagree on the
 * different MPI nodes */
int thread_data_eq(parallelizing_info_ptr pi, void *buffer, size_t sz)
{
    void ** bufs = shared_malloc(pi->m, pi->m->ncores * sizeof(void*));
    int * oks = shared_malloc(pi->m, pi->m->ncores * sizeof(int));

    bufs[pi->m->trank] = buffer;
    oks[pi->m->trank] = 1;

    /* given 5 threads, here is for example what is done at each step.
     *
     * stride=1
     *  ok0 = (x0==x1)
     *  ok2 = (x2==x3)
     * stride=2
     *  ok0 = x0==x1==x2==x3
     * stride 4
     *  ok0 = x0==x1==x2==x3==x4
     */
    for(unsigned int stride = 1 ; stride < pi->m->ncores ; stride <<= 1) {
        unsigned int i = pi->m->trank;
        serialize_threads(pi->m);       /* because of this, we don't
                                           leave the loop with break */
        if (i + stride >= pi->m->ncores) continue;
        /* many threads are doing nothing here */
        if (i & stride) continue;
        oks[i] = oks[i] && oks[i + stride] && memcmp(bufs[i], bufs[i + stride], sz) == 0;
    }
    serialize_threads(pi->m);
    int ret = oks[0];
    shared_free(pi->m, oks);
    shared_free(pi->m, bufs);
    return ret;
}

int global_data_eq(parallelizing_info_ptr pi, void *buffer, size_t sz)
{
    int ok= thread_data_eq(pi, buffer, sz);
    if (pi->m->trank == 0) {
        MPI_Allreduce(MPI_IN_PLACE, &ok, 1, MPI_INT, MPI_BAND, pi->m->pals);
    }
    thread_broadcast(pi->m, &ok, sizeof(int), 0);
    return ok;
}

