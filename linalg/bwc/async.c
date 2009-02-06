#define _POSIX_C_SOURCE 200112L

#include <signal.h>
#include <string.h>

#include "async.h"
#include "rusage.h"
#include "macros.h"
#include "bw-common.h"

// int int_caught = 0;
int hup_caught = 0;

void timing_init(struct timing_data * t, int iter)
{
    t->go_mark = iter;
    t->last_print = iter;
    t->next_print = iter + 1;
    t->next_async_check = iter + 1;
    t->async_check_period = 1;

    t->current->job = seconds();
    t->current->thread = thread_seconds();
    t->current->wct = walltime_seconds();
    memcpy(t->go, t->current, sizeof(*t->go));
}

void timing_clear(struct timing_data * t MAYBE_UNUSED)
{
}

void timing_rare_checks(pi_wiring_ptr wr, struct timing_data * t, int iter)
{
    int tcan_print = can_print && wr->trank == 0;

    /* We've decided that it was time to check for asynchronous data.
     * Since it's an expensive operation, the whole point is to avoid
     * doing this check too often. */
    t->current->job = seconds();
    t->current->thread = thread_seconds();
    t->current->wct = walltime_seconds();

    /* First, re-evaluate the async checking period */
    double av = (t->current->thread - t->go->thread) / (iter - t->go_mark);
    double good_period = PREFERRED_ASYNC_LAG / av;
    int guess = 1 + (int) good_period;
    complete_broadcast(wr, &guess, sizeof(int), 0, 0);
    if (iter == t->go_mark + 1 && tcan_print) {
        printf("Checking for asynchronous events every %d iterations\n", guess);
    }
    t->next_async_check += guess;
    t->async_check_period = guess;

    /* Now do some possibly expensive checks */

    /* This read is unsafe. We need to protect it. */
    unsigned int caught_something = hup_caught; // || int_caught;

    /* propagate the unsafe read. */
    if (wr->trank == 0) {
        MPI_Allreduce(MPI_IN_PLACE, &caught_something, 1,
                MPI_UNSIGNED, MPI_MAX, wr->pals);
    }

    /* reconcile threads */
    serialize_threads(wr);      // bug in thread_agreement
    void * ptr = &caught_something;
    thread_agreement(wr, &ptr, 0);
    caught_something = * (unsigned int *) ptr;

    /* ok, got it. */

    if (!caught_something) {
        /* And the result of all this is... nothing :-(( */
        return;
    }

    /* We have something ! */
    if (tcan_print) {
        printf("Signal caught (iteration %d), restarting timings\n", iter);
    }
    serialize(wr);
    hup_caught = 0;
    // int_caught = 0;
    timing_init(t, iter);
}


void timing_check(parallelizing_info pi, struct timing_data * timing, int iter)
{
    pi_wiring_ptr wr = pi->m;

    // printf("timing_check %d timing%d\n", iter, wr->trank);

    if (iter == timing->next_async_check) {
        // printf("timing_check async %d timing%d\n", iter, wr->trank);
        timing_rare_checks(wr, timing, iter);
    }

    if (iter < timing->next_print)
        return;

    int tcan_print = can_print && wr->trank == 0;

    const int period[] = {1,5,20,100,1000,10000,100000,1000000,0};
    int i;
    for(i = 0 ; ; i++) {
        if (period[i] == 0 || iter % period[i])
            break;
    }
    i--;
    ASSERT(iter % period[i] == 0);
    timing->next_print = iter + period[i];

    timing->current->job = seconds();
    timing->current->thread = thread_seconds();
    timing->current->wct = walltime_seconds();

    /* ceinture et bretelles */

    char buf[10];
    double dt = timing->current->thread - timing->go->thread;
    double av = dt / (iter - timing->go_mark);
    double pcpu = 100.0 * dt / (timing->current->wct - timing->go->wct);
    snprintf(buf, sizeof(buf), "%.2f %2d%%", av, (int) pcpu);

    if (tcan_print)
        printf("iteration %d\n", iter);

    grid_print(pi, buf, sizeof(buf), tcan_print);
}

void block_control_signals()
{
    /* Only the master thread receives control signals */
    sigset_t sset[1];
    sigemptyset(sset);
    sigaddset(sset, SIGHUP);
    // sigaddset(sset, SIGINT);
    pthread_sigmask(SIG_BLOCK, sset, NULL);
}

void sighandler(int sig)
{
    int caught = 0;
    if (sig == SIGHUP) hup_caught = caught = 1;
    // if (sig == SIGINT) int_caught = caught = 1;
    /* Of course, everybody is allowed to print this. */
    if (caught) printf("Signal caught, wait before it is acknowledged\n");
}

void catch_control_signals()
{
    struct sigaction sa[1];
    memset(sa, 0, sizeof(sa));
    sa->sa_handler = sighandler;
    sigaction(SIGHUP, sa, NULL);
    // sigaction(SIGINT, sa, NULL);

    sigset_t sset[1];
    sigemptyset(sset);
    sigaddset(sset, SIGHUP);
    // sigaddset(sset, SIGINT);
    pthread_sigmask(SIG_UNBLOCK, sset, NULL);
}
