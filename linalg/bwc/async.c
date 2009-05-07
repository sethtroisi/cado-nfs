#define _POSIX_C_SOURCE 200112L
#include <signal.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "bwc_config.h"
#include "async.h"
#include "rusage.h"
#include "macros.h"
#include "bw-common.h"

// int int_caught = 0;
int hup_caught = 0;

void timing_update_ticks(struct timing_data * t, int iter MAYBE_UNUSED)
{
    job_seconds(t->current->job);
    thread_seconds(t->current->thread);
    t->current->wct = walltime_seconds();
}

void timing_partial_init(struct timing_data * t, int iter)
{
    t->go_mark = iter;
    t->last_print = iter;
    t->next_print = iter + 1;
    t->next_async_check = iter + 1;
    t->async_check_period = 1;

    timing_update_ticks(t, iter);
    memcpy(t->go, t->current, sizeof(*t->go));
}

void timing_init(struct timing_data * t, int start, int end)
{
    timing_partial_init(t, start);
    memcpy(t->beginning, t->go, sizeof(*t->go));
    t->begin_mark = start;
    t->end_mark = end;
}

void timing_clear(struct timing_data * t MAYBE_UNUSED)
{
}

void timing_rare_checks(pi_wiring_ptr wr, struct timing_data * t, int iter, int print)
{
    /* We've decided that it was time to check for asynchronous data.
     * Since it's an expensive operation, the whole point is to avoid
     * doing this check too often. */
    timing_update_ticks(t, iter);

    /* First, re-evaluate the async checking period */
    double av[2];

    av[0] = (t->current->thread[0] - t->go->thread[0])/(iter - t->go_mark);
    av[1] = (t->current->thread[1] - t->go->thread[1])/(iter - t->go_mark);

    double good_period = PREFERRED_ASYNC_LAG / av[0];
    int guess = 1 + (int) good_period;
    complete_broadcast(wr, &guess, sizeof(int), 0, 0);
    /* negative stuff is most probably caused by overflows in a fast
     * context */
    if (guess <= 0) guess = 10;
    if (iter == t->go_mark + 1 && print) {
        printf("Checking for asynchronous events every %d iterations\n", guess);
    }
    t->next_async_check += guess;
    t->async_check_period = guess;

    /* Now do some possibly expensive checks */

    /* This read is unsafe. We need to protect it. */
    unsigned int caught_something = hup_caught; // || int_caught;

    /* propagate the unsafe read. */
    if (wr->trank == 0) {
        int err = MPI_Allreduce(MPI_IN_PLACE, &caught_something, 1,
                MPI_UNSIGNED, MPI_MAX, wr->pals);
        ASSERT_ALWAYS(!err);
    }

    /* reconcile threads */
    serialize_threads(wr);      // bug in thread_agreement
    void * ptr = &caught_something;
    thread_agreement(wr, &ptr, 0);
    caught_something = * (unsigned int *) ptr;

    /* ok, got it. Before we possibly leave, make sure everybody has read
     * the data from the pointer. */
    serialize(wr);

    if (!caught_something) {
        /* And the result of all this is... nothing :-(( */
        return;
    }

    /* We have something ! */
    if (print) {
        printf("Signal caught (iteration %d), restarting timings\n", iter);
    }
    serialize(wr);
    hup_caught = 0;
    // int_caught = 0;
    timing_partial_init(t, iter);
}


void timing_check(parallelizing_info pi, struct timing_data * timing, int iter, int print)
{
    pi_wiring_ptr wr = pi->m;

    // printf("timing_check %d timing%d\n", iter, wr->trank);

    if (iter == timing->next_async_check) {
        // printf("timing_check async %d timing%d\n", iter, wr->trank);
        timing_rare_checks(wr, timing, iter, print);
    }

    if (iter < timing->next_print)
        return;

    /* We're printing something, so we might as well check for signals
     * now */
    timing_rare_checks(wr, timing, iter, print);

    const int period[] = {1,5,20,100,1000,10000,100000,1000000,0};
    int i;
    for(i = 0 ; ; i++) {
        if (period[i] == 0 || iter % period[i])
            break;
    }
    i--;
    ASSERT(iter % period[i] == 0);
    timing->next_print = iter + period[i];

    timing_update_ticks(timing, iter);

    char buf[20];

    double dt_full[2];
    double dt_recent[2];
    double av[2];
    double pcpu[2];
    double wct_recent;
    double di = iter - timing->go_mark;

    dt_full[0] = timing->current->thread[0] - timing->beginning->thread[0];
    dt_full[1] = timing->current->thread[1] - timing->beginning->thread[1];
    dt_recent[0] = timing->current->thread[0] - timing->go->thread[0];
    dt_recent[1] = timing->current->thread[1] - timing->go->thread[1];
    wct_recent = timing->current->wct - timing->go->wct;

    av[0] = dt_recent[0] / di;
    av[1] = dt_recent[1] / di;

    pcpu[0] = 100.0 * dt_recent[0] / wct_recent;
    pcpu[1] = 100.0 * dt_recent[1] / wct_recent;

    snprintf(buf, sizeof(buf), "%.2f %.2f %2d%%+%2d%%",
            dt_full[0], av[0], (int) pcpu[0], (int) pcpu[1]);

    if (print)
        printf("iteration %d\n", iter);

    grid_print(pi, buf, sizeof(buf), print);
}

void timing_disp_collective_oneline(parallelizing_info pi, struct timing_data * timing, int iter, unsigned long ncoeffs, int print)
{
    double dt_full[2];
    double dt_recent[2];
    double av[2];
    double pcpu[2];
    double wct_full;
    double wct_recent;

    timing_update_ticks(timing, iter);

    double di = iter - timing->go_mark;

    dt_full[0]   = timing->current->job[0] - timing->beginning->job[0];
    dt_full[1]   = timing->current->job[1] - timing->beginning->job[1];
    dt_recent[0] = timing->current->job[0]-timing->go->job[0];
    dt_recent[1] = timing->current->job[1]-timing->go->job[1];
    wct_full   = timing->current->wct - timing->beginning->wct;
    wct_recent = timing->current->wct - timing->go->wct;


    // dt must be collected.
    int err;
    double ncoeffs_d = ncoeffs;

    SEVERAL_THREADS_PLAY_MPI_BEGIN(pi->m) {
        err = MPI_Allreduce(MPI_IN_PLACE, dt_full, 2, MPI_DOUBLE, MPI_SUM, pi->m->pals);
        ASSERT_ALWAYS(!err);

        err = MPI_Allreduce(MPI_IN_PLACE, dt_recent, 2, MPI_DOUBLE, MPI_SUM, pi->m->pals);
        ASSERT_ALWAYS(!err);

        err = MPI_Allreduce(MPI_IN_PLACE, &ncoeffs_d, 1, MPI_DOUBLE, MPI_SUM, pi->m->pals);
        ASSERT_ALWAYS(!err);
    }
    SEVERAL_THREADS_PLAY_MPI_END;

    /* dt_recent[] contains the JOB seconds. So we don't have to sum it up over
     * threads. However, we do for ncoeffs ! */
    double ncoeffs_total = 0;
    void * ptr = &ncoeffs_total;
    thread_agreement(pi->m, &ptr, 0);
    double * main_ncoeffs_total = ptr;
    my_pthread_mutex_lock(pi->m->th->m);
    * main_ncoeffs_total += ncoeffs_d;
    my_pthread_mutex_unlock(pi->m->th->m);
    serialize_threads(pi->m);
    ncoeffs_d = * main_ncoeffs_total;

    av[0] = dt_recent[0] / di;
    av[1] = dt_recent[1] / di;
    pcpu[0] = 100.0 * dt_recent[0] / wct_recent;
    pcpu[1] = 100.0 * dt_recent[1] / wct_recent;
    if (di == 0)
        pcpu[0] = pcpu[1] = NAN;


    if (print) {

        double avcpu = av[0] + av[1];
        double nsc = avcpu / ncoeffs_d * 1.0e9;

        double avwct = wct_recent / di;

        /* still to go: timing->end_mark - iter */
        time_t now[1];
        time_t eta[1];
        char eta_string[32] = "not available yet\n";
        time(now);
        eta[0] = now[0] + (timing->end_mark - iter) * avwct;
        if (di)
            ctime_r(eta, eta_string);

        char * what_cpu = "s";
        char * what_wct = "s";
        if (avcpu < 0.1) { what_cpu = "ms"; avcpu *= 1000.0; }
        if (avwct < 0.1) { what_wct = "ms"; avwct *= 1000.0; }

        printf("N=%d ; CPU"
                ": %.2f cpu + %.2f sys"
                ", %.2f %s/iter"
                ", %.2f ns/coeff"
                "\n",
                iter, dt_full[0], dt_full[1],
                avcpu, what_cpu, nsc);

        printf("N=%d ; WCT"
                ": %.2f tot"
                ", %.2f %s/iter"
                ", %.2f%% cpu + %.2f%% sys\n",
                iter, wct_full, avwct, what_wct,
                pcpu[0], pcpu[1]);

        printf("N=%d ; ETA (N=%d): %s", iter, timing->end_mark, eta_string);
    }
}


void block_control_signals()
{
    /* Only the master thread receives control signals */
    sigset_t sset[1];
    sigemptyset(sset);
    sigaddset(sset, SIGHUP);
    // sigaddset(sset, SIGINT);
    my_pthread_sigmask(SIG_BLOCK, sset, NULL);
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
    my_pthread_sigmask(SIG_UNBLOCK, sset, NULL);
}
