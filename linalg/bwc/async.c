#include "cado.h"
#include <stddef.h>
#include <signal.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <ctype.h>
#include "bwc_config.h"
#include "async.h"
#include "portability.h"
#include "macros.h"
#include "bw-common.h"
#include "utils.h"

// int int_caught = 0;
int hup_caught = 0;

static void extract_interval(timing_interval_data * since_last_reset, timing_interval_data * since_beginning, struct timing_data * t)
{
    memcpy(since_last_reset, t->since_last_reset, NTIMERS_EACH_ITERATION * sizeof(struct timing_interval_data_s));
    memcpy(since_beginning, t->since_beginning, NTIMERS_EACH_ITERATION * sizeof(struct timing_interval_data_s));
    double d[2];
    seconds_user_sys(d);
    since_last_reset[t->which]->job[0] += d[0];
    since_last_reset[t->which]->job[1] += d[1];
    since_beginning[t->which]->job[0] += d[0];
    since_beginning[t->which]->job[1] += d[1];
    thread_seconds_user_sys(d);
    since_last_reset[t->which]->thread[0] += d[0];
    since_last_reset[t->which]->thread[1] += d[1];
    since_beginning[t->which]->thread[0] += d[0];
    since_beginning[t->which]->thread[1] += d[1];
    double w = wct_seconds();
    since_last_reset[t->which]->wct += w;
    since_beginning[t->which]->wct += w;
}

void timing_next_timer(struct timing_data * t)
{
    double d[2];
    seconds_user_sys(d);
    int next = (t->which + 1) % NTIMERS_EACH_ITERATION;
    t->since_last_reset[t->which]->job[0] += d[0];
    t->since_last_reset[t->which]->job[1] += d[1];
    t->since_beginning[t->which]->job[0] += d[0];
    t->since_beginning[t->which]->job[1] += d[1];
    t->since_last_reset[next]->job[0] -= d[0];
    t->since_last_reset[next]->job[1] -= d[1];
    t->since_beginning[next]->job[0] -= d[0];
    t->since_beginning[next]->job[1] -= d[1];
    thread_seconds_user_sys(d);
    t->since_last_reset[t->which]->thread[0] += d[0];
    t->since_last_reset[t->which]->thread[1] += d[1];
    t->since_beginning[t->which]->thread[0] += d[0];
    t->since_beginning[t->which]->thread[1] += d[1];
    t->since_last_reset[next]->thread[0] -= d[0];
    t->since_last_reset[next]->thread[1] -= d[1];
    t->since_beginning[next]->thread[0] -= d[0];
    t->since_beginning[next]->thread[1] -= d[1];
    double w = wct_seconds();
    t->since_last_reset[t->which]->wct += w;
    t->since_last_reset[next]->wct -= w;
    t->since_beginning[t->which]->wct += w;
    t->since_beginning[next]->wct -= w;
    t->which = next;
}

static void timing_partial_init(struct timing_data * t, int iter)
{
    memset(t, 0, sizeof(struct timing_data));
    t->go_mark = iter;
    t->last_print = iter;
    t->next_print = iter + 1;
    t->next_async_check = iter + 1;
    t->async_check_period = 1;

    double d[2];
    seconds_user_sys(d);
    t->since_last_reset[t->which]->job[0] = -d[0];
    t->since_last_reset[t->which]->job[1] = -d[1];

    thread_seconds_user_sys(d);
    t->since_last_reset[t->which]->thread[0] = -d[0];
    t->since_last_reset[t->which]->thread[1] = -d[1];

    double w = wct_seconds();
    t->since_last_reset[t->which]->wct = -w;
}

void timing_init(struct timing_data * t, int start, int end)
{
    timing_partial_init(t, start);
    memcpy(t->since_beginning, t->since_last_reset, sizeof(t->since_last_reset));
    t->begin_mark = start;
    t->end_mark = end;
}

void timing_clear(struct timing_data * t MAYBE_UNUSED)
{
}

static void timing_rare_checks(pi_wiring_ptr wr, struct timing_data * t, int iter, int print)
{
    /* We've decided that it was time to check for asynchronous data.
     * Since it's an expensive operation, the whole point is to avoid
     * doing this check too often. */
    // timing_update_ticks(t, iter);

    timing_interval_data since_last_reset[NTIMERS_EACH_ITERATION];
    timing_interval_data since_beginning[NTIMERS_EACH_ITERATION];
    extract_interval(since_last_reset, since_beginning, t);

    /* First, re-evaluate the async checking period */
    double av;
    av = (since_last_reset[0]->wct + since_last_reset[1]->wct) / (iter - t->go_mark);

#ifndef MPI_LIBRARY_MT_CAPABLE
    /* Other threads might still be lingering in the matrix
     * multiplication routines. That's really not good, since we're also
     * going to do mpi ourselves ! */
    serialize_threads(wr);
#endif  /* MPI_LIBRARY_MT_CAPABLE */

    double good_period = PREFERRED_ASYNC_LAG / av;
    int guess = 1 + (int) good_period;
    global_broadcast(wr, &guess, sizeof(int), 0, 0);
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
    thread_broadcast(wr, &caught_something, sizeof(unsigned int), 0);

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

    if (!verbose_enabled(CADO_VERBOSE_PRINT_BWC_TIMING_GRIDS))
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

    // timing_update_ticks(timing, iter);

    char buf[20];

    timing_interval_data since_last_reset[NTIMERS_EACH_ITERATION];
    timing_interval_data since_beginning[NTIMERS_EACH_ITERATION];
    extract_interval(since_last_reset, since_beginning, timing);
    double di = iter - timing->go_mark;

#ifdef MEASURE_LINALG_JITTER_TIMINGS
    double cpu = since_last_reset[0]->wct;
    double cpuw = since_last_reset[1]->wct;
    double comm = since_last_reset[2]->wct;
    double commw = since_last_reset[3]->wct;
    double cput = cpu + cpuw;
    double commt = comm + commw;
    double thrcpu = since_last_reset[0]->thread[0];
    double thrcomm = since_last_reset[2]->thread[0];
#else
    double cput = since_last_reset[0]->wct;
    double commt = since_last_reset[1]->wct;
    double thrcpu = since_last_reset[0]->thread[0];
    double thrcomm = since_last_reset[1]->thread[0];
#endif
    // (avg wct cpu)(cpu % cpu) + (avg comm)(cpu % comm)
    snprintf(buf, sizeof(buf), "%.2f@%.0f%%+%.2f@%.0f%%",
            cput/di, 100.0 * thrcpu / cput,
            commt/di, 100.0 * thrcomm / commt);

    if (print)
        printf("iteration %d\n", iter);

    grid_print(pi, buf, sizeof(buf), print);
}

/* This adds the contents of x[0..n[ on all threads, and communicates the
 * result to everyone. This is a thread-level only function !
 */
void pi_thread_allreduce_add_double(parallelizing_info pi, double * x, int n)
{
    /* Broadcast master's x */
    double * master_x = (pi->m->trank == 0) ? x : NULL;
    thread_broadcast(pi->m, (void **) &master_x, sizeof(void*), 0);
    /* Threads other than master add their x to master's x */
    if (pi->m->trank != 0) {
        my_pthread_mutex_lock(pi->m->th->m);
        for(int i = 0 ; i < n ; i++)
            master_x[i] += x[i];
        my_pthread_mutex_unlock(pi->m->th->m);
    }
    serialize_threads(pi->m);
    /* Master has totals in x, everyone else needs to update their x */
    if (pi->m->trank != 0)
        memcpy(x, master_x, n * sizeof(int));
    /* as long as we have a read intention on master_x, we should not get
     * past this barrier */
    serialize_threads(pi->m);
}

void pi_thread_allreduce_min_double(parallelizing_info pi, double * x, int n)
{
    /* Broadcast master's x */
    double * master_x = (pi->m->trank == 0) ? x : NULL;
    thread_broadcast(pi->m, (void **) &master_x, sizeof(void*), 0);
    /* Threads other than master add their x to master's x */
    if (pi->m->trank != 0) {
        my_pthread_mutex_lock(pi->m->th->m);
        for(int i = 0 ; i < n ; i++)
            master_x[i] = MIN(master_x[i], x[i]);
        my_pthread_mutex_unlock(pi->m->th->m);
    }
    serialize_threads(pi->m);
    /* Master has totals in x, everyone else needs to update their x */
    if (pi->m->trank != 0)
        memcpy(x, master_x, n * sizeof(int));
    /* as long as we have a read intention on master_x, we should not get
     * past this barrier */
    serialize_threads(pi->m);
}

void pi_thread_allreduce_max_double(parallelizing_info pi, double * x, int n)
{
    /* Broadcast master's x */
    double * master_x = (pi->m->trank == 0) ? x : NULL;
    thread_broadcast(pi->m, (void **) &master_x, sizeof(void*), 0);
    /* Threads other than master add their x to master's x */
    if (pi->m->trank != 0) {
        my_pthread_mutex_lock(pi->m->th->m);
        for(int i = 0 ; i < n ; i++)
            master_x[i] = MAX(master_x[i], x[i]);
        my_pthread_mutex_unlock(pi->m->th->m);
    }
    serialize_threads(pi->m);
    /* Master has totals in x, everyone else needs to update their x */
    if (pi->m->trank != 0)
        memcpy(x, master_x, n * sizeof(int));
    /* as long as we have a read intention on master_x, we should not get
     * past this barrier */
    serialize_threads(pi->m);
}

static const char * timer_names[] = TIMER_NAMES;

/* stage=0 for krylov, 1 for mksol */
void timing_disp_backend(parallelizing_info pi, struct timing_data * timing, int iter, unsigned long ncoeffs, int print, int stage, int done)
{
    if (!verbose_enabled(CADO_VERBOSE_PRINT_BWC_ITERATION_TIMINGS))
        return;

    timing_interval_data since_last_reset[NTIMERS_EACH_ITERATION];
    timing_interval_data since_beginning[NTIMERS_EACH_ITERATION];
    extract_interval(since_last_reset, since_beginning, timing);

    timing_interval_data * T = done ? since_beginning : since_last_reset;

    double di = iter - timing->go_mark;
    double ncoeffs_d = ncoeffs;

    /* for each timer, we want to print:
     *
     *  - the total accumulated number of seconds spent. This really is
     *    the addition of all the per-thread WCTs.
     *  - average [min..max] for these
     *  - the CPU % for each. This is \sum job[0] / total wct.
     */

    timing_interval_data Tmin[NTIMERS_EACH_ITERATION];
    timing_interval_data Tmax[NTIMERS_EACH_ITERATION];
    timing_interval_data Tsum[NTIMERS_EACH_ITERATION];

    int ndoubles = NTIMERS_EACH_ITERATION * sizeof(timing_interval_data) / sizeof(double);

    memcpy(Tmin, T, sizeof(Tmin));
    memcpy(Tmax, T, sizeof(Tmax));
    memcpy(Tsum, T, sizeof(Tsum));

    SEVERAL_THREADS_PLAY_MPI_BEGIN(pi->m) {
        int err;
        err = MPI_Allreduce(MPI_IN_PLACE, Tsum, ndoubles, MPI_DOUBLE, MPI_SUM, pi->m->pals);
        ASSERT_ALWAYS(!err);
        err = MPI_Allreduce(MPI_IN_PLACE, Tmin, ndoubles, MPI_DOUBLE, MPI_MIN, pi->m->pals);
        ASSERT_ALWAYS(!err);
        err = MPI_Allreduce(MPI_IN_PLACE, Tmax, ndoubles, MPI_DOUBLE, MPI_MAX, pi->m->pals);
        ASSERT_ALWAYS(!err);

        err = MPI_Allreduce(MPI_IN_PLACE, &ncoeffs_d, 1, MPI_DOUBLE, MPI_SUM, pi->m->pals);
        ASSERT_ALWAYS(!err);
    }
    SEVERAL_THREADS_PLAY_MPI_END;

    pi_thread_allreduce_add_double(pi, (double*) Tsum, ndoubles);
    pi_thread_allreduce_min_double(pi, (double*) Tmin, ndoubles);
    pi_thread_allreduce_max_double(pi, (double*) Tmax, ndoubles);
    pi_thread_allreduce_add_double(pi, &ncoeffs_d, 1);

    double sum_dwct = 0;
    for(int i = 0 ; i < NTIMERS_EACH_ITERATION ; i++) {
        sum_dwct += Tsum[i]->wct;
    }
    double avdwct = sum_dwct / pi->m->totalsize / di;

    if (print) {
        for(int timer = 0 ; timer < NTIMERS_EACH_ITERATION ; timer++) {
            double avwct = Tsum[timer]->wct / pi->m->totalsize / di;

            char extra[32]={'\0'};
            if (timer == 0) {
                // nanoseconds per coefficient are computed based on the
                // total (aggregated) wall-clock time per iteration within
                // the cpu-bound part, and divided by the total number of
                // coefficients.
                double nsc = Tsum[timer]->wct / di / ncoeffs_d * 1.0e9;
                snprintf(extra, sizeof(extra), ", %.2f ns/%.1fGcoeff", nsc, ncoeffs_d*1.0e-9);
            }


            printf("%sN=%d ; %s: %.2f s/iter [%.3f..%.3f]%s\n",
                    done ? "done " : "",
                    iter,
                    timer_names[timer],
                    avwct,
                    Tmin[timer]->wct/di,
                    Tmax[timer]->wct/di,
                    extra);
        }

        if (!done) {
            /* still to go: timing->end_mark - iter */
            time_t now[1];
            time_t eta[1];
            char eta_string[32] = "not available yet\n";
            time(now);
            *eta = *now + (timing->end_mark - iter) * avdwct;
            if (di) {
#ifdef HAVE_CTIME_R
                ctime_r(eta, eta_string);
#else
                strncpy(eta_string, ctime(eta), sizeof(eta_string));
#endif
            }


            unsigned int s = strlen(eta_string);
            for( ; s && isspace((int)(unsigned char)eta_string[s-1]) ; eta_string[--s]='\0') ;

            printf("%s: N=%d ; ETA (N=%d): %s [%.3f s/iter]\n",
                   (stage == 0) ? "krylov" : "mksol",
                   iter, timing->end_mark, eta_string, avdwct);
        }
    }
    /* We're sharing via thread_broadcast data which sits on the stack of
     * one thread. So it's important that no thread exits this function
     * prematurely ! */
    serialize_threads(pi->m);
}

void timing_disp_collective_oneline(parallelizing_info pi, struct timing_data * timing, int iter, unsigned long ncoeffs, int print, int stage)
{
    timing_disp_backend(pi, timing, iter, ncoeffs, print, stage, 0);
}

void timing_final_tally(parallelizing_info pi, struct timing_data * timing, unsigned long ncoeffs, int print, int stage)
{
    timing_disp_backend(pi, timing, timing->end_mark, ncoeffs, print, stage, 1);
}
void block_control_signals()
{
#ifndef HAVE_MINGW   /* seems hopeless */
    /* Only the master thread receives control signals */
    sigset_t sset[1];
    sigemptyset(sset);
#ifdef HAVE_SIGHUP
    sigaddset(sset, SIGHUP);
#endif
    // sigaddset(sset, SIGINT);
    my_pthread_sigmask(SIG_BLOCK, sset, NULL);
#endif  /* HAVE_MINGW */
}

#ifndef HAVE_MINGW   /* seems hopeless */
void sighandler(int sig)
{
    int caught = 0;
#ifdef HAVE_SIGHUP
    if (sig == SIGHUP) hup_caught = caught = 1;
#endif
    // if (sig == SIGINT) int_caught = caught = 1;
    /* Of course, everybody is allowed to print this. */
    if (caught) printf("Signal caught, wait before it is acknowledged\n");
}
#endif  /* HAVE_MINGW */

void catch_control_signals()
{
#ifndef HAVE_MINGW   /* seems hopeless */
#ifdef HAVE_SIGACTION
    struct sigaction sa[1];
    memset(sa, 0, sizeof(sa));
    sa->sa_handler = sighandler;
#ifdef HAVE_SIGHUP
    sigaction(SIGHUP, sa, NULL);
#endif
    // sigaction(SIGINT, sa, NULL);
#endif

    sigset_t sset[1];
    sigemptyset(sset);
#ifdef HAVE_SIGHUP
    sigaddset(sset, SIGHUP);
#endif
    // sigaddset(sset, SIGINT);
    my_pthread_sigmask(SIG_UNBLOCK, sset, NULL);
#endif  /* HAVE_MINGW */
}
