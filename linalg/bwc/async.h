#ifndef ASYNC_H_
#define ASYNC_H_

#include "parallelizing_info.h"

/* How long do we typically accept to wait for asynchronous output (=
 * signals) to be checked ? */
#define PREFERRED_ASYNC_LAG     15.0

/* This structure **MUST** be made only of doubles, or disaster will
 * occur */
struct timing_interval_data_s {
    double job[2];  // user, system
    double thread[2];       // user, system
    double wct;
};
typedef struct timing_interval_data_s timing_interval_data[1];

#ifdef  MEASURE_LINALG_JITTER_TIMINGS
#define NTIMERS_EACH_ITERATION  4
#define TIMER_NAMES { "CPU", "cpu-wait", "COMM", "comm-wait", }
#else
#define NTIMERS_EACH_ITERATION  2
#define TIMER_NAMES { "CPU", "COMM", }
#endif

struct timing_data {
    int go_mark;
    int begin_mark;
    int end_mark;
    int last_print;
    int next_print;
    int next_async_check;
    int async_check_period;
    int which;
    timing_interval_data since_last_reset[NTIMERS_EACH_ITERATION];
    timing_interval_data since_beginning[NTIMERS_EACH_ITERATION];
    // since_last_reset[which^1] and since_beginning[which^1] are all
    // complete timings, while the others are running -- in the sense
    // that `running' timings must be added  the contribution of the
    // current timer to constitute something meaningful.
};

#ifdef __cplusplus
extern "C" {
#endif

void timing_init(struct timing_data * t, int start, int end);
void timing_clear(struct timing_data * t);
void timing_next_timer(struct timing_data * t);

void timing_check(parallelizing_info pi, struct timing_data * t, int iter, int print);
void timing_update_ticks(struct timing_data * t, int iter);
void timing_disp_collective_oneline(parallelizing_info pi, struct timing_data * timing, int iter, unsigned long ncoeffs, int print, int stage);
void timing_final_tally(parallelizing_info pi, struct timing_data * timing, unsigned long ncoeffs, int print, int stage);

void block_control_signals();
void catch_control_signals();

#ifdef __cplusplus
}
#endif

#endif	/* ASYNC_H_ */
