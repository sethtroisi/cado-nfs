#ifndef ASYNC_H_
#define ASYNC_H_

#include "parallelizing_info.h"

/* How long do we typically accept to wait for asynchronous output (=
 * signals) to be checked ? */
#define PREFERRED_ASYNC_LAG     15.0

struct timing_interval_data_s {
    double job[2];  // user, system
    double thread[2];       // user, system
    double wct;
};
typedef struct timing_interval_data_s timing_interval_data[1];

struct timing_data {
    int go_mark;
    int begin_mark;
    int end_mark;
    int last_print;
    int next_print;
    int next_async_check;
    int async_check_period;
    int which;  // 0: general timer ; 1: cpu bound
    timing_interval_data since_last_reset[2];
    timing_interval_data since_beginning[2];
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
void timing_flip_timer(struct timing_data * t);

void timing_check(parallelizing_info pi, struct timing_data * t, int iter, int print);
void timing_update_ticks(struct timing_data * t, int iter);
void timing_disp_collective_oneline(parallelizing_info pi, struct timing_data * timing, int iter, unsigned long ncoeffs, int print, int stage);
void timing_final_tally(const char *name, parallelizing_info pi, struct timing_data * timing, unsigned long ncoeffs, int print);

void block_control_signals();
void catch_control_signals();

#ifdef __cplusplus
}
#endif

#endif	/* ASYNC_H_ */
