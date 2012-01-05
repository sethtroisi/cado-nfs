#ifndef CADO_UTILS_TIMING_H_
#define CADO_UTILS_TIMING_H_

#include <stdint.h> /* for uint64_t */

#ifdef __cplusplus
extern "C" {
#endif

extern uint64_t microseconds (void);
extern int cputime (void);
extern double seconds (void);
extern void seconds_user_sys (double *);
extern double wct_seconds (void);
extern void print_timing_and_memory (double);

#ifdef __cplusplus
}
#endif

#endif	/* CADO_UTILS_TIMING_H_ */
