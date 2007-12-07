#ifndef TIMING_H_
#define TIMING_H_

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

extern uint64_t microseconds();

/* cputime */
static inline int cputime(void) { return (int) microseconds() / 1000; }
static inline double seconds(void) { return (double) microseconds() /1.0e6; }


#ifdef __cplusplus
}
#endif

#endif	/* TIMING_H_ */
