#ifndef CADO_UTILS_TIMING_H_
#define CADO_UTILS_TIMING_H_

#ifdef __cplusplus
extern "C" {
#endif

extern uint64_t microseconds();
extern int cputime();
extern double seconds();

#ifdef __cplusplus
}
#endif

#endif	/* CADO_UTILS_TIMING_H_ */
