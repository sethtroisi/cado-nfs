#ifndef RUSAGE_H_
#define RUSAGE_H_

#ifdef __cplusplus
extern "C" {
#endif

extern void thread_seconds_user_sys(double *);

extern double walltime_seconds();

#ifdef __cplusplus
}
#endif

#endif	/* RUSAGE_H_ */
