#include <stdlib.h>
#include <sys/types.h>		/* for cputime */
#include <sys/resource.h>	/* for cputime */
#include <sys/time.h>	/* for gettimeofday */
#include "timing.h"

uint64_t microseconds()
{
    struct rusage res[1];
    getrusage(RUSAGE_SELF, res);
    uint64_t r;
    r = (uint64_t) res->ru_utime.tv_sec;
    r *= (uint64_t) 1000000UL;
    r += (uint64_t) res->ru_utime.tv_usec;
    return r;
}


/* cputime */
int cputime()
{
    return (int) (microseconds() / (uint64_t) 1000);
}

double seconds()
{
    return (double) microseconds() / 1.0e6;
}

void seconds_user_sys(double * res)
{
    struct rusage ru[1];
    getrusage(RUSAGE_SELF, ru);
    res[0] = ru->ru_utime.tv_sec +  (double) ru->ru_utime.tv_usec / 1.0e6;
    res[1] = ru->ru_stime.tv_sec +  (double) ru->ru_stime.tv_usec / 1.0e6;
}

double wct_seconds()
{
#if 0
    clock_t c = clock();
    return (double) c / CLOCKS_PER_SEC;
#endif
    struct timeval tv[1];
    gettimeofday(tv, NULL);
    return (double) tv->tv_sec + 1.0e-6 * tv->tv_usec;
}

