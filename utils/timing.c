#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>		/* for cputime */
#include <sys/resource.h>	/* for cputime */
#include <sys/time.h>	/* for gettimeofday */
#include "timing.h"
#include "memusage.h"

uint64_t
microseconds (void)
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
int
cputime (void)
{
    return (int) (microseconds() / (uint64_t) 1000);
}

double
seconds (void)
{
    return (double) microseconds() / 1.0e6;
}

void
seconds_user_sys (double * res)
{
    struct rusage ru[1];
    getrusage(RUSAGE_SELF, ru);
    res[0] = ru->ru_utime.tv_sec +  (double) ru->ru_utime.tv_usec / 1.0e6;
    res[1] = ru->ru_stime.tv_sec +  (double) ru->ru_stime.tv_usec / 1.0e6;
}

double
wct_seconds (void)
{
    struct timeval tv[1];
    gettimeofday (tv, NULL);
    return (double) tv->tv_sec + 1.0e-6 * tv->tv_usec;
}

void
print_timing_and_memory (void)
{
  fprintf (stderr, "Memory usage: VIRT %luM (peak %luM)\n",
           Memusage () >> 10, PeakMemusage () >> 10);
  fprintf (stderr, "Total time: %1.0f seconds (cpu), %1.0f seconds (real)\n",
           seconds (), wct_seconds ());
}

