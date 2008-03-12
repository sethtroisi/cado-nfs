#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
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

double seconds()
{
    return (double) microseconds() / 1.0e6;
}

uint64_t wallclock_microseconds()
{
    struct timeval tv;
    uint64_t r;
    gettimeofday(&tv, NULL);

    r = (uint64_t) tv.tv_sec;
    r *= (uint64_t) 1000000UL;
    r += (uint64_t) tv.tv_usec;
    return r;
}


double wallclock_seconds()
{
    return (double) wallclock_microseconds() / 1.0e6;
}
