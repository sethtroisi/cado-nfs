#include <sys/types.h>    /* for cputime */
#include <sys/resource.h> /* for cputime */
#include <stdint.h>

uint64_t microseconds()
{
        struct rusage res[1];
        getrusage(RUSAGE_SELF,res);
        uint64_t r;
        r  = (uint64_t) res->ru_utime.tv_sec;
        r *= (uint64_t) 1000000UL;
        r += (uint64_t) res->ru_utime.tv_usec;
        return r;
}
