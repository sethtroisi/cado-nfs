#define _GNU_SOURCE         /* asprintf */
#define _DARWIN_C_SOURCE    /* for asprintf. _ANSI_SOURCE must be undefined */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <stdarg.h>

#include <sys/types.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <time.h>

#include "bwc_config.h"
#include "macros.h"

#if defined(__linux) && defined(__GLIBC__)
#include <linux/version.h>
#include <features.h>
#if LINUX_VERSION_CODE >= 0x20611 && LEXGE2(__GLIBC__,__GLIBC_MINOR__,2,9)
#define HAVE_GETRUSAGE_THREAD
#endif
#endif

void job_seconds(double * res)
{
    struct rusage ru[1];
    getrusage(RUSAGE_SELF, ru);
    res[0] = ru->ru_utime.tv_sec +  (double) ru->ru_utime.tv_usec / 1.0e6;
    res[1] = ru->ru_stime.tv_sec +  (double) ru->ru_stime.tv_usec / 1.0e6;
}

double walltime_seconds()
{
#if 0
    clock_t c = clock();
    return (double) c / CLOCKS_PER_SEC;
#endif
    struct timeval tv[1];
    gettimeofday(tv, NULL);
    return (double) tv->tv_sec + 1.0e-6 * tv->tv_usec;
}

/* We need some way to detect the time spent by threads. Unfortunately,
 * this is something not defined by POSIX.
 */

/* Only on linux >= 2.6.26, and IRIX, solaris. */
#if defined(HAVE_GETRUSAGE_THREAD)
void thread_seconds(double * res)
{
    struct rusage ru[1];
    getrusage(RUSAGE_THREAD, ru);
    res[0] = ru->ru_utime.tv_sec +  (double) ru->ru_utime.tv_usec / 1.0e6;
    res[1] = ru->ru_stime.tv_sec +  (double) ru->ru_stime.tv_usec / 1.0e6;
}
#elif defined(__linux)
#include <sys/syscall.h>

static inline pid_t gettid() { return syscall(SYS_gettid); }

/* On linux (well, using nptl and not linuxthreads, but we don't care
 * much), that's doable by parsing /proc/<tid>/stat
 */

/* statfields is fairly crude, for sure. ... has to be terminated by -1.
 * Arguments must come in groups of three, first the field index, then
 * its parsing format, then the adress of the pointer.
 */
static int statfields(pid_t t, ...)
{
    int nparsed = 0;
    va_list ap;
    va_start(ap, t);
    char * tmp;
    int rc = asprintf(&tmp,"/proc/%d/stat",t);
    if (rc < 0) return 0;
    char buf[1024];
    FILE * f = fopen(tmp,"r");
    char * s = fgets(buf, sizeof(buf), f);
    if (s == NULL) return 0;
    int j = va_arg(ap, int);
    for(int i = 0 ; j != -1 ; i++) {
        char * next;
        next = strchr(s, ' ');
        if (next == NULL)
            break;
        *next++='\0';
        if (i == j) {
            const char * fmt = va_arg(ap, const char *);
            void * ptr = va_arg(ap, void *);
            nparsed += sscanf(s, fmt, ptr);
            j = va_arg(ap, int);
        }
        s = next;
    }

    fclose(f);
    free(tmp);
    va_end(ap);

    return nparsed;
}

void thread_seconds(double * res)
{
    unsigned long utime = 0;
    unsigned long stime = 0;
    statfields(gettid(), 13, "%lu", &utime, 14, "%lu", &stime, -1);
    res[0] = (double) utime / sysconf(_SC_CLK_TCK);
    res[1] = (double) stime / sysconf(_SC_CLK_TCK);
}

#else   /* Otherwise we'll do something stupid */
void thread_seconds(double * res)
{
    job_seconds(res); /* really stupid */
}
#endif

