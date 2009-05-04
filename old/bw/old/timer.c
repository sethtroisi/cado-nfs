#include <stdio.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <time.h>
#include "params.h"
#include "timer.h"

#include "types.h"
/* #include "threaded.h" */

/* linuxthreads : thread == self, actually
 *
 * This is not true with nptl.
 */ 

double timer_r(struct timeval *p, int param)
{
	double res=0.;
	struct rusage now;

#if defined(ENABLE_PTHREADS) && defined(HAS_REALTIME_CLOCKS)
	/* We avoid getrusage completely. However we fill the timeval
	 * struct. */
	struct timespec tm;
	clock_gettime(tsv()->cpuclock, &tm);
	now.ru_utime.tv_sec = tm.tv_sec;
	now.ru_utime.tv_usec = tm.tv_nsec / 1000;
#elif defined(HAS_RUSAGE_THREAD)
	if (param & TIMER_MT)
		getrusage(RUSAGE_THREAD,&now);
	else
		getrusage(RUSAGE_SELF,&now);
#else	/* default, or LINUXTHREADS case */
	getrusage(RUSAGE_SELF,&now);
#endif

	if (param & TIMER_ASK)
		res=(now.ru_utime.tv_sec-p->tv_sec)+
			(now.ru_utime.tv_usec-p->tv_usec)/1000000.;

	if (param & TIMER_SET) {
		/* printf("!! TIMER RESET !!\n"); */
		*p=now.ru_utime;
	}

#if	defined(LINUXTHREADS) ||	\
	defined(HAS_RUSAGE_THREAD) ||	\
	(defined(ENABLE_PTHREADS) && defined(HAS_REALTIME_CLOCKS))
	
	/* Our timing method gives thread-specific times in all cases */
	return res;
#else
	/* This is a poor man's equivalent. The global time is averaged */
	if (param & TIMER_MT)
		return res / TIMER_NT(param);
	else
		return res;
#endif
}

/* This is bad, not reentrant and so on. Its results cannot be
 * accumulated because of the inherent inaccuracy of the measure.
 * Use timer_r instead !
 */
double timer(int param)
{
	static struct timeval last = {0, };
	return timer_r(&last,param);
}

double wall_clock_r(struct timeval *p, int param)
{
	struct timeval now;
	struct timezone tz;	/* total dummy storage... */
	double res=0.;

	gettimeofday(&now,&tz);
	if (param & TIMER_ASK)
		res=(now.tv_sec-p->tv_sec)+
			(now.tv_usec-p->tv_usec)/1000000.;
	if (param & TIMER_SET)
		*p=now;

	return res;
}

double wall_clock(int param)
{
	static struct timeval last = {0, };
	return wall_clock_r(&last,param);
}

