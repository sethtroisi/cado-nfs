#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <fstream>
#include <sstream>
#include <istream>
#include <ostream>
#include "manu.h"
#include "must_open.hpp"
#include "threads.hpp"
#include "ticks.hpp"

#ifdef	__linux__
#include <sys/types.h>
#include <unistd.h>
#include <sys/syscall.h>

unsigned long linux_tid()
{
	return syscall(SYS_gettid);
}
#endif

double oncpu_ticks()
{
	struct rusage now;
	getrusage(RUSAGE_SELF,&now);
	return now.ru_utime.tv_sec + now.ru_utime.tv_usec / 1.0e6;
}

double thread_ticks(int nt MAYBE_UNUSED)
{
#ifdef	__linux__
	std::ostringstream s;
	s << fmt("/proc/%/task/%/stat") % getpid() % linux_tid();
	s.flush();
	std::ifstream statfile;
	must_open(statfile, s.str());
	unsigned long utime;
	{
		std::string foo;
		/* go straight to the interesting field */
		for(int i = 0 ; i < 13 ; i++) statfile >> foo;
		statfile >> utime;
		BUG_ON(statfile.fail());
	}
	return utime / (double) sysconf(_SC_CLK_TCK);
#else
	return oncpu_ticks() / nt;
#endif
}


double wallclock_ticks()
{
	struct timeval now;
	struct timezone tz;	/* total dummy storage... */

	gettimeofday(&now,&tz);
	
	return now.tv_sec + now.tv_usec / 1.0e6;
}
