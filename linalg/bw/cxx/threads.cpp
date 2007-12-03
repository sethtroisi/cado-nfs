#include "manu.h"

#include <cerrno>
#ifdef ENABLE_PTHREADS
#include <pthread.h>
#endif

#include "barrier.h"

#include "constants.hpp"
#include "fmt.hpp"
#include "threads.hpp"
#include "ticks.hpp"

#include <iostream>

using namespace std;

void configure_threads(int n)
{
#ifndef	ENABLE_PTHREADS
	BUG_ON(n > 1);
#endif
}

static std::string report(const std::pair<double, double>& ticks)
{
	ostringstream s;
	s << fmt("%[F.2]s on cpu, %[F.2]s total ; ratio %[F.1]%%")
		% ticks.first
		% ticks.second
		% (ticks.second > 1.0e-3
				? (100.0 * (ticks.first / ticks.second))
				: 0.0);
	return s.str();
}

/**************************************************************************/
#ifdef ENABLE_PTHREADS

bool is_multithreaded = false;

static pthread_once_t id_alloc_once = PTHREAD_ONCE_INIT;
static pthread_key_t id_alloc_key;
static void id_alloc(void)
{
	pthread_key_create(&id_alloc_key, &free);
}
unsigned int tseqid(void)
{
	if (is_multithreaded)
		return *(unsigned int *) pthread_getspecific(id_alloc_key);
	else
		return 0;
}

barrier_t       thread_group_barrier;

struct thread_starter_info {
	unsigned int i;
	threadfunc_t f;
	void * arg;
	std::pair<double, double> ticks;
};

static void *thread_starter(struct thread_starter_info *p)
{
	pthread_once(&id_alloc_once, &id_alloc);
	pthread_setspecific(id_alloc_key, malloc(sizeof(unsigned int)));
	*(unsigned int *) pthread_getspecific(id_alloc_key) = p->i;

	cout << fmt("// thread % : tid 0x%[h]")
		% p->i % pthread_self()
#ifdef	__linux__
		<< " ltid " << linux_tid()
#endif
		<< "\n";

	barrier_wait(&thread_group_barrier, NULL, NULL);
	p->ticks.first  = -thread_ticks();
	p->ticks.second = -wallclock_ticks();

	void * res = p->f (p->arg);

	barrier_wait(&thread_group_barrier, NULL, NULL);
	p->ticks.first  += thread_ticks();
	p->ticks.second += wallclock_ticks();

	return res;
}

void start_threads(threadfunc_t fcn, std::vector<void *>& arg)
{
	int nt = arg.size();
	vector<pthread_t> tid(nt);
	vector<thread_starter_info> tsi(nt);
	double start_time = - wallclock_ticks();

	cout.flush();
	cerr.flush();

	barrier_init(&thread_group_barrier, nt + 1, NULL);

	is_multithreaded = true;
	for (int i = 0; i < nt; i++) {
		tsi[i].i = i;
		tsi[i].f = fcn;
		tsi[i].arg = arg[i];
		int retcode = pthread_create(&(tid[i]), NULL,
					 (threadfunc_t) & thread_starter,
					 (void *) &(tsi[i]));
		BUG_ON(retcode != 0);
	}

	barrier_wait(&thread_group_barrier, NULL, NULL);

	/* Here the actual computation goes. _we_ do nothing at all.
	 * Perhaps our job could be dispatching signals. */

	barrier_wait(&thread_group_barrier, NULL, NULL);
	
	std::pair<double, double> total_ticks(0,0);
	for (int i = 0; i < nt; i++) {

		/* thread functions return a malloc'ed area, containing
		 * (in my case) a pair of doubles. */
		int retcode = pthread_join(tid[i], &arg[i]);

		if (retcode == ESRCH) {
			cerr << fmt("// thread % tid 0x%[h] vanished\n")
				% i % tid[i];
			continue;
		}
		BUG_ON(retcode != 0);

		total_ticks.first += tsi[i].ticks.first;

		{
			time_t now;
			struct tm loc_now;
			char now_string[80];

			time(&now);
			localtime_r(&now, &loc_now);
			strftime(now_string, 80, "%d %b %T", &loc_now);
			cout << fmt("// % thread % done %\n")
				% now_string
				% i
				% report(tsi[i].ticks);
		}
	}
	total_ticks.second = start_time + wallclock_ticks();
	is_multithreaded = false;
	cout << fmt("// cumulative: %\n") % report(total_ticks);
}

#else

/* ENABLE_PTHREADS */

void start_threads(threadfunc_t fcn, std::vector<void *>& arg)
{
	std::pair<double, double> ticks;
	ticks.first  = -thread_ticks();
	ticks.second = -wallclock_ticks();

	BUG_ON(arg.size() != 1);

	arg[0] = (*fcn) (arg[0]);

	ticks.first  += thread_ticks();
	ticks.second += wallclock_ticks();
	cout << fmt("// spent %\n") % report(ticks);
}

#endif

