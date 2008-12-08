#ifdef ENABLE_PTHREADS
#include <sys/types.h>
#include <sys/time.h>
#include <errno.h>
#include <pthread.h>
#include <time.h>

#include "barrier.h"
#include "manu.h"

#if GNUC_VERSION_ATLEAST(4,3,0)
#pragma GCC diagnostic ignored "-Wempty-body"
/* It turns out that the newer pthread.h header from Fedora F9
 * glibc-headers-2.8-3 has a wonderful ``do; while (0)'' statement for
 * the definition of pthread_cleanup_pop
 */
#endif

int barrier_init (barrier_t *barrier, int count, int *cstate)
{
	int status;

	barrier->threshold = barrier->counter = count;
	barrier->cycle = 0;
	status = pthread_mutex_init (&barrier->mutex, NULL);
	if (status != 0)
		return status;
	status = pthread_cond_init (&barrier->cv, NULL);
	if (status != 0) {
		pthread_mutex_destroy (&barrier->mutex);
		return status;
	}
	barrier->result= 0x12345678;
	barrier->valid = BARRIER_VALID;
	if (cstate) {
		barrier->cancelstate=*cstate;
	} else {
		barrier->cancelstate=PTHREAD_CANCEL_DISABLE;
	}

	return 0;
}

int barrier_destroy (barrier_t *barrier)
{
	int status, status2;

	if (barrier->valid != BARRIER_VALID)
		return EINVAL;

	status = pthread_mutex_lock (&barrier->mutex);
	if (status != 0)
		return status;

	/*
	 * Check whether any threads are known to be waiting; report
	 * "BUSY" if so.
	 */
	if (barrier->counter != barrier->threshold) {
		pthread_mutex_unlock (&barrier->mutex);
		return EBUSY;
	}

	barrier->valid = 0;
	status = pthread_mutex_unlock (&barrier->mutex);
	if (status != 0)
		return status;

	/*
	 * If unable to destroy either 1003.1c synchronization
	 * object, return the error status.
	 */
	status = pthread_mutex_destroy (&barrier->mutex);
	status2 = pthread_cond_destroy (&barrier->cv);
	return (status != 0 ? status : status2);
}

/*
 * Once all threads waiting from inside a barrier are canceled, the
 * barrier should end up being empty.
 */
void barrier_cleanup_handler(void *b)
{
	barrier_t * bar=(barrier_t*)b;
	bar->counter++;
	pthread_mutex_unlock(&bar->mutex);
}


int barrier_wait (barrier_t *barrier, double * pwait, int (*all_in)(int))
{
	struct timeval tv;
	int result = 0, cancel, tmp;
        volatile int status;
        unsigned int cycle;
	int func_retcode=0; /* Avoid see gcc moaning because it might
			       be used unintialized (no, it won't) */

	if (pwait) gettimeofday(&tv,NULL);

	if (barrier->valid != BARRIER_VALID)
		return -EINVAL;

	status = pthread_mutex_lock (&barrier->mutex);
	if (status != 0)
		return -status;

	cycle = barrier->cycle;   /* Remember which cycle we're on */

	if (all_in)
	func_retcode = (*all_in)(barrier->threshold-barrier->counter);
	if (--barrier->counter == 0) {
		if (all_in) {
			result = barrier->result = func_retcode;
		}

		barrier->cycle++;
		barrier->counter = barrier->threshold;

		status = pthread_cond_broadcast (&barrier->cv);
		/*
		 * The last thread into the barrier will return status
		 * -1 rather than 0, so that it can be used to perform
		 * some special serial code following the barrier.
		 */
		if (status == 0)
			status = -1;
	} else {
		pthread_setcancelstate (barrier->cancelstate, &cancel);
		pthread_cleanup_push(&barrier_cleanup_handler,barrier);

		/*
		 * Wait until the barrier's cycle changes, which means
		 * that it has been broadcast, and we don't want to wait
		 * anymore.
		 */
		while (cycle == barrier->cycle) {
			status = pthread_cond_wait (
					&barrier->cv, &barrier->mutex);
			if (status != 0) break;
		}

		pthread_cleanup_pop(0);
		pthread_setcancelstate (cancel, &tmp);

		/* Now I hold the mutex, and I've been signaled. It's
		 * time to read barrier->result, since I know it's just
		 * been set by the thread who signaled me.
		 */
		if (all_in)
			result=barrier->result;
	}

	/* See how long we have been waiting */
	if (pwait) {
		struct timeval now;
		gettimeofday(&now,NULL);
		*pwait = (now.tv_sec - tv.tv_sec) +
			(now.tv_usec - tv.tv_usec) / 1.0e6;
	}
	/*
	 * Ignore an error in unlocking. It shouldn't happen, and
	 * reporting it here would be misleading -- the barrier wait
	 * completed, after all, whereas returning, for example,
	 * EINVAL would imply the wait had failed. The next attempt
	 * to use the barrier *will* return an error, or hang, due
	 * to whatever happened to the mutex.
	 */
	pthread_mutex_unlock (&barrier->mutex);
	if (all_in)
		return result;
	else
		return -status;          /* error, 1 for waker, or 0 */
}

#else	/* ENABLE_PTHREADS */

#include <stdlib.h>
#include "macros.h"
#include "barrier.h"

int barrier_wait(barrier_t * p MAYBE_UNUSED,
		double * pwait,
		int (*code)(int))
{
	if (pwait)
		*pwait=0.0;
	if (code != NULL)
		return (*code)(0);
	else
		return 1;
}

int barrier_init (barrier_t *barrier MAYBE_UNUSED,
		int count,
		int *p MAYBE_UNUSED)
{
	return (count==1);
}

int barrier_destroy (barrier_t *barrier MAYBE_UNUSED)
{
	return 0;
}


#endif	/* ENABLE_PTHREADS */
