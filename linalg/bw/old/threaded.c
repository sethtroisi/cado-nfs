#include "params.h"

#if defined(ENABLE_PTHREADS) && defined(HAS_REALTIME_CLOCKS)
/* For pthread_getcpuclockid */
#define	_POSIX_C_SOURCE	200112L
#endif

#include <sys/time.h>
#include <sys/types.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "slave_params.h"
#include <assert.h>
#include "types.h"
#include "macros.h"
#include "variables.h"
#include "bw_scalar.h"
#include "bw_lvblock.h"
#include "slave.h"
#include "timer.h"
#include "threaded.h"
#include "auxfuncs.h"

#ifdef ENABLE_PTHREADS

#include <pthread.h>
#include <unistd.h>
unsigned int	number_of_threads;
unsigned int is_multithreaded=0;
struct tsv * tsv_table = NULL;

void configure_threads(int n, coord_t cols)
{
	int i;
	number_of_threads=n;
	if (tsv_table)
		free(tsv_table);
	tsv_table=my_malloc(n * sizeof(struct tsv));
	for(i=0;i<n;i++) {
		tsv_table[i].i0=i*cols/n;
		tsv_table[i].i1=(i+1)*cols/n;
		printf("Thread %d handles indexes i, %d<=i<%d, %d in total\n",
				i,tsv_table[i].i0,tsv_table[i].i1,
				tsv_table[i].i1-tsv_table[i].i0);
	}
	for(i=0;i<n;i++) {
		tsv_table[i].x0=i*nbxs/n;
		tsv_table[i].x1=(i+1)*nbxs/n;
		printf("Thread %d handles indexes x, %d<=x<%d, %d in total\n",
				i,tsv_table[i].x0,tsv_table[i].x1,
				tsv_table[i].x1-tsv_table[i].x0);
	}
		
}

struct thread_starter_info {
	threadfunc_t	f;
	unsigned int	i;
	int		iterations;
};

static pthread_once_t id_alloc_once = PTHREAD_ONCE_INIT;

static pthread_key_t id_alloc_key;

static void id_destroy(void *buf)
{
	if (buf)
		free(buf);
	buf=NULL;
}

static void id_alloc (void)
{
	pthread_key_create (&id_alloc_key, &id_destroy);
}

static void * thread_starter(struct thread_starter_info * p)
{
	int rc;

	pthread_once(&id_alloc_once, &id_alloc);
	pthread_setspecific(id_alloc_key, my_malloc(sizeof(unsigned int)));
	*(unsigned int*)pthread_getspecific(id_alloc_key)=p->i;
	printf("Thread %d : pid %d\n",p->i,getpid());
#ifdef GROSBUG
	tsv()->dbstr=fdopen(3+p->i,"w");
#endif
#ifdef	HAS_REALTIME_CLOCKS
	rc = pthread_getcpuclockid(pthread_self(), &(tsv()->cpuclock));
	if (rc != 0) {
		fprintf(stderr, "pthread_getcpuclockid failed, rc=%d\n", rc);
		exit(1);
	}
#endif
	return (*(p->f))(&(p->iterations));
}

unsigned int tseqid(void)
{
	if (is_multithreaded)
		return *(unsigned int*)pthread_getspecific(id_alloc_key);
	else
		return 0;
}

int is_a_working_thread(void)
{
	if (is_multithreaded)
		return pthread_getspecific(id_alloc_key)!=NULL;
	else 
		return 1;
}

double start_threads(void* (*fcn) (void*), int iterations)
{
	int i;
	pthread_t * tid;
	int retcode;
	double t,t0;
	struct thread_starter_info * tsi;

	if (iterations==0) {
		printf("Are you kidding me ?\n");
		return 0.0;
	}

	/* XXX ! Don't think about using one single struct for the tsi.
	 * The data would be shared between the different invocations of
	 * thread_starter(), for a really ugly result (threads sharing
	 * their tseqid(), for instance...)
	 */

	t=t0=timer(TIMER_ASK);
	tid=my_malloc(number_of_threads*sizeof(pthread_t));
	tsi=my_malloc(number_of_threads*sizeof(struct thread_starter_info));
	is_multithreaded=1;
	for(i=0;i<number_of_threads;i++) {
		tsi[i].i=i;
		tsi[i].f=fcn;
		tsi[i].iterations=iterations;
		retcode=pthread_create(&(tid[i]), NULL,
				(threadfunc_t) &thread_starter,
				(void*) &(tsi[i]));
		switch(retcode) {
			case 0: printf("Successfully started "
						"thread %d (tid %lx)\n",
						i,tid[i]);
				break;
			case EAGAIN:
				fprintf(stderr,
					"Too many active threads (%d)\n", i);
				exit(retcode);
			default:
				fprintf(stderr,"pthread_create : %s\n",
						strerror(retcode));
				exit(retcode);
		}
	}
	t=0.0;
	for(i=0;i<number_of_threads;i++) {
		void * tret;
		time_t now;
		struct tm loc_now;
		char now_string[80];
		double dt;

		/* thread functions return a malloced area, containing
		 * (in my case) a double.
		 */
		retcode=pthread_join(tid[i],(void**)&tret);
		if (tret != NULL) {
			dt=(*(double *)tret);
		} else {
			dt=0;
		}
		t+=dt;
		free(tret);
		time(&now);
		localtime_r(&now,&loc_now);
		strftime(now_string,80,"%d %b %T", &loc_now);
		printf("Thread %d finished on %s\n",i,now_string);
		switch(retcode) {
			case 0: 
#if 0
				/* FPE happens here on MT. Compiler bug ? */
				printf("Thread %d (tid %lx) finished. "
				"Spent %.2fs, avg %.3fs\n",
				i,tid[i],(*tret),
				(*tret)/iterations);
#else
				printf("Thread %d (tid %lx) finished. "
				"Spent %.2fs, %d iterations\n",
				i,tid[i],dt,iterations);
#endif
		
				break;
			case ESRCH:
				fprintf(stderr,"Thread %d cannot be found\n",i);
				break;
			default:
				fprintf(stderr,"pthread_join : %s\n",
						strerror(retcode));
				exit(retcode);
		}
	}
	is_multithreaded=0;
	free(tid);
	free(tsi);
	return t;
}

#else	/* ENABLE_PTHREADS */

struct tsv tsv_table[1];

void configure_threads(int n, coord_t cols)
{
	if (n > 1) {
		fprintf(stderr,"This program has threads disabled !\n");
		exit(EINVAL);
	}
	tsv()->i0=0;
	tsv()->i1=cols;
	tsv()->x0=0;
	tsv()->x1=nbxs;
	tsv()->idx_off=0;
	tsv()->val_off=0;
}

double start_threads(void* (*fcn) (void*), int n UNUSED_VARIABLE)
{
	double *tret;
	double t;
	tret=(*fcn)(NULL);
	t=*tret;
	free(tret);
	return t;
}

#endif
