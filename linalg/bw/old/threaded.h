#ifndef BW_THREADED_H_
#define BW_THREADED_H_

#ifdef	__cplusplus
extern "C" {
#endif

/* The thread-specific variables */
struct tsv {
	coord_t 	i0;
	coord_t		i1;
	size_t		idx_off;
	size_t		val_off;
	unsigned int	x0;
	unsigned int	x1;
	bw_vector_block	src;
	bw_vector_block	dst;
	int		iter;
	bw_scalar	s;
	bw_scalar	t;
	int		go_mark;
	int		ncoeffs;
#ifdef	HAS_REALTIME_CLOCKS
	clockid_t	cpuclock;
#endif
	struct timeval	start_time;
	struct timeval	elapsed_time;
	int		reset_notice;	/* Only for start times... */
	FILE	     ** aux_files;	/* for coeffs in mksol */
#ifdef GROSBUG
	FILE *		dbstr;
#endif
};

typedef void * (*threadfunc_t)(void *);
extern void configure_threads(int, coord_t);
extern double start_threads(threadfunc_t, int);

#ifdef	ENABLE_PTHREADS
extern unsigned int tseqid(void);

extern unsigned int number_of_threads;
extern unsigned int is_multithreaded;
extern int is_a_working_thread(void);
extern struct tsv * tsv_table;
#else
extern struct tsv tsv_table[1];
#define	tseqid()		0
#define	number_of_threads	1
#define	is_multithreaded	0
#define	is_a_working_thread()	1
#endif	/* ENABLE_PTHREADS */

#define is_master_thread() (tseqid()==0)
#define tsv() (&(tsv_table[tseqid()]))

#ifdef	__cplusplus
}
#endif

#endif	/* BW_THREADED_H_ */
