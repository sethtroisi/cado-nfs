#ifndef SLAVE_PARAMS_H_
#define SLAVE_PARAMS_H_

#define	SMALL_CHECK

#define	REALLY_WANT_ADDMUL_INLINES

/* Here come all the specific parameters, designed only for the slave
 * thread */

/* Order vector blocks so that bytes from the same level but from
 * different blocks will be contiguous. If this option is not defined,
 * the bytes from the same big integer will be contiguous instead. The
 * first setting accelerates the matrix multiplication phase if the MMX
 * optimization is possible, while the second is probably better for
 * modular reduction. This flags is meaningless when nbys==1. */
/* #define BLOCKS_TOGETHER */

/* In order to help the compiler in his optimizations, you might want to
 * hard-code most of the parameters in here. Note that it does not waive
 * you from the task of specifying them on the command line : The program
 * wants to check them for consistency */

/*
#define HARD_bw_allocsize 3
#define HARD_bw_longsize 10
#define HARD_m_param	2
#define HARD_n_param	2
#define HARD_nbxs	2
#define HARD_nbys	1
*/
#if 0
#ifdef __alpha
#define HARD_bw_allocsize 5
#define HARD_bw_longsize 14
#define HARD_m_param	8
#define HARD_n_param	8
#define HARD_nbxs	8
#define HARD_nbys	1
#endif
#ifdef __i386__
#define HARD_bw_allocsize 10
#define HARD_bw_longsize 24
#define HARD_m_param	8
#define HARD_n_param	8
#define HARD_nbxs	8
#define HARD_nbys	1
#endif
#endif
#ifndef SMALL_CHECK
#define HARDCODE_PARAMS
#if 0
#if defined(__i386__) || defined(__sparc)
#define HARD_bw_allocsize 19
#define HARD_m_param	4
#define HARD_n_param	4
#else
#define HARD_bw_allocsize 10
#define HARD_m_param	4
#define HARD_n_param	4
#endif
#endif
#define HARD_bw_allocsize	2
#define HARD_m_param	1
#define HARD_n_param	1

#define HARD_nbxs	HARD_m_param
#define HARD_nbys	1
#define HARD_bw_longsize	(2*HARD_bw_allocsize+4)
#endif

/* Default number of threads to spawn when multithreading... */
#define NTHREADS	2

#define	DEFAULT_KEEP	10


#endif /* SLAVE_PARAMS_H_ */
