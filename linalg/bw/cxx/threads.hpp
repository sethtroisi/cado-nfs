#ifndef THREADS_HPP_
#define THREADS_HPP_

#include <boost/cstdint.hpp>
#include <gmp.h>
#include <gmpxx.h>
#include <vector>
#include <istream>
#include <ostream>
#include <iterator>

#include "constants.hpp"


#if 0
	boost::uint32_t	idx_off;
	boost::uint32_t	val_off;
	scalar		src;
	scalar		dst;
	int		iter;
	scalar		s;
	scalar		t;
	int		go_mark;
	int		ncoeffs;
	struct timeval	start_time;
	struct timeval	elapsed_time;
	int		reset_notice;

	/* For A files (dot products), and polyeval results */
	std::vector<std::ostream_iterator<mpz_class> > out;

	/* For poly coeffs */
	std::vector<std::istream_iterator<mpz_class> > in;
#endif

typedef void * (*threadfunc_t)(void *);

#ifdef	ENABLE_PTHREADS
extern unsigned int tseqid(void);
extern bool is_multithreaded;
#else
#define	is_multithreaded	false
#define	tseqid()	0
#endif

extern void configure_threads(int);
extern void start_threads(threadfunc_t, std::vector<void *>&);
inline void start_threads(threadfunc_t fcn, int nt) {
	std::vector<void *> dummy(nt, (void *) NULL);
	start_threads(fcn, dummy);
}

#define is_master_thread() (tseqid()==0)

#endif	/* THREADS_HPP_ */
