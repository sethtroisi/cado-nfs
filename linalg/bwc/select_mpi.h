#ifndef SELECT_MPI_H_
#define SELECT_MPI_H_
#include "bwc_config.h"
#if !(defined(_POSIX_C_SOURCE) && _POSIX_C_SOURCE >= 200112L)
#include "macros.h"
#if _POSIX_C_SOURCE == 199506L && defined(_GNU_SOURCE) && defined(__GLIBC__) && LEXLE2(__GLIBC__,__GLIBC_MINOR__,2,4)
/* With glibc 2.4, if _GNU_SOURCE is defined, then _POSIX_C_SOURCE is set
 * to 199506L inconditionally, and there's nothing we can do about it
 * (see /usr/include/features.h). Then the compilation will fail on
 * missing prototypes, and we have to resort to kludges.
 * Note that the update was pushed before the glibc 2.4 release, but some
 * broken pre-2.4 are floating around with the patch unapplied (ubuntu).
 * http://sourceware.org/cgi-bin/cvsweb.cgi/libc/include/features.h?rev=1.43&content-type=text/x-cvsweb-markup&cvsroot=glibc
 * */
#ifndef __USE_XOPEN2K
#warning "hacking glibc 2.6 for _POSIX_C_SOURCE. See comment"
#define __USE_XOPEN2K
#endif
#undef  _POSIX_C_SOURCE
#else
/* Otherwise, the programmer is to blame at first, so we prefer to have
 * the source file fixed. */
#error "Define _POSIX_C_SOURCE to at least 200112L on top of the translation unit"
#endif
/* Why all this ? Because we need it for the prototype of
 * pthread_barrier_wait */
#define _POSIX_C_SOURCE 200112L
#endif

// for some absurd reason, the gnu standard C++ library <string> header
// includes pthread. And pollutes the global namespace with it. So we
// can't afford having our pthread proxies, we have to resort to
// alternative naming...

#ifdef WITH_PTHREADS
#include <pthread.h>
#define my_pthread_t                       pthread_t
#define my_pthread_create                  pthread_create
#define my_pthread_join                    pthread_join
#define my_pthread_barrier_t               pthread_barrier_t
#define my_pthread_barrier_init            pthread_barrier_init
#define MY_PTHREAD_BARRIER_SERIAL_THREAD   PTHREAD_BARRIER_SERIAL_THREAD
#define my_pthread_barrier_wait            pthread_barrier_wait
#define my_pthread_barrier_destroy         pthread_barrier_destroy
#define my_pthread_mutex_t                 pthread_mutex_t
#define my_pthread_mutex_init              pthread_mutex_init
#define my_pthread_mutex_lock              pthread_mutex_lock
#define my_pthread_mutex_unlock            pthread_mutex_unlock
#define my_pthread_mutex_trylock           pthread_mutex_trylock
#define my_pthread_mutex_destroy           pthread_mutex_destroy
#define my_pthread_sigmask                 pthread_sigmask
#else
#include "fakepthread.h"
#endif

#ifdef WITH_MPI
#include <mpi.h>
#else
#include "fakempi.h"
#endif

#endif	/* SELECT_MPI_H_ */
