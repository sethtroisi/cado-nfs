#ifndef SELECT_MPI_H_
#define SELECT_MPI_H_

#if !(defined(_POSIX_C_SOURCE) && _POSIX_C_SOURCE >= 200112L)
#error "Define _POSIX_C_SOURCE to at least 200112L on top of the translation unit"
/* Why ? Because we need it for the prototype of pthread_barrier_wait */
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
#else
#include "fakepthread.h"
#endif

#ifdef WITH_MPI
#include <mpi.h>
#else
#include "fakempi.h"
#endif

#endif	/* SELECT_MPI_H_ */
