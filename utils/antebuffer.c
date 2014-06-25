/* Very efficient antebuffer to preempt the file(s) load.
   Syntax: antebuffer X file [file] | YourCommand
   with buffer size = 2^X.
   The files are written on stdout.
   Best size for the buffer is about 4MB (local disk) to
   16MB (NFS), eventually until 128 MB.
   Limitation: X <= 31; more than 2GB antebuffer has no sense anyway.
*/

/* To avoid the warning: implicit declaration of nanosleep for c99 compliant */
#include "cado.h"
#ifdef _POSIX_C_SOURCE
#undef _POSIX_C_SOURCE
#endif
#define _POSIX_C_SOURCE 200809L

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <fcntl.h>
#include <stdint.h>
#include <pthread.h>
#include <errno.h>
#include <time.h>
#include "portability.h"
#include "timing.h"

static pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
static pthread_cond_t cond = PTHREAD_COND_INITIALIZER;

#ifdef HAVE_MINGW
int _CRT_fmode = _O_BINARY; /* Binary open for stdin/out/err */
#endif

#ifndef HAVE_NANOSLEEP
  int nanosleep(const struct timespec *req, struct timespec *rem) {
    if (rem == NULL) {
      /* Dummy to shut up the warning */
    }
#ifdef HAVE_USLEEP
    unsigned long usec = req->tv_sec * 1000000UL + req->tv_nsec / 1000UL;
    usleep(usec);
#else
    sleep(req->tv_sec);
#endif
    return 0;
  }
#endif

/* These variables are used by pthread and main */
static volatile uint64_t ab_cptp = 0, ab_cptc = 0;   /* Main counters */
static volatile int ab_end = 0;                      /* 1 if end */
static char *ab_buf;                                 /* Main buffer */
static const struct timespec waiting = { 0, 1<<13 }; /* About 5 to 20 microseconds */
static int ab_in;                                    /* fd for loading */
static unsigned int ab_size, ab_sizeio;              /* size for buffer & in */

static
void ab_cons () {
  size_t cpy_ab_cptc = 0;
  int w;
  unsigned int c, t = 0;
  
  pthread_setcanceltype (PTHREAD_CANCEL_DEFERRED, NULL);
  for ( ; ; ) {
    pthread_mutex_lock (&mutex);
    while (!(c = (int) (ab_cptp - cpy_ab_cptc))) {
      if (UNLIKELY(ab_end)) {
	pthread_mutex_unlock (&mutex);
	pthread_exit(NULL);
      }
      pthread_cond_wait (&cond, &mutex);
    }
    pthread_mutex_unlock (&mutex);

    if (c > ab_sizeio) c = ab_sizeio;
    if (c > ab_size - t) c = ab_size - t;
    while ((w = (int) write(1, &(ab_buf[t]), (size_t) c)) <= 0) nanosleep(&waiting, NULL);
    cpy_ab_cptc += w;
    t = (t + w) & (ab_size - 1);
    pthread_mutex_lock (&mutex);
    ab_cptc = cpy_ab_cptc; pthread_cond_signal(&cond);
    pthread_mutex_unlock (&mutex);
  }
}

int main(int argc, char **argv) {
  pthread_t ab_tc;
  pthread_attr_t ab_attr;
  size_t cpy_ab_cptp = 0;
  int r;
  unsigned int t, c, p;

#ifdef HAVE_MINGW
  _fmode = _O_BINARY;     /* Binary open for all others files */
#endif

  if (argc < 3) {
  error:
    fprintf (stderr, "%s syntax:  %s SIZE file [file]  or  %s SIZE -\n"
	     "This command is an ante buffer which loads the files or stdin in a\n"
	     "buffer and write them on stdout. The size of this buffer is 2^SIZE.\n"
	     "16<=SIZE<=31, best for 20 to 24.\n", argv[0], argv[0], argv[0]);
    exit (1);
  }
  if ((sscanf(argv[1], "%u", &ab_size)) != 1) goto error;
  if (ab_size < 16 || ab_size > 31) goto error;
  ab_size = 1UL << ab_size;
  ab_sizeio = ab_size >> 4;
  if (!(ab_buf = malloc(ab_size))) {
    fprintf (stderr, "%s: malloc error: %s\n", argv[0], strerror(errno));
    exit(1);
  }
  pthread_attr_init(&ab_attr);
  pthread_attr_setstacksize(&ab_attr, 1<<12);
  pthread_attr_setdetachstate(&ab_attr, PTHREAD_CREATE_JOINABLE);
  if (pthread_create(&ab_tc, &ab_attr, (void *) ab_cons, NULL)) {
    fprintf (stderr, "%s: pthread_create error: %s\n", argv[0], strerror(errno));
    exit(1);
  }
  if (strcmp(argv[2], "-")) {
    p = 3;
    if ((ab_in = open(argv[2], O_RDONLY)) == -1) {
      fprintf (stderr, "%s: open or load error in file %s: %s\n", argv[0], argv[2], strerror(errno));
      exit (1);
    }
  }
  else {
    p = 0;
    ab_in = 0;
  }
  t = 0;
  for ( ; ; ) {
    pthread_mutex_lock (&mutex);
    while (!(c = (int) ((ab_cptc + ab_size) - cpy_ab_cptp))) pthread_cond_wait (&cond, &mutex);
    pthread_mutex_unlock (&mutex);
    
    if (c > ab_sizeio)   c = ab_sizeio;
    if (c > ab_size - t) c = ab_size - t;
    r = read(ab_in, &(ab_buf[t]), (size_t) c);
    if (r > 0) {
      cpy_ab_cptp += r;
      t = (t + r) & (ab_size - 1);
      pthread_mutex_lock (&mutex);
      ab_cptp = cpy_ab_cptp; pthread_cond_signal(&cond);
      pthread_mutex_unlock (&mutex);
    }
    else 
      if (!r) {
	if (!p || p == (unsigned int) argc) {
	  pthread_mutex_lock (&mutex);
	  ab_end = 1; pthread_cond_signal(&cond);
	  pthread_mutex_unlock (&mutex);
	  break;
	}
	close(ab_in);
	if ((ab_in = open(argv[p++], O_RDONLY)) == -1) {
	  fprintf (stderr, "%s: open or load error in file %s: %s\n", argv[0], argv[p-1], strerror(errno));
	  exit (1);
	}
      }
      else {
	fprintf (stderr, "%s: read error in file %s: %s\n", argv[0], argv[p-1], strerror(errno));
	  exit (1);
	}
  }
  if (p) close(ab_in);
  pthread_join(ab_tc, NULL);
  pthread_mutex_destroy (&mutex);
  pthread_cond_destroy (&cond);
  free (ab_buf);
  ab_buf = NULL;
  double tt[2];
  seconds_user_sys(tt);
  fprintf(stderr, "antebuffer exits after having spent %.2fs+%2.fs on cpu\n",
          tt[0], tt[1]);
  exit (0);
}
