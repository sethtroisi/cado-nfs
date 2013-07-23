#ifndef ROLLBUF_H_
#define ROLLBUF_H_

/* Example of a rotating, and reallocating buffer. A separate thread has
 * to fetch the data from the real source and fill the buffer. This way,
 * reads block as little as possible -- sort of a userland pipe, only we
 * get more control. Not clear we gain anything, but it's less fragile
 * than forking an external program.
 */
#include <pthread.h>

/* Getting data from a rolling buffer can be done via two interfaces. One
 * interface (rollbuf_get) provides a pointer, while the other
 * (rollbuf_get2) does not, and asks the rollbuf interface to provide a
 * pointer.
 *
 * The area pointed to by the pointer returned from rollbuf_get2 is free
 * to be read from by the calling thread until the next call to
 * rollbuf_get. This implies in particular that we assume that there is
 * exactly one threading getting data, no more.
 *
 * Note though that this is not implemented by means of a zero-copy
 * mechanism (doing so would make rollback_put break the active area of
 * data being read, in case a realloc() occurs). The data is copied to a
 * buffer exclusively dedicated to reading.
 *
 * (the reading buffer is allocated only if needed, and freed by
 * rollbuf_clear eventually)
 */
struct rollbuf_s {
    char * p;
    size_t alloc;
    size_t avail_to_read;
    size_t avail_to_write;
    const char * rhead;
    char * whead;
    char * rbuf;        /* Only for rollbuf_get2 */
    pthread_mutex_t mx[1];
    pthread_cond_t bored[1];
    int empty_count;
    int full_count;
    int done:1;
};

typedef struct rollbuf_s rollbuf[1];
typedef struct rollbuf_s * rollbuf_ptr;
typedef const struct rollbuf_s * rollbuf_srcptr;

#ifdef __cplusplus
extern "C" {
#endif

extern void rollbuf_init(rollbuf_ptr r);
extern void rollbuf_clear(rollbuf_ptr r);
extern int rollbuf_put(rollbuf_ptr r, char * p, size_t s);
extern void rollbuf_mark_done(rollbuf_ptr r);
extern int rollbuf_is_done(rollbuf_ptr r);
extern int rollbuf_get(rollbuf_ptr r, char * p, size_t s);
extern int rollbuf_get2(rollbuf_ptr r, void ** p, size_t s);

#ifdef __cplusplus
}
#endif

#endif	/* ROLLBUF_H_ */
