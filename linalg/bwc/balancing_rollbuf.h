#ifndef ROLLBUF_H_
#define ROLLBUF_H_

/* Example of a rotating, and reallocating buffer. A separate thread has
 * to fetch the data from the real source and fill the buffer. This way,
 * reads block as little as possible -- sort of a userland pipe, only we
 * get more control. Not clear we gain anything, but it's less fragile
 * than forking an external program.
 */
#include <pthread.h>

struct rollbuf_s {
    char * p;
    size_t alloc;
    size_t avail_to_read;
    size_t avail_to_write;
    const char * rhead;
    char * whead;
    char * rbuf;
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
