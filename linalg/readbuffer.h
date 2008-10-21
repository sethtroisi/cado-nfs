#ifndef READBUFFER_H_
#define READBUFFER_H_

#include <stdio.h>
#include <sys/types.h>

/* reading buffer stuff
 *
 * This interface provides convenience functions for accessing and
 * parsing files. The typical scenario goes as follows.
 *
 * rb_open
 * many times:
 *      rb_read_line
 *      the line is in b->buf.
 *      parse the line using your favorite tokenizing routine, e.g.
 *      strtoul (*). It is the caller's responsibility to maintain a
 *      index to the current still-to-parse data (as a position within
 *      the b->buf buffer).
 *      if you ever fear that the default 1k buffer won't be long enough,
 *      you may periodically use rb_feed_buffer_again_if_lowwater, with
 *      *p_pos pointing to the position mentioned above, and the lowwater
 *      argument telling how many bytes are needed.
 * rb_close
 *
 * */
struct reading_buffer_s {
    char buf[1024];
    FILE * f;
    int siz;
    off_t o;
    char * filename;
};

typedef struct reading_buffer_s reading_buffer[1];

#ifdef __cplusplus
extern "C" {
#endif

extern void rb_open(reading_buffer b, const char * filename);
extern void rb_close(reading_buffer b);
extern void rb_read_line(reading_buffer b);
extern void rb_gobble_long_line(reading_buffer b);
extern void rb_feed_buffer_again_if_lowwater(reading_buffer b,
        int * ppos, int low);

#ifdef __cplusplus
}
#endif

#endif	/* READBUFFER_H_ */
