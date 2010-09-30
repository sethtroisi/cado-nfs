#ifndef PURGEDFILE_H_
#define PURGEDFILE_H_
#include <stdio.h>
#include <stdint.h>

struct purgedfile_stream_s {
    FILE * source;
    int nrows, ncols;

    int nodup_index;    // parsed, but useless.
    int64_t a;
    uint64_t b;
    int nc;
    int * cols;
    int nc_alloc;

    // parameters. may be set by the caller after calling
    // purgedfile_stream_init.
    
    int parse_only_ab;  // used when only the (a,b) pair information
                        // is interesting, e.g. for characters
                        
    // various stats stuff. May be used by the caller.
    size_t pos;
    int rrows;
    unsigned long lnum;
    // only valid after purgedfile_stream_trigger_disp_progress
    // of purgedfile_stream_disp_progress_now_p
    double dt, mb_s, rows_s;

    // temporaries, + various stuff for internal use.
    int pipe;
    double t0, t1;
};

typedef struct purgedfile_stream_s purgedfile_stream[1];
typedef struct purgedfile_stream_s * purgedfile_stream_ptr;
typedef const struct purgedfile_stream_s * purgedfile_stream_srcptr;

#ifdef __cplusplus
extern "C" {
#endif

extern void purgedfile_stream_init(purgedfile_stream_ptr ps);
extern void purgedfile_stream_closefile(purgedfile_stream_ptr ps);
extern void purgedfile_stream_clear(purgedfile_stream_ptr ps);
extern int purgedfile_stream_disp_progress_now_p(purgedfile_stream_ptr ps);
extern void purgedfile_stream_trigger_disp_progress(purgedfile_stream_ptr ps);
extern void purgedfile_stream_openfile(purgedfile_stream_ptr ps, const char * fname);
extern void purgedfile_stream_rewind(purgedfile_stream_ptr ps);
extern int purgedfile_stream_get(purgedfile_stream_ptr ps, char * line);

#ifdef __cplusplus
}
#endif

#endif	/* PURGEDFILE_H_ */
