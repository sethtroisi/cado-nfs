/* This program takes a matrix file -- hopefully with sorted rows and
 * columns -- and prepares ready-to-work in-memory datasets for a grid of
 * programs.
 *
 * The program must be invoked with as many mpi instances as there are
 * nodes in the desired grid. Each instance will handle the preparation
 * of the memory data it feels concerned with.
 *
 * Usage: ./preparegrid h=42 v=17 f=c131.small
 */

#define __STDC_FORMAT_MACROS    /* for c++ compilers */
#define _GNU_SOURCE             /* strndup */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <inttypes.h>   /* SCNu32 etc */
// #include <time.h>
// #include <math.h> /* for log, pow, fabs */
#include "cado.h"
#include "utils/utils.h"

#include "select_mpi.h"

#include "readbuffer.h"

#define PARSE_CHECK(tst, kind, text) do {                               \
    if ((tst)) {                                                        \
        fprintf(stderr, "parse error while reading %s: %s",             \
                kind, text);                                            \
        exit(1);                                                        \
    }                                                                   \
} while (0)

void usage()
{
    fprintf(stderr, "Usage: ./preparegrid h=42 v=17 f=c131.small\n");
    exit(1);
}

uint32_t nr;
uint32_t nc;
char * header /* = NULL */;

void read_matrix(const char * filename)
{
    reading_buffer b;
    int rc;

    rb_open(b, filename);

    /* Read matrix header */
    rb_read_line(b);
    if (strncmp(b->buf, "//", 2) == 0) {
        rc = sscanf(b->buf, "// %" SCNu32 " ROWS %" SCNu32, &nr, &nc);
        PARSE_CHECK(rc < 2, "header", b->buf);
    } else {
        rc = sscanf(b->buf, "%" SCNu32 " %" SCNu32, &nr, &nc);
        PARSE_CHECK(rc < 2, "header", b->buf);
    }
    rb_gobble_long_line(b);
    // header_bytes = b->o;
    header = strndup(b->buf, b->o);



    rb_close(b);
}

int
main (int argc, char *argv[])
{
    param_list pl;
    int hslices;
    int vslices;
    int i;
    int verbose=0;
    char matrix_filename[FILENAME_MAX];

    MPI_Init(&argc,&argv);

    /* print command line */
    fprintf (stderr, "# %s.r%s", argv[0], REV);
    for (i = 1; i < argc; i++)
        fprintf (stderr, " %s", argv[i]);
    fprintf (stderr, "\n");

    param_list_init (pl);
    argv++, argc--;
    for( ; argc ; ) {
        if (strcmp(argv[0], "-v") == 0) { verbose++; argv++,argc--; continue; }
        if (param_list_update_cmdline(pl, NULL, &argc, &argv)) { continue; }
        fprintf(stderr, "Unhandled parameter %s\n", argv[0]);
        usage();
    }

    param_list_parse_string(pl, "f", matrix_filename, sizeof(matrix_filename));
    param_list_parse_int(pl, "h", &hslices);
    param_list_parse_int(pl, "v", &vslices);
    if (verbose)
        param_list_display (pl, stderr);
    if (param_list_warn_unused(pl)) {
        usage();
    }
    param_list_clear(pl);

    read_matrix(matrix_filename);


    return 0;
}
