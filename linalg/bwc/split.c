#define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <limits.h>
#include "macros.h"

const char * ifile = "Y.0.twisted";
const char * ofile_fmt = "V%d-%d.0.twisted";

int main(int argc, char * argv[])
{
    int * stops;
    FILE ** files;

    if (argc <= 1) {
        fprintf(stderr, "Usage: split <n_0> <n_1> ...\n");
        exit(1);
    }

    int nstops = argc-1;
    stops = malloc((nstops+1) * sizeof(unsigned int));
    files = malloc(nstops * sizeof(FILE *));
    stops[0] = 0;
    for(int i = 1 ; i < argc ; i++) {
        unsigned int len = atoi(argv[i]);
        ASSERT_ALWAYS(len);
        ASSERT_ALWAYS(len % CHAR_BIT == 0);
        stops[i] = stops[i-1] + len / CHAR_BIT;
    }

    struct stat sbuf[1];
    stat(ifile, sbuf);

    if ((sbuf->st_size) % stops[nstops] != 0) {
        fprintf(stderr, 
                "Size of %s (%ld bytes) is not a multiple of %d bytes\n",
                ifile, (long) sbuf->st_size, stops[nstops]);
    }

    int rc;

    FILE * f = fopen(ifile, "r");
    if (f == NULL) {
        fprintf(stderr,"%s: %s\n", ifile, strerror(errno));
        exit(1);
    }

    if (nstops == 1) {
        char * fname;
        int i = 0;
        /* Special case ; a hard link is enough */
        asprintf(&fname, ofile_fmt, stops[i], CHAR_BIT * stops[i+1]);
        rc = stat(fname, sbuf);
        if (rc == 0) {
            fprintf(stderr,"%s already exists\n", fname);
            exit(1);
        }
        if (rc < 0 && errno != ENOENT) {
            fprintf(stderr,"%s: %s\n", fname, strerror(errno));
            exit(1);
        }
        rc = link(ifile, fname);
        if (rc < 0) {
            fprintf(stderr,"%s: %s\n", fname, strerror(errno));
            exit(1);
        }
        free(fname);
        return 0;
    }

    for(int i = 0 ; i < nstops ; i++) {
        char * fname;
        asprintf(&fname, ofile_fmt, stops[i], CHAR_BIT * stops[i+1]);
        rc = stat(fname, sbuf);
        if (rc == 0) {
            fprintf(stderr,"%s already exists\n", fname);
            exit(1);
        }
        if (rc < 0 && errno != ENOENT) {
            fprintf(stderr,"%s: %s\n", fname, strerror(errno));
            exit(1);
        }
        files[i] = fopen(fname, "w");
        if (files[i] == NULL) {
            fprintf(stderr,"%s: %s\n", fname, strerror(errno));
            exit(1);
        }
        free(fname);
    }

    void * ptr = malloc(stops[nstops]);

    for(;;) {
        rc = fread(ptr, 1, stops[nstops], f);
        if (rc != stops[nstops] && rc != 0) {
            fprintf(stderr, "Unexpected short read\n");
            exit(1);
        }
        if (rc == 0)
            break;

        char * q = ptr;

        for(int i = 0 ; i < nstops ; i++) {
            rc = fwrite(q, 1, stops[i+1]-stops[i], files[i]);
            if (rc != stops[i+1]-stops[i]) {
                fprintf(stderr, "short write\n");
                exit(1);
            }
            q += stops[i+1]-stops[i];
        }
    }

    for(int i = 0 ; i < nstops ; i++) {
        fclose(files[i]);
    }
    fclose(f);
}

