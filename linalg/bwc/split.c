#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <limits.h>
#include "bwc_config.h"
#include "portability.h"
#include "macros.h"
#include "filenames.h"
#include "bw-common.h"
#include "params.h"
#include "misc.h"


/* Maximum number of possible splits */
#define MAXSPLITS 16

/* splits for the different sites */
int splits[MAXSPLITS + 1];

int force = 0;

int main(int argc, char * argv[])
{
    param_list pl;

    bw_common_init_new(bw, &argc, &argv);
    param_list_init(pl);

    bw_common_decl_usage(pl);
    /* {{{ declare local parameters and switches */
    param_list_decl_usage(pl, "force",
            "overwrite output files if already present");
    param_list_decl_usage(pl, "splits",
            "sub-ranges into which the input file should be split");
    param_list_decl_usage(pl, "binary-ratio",
            "number of bytes per field element (fractional, e.g. 1/8)");
    param_list_decl_usage(pl, "ifile",
            "input file");
    param_list_decl_usage(pl, "ofile-fmt",
            "output file pattern (expects two %u's)");
    param_list_configure_switch(pl, "--force", &force);
    /* }}} */

    bw_common_parse_cmdline(bw, pl, &argc, &argv);

    bw_common_interpret_parameters(bw, pl);

    
    FILE ** files;
    int nsplits;
    int scale[2] = {1, 8};      /* default for binary: one byte of data
                                   contains 8 elements of the field we
                                   are considering */
    int bits_per_coeff;
    char * ifile = NULL;
    char * ofile_fmt = NULL;

    /* {{{ interpret our parameters */
    {
        nsplits = param_list_parse_int_list(pl, "splits", splits, MAXSPLITS, ",");

        bits_per_coeff = 64 * iceildiv(mpz_sizeinbase(bw->p, 2), 64);
        if (mpz_cmp_ui(bw->p, 2) > 0) {
            scale[0] = bits_per_coeff/8;
            scale[1] = 1;
        }

        param_list_parse_int_and_int(pl, "binary-ratio", scale, "/");

        const char * tmp;
        if ((tmp = param_list_lookup_string(pl, "ifile")) != NULL)
            ifile = strdup(tmp);
        if ((tmp = param_list_lookup_string(pl, "ofile-fmt")) != NULL)
            ofile_fmt = strdup(tmp);

        ASSERT_ALWAYS(ifile);
        ASSERT_ALWAYS(ofile_fmt);
    }
    /* }}} */

    if (param_list_warn_unused(pl)) {
        param_list_print_usage(pl, argv[0], stderr);
        exit(EXIT_FAILURE);
    }


    param_list_clear(pl);


    if (nsplits <= 0) {
        fprintf(stderr, "Please indicate the splitting points\n");
        param_list_print_usage(pl, argv[0], stderr);
        exit(EXIT_FAILURE);
    }

    files = malloc(nsplits * sizeof(FILE *));

    for(int i = 0 ; i < nsplits ; i++) {
        ASSERT_ALWAYS((i == 0) == (splits[i] == 0));
        ASSERT_ALWAYS((i == 0) || (splits[i-1] < splits[i]));
    }
    if (splits[nsplits-1] != bw->n) {
        fprintf(stderr, "last split does not coincide with configured n\n");
        exit(EXIT_FAILURE);
    }

    for(int i = 0 ; i < nsplits ; i++) {
        ASSERT_ALWAYS(splits[i] % scale[1] == 0);
    }

    nsplits--;

    int rc;

    struct stat sbuf[1];
    if (stat(ifile, sbuf) < 0) {
        perror(ifile);
        exit(EXIT_FAILURE);
    }

    if ((sbuf->st_size) % (splits[nsplits]*scale[0]/scale[1]) != 0) {
        fprintf(stderr, 
                "Size of %s (%ld bytes) is not a multiple of %d bytes\n",
                ifile, (long) sbuf->st_size,
                splits[nsplits]*scale[0]/scale[1]);
    }

#ifndef HAVE_MINGW
    /* mingw does not have link() */
    if (nsplits == 1) {
        char * fname;
        int i = 0;
        /* Special case ; a hard link is enough */
        rc = asprintf(&fname, ofile_fmt, splits[i], splits[i+1]);
        rc = stat(fname, sbuf);
        if (rc == 0 && !force) {
            fprintf(stderr,"%s already exists\n", fname);
            exit(EXIT_FAILURE);
        }
        if (rc == 0 && force) { unlink(fname); }
        if (rc < 0 && errno != ENOENT) {
            fprintf(stderr,"%s: %s\n", fname, strerror(errno));
            exit(EXIT_FAILURE);
        }
        rc = link(ifile, fname);
        if (rc < 0) {
            fprintf(stderr,"%s: %s\n", fname, strerror(errno));
            exit(EXIT_FAILURE);
        }
        free(fname);
    } else
#endif
    /* mingw always takes this branch */
    {
        FILE * f = fopen(ifile, "rb");
        if (f == NULL) {
            fprintf(stderr,"%s: %s\n", ifile, strerror(errno));
            exit(EXIT_FAILURE);
        }

        for(int i = 0 ; i < nsplits ; i++) {
            char * fname;
            rc = asprintf(&fname, ofile_fmt, splits[i], splits[i+1]);
            rc = stat(fname, sbuf);
            if (rc == 0 && !force) {
                fprintf(stderr,"%s already exists\n", fname);
                exit(EXIT_FAILURE);
            }
            if (rc == 0 && force) { unlink(fname); }
            if (rc < 0 && errno != ENOENT) {
                fprintf(stderr,"%s: %s\n", fname, strerror(errno));
                exit(EXIT_FAILURE);
            }
            files[i] = fopen(fname, "wb");
            if (files[i] == NULL) {
                fprintf(stderr,"%s: %s\n", fname, strerror(errno));
                exit(EXIT_FAILURE);
            }
            free(fname);
        }

        void * ptr = malloc(splits[nsplits]*scale[0]/scale[1]);

        for(;;) {
            rc = fread(ptr, 1, splits[nsplits]*scale[0]/scale[1], f);
            if (rc != splits[nsplits]*scale[0]/scale[1] && rc != 0) {
                fprintf(stderr, "Unexpected short read\n");
                exit(EXIT_FAILURE);
            }
            if (rc == 0)
                break;

            char * q = ptr;

            for(int i = 0 ; i < nsplits ; i++) {
                int d = splits[i+1]*scale[0]/scale[1] - splits[i]*scale[0]/scale[1];
                rc = fwrite(q, 1, d, files[i]);
                if (rc != d) {
                    fprintf(stderr, "short write\n");
                    exit(EXIT_FAILURE);
                }
                q += d;
            }
        }
        for(int i = 0 ; i < nsplits ; i++) {
            fclose(files[i]);
        }
        fclose(f);
    }

    free(files);
    free(ifile);
    free(ofile_fmt);

    bw_common_clear_new(bw);

    return 0;
}

