/* dup1: 1st duplicate pass, split relation files into 'nslices'
         slices (adapted from check).

   Usage:
   dup1 [-bz] [-n nslices_log] -out <dir> file1 ... filen
   by default nslices_log = 1 (nslices = 2).

   Files file1 ... filen are split into 'nslices' slices in
   <dir>/0/filej ... <dir>/31/filej.

   If option -bz is given, then the output is compressed with bzip2
   instead of gzip.
   Input can be in gzipped or bzipped format.
*/

#include "cado.h"

#define MAX_NSLICES 32

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <inttypes.h>
#include <ctype.h>
#include <time.h>
#include <limits.h> /* for CHAR_BIT */
#include <unistd.h>
#include <assert.h>
#include <fcntl.h>   /* for _O_BINARY */

#include "portability.h"
#include "macros.h"
#include "utils.h"
#include "relation.h"

#include "filter_utils.h"

/* Only (a,b) are parsed on input. This flags control whether we copy the
 * rest of the relation data to the output file, or if we content
 * ourselves with smaller .ab files */
static int only_ab = 0;

/* output relations in dirname/0/name, ..., dirname/31/name */
/* Adds the number of relation added to each slice to *nr_rels */
int split_relfile (relation_stream_ptr rs, const char *name,
                   const char *dirname, const char * outfmt,
                   int *do_slice, int nslices_log, int nslices, 
                   unsigned int ab_base, index_t *nr_rels)
{
    FILE * f_in;
    int p_in;
    const char * suffix_in;

    int ok MAYBE_UNUSED;
    char * oname[MAX_NSLICES];
    FILE * ofile[MAX_NSLICES];
    int p_out[MAX_NSLICES];

    uint64_t h;

    f_in = fopen_maybe_compressed2(name, "r", &p_in, &suffix_in);
    ASSERT_ALWAYS(f_in != NULL);

    char * newname = strdup(name);
    const char * suffix_out = outfmt;
    if (!suffix_out) suffix_out = suffix_in;
    // remove recognized suffix.
    ASSERT_ALWAYS(strlen(suffix_in) <= strlen(newname));
    newname[strlen(newname)-strlen(suffix_in)]='\0';
    for(int i = 0 ; i < nslices ; i++) {
        int rc = asprintf(&oname[i],
                only_ab ? "%s/%d/%s.ab%s" : "%s/%d/%s%s",
                dirname, i, path_basename(newname), suffix_out);
        ASSERT_ALWAYS(rc >= 0);
        ofile[i] = fopen_maybe_compressed2(oname[i], "w", &p_out[i], NULL);
        ASSERT_ALWAYS(ofile[i] != NULL);
    }
    free(newname);

    relation_stream_bind(rs, f_in);
    rs->parse_only_ab = 1;

    size_t pos0 = rs->pos;
    for (int lnum = 0 ; ; lnum++) {
        char line[RELATION_MAX_BYTES];
        size_t rpos = rs->pos - pos0;

        if (relation_stream_get(rs, line, 0, ab_base) < 0)
            break;

	ok = 1;

        h = CA_DUP1 * (uint64_t) rs->rel.a + CB_DUP1 * rs->rel.b;
        /* Using the low bit of h is not a good idea, since then
           odd values of i are twice more likely. The second low bit
           also gives a small bias with RSA768 (but not for random
           coprime a, b). We use here the nslices_log high bits.
        */
        int i = h >> (64 - nslices_log);
        nr_rels[i]++;
        /* print relation */
        if (do_slice[i]) {
            if (only_ab) {
                struct {
                    int64_t a;
                    uint64_t b;
                    size_t pos;
                } foo = { .a = rs->rel.a, .b=rs->rel.b, .pos=rpos, };
                fwrite(&foo, sizeof(foo), 1, ofile[i]);
            } else {
                if (fputs(line, ofile[i]) < 0) {
                    fprintf (stderr, "Problem writing to file %d\n", i);
                    exit(1);
                }
            }
        }

        if (relation_stream_disp_progress_now_p(rs)) {
            fprintf (stderr,
                    "# split %" PRid " relations in %.1fs"
                    " -- %.1f MB/s -- %.1f rels/s\n",
                    rs->nrels, rs->dt, rs->mb_s, rs->rels_s);
        }
    }

    relation_stream_unbind(rs);

    for (int i = 0; i < nslices; i++) {
      if (do_slice[i])
        {
          if (p_out[i]) pclose (ofile[i]); else fclose(ofile[i]);
        }
      free(oname[i]);
    }
    fprintf (stderr, "# done %s\n", name);

    if (p_in) pclose(f_in); else fclose(f_in);
    return 0;
}

/*
void usage()
{
    fprintf(stderr, "Usage: ./dup1 [-only <i>] [-n <nslices_log>] -out <output_dir> [files...]\n");
}
*/
static void
usage(const char *argv0)
{
    fprintf (stderr, "Usage: %s [options] ", argv0);
    fprintf (stderr, "[ -filelist <fl> [-basepath <dir>] | file1 ... filen ]\n");
    fprintf (stderr, "Mandatory command line options:\n");
    fprintf (stderr, "     -out <dir>  - output directory\n");
    fprintf (stderr, "\nOther command line options:\n");
    fprintf (stderr, "    -n           - log of number of slices (default 1)\n");
    fprintf (stderr, "    -only        - do only slice i (default all)\n");
    fprintf (stderr, "    -ab          - only print a and b in the ouput\n");
    fprintf (stderr, "    -outfmt .ext - output is written in .ext files\n");
    fprintf (stderr, "    -abhexa      - read a and b as hexa no decimal\n");
    exit (1);
}

int
main (int argc, char * argv[])
{
    char * argv0 = argv[0];
    int had_error = 0;

    param_list pl;
    param_list_init(pl);
    argv++,argc--;

    int ab_hexa = 0;
    param_list_configure_switch(pl, "ab", &only_ab);
    param_list_configure_switch(pl, "abhexa", &ab_hexa);

#ifdef HAVE_MINGW
    _fmode = _O_BINARY;     /* Binary open for all files */
#endif

    for( ; argc ; ) {
        if (param_list_update_cmdline(pl, &argc, &argv)) { continue; }
        /* Since we accept file names freeform, we decide to never abort
         * on unrecognized options */
        break;
        // fprintf (stderr, "Unknown option: %s\n", argv[0]);
        // abort();
    }

    int nslices_log = 1;
    param_list_parse_int(pl, "n", &nslices_log);
    int nslices = 1 << nslices_log;
    const char * dirname = param_list_lookup_string(pl, "out");
    int only_slice = -1;
    param_list_parse_int(pl, "only", &only_slice);
    const char * outfmt = param_list_lookup_string(pl, "outfmt");
    const char * filelist = param_list_lookup_string(pl, "filelist");
    const char * basepath = param_list_lookup_string(pl, "basepath");
    const char * path_antebuffer = param_list_lookup_string(pl, "path_antebuffer");
    set_antebuffer_path (argv0, path_antebuffer);

    if (param_list_warn_unused(pl)) {
        exit(1);
    }

    if (basepath && !filelist) {
        fprintf(stderr, "-basepath only valid with -filelist\n");
        exit(1);
    }

    if (!dirname)
        usage(argv0);
    if (outfmt && !is_supported_compression_format(outfmt)) {
        fprintf(stderr, "output compression format unsupported\n");
        usage(argv0);
    }

    int do_slice[MAX_NSLICES];

    if (only_slice == -1) { /* split all slices */
        for (int i = 0; i < nslices; i++)
          do_slice[i] = 1;
    } else { /* split only slide i */
        for (int i = 0; i < nslices; i++)
          do_slice[i] = i == only_slice;
    }

    if ((filelist != NULL) + (argc != 0) != 1) {
        fprintf(stderr, "Provide either -filelist or freeform file names\n");
        usage(argv0);
    }

    char ** files = filelist ? filelist_from_file(basepath, filelist, 0) : argv;
    index_t nr_rels[MAX_NSLICES];
    memset (nr_rels, 0, sizeof(index_t) * nslices);

    //process_rels (files, &thread_dup1, NULL, 0, NULL, NULL, STEP_DUP1);

    relation_stream rs;
    relation_stream_init(rs);
    for (char ** fp = files ; *fp ; fp++) {
        had_error |= split_relfile (rs, *fp, dirname, outfmt, do_slice,
                                    nslices_log, nslices, (ab_hexa)?16:10, 
                                    nr_rels);
    }
    relation_stream_trigger_disp_progress(rs);
    fprintf (stderr,
            "# split %" PRid " relations in %.1fs"
            " -- %.1f MB/s -- %.1f rels/s\n",
            rs->nrels, rs->dt, rs->mb_s, rs->rels_s);

    for (int i = 0; i < nslices; i++) {
        fprintf (stderr, "# slice %d received %"PRid" relations\n", i,
                                                                    nr_rels[i]);
    }
    //relation_stream_clear(rs);

    if (filelist) filelist_clear(files);

    param_list_clear(pl);

    if (had_error)
      return 1;
    else
      return 0;
}
