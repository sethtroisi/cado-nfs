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

#define DEFAULT_LOG_MAX_NRELS_PER_FILES 25

/* Only (a,b) are parsed on input. This flags control whether we copy the
 * rest of the relation data to the output file, or if we content
 * ourselves with smaller .ab files */
static int only_ab = 0;

static index_t nr_rels_tot[MAX_NSLICES], nr_rels_file[MAX_NSLICES];
static int nslices_log = 1, nslices, do_slice[MAX_NSLICES];

static unsigned int log_max_nrels_per_files = DEFAULT_LOG_MAX_NRELS_PER_FILES;

static unsigned int next_num_files = 0;
static const char * prefix_files = NULL;
static FILE * outfiles[MAX_NSLICES];
static const char * outfmt;
static const char * outdir;
static char * filename[MAX_NSLICES];

static void
init_outfiles ()
{
  for(int i = 0 ; i < nslices ; i++)
  {
    int rc = asprintf(&filename[i],
                      only_ab ? "%s/%d/%s.%04x.ab%s" : "%s/%d/%s.%04x%s",
                      outdir, i, prefix_files, next_num_files, outfmt);

    ASSERT_ALWAYS(rc >= 0);
    outfiles[i] = fopen_maybe_compressed(filename[i], "w");
    ASSERT_ALWAYS(outfiles[i] != NULL);
    nr_rels_file[i] = 0;
  }
}

static void
close_outfiles ()
{
  for(int i = 0 ; i < nslices ; i++)
  {
    fclose_maybe_compressed(outfiles[i], filename[i]);
    free(filename[i]);
    nr_rels_tot[i] += nr_rels_file[i];
  }
}

static void
open_next_outfile ()
{
  close_outfiles ();
  next_num_files++;
  init_outfiles ();
}


static inline unsigned int
compute_slice (int64_t a, uint64_t b)
{
  uint64_t h = CA_DUP1 * (uint64_t) a + CB_DUP1 * b;
  /* Using the low bit of h is not a good idea, since then
     odd values of i are twice more likely. The second low bit
     also gives a small bias with RSA768 (but not for random
     coprime a, b). We use here the nslices_log high bits.
  */
  return (unsigned int) h >> (64 - nslices_log);
}

/* Callback function called by prempt_scan_relations */

static void *
thread_dup1(buf_arg_t *arg)
{
  unsigned int slice, j;
  unsigned long cpy_cpt_rel_b;
  buf_rel_t *my_rel;

  cpy_cpt_rel_b = cpt_rel_b;
  for ( ; ; )
  {
    while (cpt_rel_a == cpy_cpt_rel_b)
      if (!is_finish())
        NANOSLEEP();
      else if (cpt_rel_a == cpy_cpt_rel_b)
          pthread_exit(NULL);

    if (cpt_rel_a == cpy_cpt_rel_b + 1)
      NANOSLEEP();

    j = (unsigned int) (cpy_cpt_rel_b & (SIZE_BUF_REL - 1));
    my_rel = &(arg->rels[j]);

    slice = compute_slice (my_rel->a, my_rel->b);

    if (do_slice[slice])
    {
      if (only_ab)
      {
        char *p = my_rel->line;
        while (*p != ':')
          p++;
        *p = '\n';
      }

      fputs (my_rel->line, arg->fd[slice]);
    }

    nr_rels_file[slice]++;
    if (nr_rels_file[slice] >> log_max_nrels_per_files)
      open_next_outfile();

    test_and_print_progress_now ();
    cpy_cpt_rel_b++;
    cpt_rel_b = cpy_cpt_rel_b;
  }
}

static void
usage(const char *argv0)
{
    fprintf (stderr, "Usage: %s [options] ", argv0);
    fprintf (stderr, "[ -filelist <fl> [-basepath <dir>] | file1 ... filen ]\n");
    fprintf (stderr, "Mandatory command line options:\n");
    fprintf (stderr, "     -out <dir>  - output directory\n");
    fprintf (stderr, "     -prefix <s> - prefix for output files\n");
    fprintf (stderr, "\nOther command line options:\n");
    fprintf (stderr, "    -lognrelsoutfile n - log of number of rels in output"
                     " files (default %d)\n", DEFAULT_LOG_MAX_NRELS_PER_FILES);
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

    /* Update parameter list at least once to register argc/argv pointers. */
    param_list_update_cmdline (pl, &argc, &argv);
    /* print command-line arguments */
    param_list_print_command_line (stdout, pl);
    fflush(stdout);

    param_list_parse_int(pl, "n", &nslices_log);
    nslices = 1 << nslices_log;
    outdir = param_list_lookup_string(pl, "out");
    int only_slice = -1;
    param_list_parse_int(pl, "only", &only_slice);
    param_list_parse_uint(pl, "lognrelsoutfile", &log_max_nrels_per_files);
    outfmt = param_list_lookup_string(pl, "outfmt");
    const char * filelist = param_list_lookup_string(pl, "filelist");
    const char * basepath = param_list_lookup_string(pl, "basepath");
    const char * path_antebuffer = param_list_lookup_string(pl, "path_antebuffer");
    set_antebuffer_path (argv0, path_antebuffer);
    prefix_files = param_list_lookup_string(pl, "prefix");

    if (param_list_warn_unused(pl)) {
        exit(1);
    }

    if (basepath && !filelist) {
        fprintf(stderr, "-basepath only valid with -filelist\n");
        exit(1);
    }

    if (!prefix_files)
        usage(argv0);

    if (!outdir)
        usage(argv0);
    if (outfmt && !is_supported_compression_format(outfmt)) {
        fprintf(stderr, "output compression format unsupported\n");
        usage(argv0);
    }

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

    if (!outfmt) // In that case, the suffix_in of first input file is used
      get_suffix_from_filename (files[0], &outfmt);

    memset (nr_rels_tot, 0, sizeof(index_t) * nslices);
    init_outfiles();

    unsigned int step;
    if (ab_hexa)
      step = STEP_DUP1_HEXA;
    else
      step = STEP_DUP1_DECIMAL;
    
    process_rels (files, &thread_dup1, NULL, 0, outfiles, NULL, step);

    close_outfiles();

    for (int i = 0; i < nslices; i++)
        fprintf (stderr, "# slice %d received %" PRid " relations\n", i,
                                                                nr_rels_tot[i]);

    if (filelist) filelist_clear(files);

    param_list_clear(pl);

    return 0;
}
