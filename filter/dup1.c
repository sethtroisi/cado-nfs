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

static index_t nr_rels_tot[MAX_NSLICES];
static int nslices_log = 1, do_slice[MAX_NSLICES];


typedef struct {
  const char *prefix, *suffix;
  char *filename;
  FILE *file;
  const char *msg;
  unsigned int next_idx;
  size_t lines_per_file, lines_left;
} split_output_iter_t;

static split_output_iter_t *
split_iter_init(const char *prefix, const char *suffix, 
                const size_t lines_per_file, const char *msg)
{
  split_output_iter_t *iter = malloc(sizeof(split_output_iter_t));
  ASSERT_ALWAYS(iter != NULL);
  iter->prefix = strdup(prefix);
  iter->suffix = strdup(suffix);
  iter->next_idx = 0;
  iter->filename = NULL;
  iter->file = NULL;
  if (msg)
    iter->msg = strdup(msg);
  else
    iter->msg = NULL;
  ASSERT_ALWAYS(lines_per_file > 0);
  iter->lines_per_file = lines_per_file;
  iter->lines_left = 0; /* Force opening of file on next write */
  return iter; 
}

/* used for counting time in different processes */
timingstats_dict_t stats;


static void
split_iter_end(split_output_iter_t *iter)
{
  if (iter->file != NULL)
    fclose(iter->file);
  free(iter->filename);
  free((void *) iter->prefix);
  free((void *) iter->suffix);
  free((void *) iter->msg);
  free(iter);
}

/* Closes the currently open file, if any, and opens the next one */
void 
split_iter_open_next_file(split_output_iter_t *iter)
{
  if (iter->file != NULL) {
    struct rusage r[1];
    fclose_maybe_compressed2(iter->file, iter->filename, r);
    timingstats_dict_add(stats, iter->prefix, r);
  }

  free (iter->filename);
  int rc = asprintf(&(iter->filename), "%s%04x%s", 
                    iter->prefix, iter->next_idx++, iter->suffix);
  ASSERT_ALWAYS (rc >= 0);
  if (iter->msg != NULL)
    fprintf (stderr, "%s%s\n", iter->msg, iter->filename);
  iter->file = fopen_maybe_compressed(iter->filename, "w");
  iter->lines_left = iter->lines_per_file;
}

static void
split_iter_write_next(split_output_iter_t *iter, const char *line)
{
  if (iter->lines_left == 0)
    split_iter_open_next_file(iter);
  fputs (line, iter->file);
  iter->lines_left--;
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
thread_dup1 (void * context_data, earlyparsed_relation_ptr rel)
{
    unsigned int slice = compute_slice (rel->a, rel->b);
    split_output_iter_t **outiters = (split_output_iter_t**)context_data;

    if (do_slice[slice])
    {
      if (only_ab)
      {
        char *p = rel->line;
        while (*p != ':')
          p++;
        *p = '\n';
      }

      split_output_iter_t *iter = outiters[slice];
      split_iter_write_next(iter, rel->line);
      nr_rels_tot[slice]++;
    }
    return NULL;
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

    int abhexa = 0;
    param_list_configure_switch(pl, "ab", &only_ab);
    param_list_configure_switch(pl, "abhexa", &abhexa);

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
    int nslices = 1 << nslices_log;
    const char *outdir = param_list_lookup_string(pl, "out");
    int only_slice = -1;
    param_list_parse_int(pl, "only", &only_slice);
    unsigned int log_max_nrels_per_files = DEFAULT_LOG_MAX_NRELS_PER_FILES;
    param_list_parse_uint(pl, "lognrelsoutfile", &log_max_nrels_per_files);
    const char *outfmt = param_list_lookup_string(pl, "outfmt");
    const char * filelist = param_list_lookup_string(pl, "filelist");
    const char * basepath = param_list_lookup_string(pl, "basepath");
    const char * path_antebuffer = param_list_lookup_string(pl, "path_antebuffer");
    set_antebuffer_path (argv0, path_antebuffer);
    const char *prefix_files = param_list_lookup_string(pl, "prefix");

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

    // If not output suffix is specified, use suffix of first input file
    if (!outfmt && files[0] != NULL)
      get_suffix_from_filename (files[0], &outfmt);

    memset (nr_rels_tot, 0, sizeof(index_t) * nslices);

    split_output_iter_t **outiters;
    outiters = malloc(sizeof(split_output_iter_t *) * nslices);
    ASSERT_ALWAYS(outiters != NULL);
    for(int i = 0 ; i < nslices ; i++)
    {
      char *prefix, *suffix, *msg;
      int rc = asprintf(&prefix, "%s/%d/%s.",
                        outdir, i, prefix_files);
      ASSERT_ALWAYS(rc >= 0);
      rc = asprintf(&suffix, only_ab ? ".ab%s" : "%s", outfmt);
      ASSERT_ALWAYS(rc >= 0);
      rc = asprintf (&msg, "# Opening output file for slice %d : ", i);
      ASSERT_ALWAYS(rc >= 0);
      outiters[i] = split_iter_init(prefix, suffix, 1UL<<log_max_nrels_per_files, msg);
      free(prefix);
      free(suffix);
      free(msg);
    }

    timingstats_dict_init(stats);
    filter_rels(files, (filter_rels_callback_t) &thread_dup1, (void*)outiters,
            EARLYPARSE_NEED_LINE |
            (abhexa ? EARLYPARSE_NEED_AB_HEXA : EARLYPARSE_NEED_AB_DECIMAL),
            NULL, stats);

    for(int i = 0 ; i < nslices ; i++)
      split_iter_end(outiters[i]);

    for (int i = 0; i < nslices; i++)
        fprintf (stderr, "# slice %d received %" PRid " relations\n", i,
                                                                nr_rels_tot[i]);

    if (filelist) filelist_clear(files);

    free(outiters);

    param_list_clear(pl);

    // double thread_times[2];
    // thread_seconds_user_sys(thread_times);
    timingstats_dict_add_mythread(stats, "main");
    // fprintf(stderr, "Main thread ends after having spent %.2fs+%.2fs on cpu \n", thread_times[0], thread_times[1]);
    timingstats_dict_disp(stats);
    timingstats_dict_clear(stats);

    return 0;
}
