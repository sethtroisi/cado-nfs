/* merge --- main program to merge relations into relation-sets (cycles)

Copyright 2008-2009 Francois Morain.
Reviewed by Paul Zimmermann, February 2009.

This file is part of CADO-NFS.

CADO-NFS is free software; you can redistribute it and/or modify it under the
terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 2.1 of the License, or (at your option)
any later version.

CADO-NFS is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
details.

You should have received a copy of the GNU Lesser General Public License
along with CADO-NFS; see the file COPYING.  If not, write to the Free Software
Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA.
*/

/* Algorithm: digest and interpolation from Cavallar02.

  Stefania Cavallar, On the Number Field Sieve Integer Factorisation Algorithm,
  PhD thesis, University of Leiden, 2002, 108 pages.
 */

#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>   /* for _O_BINARY */
#include <string.h> /* for strcmp */

#include "portability.h"
#include "utils.h" /* for fopen_maybe_compressed */

#include "filter_utils.h"
#include "filter_matrix.h" /* for filter_matrix_t */
#include "report.h"     /* for report_t */
#include "markowitz.h" /* for MkzInit */
#include "merge_mono.h" /* for mergeOneByOne */
#include "sparse.h"

static void declare_usage(param_list pl)
{
  param_list_decl_usage(pl, "mat", "input purged file");
  param_list_decl_usage(pl, "out", "output history file");
  param_list_decl_usage(pl, "keep", "excess to keep (default "
                                    STR(DEFAULT_FILTER_EXCESS) ")");
  param_list_decl_usage(pl, "skip", "number of heavy columns to bury (default "
                                    STR(DEFAULT_MERGE_SKIP) ")");
  param_list_decl_usage(pl, "maxlevel", "maximum number of rows in a merge "
                            "(default " STR(DEFAULT_MERGE_MAXLEVEL) ")");
  param_list_decl_usage(pl, "forbw", "controls the optimization function "
                            "(see below, default " STR(DEFAULT_MERGE_FORBW) ")");
  param_list_decl_usage(pl, "ratio", "maximal ration cN(final)/cN(min) with "
                            "-forbw 0 (default " STR(DEFAULT_MERGE_RATIO) ")");
  param_list_decl_usage(pl, "coverNmax", "stop when c/N exceeds this value with"
                            " -forbw 3 (default " STR(DEFAULT_MERGE_COVERNMAX) ")");
  param_list_decl_usage(pl, "itermax", "maximum number of columns that can be "
                                       "removed (0 means no maximum)");
  param_list_decl_usage(pl, "resume", "resume from history file (cf -itermax)");
  param_list_decl_usage(pl, "mkztype", "controls how the weight of a merge is "
                            "approximated (default " STR(DEFAULT_MERGE_MKZTYPE) ")");
  param_list_decl_usage(pl, "wmstmax", "controls until when a mst is used with "
                            "-mkztype 2 (default " STR(DEFAULT_MERGE_WMSTMAX) ")");
  param_list_decl_usage(pl, "force-posix-threads", "(switch)");
  param_list_decl_usage(pl, "path_antebuffer", "path to antebuffer program");
}

static void
usage (param_list pl, char *argv0)
{
    param_list_print_usage(pl, argv0, stderr);
    fprintf (stderr, "\nThe different optimization functions are, where c is "
                     "the total matrix weight and N \nthe number of rows "
                     "(relation-sets):\n");
  fprintf (stderr, "   -forbw 0 - stop when cN exceeds ratio*min(cN)\n");
  fprintf (stderr, "   -forbw 1 - stop when the product cN is minimal\n");
  fprintf (stderr, "   -forbw 3 - stop when the ratio c/N exceeds coverNmax\n");
    exit(EXIT_FAILURE);
}


int
main (int argc, char *argv[])
{
    char *argv0 = argv[0];

    filter_matrix_t mat[1];
    report_t rep[1];

    int maxlevel = DEFAULT_MERGE_MAXLEVEL;
    uint32_t keep = DEFAULT_FILTER_EXCESS;
    uint32_t skip = DEFAULT_MERGE_SKIP;
    double ratio = DEFAULT_MERGE_RATIO; /* bound on cN_new/cN to stop the merge */
    uint32_t forbw = DEFAULT_MERGE_FORBW;
    double coverNmax = DEFAULT_MERGE_COVERNMAX;
    uint32_t mkztype = DEFAULT_MERGE_MKZTYPE;
    uint32_t wmstmax = DEFAULT_MERGE_WMSTMAX;
                               /* use real MST minimum for wt[j] <= wmstmax*/
    int64_t nbmergemax = -1; /* Negative value means no maximum */

#ifdef HAVE_MINGW
    _fmode = _O_BINARY;     /* Binary open for all files */
#endif

    double tt;
    double wct0 = wct_seconds ();
    param_list pl;
    param_list_init (pl);
    declare_usage(pl);
    argv++,argc--;

    param_list_configure_switch(pl, "force-posix-threads", &filter_rels_force_posix_threads);

#ifdef HAVE_MINGW
    _fmode = _O_BINARY;     /* Binary open for all files */
#endif

    if (argc == 0)
      usage (pl, argv0);

    for( ; argc ; ) {
      if (param_list_update_cmdline(pl, &argc, &argv)) continue;
      fprintf (stderr, "Unknown option: %s\n", argv[0]);
      usage (pl, argv0);
    }
    /* print command-line arguments */
    param_list_print_command_line (stdout, pl);
    fflush(stdout);

    const char * purgedname = param_list_lookup_string (pl, "mat");
    const char * outname = param_list_lookup_string (pl, "out");
    /* -resume can be useful to continue a merge stopped due  */
    /* to a too small value of -maxlevel                      */
    const char * resumename = param_list_lookup_string (pl, "resume");
    const char *path_antebuffer = param_list_lookup_string(pl, "path_antebuffer");

    param_list_parse_int (pl, "maxlevel", &maxlevel);
    param_list_parse_uint (pl, "keep", &keep);
    param_list_parse_uint (pl, "skip", &skip);
    param_list_parse_uint (pl, "forbw", &forbw);

    param_list_parse_double (pl, "ratio", &ratio);
    param_list_parse_double (pl, "coverNmax", &coverNmax);

    param_list_parse_uint (pl, "mkztype", &mkztype);
    param_list_parse_uint (pl, "wmstmax", &wmstmax);

    param_list_parse_int64 (pl, "nbmergemax", &nbmergemax);

    /* Some checks on command line arguments */
    if (param_list_warn_unused(pl))
    {
      fprintf(stderr, "Error, unused parameters are given\n");
      usage(pl, argv0);
    }

    if (purgedname == NULL)
    {
      fprintf(stderr, "Error, missing -mat command line argument\n");
      usage (pl, argv0);
    }
    if (outname == NULL)
    {
      fprintf(stderr, "Error, missing -out command line argument\n");
      usage (pl, argv0);
    }
    if (maxlevel <= 0 || maxlevel > MERGE_LEVEL_MAX)
    {
      fprintf (stderr, "Error: maxlevel should be positive and less than %d\n",
                       MERGE_LEVEL_MAX);
      usage (pl, argv0);
    }
    if (forbw > 3 || forbw == 2)
    {
      fprintf (stderr, "Error: -forbw should be 0, 1 or 3.\n");
      usage (pl, argv0);
    }
    if (mkztype > 2)
    {
      fprintf (stderr, "Error: -mkztype should be 0, 1, or 2.\n");
      usage (pl, argv0);
    }

    set_antebuffer_path (argv0, path_antebuffer);

    /* Read number of rows and cols on first line of purged file */
    purgedfile_read_firstline (purgedname, &(mat->nrows), &(mat->ncols));

    /* initialize rep (i.e., mostly opens outname) and write matrix dimension */
    rep->type = 0;
    rep->outfile = fopen_maybe_compressed (outname, "w");
    ASSERT_ALWAYS(rep->outfile != NULL);

    /* Init structure containing the matrix and the heap of potential merges */
    initMat (mat, maxlevel, keep, skip);

    /* Read all rels and fill-in the mat structure */
    tt = seconds ();
    filter_matrix_read (mat, purgedname);
    printf ("Time for filter_matrix_read: %2.2lfs\n", seconds () - tt);

    /* resume from given history file if needed */
    if (resumename != NULL)
      resume (rep, mat, resumename);

    mat->wmstmax = wmstmax;
    mat->mkztype = mkztype;
    tt = seconds();
    MkzInit (mat);
    printf ("Time for MkzInit: %2.2lfs\n", seconds()-tt);

    mergeOneByOne (rep, mat, maxlevel, forbw, ratio, coverNmax, nbmergemax);

    fclose_maybe_compressed (rep->outfile, outname);
    printf ("Final matrix has N=%" PRIu64 " nc=%" PRIu64 " (%" PRIu64 ") "
            "w(M)=%" PRIu64 " N*w(M)=%" PRIu64 "\n ", mat->rem_nrows,
            mat->rem_ncols, mat->rem_nrows - mat->rem_ncols, mat->weight,
	          mat->rem_nrows * mat->weight);
    fflush (stdout);
    MkzClose (mat);
    clearMat (mat);

    param_list_clear (pl);

    printf ("Total merge time: %1.0f seconds\n", seconds ());

    print_timing_and_memory (wct0);

    return 0;
}
