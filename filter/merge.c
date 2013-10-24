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

#define MAXLEVEL_DEFAULT 10
#define KEEP_DEFAULT 160
#define FORBW_DEFAULT 0
#define RATIO_DEFAULT 1.1
#define COVERNMAX_DEFAULT 100.0
#define MKZTYPE_DEFAULT 1 /* pure Markowitz */
#define WMSTMAX_DEFAULT 7 /* relevant only if mkztype == 2 */

static void
usage (const char *argv0)
{
  fprintf (stderr, "Usage: %s [options]\n", argv0);
  fprintf (stderr, "\nMandatory command line options: \n");
  fprintf (stderr, "   -mat   xxx     - input (purged) file is xxx\n");
  fprintf (stderr, "   -out   xxx     - output (history) file is xxx\n");
  fprintf (stderr, "\nOther command line options: \n");
  fprintf (stderr, "   -maxlevel nnn  - merge up to nnn rows (default %u)\n",
	   MAXLEVEL_DEFAULT);
  fprintf (stderr, "   -keep nnn      - keep an excess of nnn (default %u)\n",
	   KEEP_DEFAULT);
  fprintf (stderr, "   -skip nnn      - bury the nnn heaviest columns (default %u)\n",
	   SKIP_DEFAULT);
  fprintf (stderr, "   -forbw nnn     - controls the optimization function (default %u, see below)\n",
	   FORBW_DEFAULT);
  fprintf (stderr, "   -ratio rrr     - maximal ratio cN(final)/cN(min) with forbw=0 (default %1.1f)\n",
	   RATIO_DEFAULT);
  fprintf (stderr, "   -coverNmax nnn - with forbw=3, stop when c/N exceeds nnn (default %1.2f)\n", COVERNMAX_DEFAULT);
  fprintf (stderr, "   -itermax nnn   - if non-zero, stop when nnn columns have been removed (cf -resume)\n");
  fprintf (stderr, "   -resume xxx    - resume from history file xxx (cf -itermax)\n");
  fprintf (stderr, "   -mkztype nnn   - controls how the weight of a merge is approximated (default %d)\n", MKZTYPE_DEFAULT);
  fprintf (stderr, "   -wmstmax nnn   - if mkztype = 2, controls until when a mst is used (default %d)\n", WMSTMAX_DEFAULT);
  fprintf (stderr, "   -path_antebuffer <dir> - where is antebuffer\n");
  fprintf (stderr, "\nThe different optimization functions are, where c is the total matrix weight\n");
  fprintf (stderr, "and N the number of rows (relation-sets):\n");
  fprintf (stderr, "   -forbw 0 - optimize the matrix size N (cf -ratio)\n");
  fprintf (stderr, "   -forbw 1 - stop when the product cN is minimal\n");
  fprintf (stderr, "   -forbw 3 - stop when the ratio c/N exceeds coverNmax\n");
  exit (1);
}

int
main (int argc, char *argv[])
{
    char *argv0 = argv[0];

    filter_matrix_t mat[1];
    report_t rep[1];

    uint32_t maxlevel = MAXLEVEL_DEFAULT;
    uint32_t keep = KEEP_DEFAULT;
    uint32_t skip = SKIP_DEFAULT;
    double ratio = RATIO_DEFAULT; /* bound on cN_new/cN to stop the merge */
    uint32_t forbw = FORBW_DEFAULT;
    double coverNmax = COVERNMAX_DEFAULT;
    uint32_t mkztype = MKZTYPE_DEFAULT;
    uint32_t wmstmax = WMSTMAX_DEFAULT;
                               /* use real MST minimum for wt[j] <= wmstmax*/
    uint32_t itermax = 0;

#ifdef HAVE_MINGW
    _fmode = _O_BINARY;     /* Binary open for all files */
#endif

    double tt;
    double wct0 = wct_seconds ();
    param_list pl;
    param_list_init (pl);

    param_list_configure_switch(pl, "--force-posix-threads", &filter_rels_force_posix_threads);
    argv++, argc--;

    for( ; argc ; ) {
      if (param_list_update_cmdline(pl, &argc, &argv)) continue;
      fprintf (stderr, "Unknown option: %s\n", argv[0]);
      usage (argv0);
    }

    /* Update parameter list at least once to register argc/argv pointers. */
    param_list_update_cmdline (pl, &argc, &argv);
    /* print command-line arguments */
    param_list_print_command_line (stdout, pl);
    fflush(stdout);

    const char * purgedname = param_list_lookup_string (pl, "mat");
    const char * outname = param_list_lookup_string (pl, "out");
    /* -resume can be useful to continue a merge stopped due  */
    /* to a too small value of -maxlevel                      */
    const char * resumename = param_list_lookup_string (pl, "resume");
  const char * path_antebuffer = param_list_lookup_string(pl, "path_antebuffer");

    set_antebuffer_path (argv0, path_antebuffer);
    
    param_list_parse_uint (pl, "maxlevel", &maxlevel);
    param_list_parse_uint (pl, "keep", &keep);
    param_list_parse_uint (pl, "skip", &skip);
    param_list_parse_uint (pl, "forbw", &forbw);

    param_list_parse_double (pl, "ratio", &ratio);
    param_list_parse_double (pl, "coverNmax", &coverNmax);

    param_list_parse_uint (pl, "mkztype", &mkztype);
    param_list_parse_uint (pl, "wmstmax", &wmstmax);

    param_list_parse_uint (pl, "itermax", &itermax);

    /* Some checks on command line arguments */
    if (purgedname == NULL || outname == NULL)
    {
      fprintf (stderr, "Error: -mat and -out are mandatory.\n");
      usage (argv0);
    }

    if (maxlevel == 0 || maxlevel > MERGE_LEVEL_MAX)
    {
      fprintf (stderr, "Error: maxlevel should be positive and less than %d\n",
                       MERGE_LEVEL_MAX);
      exit (1);
    }

    if (forbw > 3)
    {
      fprintf (stderr, "Error: -forbw should be 0, 1, 2 or 3.\n");
      exit (1);
    }

    if (mkztype > 2)
    {
      fprintf (stderr, "Error: -mkztype should be 0, 1, or 2.\n");
      exit (1);
    }

    if (param_list_warn_unused (pl))
      usage (argv0);

    /* Read number of rows and cols on first line of purged file */
    purgedfile_read_firstline (purgedname, &(mat->nrows), &(mat->ncols));

    /* initialize rep (i.e., mostly opens outname) and write matrix dimension */
    rep->type = 0;
    rep->outfile = fopen_maybe_compressed (outname, "w");
    ASSERT_ALWAYS(rep->outfile != NULL);

    /* Init structure containing the matrix and the heap of potential merges */
    initMat (mat, maxlevel, keep, skip, itermax);

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

    mergeOneByOne (rep, mat, maxlevel, forbw, ratio, coverNmax);

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
