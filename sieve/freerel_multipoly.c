/*
 * Program: free relations
 * Original author : F. Morain
 * Purpose: creating free relations in a suitable format
 * Modified / rewritten by C. Bouvier (and others)

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

#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>   /* for _O_BINARY */
#include <string.h>
#include <math.h>

#include <gmp.h>
#include "portability.h"
#include "mod_ul.c"
#include "utils.h"
#include "typedefs.h"

/* generate all free relations up to the large prime bound */
/* generate the renumbering table */


static unsigned long MAYBE_UNUSED
allFreeRelations (cado_poly pol, unsigned long pmin, unsigned long pmax,
                  renumber_t renumber_table, const char *outfilename)
{
  unsigned long lpb[NB_POLYS_MAX], max_lpb;
  unsigned long p, *roots[NB_POLYS_MAX], nfree = 0;
  int d[NB_POLYS_MAX], k[NB_POLYS_MAX];
  index_t old_table_size = renumber_table->size;
  FILE *fpout = NULL;
  mpz_srcptr lc[NB_POLYS_MAX];

  fpout = fopen_maybe_compressed (outfilename, "w");
  ASSERT_ALWAYS (fpout != NULL);

  /* we generate all free relations up to the *maximum* of the two large
     prime bounds (the optimal would be up to the second maximum). */
  /* we generate the renumbering table up to the *maximun* of the two large
     prime bounds */
  for (unsigned int i = 0; i < renumber_table->nb_polys; i++)
  {
    ASSERT_ALWAYS(renumber_table->lpb[i] < sizeof(unsigned long) * CHAR_BIT);
    lpb[i] = 1UL << renumber_table->lpb[i];
    d[i] = pol->pols[i]->deg;
    lc[i] = pol->pols[i]->coeff[d[i]];
    roots[i] = (unsigned long *) malloc (d[i] * sizeof(unsigned long));
    ASSERT_ALWAYS (roots[i] != NULL);
  }
  max_lpb = 1UL << renumber_table->max_lpb;

  if (pmax == 0)
    pmax = max_lpb;
  ASSERT_ALWAYS (pmax <= max_lpb);

  printf ("Generating freerels for %lu <= p <= %lu\n", pmin, pmax);
  printf ("Generating renumber table for 2 <= p <= %lu\n", max_lpb);
  fflush (stdout);

  stats_data_t stats; /* struct for printing progress */
  uint64_t nb_p = 0;
  /* will print report at 2^10, 2^11, ... 2^23 computed primes and every
   * 2^23 primes after that */
  stats_init (stats, stdout, &nb_p, 23, "Looked into", "primes", "", "p");
  for (p = 2; p <= max_lpb; p = getprime (p))
  {
    /* first compute the roots */
    for (unsigned int side = 0; side < renumber_table->nb_polys; side++)
    {
      if (p >= lpb[side])
      {
        k[side] = 0;
        continue;
      }
      if (d[side] == 1)
      {
        k[side] = 1;
        continue;
      }
      k[side] = mpz_poly_roots_ulong (roots[side], pol->pols[side], p);
      // Check for a projective root
      if (mpz_divisible_ui_p (lc[side], p))
        roots[side][k[side]++] = p;
    }

    renumber_write_p (renumber_table, p, roots, k);

    if (pmin <= p && p <= pmax)
    {
      index_t begin_i = old_table_size;
      for (unsigned int i = 0; i < renumber_table->nb_polys; i++)
      {
        if (i > 0)
          begin_i += k[i-1];
        index_t begin_j = begin_i;
        for (unsigned int j = i+1; j < renumber_table->nb_polys; j++)
        {
          begin_j += k[j-1];
          if (k[i] == d[i] && k[j] == d[j]) /* freerel between side i and j */
          {
            fprintf (fpout, "# Freerel for p = %lu between side %u and %u\n",
                            p, i, j);
            fprintf (fpout, "%lx,0:", p);
            for (int l = 0; l < k[i]; l++)
              fprintf (fpout, "%" PRid ",", (index_t) begin_i + l);
            for (int l = 0; l < k[j]; l++)
              fprintf (fpout, "%" PRid "%s", (index_t) begin_j + l,
                                             (l == k[j]-1) ? "\n" : ",");
            nfree++;
          }
        }
      }
    }

    old_table_size = renumber_table->size;
    nb_p++;
    if (stats_test_progress(stats))
      stats_print_progress (stats, nb_p, 0, 0, 0);
  }
  stats_print_progress (stats, nb_p, 0, 0, 1);

  getprime (0);
  for (unsigned int i = 0; i < renumber_table->nb_polys; i++)
    free (roots[i]);

  fclose_maybe_compressed (fpout, outfilename);

  return nfree;
}

static void declare_usage(param_list pl)
{
  param_list_decl_usage(pl, "poly", "input polynomial file");
  param_list_decl_usage(pl, "renumber", "output file for renumbering table");
  param_list_decl_usage(pl, "out", "output file for free relations");
  for (unsigned int i = 0; i < NB_POLYS_MAX; i++)
  {
    char desc[64], name[8];
    snprintf (desc, 64, "large prime bound on side %u", i);
    snprintf (name, 8, "lpb%u", i);
    param_list_decl_usage(pl, name, desc);
  }
  param_list_decl_usage(pl, "pmin", "do not create freerel below this bound");
  param_list_decl_usage(pl, "pmax", "do not create freerel beyond this bound");
  param_list_decl_usage(pl, "badideals", "file describing bad ideals (for DL)");
  param_list_decl_usage(pl, "addfullcol", "(switch) add a column of 1 in the matrix (for DL)");
}

static void
usage (param_list pl, char *argv0)
{
    param_list_print_usage(pl, argv0, stderr);
    exit(EXIT_FAILURE);
}

int
main (int argc, char *argv[])
{
  char *argv0 = argv[0];
  cado_poly cpoly;
  unsigned long pmin = 2, pmax = 0, nfree;
  renumber_t renumber_table;
  int add_full_col = 0;
  unsigned long lpb[NB_POLYS_MAX] = { 0 };

  param_list pl;
  param_list_init(pl);
  declare_usage(pl);
  param_list_configure_switch(pl, "-addfullcol", &add_full_col);

#ifdef HAVE_MINGW
  _fmode = _O_BINARY;     /* Binary open for all files */
#endif

  argv++, argc--;
  if (argc == 0)
    usage (pl, argv0);

  for( ; argc ; )
  {
    if (param_list_update_cmdline(pl, &argc, &argv))
      continue;
    FILE *f;
    if ((f = fopen(argv[0], "r")) != NULL)
    {
      param_list_read_stream(pl, f, 0);
      fclose(f);
      argv++,argc--;
      continue;
    }
    fprintf(stderr, "Unhandled parameter %s\n", argv[0]);
    usage (pl, argv0);
  }
  /* print command-line arguments */
  param_list_print_command_line (stdout, pl);
  printf ("\n");
  fflush(stdout);

  const char * polyfilename = param_list_lookup_string(pl, "poly");
  const char * outfilename = param_list_lookup_string(pl, "out");
  const char * badidealsfilename = param_list_lookup_string(pl, "badideals");
  const char * renumberfilename = param_list_lookup_string(pl, "renumber");
  param_list_parse_ulong(pl, "pmin", &pmin);
  param_list_parse_ulong(pl, "pmax", &pmax);
  for (unsigned int i = 0; i < NB_POLYS_MAX; i++)
  {
    char name[8];
    snprintf (name, 8, "lpb%u", i);
    param_list_parse_ulong(pl, name, &lpb[i]);
  }

  if (polyfilename == NULL)
  {
    fprintf (stderr, "Error, missing -poly command line argument\n");
    usage (pl, argv0);
  }
  if (renumberfilename == NULL)
  {
    fprintf (stderr, "Error, missing -renumber command line argument\n");
    usage (pl, argv0);
  }
  if (outfilename == NULL)
  {
    fprintf (stderr, "Error, missing -out command line argument\n");
    usage (pl, argv0);
  }

  cado_poly_init(cpoly);
  if (!cado_poly_read (cpoly, polyfilename))
  {
    fprintf (stderr, "Error reading polynomial file\n");
    exit (EXIT_FAILURE);
  }

  for (unsigned int i = 0; i < cpoly->nb_polys; i++)
  {
    if (lpb[i] == 0)
    {
      fprintf (stderr, "Error, missing -lpb%u command line argument\n", i);
      usage (pl, argv0);
    }
  }

  if (param_list_warn_unused(pl))
    usage (pl, argv0);

  int ratside = cado_poly_get_ratside (cpoly);
  renumber_init_for_writing (renumber_table, cpoly->nb_polys, ratside,
                                                           add_full_col, lpb);
  renumber_write_open (renumber_table, renumberfilename, badidealsfilename,
                       cpoly);

  nfree = allFreeRelations (cpoly, pmin, pmax, renumber_table, outfilename);

  /* /!\ Needed by the Python script. /!\ */
  fprintf (stderr, "# Free relations: %lu\n", nfree);
  fprintf (stderr, "Renumbering struct: nprimes=%" PRIu64 "\n",
                   renumber_table->size);

  renumber_write_close (renumber_table, renumberfilename);
  renumber_clear (renumber_table);
  cado_poly_clear (cpoly);
  param_list_clear(pl);

  return 0;
}
