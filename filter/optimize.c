/* optimize --- optimize cycles in the merged matrix

Copyright 2008, 2009, 2010, 2011, 2012 Francois Morain, Paul Zimmermann

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

/*
Reference: Jason S. Papadopoulos, A Self-Tuning Filtering Implementation for
the Number Field Sieve, slides presented at the CADO Workshop on Integer
Factorization, October 2008,
http://cado.gforge.inria.fr/workshop/slides/papadopoulos.pdf
*/

#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "utils.h"

#include "sparse.h"
#include "gzip.h"

#define MERGE_MAX 64 /* maximal number of entries in a merge row */

uint64_t total_weight;

/* return the length of the merge of row[i] and row[j] */
static uint32_t
len_merge (uint32_t **row, uint8_t *len_row, uint32_t i, uint32_t j)
{
  uint32_t li, lj, ki, kj, ltmp;

  li = len_row[i];
  lj = len_row[j];
  for (ki = 0, kj = 0, ltmp = 0; ki < li && kj < lj;)
    {
      if (row[i][ki] < row[j][kj])
        ltmp++, ki++;
      else if (row[i][ki] > row[j][kj])
        ltmp++, kj++;
      else
        ki++, kj++;
    }
  /* only one of the following loops is non-empty */
  return ltmp + (li - ki) + (lj - kj);
}

static void
print_row (uint32_t **row, uint8_t *len_row, uint32_t i)
{
  uint32_t j;

  for (j = 0; j < len_row[i]; j++)
    printf ("%u ", row[i][j]);
  printf ("\n");
}

/* add row j in column r */
void
add_relation (uint32_t **transpose, uint8_t *len_transpose, uint32_t j,
              uint32_t r)
{
  if (transpose == NULL)
    return;
  ASSERT_ALWAYS(len_transpose[r] < 255);
  len_transpose[r] ++;
  transpose[r] = (uint32_t*) realloc (transpose[r],
                                      len_transpose[r] * sizeof(uint32_t));
  transpose[r][len_transpose[r] - 1] = j;
}

/* subtract row j from column r */
void
sub_relation (uint32_t **transpose, uint8_t *len_transpose, uint32_t j,
              uint32_t r)
{
  uint32_t k;

  if (transpose == NULL)
    return;
  for (k = 0; k < len_transpose[r]; k++)
    if (transpose[r][k] == j)
      break;
  if (k >= len_transpose[r])
    {
      fprintf (stderr, "Not found row %u in column %u:\n", j, r);
      print_row (transpose, len_transpose, r);
      exit (1);
    }
  ASSERT_ALWAYS(k < len_transpose[r]);
  ASSERT_ALWAYS(len_transpose[r] > 0);
  len_transpose[r] --;
  transpose[r][k] = transpose[r][len_transpose[r]];
  transpose[r] = (uint32_t*) realloc (transpose[r],
                                      len_transpose[r] * sizeof(uint32_t));
}

/* row[j] += row[i].
   Update row and len_row.
   If transpose <> NULL, update transpose and len_transpose:
   (a) each relation only in i appears once more in j
   (b) each relation only in j appears with the same count
   (c) each relation in both i and j appears once less in j */
static void
add_row (uint32_t **row, uint8_t *len_row, uint32_t j, uint32_t i,
         uint32_t **transpose, uint8_t *len_transpose)
{
  uint32_t *tmp, li, lj, ki, kj, ltmp;

  ASSERT_ALWAYS(row[i] != NULL);
  ASSERT_ALWAYS(row[j] != NULL);
  li = len_row[i];
  lj = len_row[j];
  tmp = (uint32_t*) malloc ((li + lj) * sizeof (uint32_t));
  for (ki = 0, kj = 0, ltmp = 0; ki < li && kj < lj;)
    {
      if (row[i][ki] < row[j][kj])
        {
          add_relation (transpose, len_transpose, j, row[i][ki]);
          tmp[ltmp++] = row[i][ki++];
        }
      else if (row[i][ki] > row[j][kj])
        tmp[ltmp++] = row[j][kj++];
      /* if row[i][ki] = row[j][kj], we add twice the same relation, and all
         exponents vanish modulo 2 */
      else
        {
          sub_relation (transpose, len_transpose, j, row[i][ki]);
          ki++, kj++;
        }
    }
  /* only one of the following loops is non-empty */
  while (ki < li)
    {
      add_relation (transpose, len_transpose, j, row[i][ki]);
      tmp[ltmp++] = row[i][ki++];
    }
  while (kj < lj)
    tmp[ltmp++] = row[j][kj++];
  if (ltmp < li + lj)
    tmp = (uint32_t*) realloc (tmp, ltmp * sizeof (uint32_t));
  free (row[j]);
  row[j] = tmp;
  ASSERT_ALWAYS(ltmp <= 255);
  len_row[j] = ltmp;
}

/* try a merge from the relation-sets containing relation k */
static void
try_merge (uint32_t **row, uint8_t *len_row, uint32_t **transpose,
           uint8_t *len_transpose, uint32_t k, FILE *fp)
{
  uint32_t i, j, li, lj, lij, imax, jmax, *tk;
  int32_t gain, gain_max = 0;

  tk = transpose[k];
  for (i = 0; i < len_transpose[k]; i++)
    {
      li = len_row[tk[i]];
      for (j = i + 1; j < len_transpose[k]; j++)
        {
          lj = len_row[tk[j]];
          lij = len_merge (row, len_row, tk[i], tk[j]);
          gain = (int32_t) ((li > lj) ? li : lj) - (int32_t) lij;
          if (gain > gain_max)
            {
              gain_max = gain;
              if (li > lj)
                {
                  imax = i;
                  jmax = j;
                }
              else
                {
                  imax = j;
                  jmax = i;
                }
            }
        }
    }

  if (gain_max > 0) /* perform the best merge, len_row[imax]>=len_row[jmax] */
    {
      fprintf (fp, "-%u %u\n", tk[jmax] + 1, tk[imax]);
      add_row (row, len_row, tk[imax], tk[jmax], transpose, len_transpose);
      /* update the transpose matrix */
      total_weight -= gain_max;
    }
}

#if 0
static void
merge2 (uint32_t **row, uint8_t *len_row, uint32_t i, uint32_t j, uint32_t r)
{
  uint32_t l;
  static int count = 0;

  l = len_merge (row, len_row, i, j);
  if (l < len_row[i] || l < len_row[j])
    {
      printf ("Merge row %u and %u for relation %u and save %u\n", i, j, r,
              ((len_row[i] > len_row[j]) ? len_row[i] : len_row[j]) - l);
      print_row (row, len_row, i);
      print_row (row, len_row, j);
      if (count++ == 10)
        abort();
    }
}
#endif

static void
optimize (FILE *hisfile, FILE *optfile)
{
  uint32_t nrows, ncols, **row, i, j, k, merge[MERGE_MAX], line_number = 0;
  uint32_t tot_count[256], small_nrows, small_ncols, **transpose, pass = 0;
  uint8_t *len_row, *count, *len_transpose;
  int res, neg;
  char *lineptr[1], *str;
  size_t linesize;
  uint64_t old_total_weight;

  linesize = 1024;
  *lineptr = malloc (linesize * sizeof(char));

  /* read first line: number of relations-sets and number of relations */
  res = getline (lineptr, &linesize, hisfile);
  ASSERT_ALWAYS(res != -1);

  res = sscanf (*lineptr, "%u %u", &nrows, &ncols);
  ASSERT_ALWAYS(res == 2);

  fprintf (stderr, "Optimizing %u relation-sets over %u relations\n",
           nrows, ncols);

  /* initialize the ith relation-set to the single relation i */
  len_row = (uint8_t*) malloc (nrows * sizeof (uint8_t));
  row = (uint32_t**) malloc (nrows * sizeof (uint32_t*));
  for (i = 0; i < nrows; i++)
    {
      len_row[i] = 1;
      row[i] = (uint32_t*) malloc (1 * sizeof (uint32_t));
      row[i][0] = i;
    }

  /* read the history file and perform merges:
     each line (except the first one) is of the form "i i1 ... ik":
     If i >= 0 then row[i] is to be added to rows i1...ik and destroyed at
               the end of the process. Works also is i is alone (hence:
               destroyed row).
     If i < 0 then row[-i-1] is to be added to rows i1...ik and NOT destroyed.
  */

  line_number = 1;
  while (!feof (hisfile))
    {
      res = getline (lineptr, &linesize, hisfile);
      if (res == -1) /* end of file */
        break;
      line_number ++;
      str = *lineptr; /* start of the line */
      if (str[0] == '-')
        {
          neg = 1;
          str ++;
        }
      else
        neg = 0;
      /* read first integer i */
      res = sscanf (str, "%u", &i);
      ASSERT_ALWAYS(res == 1);
      while (isdigit (str[0]))
        str ++;
      k = 0;
      while (str[0] == ' ') /* there is another integer */
        {
          res = sscanf (++str, "%u", merge + k);
          ASSERT_ALWAYS(res == 1);
          k ++;
          ASSERT_ALWAYS(k <= MERGE_MAX);
          while (isdigit (str[0]))
            str ++;
        }
      ASSERT_ALWAYS(str[0] == '\n');

      /* perform the merge */
      i -= neg;
      ASSERT_ALWAYS(row[i] != NULL); /* check row i is still active */
      for (j = 0; j < k; j++)
        /* add row i to row merge[j] */
        add_row (row, len_row, merge[j], i, NULL, NULL);
      if (neg == 0) /* destroy row i */
        {
          free (row[i]);
          row[i] = NULL;
          len_row[i] = 0;
        }
    }
  free (*lineptr);

  printf ("Read %u merges\n", line_number - 1);

  /* now count frequency of each relation */
  count = (uint8_t*) malloc ((nrows + 1) * sizeof (uint8_t));
  for (i = 0; i <= nrows; i++)
    count[i] = 0;
  small_nrows = 0;
  total_weight = 0;
  for (i = 0; i < nrows; i++)
    {
      if (row[i] != NULL)
        {
          small_nrows ++;
          total_weight += len_row[i];
          for (j = 0; j < len_row[i]; j++)
            {
              k = row[i][j];
              ASSERT_ALWAYS(k <= nrows);
              ASSERT_ALWAYS(count[k] < 255);
              count[k] ++;
            }
        }
    }

  /* compute number of relations with given frequency */
  printf ("Total relation weight: %lu (aver. %1.2f)\n", total_weight,
          (double) total_weight / (double) small_nrows);
  for (i = 0; i < 256; i++)
    tot_count[i] = 0;
  small_ncols = 0;
  for (i = 0; i <= nrows; i++)
    {
      tot_count[count[i]] ++;
      small_ncols += (count[i] > 0);
    }
#if 0
  for (i = 0; i < 256; i++)
    {
      if (tot_count[i] != 0)
        printf ("# of relations appearing %u times: %u\n", i, tot_count[i]);
    }
  printf ("Remaining %u relations-sets on %u relations\n", small_nrows,
          small_ncols);
#endif

  ASSERT_ALWAYS(nrows >= ncols);
  transpose = (uint32_t**) malloc ((nrows + 1) * sizeof(uint32_t*));
  len_transpose = (uint8_t*) malloc ((nrows + 1) * sizeof(uint8_t));
  for (i = 0; i <= nrows; i++)
    {
      transpose[i] = NULL;
      len_transpose[i] = 0;
      if (count[i] > 0)
        transpose[i] = (uint32_t*) malloc (count[i] * sizeof(uint32_t));
    }

  for (j = 0; j < nrows; j++)
    if (row[j] != NULL)
      {
        for (k = 0; k < len_row[j]; k++)
          {
            i = row[j][k];
            ASSERT_ALWAYS(i <= nrows);
            ASSERT_ALWAYS(len_transpose[i] < count[i]);
            transpose[i][len_transpose[i]] = j;
            ASSERT_ALWAYS(len_transpose[i] < 255);
            len_transpose[i] ++;
          }
      }

  /* check len_transpose is correct */
  for (i = 0; i <= nrows; i++)
    ASSERT_ALWAYS(len_transpose[i] == count[i]);

  do
    {
      old_total_weight = total_weight;
      for (j = 0; j <= nrows; j++)
        {
          if (count[j] >= 2)
            try_merge (row, len_row, transpose, len_transpose, j, optfile);
        }
      printf ("Pass %u: total weight: %lu (aver. %1.2f)\n", ++pass,
              total_weight, (double) total_weight / (double) small_nrows);
    }
  while (total_weight < old_total_weight);

  for (i = 0; i <= nrows; i++)
    free (transpose[i]);
  free (transpose);
  free (len_transpose);
  free (count);
  for (i = 0; i < nrows; i++)
    free (row[i]);
  free (len_row);
  free (row);
}

int
main (int argc, char *argv[])
{
  FILE *hisfile, *optfile;
  int verbose = 0;
  double wct0 = wct_seconds ();
  param_list pl;
  const char *hisname, *optname;

  setbuf (stdout, NULL);
  setbuf (stderr, NULL);

  /* printing the arguments as everybody does these days */
  fprintf (stderr, "%s.r%s", argv[0], CADO_REV);
  for (int i = 1; i < argc; i++)
    fprintf (stderr, " %s", argv[i]);
  fprintf (stderr, "\n");

  param_list_init (pl);
  argv++,argc--;
  param_list_configure_knob (pl, "--verbose", &verbose);
  param_list_configure_alias(pl, "--verbose", "-v");

  for( ; argc ; )
    {
      if (param_list_update_cmdline (pl, &argc, &argv))
        continue;
      fprintf (stderr, "Unknown option: %s\n", argv[0]);
      exit (1);
    }

  hisname = param_list_lookup_string (pl, "his");
  optname = param_list_lookup_string (pl, "opt");
  
  hisfile = fopen (hisname, "r");
  ASSERT_ALWAYS(hisfile != NULL);

  optfile = fopen (optname, "w");
  ASSERT_ALWAYS(optfile != NULL);

  optimize (hisfile, optfile);

  fclose (hisfile);
  fclose (optfile);

  param_list_clear (pl);

  print_timing_and_memory (wct0);

  return 0;
}
