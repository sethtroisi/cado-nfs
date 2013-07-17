#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mod_ul.c"
#include "portability.h"
#include "utils.h"

#ifdef FOR_FFS
#include "fppol.h"
#include "fq.h"
#include "utils_ffs.h"
#endif

static void
usage (char *argv0)
{
  fprintf (stderr, "Usage: %s [-check] -poly xxx.poly -renumberfile outfile\n",
                   argv0);
  exit (1);
}

int
main (int argc, char *argv[])
{
    char *polyfilename = NULL;
    int k, check = 0;
    char *renumberfilename = NULL;
    char *argv0 = argv[0];
    cado_poly cpoly;
    renumber_t tab;

    fprintf (stderr, "%s.r%s", argv[0], CADO_REV);
    for (k = 1; k < argc; k++)
      fprintf (stderr, " %s", argv[k]);
    fprintf (stderr, "\n");

    while (argc > 1 && argv[1][0] == '-')
      {
        if (argc > 2 && strcmp (argv[1], "-poly") == 0)
          {
            polyfilename = argv[2];
            argc -= 2;
            argv += 2;
          }
        else if (argc > 2 && strcmp (argv[1], "-renumber") == 0)
          {
            renumberfilename = argv[2];
            argc -= 2;
            argv += 2;
          }
        else if (argc > 1 && strcmp (argv[1], "-check") == 0)
          {
            check = 1;
            argc -= 1;
            argv += 1;
          }
        else
          usage (argv0);
      }

    if (polyfilename == NULL)
      usage (argv0);
    if (renumberfilename == NULL)
      usage (argv0);

    cado_poly_init(cpoly);
#ifndef FOR_FFS
    if (!cado_poly_read (cpoly, polyfilename))
#else
    if (!ffs_poly_read (cpoly, polyfilename))
#endif
    {
      fprintf (stderr, "Error reading polynomial file\n");
      exit (EXIT_FAILURE);
    }

    renumber_debug_print_tab(stdout, renumberfilename, cpoly);

  /* Check for all index i if it can retrieved p and r from it and if it can
   * retrieved the same index from this p and r
   */
  if (check)
  {
    renumber_init (tab, cpoly);
    renumber_read_table (tab, renumberfilename);

    index_t i, j;
    p_r_values_t p, r;
    int side;
    for (i = 0; i < tab->size; i++)
    {
      if (tab->table[i] != RENUMBER_SPECIAL_VALUE)
      {
        renumber_get_p_r_from_index(tab, &p, &r, &side, i, cpoly);
        j = renumber_get_index_from_p_r (tab, p, r, side);
        if (i == j)
          fprintf (stderr, "## %"PRid":Ok\n", i);
        else
          fprintf (stderr, "#### %"PRid":Error:%"PRid"\n", i, j);
      }
    }

    renumber_free (tab);
  }

    cado_poly_clear (cpoly);
    return 0;
}
