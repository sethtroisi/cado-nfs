#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mod_ul.c"
#include "portability.h"
#include "utils.h"

static void
usage (char *argv0)
{
  fprintf (stderr, "Usage: %s -poly xxx.poly -renumberfile outfile\n", argv0);
  exit (1);
}

int
main (int argc, char *argv[])
{
    char *polyfilename = NULL;
    int k;
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
        else if (argc > 2 && strcmp (argv[1], "-renumberfile") == 0)
          {
            renumberfilename = argv[2];
            argc -= 2;
            argv += 2;
          }
        else
          usage (argv0);
      }

    if (polyfilename == NULL)
      usage (argv0);
    if (renumberfilename == NULL)
      usage (argv0);

    cado_poly_init(cpoly);
    if (!cado_poly_read (cpoly, polyfilename))
      {
        fprintf (stderr, "Error reading polynomial file\n");
        exit (EXIT_FAILURE);
      }

    /* check that n divides Res(f,g) [might be useful to factor n...] */
    cado_poly_check (cpoly);

    renumber_debug_print_tab(stdout, renumberfilename, cpoly);

/*
    renumber_init (tab, cpoly);
    renumber_read_table (tab, renumberfilename);

    index_t i, j;
    p_r_values_t p, r;
    int side;
    for (i = 0; i < tab->size; i++)
    {
      renumber_get_p_r_from_index(tab, &p, &r, i, cpoly); 
      side = (r > p) ? 0 : 1;
      j = renumber_get_index_from_p_r (tab, p, r, side);
      if (i == j)
        fprintf (stderr, "%"PRid":Ok\n", i);
      else
        fprintf (stderr, "%"PRid":Error:%"PRid"\n", i, j);
    }

    renumber_free (tab);
*/
    cado_poly_clear (cpoly);
    return 0;
}
