#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include "macros.h"

#include "types.h"
#include "fppol.h"
#include "ffspol.h"
#include "timing.h"
#include "params.h"
#include "smoothness.h"
#include "polyfactor.h"
#include "fq.h"
#include "fqpol.h"

#include "renumber.h"

#define MAX_FFS_DEG 20

int sq_is_split(sq_t * roots, sq_srcptr q, ffspol_srcptr F) {
    fq_info_t Fq;
    fq_info_init(Fq, q);

    fqpol_t f;
    fqpol_init(f);
    fqpol_set_ffspol(f, F, Fq);

    int ret = fqpol_is_split(roots, f, Fq);
    fqpol_clear(f);
    fq_info_clear(Fq);
    return ret;
}


int sq_roots(unsigned long *roots, sq_srcptr q, ffspol_srcptr F)
{
  fq_info_t Fq;
  fq_info_init(Fq, q);
  fq_t r[MAX_FFS_DEG];

  fqpol_t f;
  fqpol_init(f);
  fqpol_set_ffspol(f, F, Fq);

  int nb = fqpol_roots(r, f, Fq);
  fqpol_clear(f);
  fq_info_clear(Fq);

  for (int i = 0; i < nb; i++)
    roots[i] = fppol64_get_ui_sparse (r[i]);
  if (f->deg != F->deg) // There is a projective root
  {
    roots[nb] = fppol64_get_ui_sparse (q);
    nb++;
  }
  return nb;
}

int sq_is_irreducible(sq_srcptr p) {
    fppol_t P;
    fppol_init(P);
    fppol_set_sq(P, p);
    int ret = fppol_is_irreducible(P);
    fppol_clear(P);
    return ret;
}

void usage(const char *argv0, const char * missing)
{
    fprintf(stderr, "Usage: %s [options | optionfile] \n", argv0);
    fprintf(stderr, "  Command line options have the form '-optname optval'\n");
    fprintf(stderr, "  and take the form 'optname=optval' in the optionfile\n");
    fprintf(stderr, "List of options (a * means mandatory, default value in []):\n");
    fprintf(stderr, "  pol0 *           function field polynomial on side 0\n");
    fprintf(stderr, "  pol1 *           function field polynomial on side 1\n");
    fprintf(stderr, "  lpb0 *           large prime bound on side 0\n");
    fprintf(stderr, "  lpb1 *           large prime bound on side 1\n");
    fprintf(stderr, "  gf               indicate the base field for sanity check\n");
 
    if (missing != NULL)
        fprintf(stderr, "Missing parameter: %s\n", missing);
    exit(1);
}

int main(int argc, char **argv)
{
    ffspol_t ffspol[2]; 
    int lpb[2] = {0, 0}, add_full_col = 0;
    char *argv0 = argv[0];
    int gf = 0;
    const char * renumberfilename = NULL;
    const char * badidealsfilename = NULL;
    
    param_list pl;
    param_list_init(pl);
    param_list_configure_knob(pl, "-addfullcol", &add_full_col);
    argv++, argc--;
    for (; argc;) {
        if (param_list_update_cmdline(pl, &argc, &argv)) {
            continue;
        }
        /* Could also be a parameter file */
        FILE *f;
        if ((f = fopen(argv[0], "r")) != NULL) {
            param_list_read_stream(pl, f);
            fclose(f);
            argv++, argc--;
            continue;
        }
        fprintf(stderr, "Unhandled parameter %s\n", argv[0]);
        usage(argv0, NULL);
    }
    // Update parameter list at least once to register argc/argv pointers.
    param_list_update_cmdline(pl, &argc, &argv);
    param_list_print_command_line(stdout, pl);

    param_list_parse_int(pl, "gf", &gf);
    if (gf) {
        if (gf != FP_SIZE) {
            fprintf(stderr, "Error: base field mismatch.\n");
            fprintf(stderr, "  The binary is compiled for GF(%d)\n", FP_SIZE);
            fprintf(stderr, "  The parameters are for GF(%d)\n", gf);
            exit(EXIT_FAILURE);
        }
    }

    // read function field polynomials
    {
        const char * polstr;
        ffspol_init(ffspol[0]);
        ffspol_init(ffspol[1]);
        polstr = param_list_lookup_string(pl, "pol0");
        if (polstr == NULL) usage(argv0, "pol0");
        ffspol_set_str(ffspol[0], polstr);
        polstr = param_list_lookup_string(pl, "pol1");
        if (polstr == NULL) usage(argv0, "pol1");
        ffspol_set_str(ffspol[1], polstr);
    }
    // read various bounds
    param_list_parse_int(pl, "lpb0", &lpb[0]);
    param_list_parse_int(pl, "lpb1", &lpb[1]);
    if (lpb[0] == 0) usage(argv0, "lpb0");
    if (lpb[1] == 0) usage(argv0, "lpb1");

    //read filename
    renumberfilename = param_list_lookup_string(pl, "renumber");
    if (renumberfilename == NULL)
      usage (argv0, "renumber");
    badidealsfilename = param_list_lookup_string(pl, "badideals");


    cado_poly dummy_poly;
    renumber_t tab;
    cado_poly_init(dummy_poly);

    dummy_poly->pols[0]->degree = ffspol[0]->deg; 
    dummy_poly->pols[0]->lpb = __FP_BITS + __FP_BITS * lpb[0]; 
    dummy_poly->pols[1]->degree = ffspol[1]->deg; 
    dummy_poly->pols[1]->lpb = __FP_BITS + __FP_BITS * lpb[1]; 

    if (dummy_poly->pols[1]->degree == 1)
    {
      dummy_poly->rat  = dummy_poly->pols[1];
      dummy_poly->alg  = dummy_poly->pols[0];
    }
    else if (dummy_poly->pols[0]->degree == 1)
    {
      dummy_poly->rat  = dummy_poly->pols[0];
      dummy_poly->alg  = dummy_poly->pols[1];
    }

    renumber_init(tab, dummy_poly);
    renumber_init_write (tab, renumberfilename, badidealsfilename, add_full_col);

    uint64_t nrel = 0;
    sq_t p;
    sq_set_ti(p, 1);   // start with the polynomial t, always irreducible
    unsigned long *roots[2];
    int max_deg = MAX(lpb[0], lpb[1]); 
    int min_deg = MIN(lpb[0], lpb[1]); 
    int k[2], d[2];
    index_t old_table_size = tab->size;

    d[0] = ffspol[0]->deg;
    d[1] = ffspol[1]->deg;
    roots[0] = (unsigned long*) malloc (d[0] * sizeof (unsigned long));
    roots[1] = (unsigned long*) malloc (d[1] * sizeof (unsigned long));

    fprintf (stderr, "Generating freerels up to degree %d\n", min_deg);
    fprintf (stderr, "Generating renumber up to degree %d\n", max_deg);

    while (sq_deg(p) <= max_deg) 
    {
      /* first compute the roots */
      for (int i = 0; i < 2; i++)
      {
        if (sq_deg(p) <= lpb[i])
          k[i] = sq_roots(roots[i], p, ffspol[i]);
        else
          k[i] = 0;
      }

      unsigned long p_int = fppol64_get_ui_sparse (p);
      renumber_write_p (tab, p_int, roots, k);

      if (sq_deg(p) <= min_deg && k[0] == d[0] && k[1] == d[1])
      { 
        //print the free rels
        index_t l;
        printf ("%lx,0:%lx", p_int, (unsigned long) old_table_size);
        for (l = old_table_size + 1; l < tab->size; l++)
          printf (",%lx", (unsigned long) l);
        printf ("\n");
        nrel++;
      }
      old_table_size = tab->size;

      do 
      {
        sq_monic_set_next(p, p, 64);
      } while (!sq_is_irreducible(p));
    } 

    fprintf(stdin, "# Computed %" PRIu64 " free relations\n", nrel);
    fprintf(stderr, "# Computed %" PRIu64 " free relations\n", nrel);
    param_list_clear(pl);
    cado_poly_clear(dummy_poly);
    renumber_close_write(tab);
    ffspol_clear(ffspol[0]);
    ffspol_clear(ffspol[1]);

    return EXIT_SUCCESS;
}
