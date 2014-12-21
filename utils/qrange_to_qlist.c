/* Usage: cut_n_roots -poly xxx.poly -q0 nnn -q1 nnn -N nnn 
     split [q0, q1[ into maximal intervals of
     N special q's each. */

#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include "modul_poly.h"
#include <gmp.h>
#include "portability.h"
#include "utils.h"
#include "macros.h"

void usage(const char *argv0)
{
    fprintf(stderr, "Usage: %s -poly xxx.poly -q0 nnn -q1 nnn\n", argv0);
    fprintf(stderr, "  prints all (q,r) in an interval [q0, q1[.\n");
    exit(1);
}

int
main (int argc0, char *argv0[])
{
  const char *polyfilename = NULL;
  cado_poly pol;
  param_list pl;
  mpz_t q0, q1;
  mpz_t P;
  int argc = argc0;
  char **argv = argv0;
  int ratq = 0;
  int side = ALGEBRAIC_SIDE;

  param_list_init(pl);
  param_list_configure_switch(pl, "-ratq", &ratq);
  argv++, argc--;
  for( ; argc ; ) {
    if (param_list_update_cmdline(pl, &argc, &argv)) { continue; }
    /* Could also be a file */
    FILE * f;
    if ((f = fopen(argv[0], "r")) != NULL) {
      param_list_read_stream(pl, f, 0);
      fclose(f);
      argv++,argc--;
      continue;
    }
    fprintf(stderr, "Unhandled parameter %s\n", argv[0]);
    usage(argv0[0]);
  }
  if (ratq) side = RATIONAL_SIDE;

  mpz_init_set_ui(q0,0);
  mpz_init_set_ui(q1,0);
  param_list_parse_mpz(pl, "q0", q0);
  param_list_parse_mpz(pl, "q1", q1);
  if (mpz_cmp_ui(q0,0) == 0 || mpz_cmp_ui(q1,0) ==0) usage(argv0[0]);

  /* check that q1 fits into an unsigned long */
  if (!mpz_fits_ulong_p(q1)) {
    gmp_fprintf (stderr, "Error, q1=%Zd exceeds ULONG_MAX\n", q1);
    exit (EXIT_FAILURE);
  }

  cado_poly_init (pol);
  ASSERT_ALWAYS((polyfilename = param_list_lookup_string(pl, "poly")) != NULL);
  if (!cado_poly_read (pol, polyfilename)) {
    fprintf (stderr, "Error reading polynomial file\n");
    usage(argv[0]);
    exit (EXIT_FAILURE);
  }
  if (pol->skew <= 0.0) {
    fprintf (stderr, "Error, please provide a positive skewness\n");
    exit (EXIT_FAILURE);
  }

  mpz_sub_ui(q0,q0,1);
  mpz_init_set (P, q0);
  mpz_nextprime (P, P);
  mpz_poly_ptr ps = pol->pols[ALGEBRAIC_SIDE];
  mpz_t * roots;
  roots = malloc(ps->deg * sizeof(mpz_t));
  for(int i = 0 ; i < ps->deg ; i++) {
      mpz_init(roots[i]);
  }
  while (mpz_cmp (P, q1) < 0) {
    int nr = mpz_poly_roots (roots, ps, P);
    for(int i = 0 ; i < nr ; i++) {
        gmp_printf("%s %Zd %Zd\n", sidenames[side], P, roots[i]);
    }
    mpz_nextprime(P,P);
  }
  for(int i = 0 ; i < ps->deg ; i++) {
      mpz_clear(roots[i]);
  }
  free(roots);

  mpz_clear(q0);
  mpz_clear(q1);
  mpz_clear (P);
  cado_poly_clear (pol);
  param_list_clear(pl);
  return 0;
}
