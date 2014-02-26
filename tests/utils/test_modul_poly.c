#include "cado.h"
#include "macros.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "modul_poly.h"
#include "gmp_aux.h"
#include "tests_common.h"
#include "cado_poly.h"

void
test_modul_poly_is_irreducible (unsigned long iter)
{
  modul_poly_t f;
  modulusul_t p;
  int d, i, irred, n;
  unsigned long q;
  residueul_t r[MAXDEGREE];

  while (iter--)
    {
      d = 1 + lrand48 () % (MAXDEGREE - 1);
      modul_poly_init (f, d);
      q = lrand48 () >> 16;
      q = ulong_nextprime (q);
      modul_initmod_ul (p, q);
      for (i = 0; i <= d; i++)
        modul_set_ul (f->coeff[i], lrand48 (), p);
      while (modul_is0 (f->coeff[d], p))
        modul_set_ul (f->coeff[d], lrand48 (), p);
      f->degree = d;
      irred = modul_poly_is_irreducible (f, p);
      if (d == 1)
        ASSERT_ALWAYS(irred != 0);
      if (d == 1 || (d == 2 && irred == 0))
        {
          n = modul_poly_cantor_zassenhaus (r, f, p);
          ASSERT_ALWAYS(n == d);
        }
      modul_clearmod (p);
      modul_poly_clear (f);
    }
}

int
main (int argc, const char *argv[])
{
  unsigned long iter = 1000;

  tests_common_cmdline (&argc, &argv, PARSE_SEED | PARSE_ITER);
  tests_common_get_iter (&iter);
  test_modul_poly_is_irreducible (iter);
  tests_common_clear ();
  exit (EXIT_SUCCESS);
}

