#include "cado.h"
#include <time.h>
#include "cado_poly.h"

void
test_cado_poly_set ()
{
  cado_poly p, q;

  cado_poly_init (p);
  cado_poly_init (q);

  mpz_set_ui (q->n, 123456789);
  q->skew = 3.1415;
  q->pols[0]->deg = 1;
  mpz_set_ui (q->pols[0]->coeff[0], 17);
  mpz_set_ui (q->pols[0]->coeff[1], 42);
  q->pols[1]->deg = 0;
  mpz_set_ui (q->pols[1]->coeff[0], 59);
  mpz_set_ui (q->m, 3613);

  cado_poly_set (p, q);

  assert (mpz_cmp_ui (p->n, 123456789) == 0);
  assert (p->skew == 3.1415);
  assert (mpz_cmp_ui (p->pols[0]->coeff[0], 17) == 0);
  assert (mpz_cmp_ui (p->pols[0]->coeff[1], 42) == 0);
  assert (mpz_cmp_ui (p->pols[1]->coeff[0], 59) == 0);
  assert (mpz_cmp_ui (p->m, 3613) == 0);

  cado_poly_clear (p);
  cado_poly_clear (q);
}

int
main (int argc, char *argv[])
{
  long int seed;

  seed = (argc >= 2) ? atoi (argv[1]) : time (NULL);
  fprintf (stderr, "Using random seed=%ld\n", seed);
  srand48 (seed);
  test_cado_poly_set ();
  exit (EXIT_SUCCESS);
}
