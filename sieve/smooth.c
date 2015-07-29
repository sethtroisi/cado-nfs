/* Bernstein's smoothness test */

/* This is a work-in-progress implementation. To use it:
   1) create a file (say cofac) containing lines of the following form:
      a b cofac_0 cofac_1
      where cofac_0 and cofac_1 are cofactors on sides 0 and 1 respectively
   2) run ./smooth -poly ... -cofac cofac -lpb0 ... -lpb1 ... -lim0 ...
      -lim1 ... [-v] [-split ...]

   The algorithm is the following:
   (a) compute the product P of all primes in [B, L]
   (b) compute the product tree T of all cofactors
   (c) compute the remainder tree of P mod T: at the leaves we have 0 for
       smooth cofactors (assuming all primes in [B,L] have multiplicity 1)
   In practice we split P into smaller chunks P_j of bit-size near that of all
   cofactors, and compute a gcd at leaves to reveal primes in P_j.
   If the number of cofactors doubles, the cost goes from M(n)*log(n) to
   M(2n)*log(2n), thus essentially doubles.
*/

#include "cado.h"
#include <stdio.h>
#include <gmp.h>
#include "facul.h"
#include "facul_doit.h"
#include "ecm/batch.h"
#include "utils.h"
#include "portability.h"

static void
declare_usage (param_list pl)
{
  param_list_decl_usage (pl, "poly", "Polynomial file");
  param_list_decl_usage (pl, "cofac", "Cofactor file");
  param_list_decl_usage (pl, "lpb0", "side-0 large prime bound is 2^lpb0");
  param_list_decl_usage (pl, "lpb1", "side-1 large prime bound is 2^lpb1");
  param_list_decl_usage (pl, "lim0", "side-0 factor base bound");
  param_list_decl_usage (pl, "lim1", "side-1 factor base bound");
  param_list_decl_usage (pl, "v",    "(switch) verbose mode");
  param_list_decl_usage (pl, "split", "number of splits");
  param_list_decl_usage (pl, "batch0", "side-0 batch file");
  param_list_decl_usage (pl, "batch1", "side-1 batch file");
}

int
main (int argc, char* argv[])
{
  cofac_list L;
  double start;
  FILE *cofac;
  FILE *batch[2] = {NULL, NULL};
  int64_t a;
  uint64_t b;
  mpz_t R, A;
  param_list pl;
  const char *poly_file, *cofac_file, *batch0_file, *batch1_file;
  int lpb[2], verbose, split = 10;
  unsigned long lim[2];

  start = seconds ();

  param_list_init (pl);
  declare_usage (pl);

  param_list_configure_switch (pl, "-v", NULL);

  argc--, argv++;
  for( ; argc ; )
    {
      if (param_list_update_cmdline (pl, &argc, &argv))
        continue;
      else
        break;
    }

  poly_file = param_list_lookup_string (pl, "poly");
  cofac_file = param_list_lookup_string (pl, "cofac");
  batch0_file = param_list_lookup_string (pl, "batch0");
  batch1_file = param_list_lookup_string (pl, "batch1");
  ASSERT_ALWAYS(param_list_parse_int (pl, "lpb1", &(lpb[1])));
  ASSERT_ALWAYS(param_list_parse_int (pl, "lpb0", &(lpb[0])));
  ASSERT_ALWAYS(param_list_parse_ulong (pl, "lim1", &(lim[1])));
  ASSERT_ALWAYS(param_list_parse_ulong (pl, "lim0", &(lim[0])));
  param_list_parse_int (pl, "split", &split);
  verbose = param_list_parse_switch (pl, "-v");

  if (poly_file == NULL)
    {
      fprintf (stderr, "Error, missing -poly <file>\n");
      exit (1);
    }
  cado_poly pol;
  cado_poly_init (pol);
  if (cado_poly_read (pol, poly_file) == 0)
    {
      fprintf (stderr, "Could not read polynomial file\n");
      exit (1);
    }

  if (cofac_file == NULL)
    {
      fprintf (stderr, "Error, missing -cofac <file>\n");
      exit (1);
    }

  if (batch0_file == NULL)
    {
      fprintf (stderr, "Error, missing -batch0 <file>\n");
      exit (1);
    }
  batch[0] = fopen (batch0_file, "r");
  if (batch[0] == NULL)
    {
      create_batch_file (batch0_file, lim[0], 1UL << lpb[0], pol->pols[0],
                         split);
      batch[0] = fopen (batch0_file, "r");
    }
  ASSERT_ALWAYS(batch[0] != NULL);

  if (batch1_file == NULL)
    {
      fprintf (stderr, "Error, missing -batch1 <file>\n");
      exit (1);
    }
  batch[1] = fopen (batch1_file, "r");
  if (batch[1] == NULL)
    {
      create_batch_file (batch1_file, lim[1], 1UL << lpb[1], pol->pols[1],
                         split);
      batch[1] = fopen (batch1_file, "r");
    }
  ASSERT_ALWAYS(batch[1] != NULL);

  /* Initialization */
  cofac_list_init (L);

  cofac = fopen (cofac_file, "r");

  mpz_init (R);
  mpz_init (A);
  while (1)
  {
    char str[1024];
    int ret;

    if (fgets (str, 1024, cofac) == NULL)
      break;
    ret = gmp_sscanf (str, "%ld %lu %Zd %Zd\n", &a, &b, R, A);
    if (ret != 4)
      break;
    cofac_list_add (L, a, b, R, A);
  }
  fclose (cofac);
  cofac_list_realloc (L, L->size);
  fprintf (stderr, "Read %zu cofactor pairs in %.0fs\n", L->size,
           seconds () - start);
  fflush (stderr);

  double start0 = seconds ();
  start = seconds ();
  find_smooth (L, lpb, lim, batch, verbose);
  fprintf (stderr, "Detecting %zu smooth cofactors took %.1f s\n", L->size,
           seconds() - start);

  start = seconds ();
  factor (L, pol, lpb[0], lpb[1], verbose);
  fprintf (stderr, "Factoring %zu smooth cofactors took %.1f s\n", L->size,
           seconds() - start);
  fprintf (stderr, "Detecting + factoring: %.1f s\n", seconds () - start0);

  cado_poly_clear (pol);
  mpz_clear (R);
  mpz_clear (A);
  cofac_list_clear (L);

  if (batch0_file != NULL)
    fclose (batch[0]);

  if (batch1_file != NULL)
    fclose (batch[1]);

  /* Clear the ecm addition chains that are stored as global variables
     to avoid re-computations. */
  free_saved_chains ();

  param_list_clear (pl);

  return 0;
}
