#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <ctype.h>

#include "macros.h"
#include "fb.h"



// Return the last character read, or EOF.
static
int skip_spaces(FILE *file)
{
  int c;
  while (isspace(c = getc(file)));
  return ungetc(c, file);
}


// Initialize a factor base, reading the ideals from a file and computing
// the corresponding lambda using the basis of the given q-lattice.
// Return 1 if successful.
int factor_base_init(factor_base_ptr FB, const char *filename, unsigned maxdeg)
{
  FILE *file;
  file = fopen(filename, "r");
  if (file == NULL) {
    perror("Error reading factor base");
    return 0;
  }

  FB->n    = 0;
  FB->elts = NULL;
  unsigned alloc = 0;
  for (fbideal_ptr gothp = NULL; skip_spaces(file) != EOF; ++FB->n, ++gothp) {
    // Need realloc?
    if (alloc <= FB->n) {
      alloc = alloc ? alloc * 2 : 256;
      FB->elts = (fbideal_t *)realloc(FB->elts, alloc*sizeof(fbideal_t));
      ASSERT_ALWAYS(FB->elts != NULL);
      gothp = FB->elts[FB->n];
    }

    // Get p.
    if (!fbprime_inp(gothp->p, file)) {
      fprintf(stderr, "Error parsing factor base %s.\n", filename);
      fprintf(stderr, "  The error occured at the %u-th ideal.\n", FB->n);
      fclose(file);
      return 0;
    }
    gothp->degp = fbprime_deg(gothp->p);
    if (maxdeg && gothp->degp > maxdeg) {
      --FB->n;
      break;
    }

    // Remove spaces.
    if (skip_spaces(file) == EOF) {
      fprintf(stderr, "Error parsing factor base %s.\n", filename);
      fprintf(stderr, "  The error occured at the %u-th ideal.\n", FB->n);
      fclose(file);
      return 0;
    }

    // Read ":".
    if (getc(file) != ':') {
      fprintf(stderr, "Error parsing ':' in factor base %s.\n", filename);
      fprintf(stderr, "  The error occured at the %u-th ideal.\n", FB->n);
      fclose(file);
      return 0;
    }

    // Get r.
    if (!fbprime_inp(gothp->r, file)) {
      fprintf(stderr, "Error parsing factor base %s.\n", filename);
      fprintf(stderr, "  The error occured at the %u-th ideal.\n", FB->n);
      fclose(file);
      return 0;
    }
  }

  FB->elts = (fbideal_t *)realloc(FB->elts, FB->n*sizeof(fbideal_t));
  ASSERT_ALWAYS(FB->elts != NULL);
  fclose(file);
  return 1;
}


// Precompute lambda for each element of the factor base.
void factor_base_precomp_lambda(factor_base_ptr FB, qlat_srcptr qlat,
        sublat_ptr sublat)
{
  ij_t ijmod;
  if (use_sublat(sublat))
    ij_set_16(ijmod, sublat->modulus);

  for (fbideal_ptr gothp = *FB->elts; gothp != FB->elts[FB->n]; ++gothp) {
    fbprime_t t0, t1;
    fbprime_mulmod(t0, qlat->b0, gothp->r, gothp->p);
    fbprime_sub   (t0, qlat->a0, t0);
    fbprime_rem   (t0, t0, gothp->p);
    if ((gothp->proj = fbprime_is_zero(t0))) {
      // This is a projective root! Yurk!
      // TODO: We'll have to handle these, someday.
      continue;
    }
    fbprime_invmod(t0, t0, gothp->p);
    fbprime_mulmod(t1, qlat->b1, gothp->r, gothp->p);
    fbprime_sub   (t1, t1, qlat->a1);
    fbprime_rem   (t1, t1, gothp->p);
    fbprime_mulmod(gothp->lambda, t0, t1, gothp->p);

    // In case we are using sublattices, precompute 1/p mod modulus
    if (use_sublat(sublat)) {
      int ret;
      ij_t tmp;
      ij_set_fbprime(tmp, gothp->p);
      ij_rem(tmp, tmp, ijmod);
      ret = ij_invmod(gothp->tildep, tmp, ijmod);
      if (!ret)
        ij_set_zero(gothp->tildep);
    }
  }
}


// Clean up memory.
void factor_base_clear(factor_base_ptr FB)
{
  free(FB->elts);
}
