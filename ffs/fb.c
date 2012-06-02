#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <ctype.h>

#include "macros.h"
#include "fb.h"


// Return the last character read, or EOF.
static
int skip_spaces_and_comments(FILE *file)
{
  int c;
  
  for (;;) {
    while (isspace(c = getc(file)));
    if (c == '#') {
      // skip the end of the line
      do {
        c = fgetc(file);
      } while ((c != '\n') && (c != EOF));
    } else
      break;
  }
  return ungetc(c, file);
}


// Initialize a factor base, reading the ideals from a file.
// The factor base should be sorted by increasing degree, starting at degree
// sorted_min_degp; ideals of lower degree don't have to be sorted.
// If max_degp is nonzero, the factor base is read until reaching ideals of
// degree max_degp.
// FIXME: fix max_degp so that this bound is exclusive, not inclusive!
// Return 1 if successful.
int factor_base_init(factor_base_ptr FB, const char *filename,
                     unsigned sorted_min_degp, unsigned max_degp)
{
  FILE *file;
  file = fopen(filename, "r");
  if (file == NULL) {
    perror("Error reading factor base");
    return 0;
  }

  FB->n    = 0;
  FB->elts = NULL;
  unsigned last_degp = 0;
  unsigned alloc     = 0;
  int previous_prime_degp = -1;
  for (fbideal_ptr gothp = NULL;
      skip_spaces_and_comments(file) != EOF;
      ++FB->n, ++gothp) {
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
    if (max_degp && gothp->degp > max_degp) {
      --FB->n;
      break;
    }
    ASSERT_ALWAYS(last_degp < sorted_min_degp || last_degp <= gothp->degp);
    last_degp = gothp->degp;

    // Read ":".
    if (getc(file) != ':') {
      fprintf(stderr, "Error parsing ':' in factor base %s.\n", filename);
      fprintf(stderr, "  The error occured at the %u-th ideal.\n", FB->n);
      fclose(file);
      return 0;
    }

    int longversion;
    int c;
    if ((c = getc(file)) == ' ') {
      longversion = 0;
    } else {
      longversion = 1;
      ungetc(c, file);
    }

    gothp->power = 0;
    // Read the number of deg p to subtract (longversion)
    if (longversion) { 
      int ret, n0, n1;
      ret = fscanf(file, "%d,%d", &n1, &n0);
      if (ret != 2) {
        fprintf(stderr, "Error parsing n1,n0 in fb %s.\n", filename);
        fprintf(stderr, "  The error occured at the %u-th ideal.\n", FB->n);
        fclose(file);
        return 0;
      }
      if (getc(file) != ':') {
        fprintf(stderr, "Error parsing ':' in factor base %s.\n", filename);
        fprintf(stderr, "  The error occured at the %u-th ideal.\n", FB->n);
        fclose(file);
        return 0;
      }
      // n0 != 0  iff  we have a power. In that case, we should use the
      // degree of the corresponding prime. We assume that it is the last
      // non-power that we met.
      if (n0 == 0) {
        gothp->degp = (n1-n0)*fbprime_deg(gothp->p);
        previous_prime_degp = fbprime_deg(gothp->p);
      } else {
        gothp->degp = (n1-n0)*previous_prime_degp;
        gothp->power = 1;
        ASSERT_ALWAYS (previous_prime_degp != -1);
        ASSERT_ALWAYS (fbprime_deg(gothp->p) % previous_prime_degp == 0);
      }
    } else { 
      // short version is always prime
      previous_prime_degp = fbprime_deg(gothp->p);
    }


    // Remove spaces.
    if (skip_spaces_and_comments(file) == EOF) {
      fprintf(stderr, "Error parsing factor base %s.\n", filename);
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

    // Other r's ?
    while ((c = getc(file)) == ',') {
      fbprime_t r;
      if (!fbprime_inp(r, file)) {
        fprintf(stderr, "Error parsing factor base %s.\n", filename);
        fprintf(stderr, "  The error occured at the %u-th ideal.\n", FB->n);
        fclose(file);
        return 0;
      }
      fbideal_ptr oldgothp = gothp;
      FB->n++;
      gothp++;
      // Need realloc?
      if (alloc <= FB->n) {
        alloc = alloc ? alloc * 2 : 256;
        FB->elts = (fbideal_t *)realloc(FB->elts, alloc*sizeof(fbideal_t));
        ASSERT_ALWAYS(FB->elts != NULL);
        gothp = FB->elts[FB->n];
        oldgothp = FB->elts[FB->n-1];
      }
      fbprime_set(gothp->p, oldgothp->p);
      gothp->degp = oldgothp->degp;
      fbprime_set(gothp->r, r);
    }
    ungetc(c, file);
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
    fbprime_t t0, t1, r;
    int was_proj = fbprime_deg(gothp->r) == fbprime_deg(gothp->p);

    if (was_proj) {
        fbprime_sub(r, gothp->r, gothp->p);
        fbprime_mulmod(t0, qlat->a0, r, gothp->p);
        fbprime_sub   (t0, t0, qlat->b0);
        fbprime_rem   (t0, t0, gothp->p);
    } else {
        fbprime_mulmod(t0, qlat->b0, gothp->r, gothp->p);
        fbprime_sub   (t0, qlat->a0, t0);
        fbprime_rem   (t0, t0, gothp->p);
    }
    if ((gothp->proj = fbprime_is_zero(t0))) {
      // This is a projective root in the q-lattice.
      continue;
    }
    int ret = fbprime_invmod(t0, t0, gothp->p);
    if (!ret) { // this can happen for powers
        // FIXME: is this condition correct ?
        gothp->proj = 1;
        continue;
    }
    if (was_proj) {
        fbprime_mulmod(t1, qlat->a1, r, gothp->p);
        fbprime_sub   (t1, qlat->b1, t1);
        fbprime_rem   (t1, t1, gothp->p);
    } else {
        fbprime_mulmod(t1, qlat->b1, gothp->r, gothp->p);
        fbprime_sub   (t1, t1, qlat->a1);
        fbprime_rem   (t1, t1, gothp->p);
    }
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


// Return the largest degree of the ideals in the factor base.
// /!\ Assume that the factor base is sorted, at least for the ideals higher
//     degree. The highest degree should thus be the degree of the last
//     ideal in the factor base.
// FIXME: should this be inclusive or exclusive?
unsigned factor_base_max_degp(factor_base_srcptr FB)
{
  return FB->elts[FB->n-1]->degp;
}


// Clean up memory.
void factor_base_clear(factor_base_ptr FB)
{
  free(FB->elts);
}
