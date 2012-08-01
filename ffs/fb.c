#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <ctype.h>

#include "macros.h"
#include "fb.h"
#include "qlat.h"


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


/* 
 * Compute a basis of the vector space of polynomials that are multiples
 * of p and of degree less than J.
 * The basis is echelon on the monomial basis, starting with high
 * degree. This echelon form is normalized in the sense that it is lower
 * triangular, full of ones.
 *
 * p must be monic.
 *
 * The "basis" parameter must be preallocated to contain at least
 * J-degp elements.
 */
static 
void normalized_echelon_multiples(ij_t *basis, fbprime_t p, int degp, int J)
{
    ASSERT(fbprime_deg(p) == degp);
    ASSERT(fbprime_is_monic(p));
    if (degp >= J)
        return;

    // Make it diagonal (classical reduced-echelon form).
    ij_set_fbprime(basis[0], p);
    for (int i = 1; i < J-degp; ++i) {
        fbprime_t tip;
        fbprime_mul_ti(tip, p, i);
        ij_set_fbprime(basis[i], tip);
        for (int j = i-1; j >= 0; --j) {
            fp_t c;
            ij_t aux;
            ij_get_coeff(c, basis[i], degp + j);
            ij_smul(aux, basis[j], c);
            ij_sub(basis[i], basis[i], aux);
        }
    }

    // Put 1's in the lower triangle.
    for (int i = 1; i < J-degp; ++i)
        for (int j = 0; j < i; ++j)
            ij_add(basis[i], basis[i], basis[j]);
}


// Store a factor base element corresponding to a small prime ideal.
// Precompute what we can at this early stage.
static
void push_small_ideal(small_factor_base_ptr FB, fbprime_t p, fbprime_t r,
    unsigned degp, int power, unsigned I, unsigned J, sublat_ptr sublat)
{
  unsigned alloc = FB->alloc;
  if (alloc <= FB->n) {
    alloc = alloc ? alloc * 2 : 256;
    FB->elts = (small_fbideal_t *)realloc(FB->elts, alloc*sizeof(small_fbideal_t));
    ASSERT_ALWAYS(FB->elts != NULL);
    FB->alloc=alloc;
  }

  fbprime_set(FB->elts[FB->n]->q, p);
  fbprime_set(FB->elts[FB->n]->r, r);
  FB->elts[FB->n]->degp = degp;
  FB->elts[FB->n]->degq = fbprime_deg(p);
  FB->elts[FB->n]->power = power;
  // In case we are using sublattices, precompute 1/p mod modulus
  if (use_sublat(sublat)) {
    int ret;
    ij_t tmp;
    ij_t ijmod;
    if (use_sublat(sublat))
      ij_set_16(ijmod, sublat->modulus);
    ij_set_fbprime(tmp, FB->elts[FB->n]->q);
    ij_rem(tmp, tmp, ijmod);
    ret = ij_invmod(FB->elts[FB->n]->tildep, tmp, ijmod);
    if (!ret)
      ij_set_zero(FB->elts[FB->n]->tildep);
  }

  FB->elts[FB->n]->basis            = (ij_t *)malloc((I+J)*sizeof(ij_t));
  FB->elts[FB->n]->adjustment_basis = (ij_t *)malloc(   J *sizeof(ij_t));
  FB->elts[FB->n]->projective_basis = (ij_t *)malloc(   J *sizeof(ij_t));
  ASSERT_ALWAYS(FB->elts[FB->n]->basis            != NULL);
  ASSERT_ALWAYS(FB->elts[FB->n]->adjustment_basis != NULL);
  ASSERT_ALWAYS(FB->elts[FB->n]->projective_basis != NULL);

  normalized_echelon_multiples(FB->elts[FB->n]->projective_basis, p,
      fbprime_deg(p), J);

  // other fields are recomputed for each special-q.
 
  FB->n++;
}

static
void push_large_ideal(large_factor_base_ptr FB, fbprime_t p, fbprime_t r,
    unsigned degp, MAYBE_UNUSED sublat_ptr sublat)
{
  unsigned alloc = FB->alloc;
  if (alloc <= FB->n) {
    alloc = alloc ? alloc * 2 : 256;
    FB->elts = (large_fbideal_t *)realloc(FB->elts, alloc*sizeof(large_fbideal_t));
    ASSERT_ALWAYS(FB->elts != NULL);
    FB->alloc=alloc;
  }

  fbprime_set(FB->elts[FB->n]->p, p);
  fbprime_set(FB->elts[FB->n]->r, r);
  FB->elts[FB->n]->data = degp;

  FB->n++;
}

static
void push_ideal(large_factor_base_ptr LFB, small_factor_base_ptr SFB,
    fbprime_t p, fbprime_t r, unsigned degp, int power, unsigned min_degp,
    unsigned I, unsigned J, sublat_ptr sublat)
{
  if (degp < min_degp)
    push_small_ideal(SFB, p, r, degp, power, I, J, sublat);
  else {
    if (power) {
      fprintf(stderr, "Warning: large power in factor base. Ignoring...\n");
      return;
    }
    push_large_ideal(LFB, p, r, degp, sublat);
  }
}



// Initialize a factor base, reading the ideals from a file.
// The factor base should be sorted by increasing degree, starting at degree
// sorted_min_degp; ideals of lower degree don't have to be sorted.
// If max_degp is nonzero, the factor base is read until reaching ideals of
// degree max_degp.
// FIXME: fix max_degp so that this bound is exclusive, not inclusive!
// Return 1 if successful.
int factor_base_init(large_factor_base_ptr LFB, small_factor_base_ptr SFB,
    const char *filename, unsigned sorted_min_degp, unsigned max_degp,
    unsigned I, unsigned J, sublat_ptr sublat)
{
  FILE *file;
  file = fopen(filename, "r");
  if (file == NULL) {
    perror("Error reading factor base");
    return 0;
  }

  LFB->n=0;
  LFB->alloc=0;
  LFB->elts=NULL;
  
  SFB->n=0;
  SFB->alloc=0;
  SFB->elts=NULL;
  
  unsigned last_degp = 0;
  int previous_prime_degp = -1;
  int cpt = 0;
  int fbb_ok = 0;

  while(skip_spaces_and_comments(file) != EOF) {
    fbprime_t p, r;
    unsigned degp;
    int power;

    cpt++;
    // Get p.
    if (!fbprime_inp(p, file)) {
      fprintf(stderr, "Error parsing factor base %s.\n", filename);
      fprintf(stderr, "  The error occured at the %u-th ideal.\n", cpt);
      fclose(file);
      return 0;
    }
    degp = fbprime_deg(p);
    if (degp == max_degp)
      fbb_ok = 1;
    if (max_degp && degp > max_degp) {
      break;
    }
    ASSERT_ALWAYS(last_degp < sorted_min_degp || last_degp <= degp);
    last_degp = degp;

    // Read ":".
    if (getc(file) != ':') {
      fprintf(stderr, "Error parsing ':' in factor base %s.\n", filename);
      fprintf(stderr, "  The error occured at the %u-th ideal.\n", cpt);
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

    power = 0;
    // Read the number of deg p to subtract (longversion)
    if (longversion) { 
      int ret, n0, n1;
      ret = fscanf(file, "%d,%d", &n1, &n0);
      if (ret != 2) {
        fprintf(stderr, "Error parsing n1,n0 in fb %s.\n", filename);
        fprintf(stderr, "  The error occured at the %u-th ideal.\n", cpt);
        fclose(file);
        return 0;
      }
      if (getc(file) != ':') {
        fprintf(stderr, "Error parsing ':' in factor base %s.\n", filename);
        fprintf(stderr, "  The error occured at the %u-th ideal.\n", cpt);
        fclose(file);
        return 0;
      }
      // n0 != 0  iff  we have a power. In that case, we should use the
      // degree of the corresponding prime. We assume that it is the last
      // non-power that we met.
      if (n0 == 0) {
        degp = (n1-n0)*fbprime_deg(p);
        previous_prime_degp = fbprime_deg(p);
      } else {
        degp = (n1-n0)*previous_prime_degp;
        power = 1;
        ASSERT_ALWAYS (previous_prime_degp != -1);
        ASSERT_ALWAYS (fbprime_deg(p) % previous_prime_degp == 0);
      }
    } else { 
      // short version is always prime
      degp = fbprime_deg(p);
      previous_prime_degp = degp;
    }


    // Remove spaces.
    if (skip_spaces_and_comments(file) == EOF) {
      fprintf(stderr, "Error parsing factor base %s.\n", filename);
      fprintf(stderr, "  The error occured at the %u-th ideal.\n", cpt);
      fclose(file);
      return 0;
    }

    // Get r.
    if (!fbprime_inp(r, file)) {
      fprintf(stderr, "Error parsing factor base %s.\n", filename);
      fprintf(stderr, "  The error occured at the %u-th ideal.\n", cpt);
      fclose(file);
      return 0;
    }

    push_ideal(LFB, SFB, p, r, degp, power, sorted_min_degp, I, J, sublat);

    // Other r's ?
    while ((c = getc(file)) == ',') {
      if (!fbprime_inp(r, file)) {
        fprintf(stderr, "Error parsing factor base %s.\n", filename);
        fprintf(stderr, "  The error occured at the %u-th ideal.\n", cpt);
        fclose(file);
        return 0;
      }
    
      push_ideal(LFB, SFB, p, r, degp, power, sorted_min_degp, I, J, sublat);
    }
    ungetc(c, file);
  }

  // The allocated FB might be twice too large. Let's realloc.
  LFB->elts = (large_fbideal_t *)realloc(LFB->elts,
      LFB->n*sizeof(large_fbideal_t));
  ASSERT_ALWAYS(LFB->elts != NULL);
  LFB->alloc = LFB->n;
  
  SFB->elts = (small_fbideal_t *)realloc(SFB->elts,
      SFB->n*sizeof(small_fbideal_t));
  ASSERT_ALWAYS(SFB->elts != NULL);
  SFB->alloc = SFB->n;

  fclose(file);
  
  // Sanity check: if we did not read any ideal of the maximum degree,
  // there is probably a problem with the factor base file.
  if (!fbb_ok) {
    fprintf(stderr, "Warning: the factor base file %s does not contain any ideal of maximu allowed\n"
        "  degree. The file might be corrupted, or was prepared for a smaller factor\n"
        "  base bound. Running again 'makefb' is recommended.\n",
        filename);
  }

  return 1;
}


void small_factor_base_precomp(small_factor_base_ptr FB,
                               unsigned I, unsigned J, qlat_srcptr qlat)
{
  for (unsigned i = 0; i < FB->n; ++i) {
    small_fbideal_ptr gothp = FB->elts[i];
    compute_lambda(gothp->lambda, gothp->q, gothp->r, qlat);
    if (fbprime_eq(gothp->lambda, gothp->q)) {
      gothp->proj = 1;
    } else {
      gothp->proj = 0;
      ijbasis_compute_small(gothp->basis, gothp->adjustment_basis,
          gothp, gothp->lambda, I, J);
    }
  }
}


// Return the largest degree of the ideals in the factor base.
// /!\ Assume that the factor base is sorted, at least for the ideals higher
//     degree. The highest degree should thus be the degree of the last
//     ideal in the factor base.
// FIXME: should this be inclusive or exclusive?
unsigned factor_base_max_degp(large_factor_base_srcptr FB)
{
  return fbideal_deg(FB->elts[FB->n-1]);
}


void factor_base_clear(large_factor_base_ptr LFB, small_factor_base_ptr SFB)
{
  for (unsigned i = 0; i < SFB->n; ++i) {
    free(SFB->elts[i]->basis);
    free(SFB->elts[i]->adjustment_basis);
    free(SFB->elts[i]->projective_basis);
  }
  free(LFB->elts);
  free(SFB->elts);
}

double expected_hit_number(large_factor_base_srcptr LFB,
    unsigned I, unsigned J)
{
  double powers[I+J];
  powers[0] = 1.0;
  for (unsigned i = 1; i < I+J; ++i)
    powers[i] = powers[i-1]*FP_SIZE;
  double nb = 0;
  large_fbideal_t * ptr = LFB->elts; 
  for (unsigned i = 0; i < LFB->n; ++i, ptr++) {
    unsigned L = fbideal_deg(ptr[0]);
    if (L <= I+J)
      nb += powers[I+J-L];
    else
      nb += 1/powers[L-(I+J)];
  }
  return nb;
}
