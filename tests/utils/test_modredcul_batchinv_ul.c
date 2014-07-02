#include "cado.h"
#include <stdlib.h>
#include "tests_common.h"
#include "gcd.h"
#include "modredc_ul.h"

static int verbose = 0;

void
print_residues_ul(const char *prefix, const unsigned long *a, size_t len,
                  const modulusredcul_t m)
{
  if (!verbose)
    return;
  printf("%s: ", prefix);
  for(size_t i = 0; i < len; i++) {
    printf("Mod(%lu,%lu)%s", a[i], modredcul_getmod_ul(m),
           i + 1 < len ? ", " : "");
  }
  printf("\n");
}

/* Returns 1 if a1[i]*a2[i] == 1, for 0 <= i < len,
   returns 0 otherwise. */
int
check_product_c_ul(unsigned long *a1, unsigned long *a2,
                   const unsigned long orig_c, const size_t len,
                   const modulusredcul_t m)
{
   int ok = 1;
   mpz_t prod;
   const unsigned long c = orig_c % modredcul_getmod_ul(m);

   mpz_init(prod);
   for (size_t i = 0; i < len; i++) {
     mpz_set_ui(prod, a1[i]);
     mpz_mul_ui(prod, prod, a2[i]);
     unsigned long r = mpz_tdiv_ui(prod, modredcul_getmod_ul(m));

     if(r != c) {
       ok = 0;
       fprintf (stderr, "Inverse is incorrect: %lu * %lu = %lu != %lu (mod %lu)\n",
                a1[i], a2[i], r, c, modredcul_getmod_ul(m));
       break;
     }
   }
   mpz_clear(prod);
   return ok;
}

int
test_modredc_batchinv_ul (const size_t len, const unsigned long c)
{
  unsigned long *a, *r;
  modulusredcul_t m;
  int ok = 1;
  
  /* Random, odd modulus */
  modredcul_initmod_ul(m, random_uint64() | 1);

  a = (unsigned long *) malloc(len * sizeof(unsigned long));
  r = (unsigned long *) malloc(len * sizeof(unsigned long));
  
  for (size_t i = 0; i < len; i++) {
    a[i] = random_uint64() % modredcul_getmod_ul(m);
    if (a[i] == 0)
      a[i] = 1;
  }
  
  print_residues_ul("Input numbers", a, len, m);
  int rc = modredcul_batchinv_ul(r, a, c, len, m);
  if (rc != 0)
    print_residues_ul("Inverses", r, len, m);
  
  /* Check that product of input and output is 1 */
  if (rc == 1) {
    if (!check_product_c_ul(r, a, c, len, m)) {
      fprintf (stderr, "Product of input and its inverse it not %lu\n", c);
      ok = 0;
      goto finish;
    }
  }
  
  int had_common_factor = 0;
  unsigned long g;
  for (size_t i = 0; i < len; i++) {
    g = gcd_ul(a[i], modredcul_getmod_ul(m));
    /* May have to divide out factors with greater multiplicity in a than
       they have in m */
    while (g != 1) {
      had_common_factor = 1;
      /* Divide out gcd so we can try again with a set of residues for which
         the inverses exist */
      a[i] /= g;
      g = gcd_ul(a[i], modredcul_getmod_ul(m));
    }
  }

  if (had_common_factor && rc != 0) {
    fprintf(stderr, "Input numbers had common factor but modredcul_batchinv()"
            " returned %d", rc);
    ok = 0;
    goto finish;
  }

  print_residues_ul("Input numbers after dividing by gcd", a, len, m);
  rc = modredcul_batchinv_ul(r, a, c, len, m);
  if (rc == 0) {
    fprintf (stderr, "After dividing out common factors, modredcul_batchinv()"
             " still returned 0\n");
    ok = 0;
    goto finish;
  }
  print_residues_ul("Inverses", r, len, m);
  
  /* Check that product of input and output is 1 */
  if (!check_product_c_ul(r, a, c, len, m)) {
    fprintf (stderr, "Product of input and its inverse it not %lu\n", c);
    ok = 0;
    goto finish;
  }
  
  finish:
  free(a);
  free(r);
  return ok;
}

int
main (int argc, const char *argv[])
{
  unsigned long iter = 20;
  int ok = 1;
  tests_common_cmdline (&argc, &argv, PARSE_SEED | PARSE_ITER | PARSE_VERBOSE);
  tests_common_get_iter (&iter);
  verbose = tests_common_get_verbose();
  
  for (unsigned long i = 0; ok && i < iter; i++) {
    ok &= test_modredc_batchinv_ul(i, 0);
    ok &= test_modredc_batchinv_ul(i, 1);
    ok &= test_modredc_batchinv_ul(i, random_uint64());
  }
  
  tests_common_clear();
  if (ok)
    exit (EXIT_SUCCESS);
  else
    exit (EXIT_FAILURE);
}
