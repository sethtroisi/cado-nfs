#include "cado.h"
#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>
#include "tests_common.h"
#include "gcd.h"
#include "modredc_ul.h"

static int verbose = 0;

static void
print_residues_ul(const char *prefix, const unsigned long *a, size_t len,
                  const unsigned long m)
{
  if (!verbose)
    return;
  printf("%s: ", prefix);
  for(size_t i = 0; i < len; i++) {
    printf("Mod(%lu,%lu)%s", a[i], m, i + 1 < len ? ", " : "");
  }
  printf("\n");
}

/* Returns 1 if a1[i]*a2[i] == 1, for 0 <= i < len,
   returns 0 otherwise. */
static int
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

static int
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
  
  print_residues_ul("Input numbers", a, len, modredcul_getmod_ul(m));
  int rc = modredcul_batchinv_ul(r, a, c, len, m);
  if (rc != 0)
    print_residues_ul("Inverses", r, len, modredcul_getmod_ul(m));
  
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

  print_residues_ul("Input numbers after dividing by gcd", a, len, modredcul_getmod_ul(m));
  rc = modredcul_batchinv_ul(r, a, c, len, m);
  if (rc == 0) {
    fprintf (stderr, "After dividing out common factors, modredcul_batchinv()"
             " still returned 0\n");
    ok = 0;
    goto finish;
  }
  print_residues_ul("Inverses", r, len, modredcul_getmod_ul(m));
  
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

static int
check_Q_to_Fp(const unsigned long q,
              const unsigned long num, const unsigned long den,
              const unsigned long k, const unsigned long p)
{
  mpz_t r, m;

  if (q >= p) {
    fprintf (stderr, "Error, %lu >= %lu\n", q, p);
    return 0;
  }

  mpz_init(r);
  mpz_init(m);

  mpz_set_ui(m, p);
  mpz_set_ui(r, den);
  mpz_neg(r, r);
  mpz_mul_2exp(r, r, k);
  mpz_invert(r, r, m);
  mpz_mul_ui(r, r, num);
  mpz_mod(r, r, m);
  
  int ok = mpz_cmp_ui(r, q) == 0;
  if (!ok) {
    gmp_fprintf (stderr, "Error, %lu != %Zd for num=%lu, den=%lu, k=%lu, p=%lu\n",
                 q, r, num, den, k, p);
  }
  
  mpz_clear(m);
  mpz_clear(r);
  return ok;
}

static int
check_modredcul_batch_Q_to_Fp(unsigned long *r, const unsigned long num,
                              const unsigned long den, const unsigned long k,
                              const unsigned long *p, const size_t n)
{
  print_residues_ul("Trying modredcul_batch_Q_to_Fp() with p=", p, n, den);
  int ok = modredcul_batch_Q_to_Fp(r, num, den, k, p, n);
  
  if (!ok) {
    fprintf(stderr, "Error, modredcul_batch_Q_to_Fp() returned error\n");
  }
  
  for (size_t i = 0; i < n; i++) {
    ok = check_Q_to_Fp(r[i], num, den, k, p[i]);
  }
  return ok;
}

static int
test_modredcul_batch_Q_to_Fp(const unsigned long num, const unsigned long den,
                             const unsigned long k, const size_t n)
{
  unsigned long r[n], p[n];

  /* Try small p values */
  unsigned long m = 1;
  for (size_t i = 0; i < n; i++) {
    do {
      m += 2;
    } while (gcd_ul(m, den) > 1);
    p[i] = m;
  }
  int ok = check_modredcul_batch_Q_to_Fp(r, num, den, k, p, n);

  /* Try random p values */
  for (size_t i = 0; i < n; i++) {
    unsigned long m;
    do {
      m = random_uint64();
    } while (m == 0 || (k > 0 && m % 2 == 0) || gcd_ul(m, den) > 1);
    p[i] = m;
  }
  ok &= check_modredcul_batch_Q_to_Fp(r, num, den, k, p, n);

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
  
  for (unsigned long i = 1; ok && i < iter; i++) {
    const unsigned long num = (i < 20) ? i : random_uint64();
    for (unsigned long j = 3; ok && j < 2*iter; j += 2) {
      const unsigned long den = (j < 20) ? j : (random_uint64() | 1UL);
      ok &= test_modredcul_batch_Q_to_Fp(num, den, 0, iter);
      if (num % 2 == 1) {
        ok &= test_modredcul_batch_Q_to_Fp(num, den, 1, iter);
        ok &= test_modredcul_batch_Q_to_Fp(num, den, 2, iter);
      }
    }
  }
  
  tests_common_clear();
  if (ok)
    exit (EXIT_SUCCESS);
  else
    exit (EXIT_FAILURE);
}
