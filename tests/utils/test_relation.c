#include "cado.h"
#include "utils.h"
#include "tests_common.h"

unsigned long
mpz_compute_r (mpz_t a, mpz_t b, mpz_t p)
{
  if (mpz_invert (b, b, p) == 0)
    return mpz_get_ui(p);
  else
  {
    mpz_mul (b, a, b);
    mpz_mod (b, b, p);

    return mpz_get_ui (b);
  }
 }

int
test_findroot (unsigned int nb)
{
  int err = 0;
  mpz_t tp, ta, tb;
  mpz_init (tp);
  mpz_init (ta);
  mpz_init (tb);

  for (unsigned int i = 0; i < nb; i++)
  {
    uint64_t a;
    uint64_t b;
    unsigned long p;

    mpz_set_ui (tp, lrand48());
    mpz_nextprime(tp, tp);
    p = mpz_get_ui (tp);

    a = random_int64 ();
    /* 5% of tests are for the case where b = 0 mod p (with b > 0)
     * We do not need to test for free relations as they never go throught
     * findroot
     */
    if (i < (nb / 20))
      b = lrand48() * p;
    else
      b = random_uint64 () + 1; /* b > 0 */
    mpz_set_int64 (ta, a);
    mpz_set_uint64 (tb, b);

    unsigned long r = findroot (a, b, p);

    unsigned long r2 = mpz_compute_r (ta, tb, tp);
    if (r != r2)
    {
      gmp_fprintf (stderr, "ERROR: a=%" PRId64 " b=%" PRIu64" p=%" PRpr "\n"
                   "Got r=%" PRpr " instead of %" PRpr "\n", a, b, p, r, r2);
      err++;
    }
  }
  mpz_clear (tp);
  mpz_clear (ta);
  mpz_clear (tb);
  return err;
}

int test_computeroots (unsigned int nb)
{
  int err = 0;
  relation_t t1[1], t2[1];
  mpz_t tp, ta, tb;
  mpz_init (tp);
  mpz_init (ta);
  mpz_init (tb);
  unsigned long p;

  for (unsigned int i = 0; i < nb; i++)
  {
    relation_init (t1);
    relation_init (t2);
    t1->a = random_int64 ();
    t1->b = random_uint64 () + 1; /* b > 0 */

    for (uint8_t k = 0; k <= lrand48() % 5; k++)
    {
      mpz_set_ui (tp, lrand48());
      mpz_nextprime(tp, tp);
      p = mpz_get_ui (tp);
      relation_add_prime(t1, RATIONAL_SIDE, p);
    }
    for (uint8_t k = 0; k <= lrand48() % 5; k++)
    {
      mpz_set_ui (tp, lrand48());
      mpz_nextprime(tp, tp);
      p = mpz_get_ui (tp);
      relation_add_prime(t1, ALGEBRAIC_SIDE, p);
    }

    relation_copy (t2, t1);
    computeroots (t1);

    for (uint8_t k = 0; k < t2->nb_ap; k++)
    {
      mpz_set_int64 (ta, t2->a);
      mpz_set_uint64 (tb, t2->b);
      mpz_set_ui (tp, t2->ap[k].p);
      unsigned long r = mpz_compute_r (ta, tb, tp);
      if (r != t1->ap[k].r)
      {
        gmp_fprintf (stderr, "ERROR: a=%" PRId64 " b=%" PRIu64" p=%" PRpr "\n"
                     "Got r=%" PRpr " instead of %" PRpr "\n", t2->a, t2->b,
                     t2->ap[k].p, t1->ap[k].r, r);
        err++;
      }
    }

    relation_clear (t1);
    relation_clear (t2);
  }

  mpz_clear (ta);
  mpz_clear (tb);
  mpz_clear (tp);

  return err;
}

int
check_str_err (const char * s1, const char * s2, mpz_t t)
{
  if (strcmp(s1, s2) != 0)
  {
    gmp_fprintf(stderr, "ERROR with integer %Zd: got \"%s\" instead of "
                        "\"%s\"\n", t, s2, s1);
    return 1;
  }
  else
    return 0;

}

int
test_conversion (unsigned int nb)
{
  int err = 0;
  char *s1, *s2, *tmp;
  s1 = (char *) malloc (25 * sizeof(char));
  s2 = (char *) malloc (25 * sizeof(char));
  mpz_t t;
  mpz_init(t);


  for (unsigned int i = 0; i < nb; i++)
  {
    uint64_t a = random_uint64 ();
    uint64_t b = random_int64 ();

    mpz_set_uint64 (t, a);

    mpz_get_str(s1, 10, t);
    tmp = u64toa10 (s2, a);
    *tmp = '\0';
    err += check_str_err (s1, s2, t);

    mpz_get_str(s1, 16, t);
    tmp = u64toa16 (s2, a);
    *tmp = '\0';
    err += check_str_err (s1, s2, t);


    mpz_set_int64 (t, b);

    mpz_get_str(s1, 10, t);
    tmp = d64toa10 (s2, b);
    *tmp = '\0';
    err += check_str_err (s1, s2, t);

    mpz_get_str(s1, 16, t);
    tmp = d64toa16 (s2, b);
    *tmp = '\0';
    err += check_str_err (s1, s2, t);
  }
  mpz_clear(t);
  free(s1);
  free(s2);

  return err;
}

int
main (int argc, const char *argv[])
{
  int err = 0;
  unsigned long iter = 10000;

  tests_common_cmdline(&argc, &argv, PARSE_SEED | PARSE_ITER);
  tests_common_get_iter(&iter);

  err += test_findroot (iter);
  err += test_computeroots (iter / 10);
  err += test_conversion (iter);

  if (err)
    fprintf (stderr, "# %d erro%s found\n", err, (err == 1) ? "r" : "rs");
  tests_common_clear();
  return (err) ? EXIT_FAILURE : EXIT_SUCCESS;
}
