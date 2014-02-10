#include "cado.h"
#include "utils.h"
#include "tests_common.h"

int
test_findroot (unsigned int nb)
{
  int err = 0;
  uint64_t a;
  uint64_t b;
  unsigned long p, r;
  mpz_t tp, ta, tb;
  mpz_init (tp);
  mpz_init (ta);
  mpz_init (tb);

  for (unsigned int i = 0; i < nb; i++)
  {
    a = random_int64 ();
    b = random_uint64 ();
    mpz_set_int64 (ta, a);
    mpz_set_uint64 (tb, b);

    mpz_set_ui (tp, lrand48());
    mpz_nextprime(tp, tp);
    p = mpz_get_ui (tp);

    r = findroot (a, b, p);
    
    if (mpz_invert (tb, tb, tp) == 0)
    {
      if (r != p)
      {
        gmp_fprintf (stderr, "ERROR: a=%" PRId64 " b=%" PRIu64" p=%" PRpr "\n"
                     "Got r=%" PRpr " instead of %" PRpr "\n", a, b, p, r, p);
        err++;
      }
    }
    else
    {
      mpz_mul (tb, ta, tb);
      mpz_mod (tb, tb, tp);

      mp_limb_t r2 = mpz_get_ui (tb);
      if (r2 != r)
      {
        gmp_fprintf (stderr, "ERROR: a=%" PRId64 " b=%" PRIu64" p=%" PRpr "\n"
                     "Got r=%" PRpr " instead of %Mx\n", a, b, p, r, r2);
        err++;
      }
    }
  }
  mpz_clear (tp);
  mpz_clear (ta);
  mpz_clear (tb);
  return err;
}

int
main (int argc, const char *argv[])
{
  int err;
  unsigned long iter = 10000;

  tests_common_cmdline(&argc, &argv, PARSE_SEED | PARSE_ITER);
  tests_common_get_iter(&iter);

  err = test_findroot (iter);
  
  if (err)
    fprintf (stderr, "# %d erro%s found\n", err, (err == 1) ? "r" : "rs");
  tests_common_clear();
  return (err) ? EXIT_FAILURE : EXIT_SUCCESS;
}
