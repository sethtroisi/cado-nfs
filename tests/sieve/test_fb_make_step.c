#include "cado.h"
#include "tests_common.h"
#include "sieve/fb.h"

int main(int argc, const char **argv)
{
  fbprime_t fbb;
  fbprime_t steps[256];
  unsigned long i, iter = 1000;

  tests_common_cmdline(&argc, &argv, PARSE_SEED | PARSE_ITER);
  tests_common_get_iter(&iter);

  /* Try small FBB */
  for (fbb = 0; fbb < 20; fbb++) {
    double scale = 0.5 + 5. * drand48();
    fb_make_steps(steps, fbb, scale);
  }
  
  /* Try random FBB */
  for (i = 0; i < iter; i++) {
    fbb = lrand48();
    double scale = 0.5 + 5. * drand48();
    fb_make_steps(steps, fbb, scale);
  }

  tests_common_clear ();
  
  return 0;
}
