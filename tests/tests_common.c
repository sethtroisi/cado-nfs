#include "cado.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "tests_common.h"

static int rng_state_inited = 0;
gmp_randstate_t state;

/* Return non-zero iff |d2| is in the interval |d1| * (1 +- err_margin) */
int
cmp_double(const double d1, const double d2, const double err_margin)
{
  return fabs(d1) * (1. - err_margin) <= fabs(d2) && fabs(d2) <= fabs(d1) * (1. + err_margin);
}

int64_t
random_int64 ()
{
  return ((int64_t) mrand48 () << 32) + (int64_t) mrand48 ();
}

uint64_t
random_uint64 ()
{
  return ((uint64_t) lrand48 () << 33) + ((uint64_t) lrand48 () << 2) + ((uint64_t) lrand48 () & 3);
}

void
tests_common_cmdline(int *argc, const char ***argv, const uint64_t flags)
{
  long seed = time (NULL);
  while (1) {
    if ((flags & PARSE_SEED) != 0 && (*argc) > 1 && 
        strcmp("-seed", (*argv)[1]) == 0) {
      if ((*argc) <= 2) {
        fprintf (stderr, "No value given for -seed parameter\n");
        exit(EXIT_FAILURE);
      }
      char * endptr;
      seed = strtol((*argv)[2], &endptr, 10);
      if ((*argv)[2][0] == '\0' || endptr[0] != '\0') {
        fprintf (stderr, "Invalid value \"%s\" given for -seed parameter\n",
                 (*argv)[2]);
        exit(EXIT_FAILURE);
      }
      *argc -= 2;
      *argv += 2;
    } else {
      break;
    }

  }
  if ((flags & PARSE_SEED) != 0) {
    printf ("Using random seed=%ld\n", seed);
    srand48 (seed);
    gmp_randinit_default (state);
    gmp_randseed_ui (state, seed);
    rng_state_inited = 1;
  }
}

/* Clean up, free any memory that may have been allocated */
void
tests_common_clear()
{
  if (rng_state_inited) {
    gmp_randclear (state);
    rng_state_inited = 0;
  }
}
