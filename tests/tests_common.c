#include "cado.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <errno.h>
#include "tests_common.h"
#include "macros.h"

static int rng_state_inited = 0;
gmp_randstate_t state;
static int parsed_iter = 0;
unsigned long iter;

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


/* Set *output to the value from -iter, if -iter was given,
   otherwise do nothing */
void
tests_common_get_iter(unsigned long *output)
{
  if (parsed_iter)
    *output = iter;
}

static int
parse_l(long *result, const char *s)
{
  char *endptr;
  long r = strtol(s, &endptr, 10);
  if (s[0] == '\0' || endptr[0] != '\0' || errno == ERANGE)
    return 0;
  *result = r;
  return 1;
}

static int
parse_ul(unsigned long *result, const char *s)
{
  char *endptr;
  unsigned long r = strtoul(s, &endptr, 10);
  if (s[0] == '\0' || endptr[0] != '\0' || errno == ERANGE)
    return 0;
  *result = r;
  return 1;
}

void
tests_common_cmdline(int *argc, const char ***argv, const uint64_t flags)
{
  long seed;

  if ((flags & PARSE_SEED) != 0)
    seed = time (NULL);

  while (1) {
    const char *name = "-seed";
    if ((flags & PARSE_SEED) != 0 && (*argc) > 1 && 
        strcmp(name, (*argv)[1]) == 0) {
      if ((*argc) <= 2) {
        fprintf (stderr, "No value given for %s parameter\n", name);
        exit(EXIT_FAILURE);
      }
      if (!parse_l(&seed, (*argv)[2])) {
        fprintf (stderr, "Invalid value \"%s\" given for %s parameter\n",
                 (*argv)[2], name);
        exit(EXIT_FAILURE);
      }
      *argc -= 2;
      *argv += 2;
      continue;
    } 

    name = "-iter";
    if ((flags & PARSE_ITER) != 0 && (*argc) > 1 && 
        strcmp(name, (*argv)[1]) == 0) {
      if ((*argc) <= 2) {
        fprintf (stderr, "No value given for %s parameter\n", name);
        exit(EXIT_FAILURE);
      }
      if (!parse_ul(&iter, (*argv)[2])) {
        fprintf (stderr, "Invalid value \"%s\" given for %s parameter\n",
                 (*argv)[2], name);
        exit(EXIT_FAILURE);
      }
      *argc -= 2;
      *argv += 2;
      parsed_iter = 1;
      printf ("Using %lu iterations\n", iter);
      continue;
    }
    break;
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
