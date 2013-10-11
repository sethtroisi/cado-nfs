#include "cado.h"

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <inttypes.h>
#include <ctype.h>
#include <time.h>
#include <limits.h> /* for CHAR_BIT */
#include <unistd.h>
#include <assert.h>
#include <fcntl.h>   /* for _O_BINARY */

#include "portability.h"
#include "macros.h"
#include "utils.h"
#include "relation.h"

#include "mpz_poly.h"
#include "filter_utils.h"

static index_t nrels_read, nrels_err, nrels_noprime;
cado_poly cpoly;
int verbose = 0;
int check_primality = 0; /* By default no primality check */

/* used for counting time in different processes */
timingstats_dict_t stats;

/* return 0 if an error exists in the factorization of a norm (rat or alg)
 * return -1 if the factorization are correct but al least one "ideal" is not
 *           prime
 * return 1 if the factorization are correct and all ideals are primes (or no
 *          primality check was asked)
 */
int
check_relation (earlyparsed_relation_ptr rel)
{
  mpz_t norm[2];
  mpz_init (norm[0]);
  mpz_init (norm[1]);
  int ret = 1;

  // compute the norm on alg and rat sides
  for(unsigned int side = 0 ; side < 2 ; side++)
  {
    cado_poly_side_ptr ps = cpoly->pols[side];
    mp_poly_homogeneous_eval_siui(norm[side], ps->f, ps->degree, rel->a, rel->b);
  }

  // check for correctness of the factorization of the norms
  for(weight_t i = 0; i < rel->nb ; i++)
  {
    unsigned int side = rel->primes[i].h;
    p_r_values_t p = rel->primes[i].p;
    exponent_t e = rel->primes[i].e;
    for (int j = 0; j < e; ++j)
    {
      if (mpz_fdiv_q_ui (norm[side], norm[side], p) != 0)
      {
        if (verbose != 0)
        {
          fprintf (stderr, "Given factor %" PRpr " with exponent %u does not "
                           "divide the norm on side %u for (%" PRId64 ", "
                           "%" PRIu64 ")\n", p, e, side, rel->a, rel->b);
        }
        ret = 0;
      }
    }
  }

  mpz_clear(norm[0]);
  mpz_clear(norm[1]);
  if (ret == 0)
    return ret;

  // check primality
  if (check_primality != 0)
  {
    for(weight_t i = 0; i < rel->nb ; i++)
    {
      p_r_values_t p = rel->primes[i].p;
      if (!modul_isprime((const long unsigned int *)&p))
      {
        unsigned int side = rel->primes[i].h;
        if (verbose != 0)
        {
          fprintf (stderr, "Given factor %" PRpr " is not prime on side %u for "
                           "(%" PRId64 ", %" PRIu64 ")\n", p, side, rel->a,
                           rel->b);
        }
        ret=-1;
      }
    }
  }

  return ret;
}

/* Callback function called by prempt_scan_relations */

static void *
thread_check (MAYBE_UNUSED void * context_data, earlyparsed_relation_ptr rel)
{
  nrels_read++;
  int ret = check_relation (rel);
  if (ret == 0)
  {
    nrels_err++;
    fprintf (stderr, "Line with a=%" PRId64 " b=%" PRIu64 " failed due to an "
                     "error on the factorization of a norm\n", rel->a, rel->b);
  }
  if (ret == -1)
  {
    nrels_noprime++;
    fprintf (stderr, "Line with a=%" PRId64 " b=%" PRIu64 " failed due to a "
                     "non-prime ideal\n", rel->a, rel->b);
  }

  return NULL;
}

static void
usage(const char *argv0)
{
    fprintf (stderr, "Usage: %s [options] ", argv0);
    fprintf (stderr, "[ -filelist <fl> [-basepath <dir>] | file1 ... filen ]\n");
    fprintf (stderr, "Mandatory command line options:\n");
    fprintf (stderr, "     -poly <file>    - polynomials file\n");
    fprintf (stderr, "\nOther command line options:\n");
    fprintf (stderr, "    -abhexa          - read a and b as hexa no decimal\n");
    fprintf (stderr, "    -check_primality - by default no primality check\n");
    fprintf (stderr, "    -v               - more verbose output\n");
    exit (1);
}

int
main (int argc, char * argv[])
{
    char * argv0 = argv[0];

    param_list pl;
    param_list_init(pl);
    argv++,argc--;

    int abhexa = 0;
    param_list_configure_switch(pl, "abhexa", &abhexa);
    param_list_configure_switch(pl, "v", &verbose);
    param_list_configure_switch(pl, "check_primality", &check_primality);

#ifdef HAVE_MINGW
    _fmode = _O_BINARY;     /* Binary open for all files */
#endif

    for( ; argc ; ) {
        if (param_list_update_cmdline(pl, &argc, &argv)) { continue; }
        /* Since we accept file names freeform, we decide to never abort
         * on unrecognized options */
        break;
        // fprintf (stderr, "Unknown option: %s\n", argv[0]);
        // abort();
    }

    /* Update parameter list at least once to register argc/argv pointers. */
    param_list_update_cmdline (pl, &argc, &argv);
    /* print command-line arguments */
    param_list_print_command_line (stdout, pl);
    fflush(stdout);

    const char * polyfilename = param_list_lookup_string(pl, "poly");
    const char * filelist = param_list_lookup_string(pl, "filelist");
    const char * basepath = param_list_lookup_string(pl, "basepath");
    const char * path_antebuffer = param_list_lookup_string(pl, "path_antebuffer");
    set_antebuffer_path (argv0, path_antebuffer);

    if (param_list_warn_unused(pl))
      usage(argv0);

    if (basepath && !filelist)
    {
      fprintf(stderr, "-basepath only valid with -filelist\n");
      usage(argv0);
    }

    if ((filelist != NULL) + (argc != 0) != 1)
    {
      fprintf(stderr, "Provide either -filelist or freeform file names\n");
      usage(argv0);
    }

    cado_poly_init (cpoly);
    if (!cado_poly_read (cpoly, polyfilename))
    {
      fprintf (stderr, "Error reading polynomial file\n");
      exit (EXIT_FAILURE);
    }

    char ** files = filelist ? filelist_from_file(basepath, filelist, 0) : argv;

    nrels_read = nrels_err = nrels_noprime = 0;

    timingstats_dict_init(stats);
    filter_rels(files, (filter_rels_callback_t) &thread_check, NULL,
                EARLYPARSE_NEED_PRIMES |
                (abhexa ? EARLYPARSE_NEED_AB_HEXA : EARLYPARSE_NEED_AB_DECIMAL),
                NULL, stats);

    fprintf(stderr, "Number of read relations: %" PRid "\n", nrels_read);
    fprintf(stderr, "    among which %" PRid " has wrong norm factorization\n"
                    "            and %" PRid " contained at least one "
                    "non-primes \"ideal\"\n", nrels_err, nrels_noprime);
    if (filelist)
      filelist_clear(files);

    param_list_clear(pl);
    cado_poly_clear (cpoly);

    timingstats_dict_add_mythread(stats, "main");
    timingstats_dict_disp(stats);
    timingstats_dict_clear(stats);

    return (nrels_err + nrels_noprime == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
