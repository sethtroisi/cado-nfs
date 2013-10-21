/* complete_rels: same as check_rels, but completes small factors if omitted
 * and can check the primality of the factors in the norm factorization.     */

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

static index_t nrels_read, nrels_ok, nrels_err, nrels_completed, nrels_noprime;
cado_poly cpoly;
int verbose = 0;
int check_primality = 0; /* By default no primality check */

/* used for counting time in different processes */
timingstats_dict_t stats;

void
rel_add_prime (earlyparsed_relation_ptr rel, unsigned int side, p_r_values_t p,
               exponent_t e)
{
  for(weight_t i = 0; i < rel->nb ; i++)
  {
    if (rel->primes[i].p == p && rel->primes[i].h == side)
    {
      rel->primes[i].e += e;
      return;
    }
  }

  if (rel->nb == rel->nb_alloc)
    realloc_buffer_primes_c (rel);
  rel->primes[rel->nb] = (prime_t) {.h = side, .p = p, .e = e};
  rel->nb++;
}

/* return 0 if an error exists in the factorization of a norm (rat or alg)
 * return r > 0 if the factorizations were correct. The possible value of r are:
 *   1 if the factorizations were complete and contained only primes
 *   2 if the factorizations were uncomplete but contained only primes
 *   3 if the factorizations were complete but contained only a non-prime "ideal"
 *   4 if the factorizations were uncomplete and contained at least a non-prime
       "ideal"
 * At the end, if r > 0, the relation is always correct, complete and contains
 * only primes
 */
int
complete_relation (earlyparsed_relation_ptr rel)
{
  mpz_t norm[2];
  mpz_init (norm[0]);
  mpz_init (norm[1]);
  unsigned int completed = 0, nonprime = 0;

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
        mpz_clear(norm[0]);
        mpz_clear(norm[1]);
        return 0;
      }
    }
  }

  // check primality
  if (check_primality != 0)
  {
    weight_t len = rel->nb; // no need to check the prime that can be added
    for(weight_t i = 0; i < len ; i++)
    {
      p_r_values_t p = rel->primes[i].p;
      exponent_t e = rel->primes[i].e;
      if (!modul_isprime((const long unsigned int *)&p))
      {
        unsigned int side = rel->primes[i].h;
        nonprime = 1;
        if (verbose != 0)
        {
          fprintf (stderr, "Given factor %" PRpr " is not prime on side %u for "
                           "(%" PRId64 ", %" PRIu64 ")\n", p, side, rel->a,
                           rel->b);
        }

        unsigned long pr = 2;
        do
        {
          exponent_t e_pr_in_p = 0;
          while (p % pr == 0)
          {
            p = p / pr;
            e_pr_in_p++;
          }

          if (p == 1)
          {
            rel->primes[i] = (prime_t) {.h = side, .p = pr, .e = e * e_pr_in_p};
            break;
          }
          else if (e_pr_in_p != 0)
            rel_add_prime (rel, side, pr, e * e_pr_in_p);

          pr = getprime(pr);
        } while (!modul_isprime((const long unsigned int *)&p));
        getprime(0);

        if (p != 1) //means remaining p is prime
          rel->primes[i] = (prime_t) {.h = side, .p = p, .e = e};
      }
    }
  }


  // complete relation if necessary.
  for (unsigned long p = 2; mpz_cmp_ui (norm[0], 1) != 0 ||
                            mpz_cmp_ui (norm[1], 1) != 0 ; p = getprime(p))
  {
    for(unsigned int side = 0 ; side < 2 ; side++)
    {
      exponent_t e = 0;
      while (mpz_divisible_ui_p (norm[side], p))
      {
        completed = 1;
        e++;
        mpz_divexact_ui (norm[side], norm[side], p);
      }
      if (e != 0)
        rel_add_prime (rel, side, p, e);
    }
  }

  getprime(0);
  mpz_clear(norm[0]);
  mpz_clear(norm[1]);
  return 1 + completed + 2*nonprime;
}


/* Callback function called by filter_rels */

static void *
thread_complete (void * context_data, earlyparsed_relation_ptr rel)
{
  FILE *outfile = (FILE *) context_data;

  nrels_read++;
  int ret = complete_relation (rel);
  if (ret == 0)
  {
    nrels_err++;
  }
  else
  {
    nrels_ok++;
    if (ret % 2 == 0)
      nrels_completed++;
    if (ret >= 3)
      nrels_noprime++;

    //print relation
    char buf[1 << 12], *p, *op;
    size_t t;
    unsigned int i, j;

    p = d64toa16(buf, rel->a);
    *p++ = ',';
    p = u64toa16(p, rel->b);
    *p++ = ':';

    for(unsigned int side = 0 ; side < 2 ; side++)
    {
      *(--p) = ':';
      p++;
      for (i = 0; i < rel->nb; i++)
      {
        if (rel->primes[i].h == side && rel->primes[i].e > 0)
        {
          op = p;
          p = u64toa16(p, (uint64_t) rel->primes[i].p);
          *p++ = ',';
          t = p - op;
          for (j = (unsigned int) ((rel->primes[i].e) - 1); j--; p += t)
            memcpy(p, op, t);
        }
      }
    }
    *(--p) = '\n';
    p[1] = 0;
    fputs(buf, outfile);
  }

  return NULL;
}

static void
usage(const char *argv0)
{
    fprintf (stderr, "Usage: %s [options] ", argv0);
    fprintf (stderr, "[ -filelist <fl> [-basepath <dir>] | file1 ... filen ]\n");
    fprintf (stderr, "Mandatory command line options:\n");
    fprintf (stderr, "     -out <file>     - output file\n");
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

    FILE *outfile = NULL;
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
    const char *outfilename = param_list_lookup_string(pl, "out");
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

    if (!outfilename)
        usage(argv0);

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

    nrels_read = nrels_err = nrels_ok = nrels_completed = nrels_noprime = 0;
    outfile =  fopen_maybe_compressed (outfilename, "w");

    timingstats_dict_init(stats);
    filter_rels(files, (filter_rels_callback_t) &thread_complete, (void*)outfile,
                EARLYPARSE_NEED_PRIMES |
                (abhexa ? EARLYPARSE_NEED_AB_HEXA : EARLYPARSE_NEED_AB_DECIMAL),
                NULL, stats);

    fprintf(stderr, "Number of read relations: %" PRid "\n", nrels_read);
    fprintf(stderr, "Number of deleted relations: %" PRid "\n", nrels_err);
    fprintf(stderr, "Number of keeped relations: %" PRid "\n", nrels_ok);
    fprintf(stderr, "    among which %" PRid " were completed\n"
                    "            and %" PRid " contained at least one "
                    "non-primes \"ideal\"\n", nrels_completed, nrels_noprime);
    if (filelist)
      filelist_clear(files);

    fclose_maybe_compressed (outfile, outfilename);

    param_list_clear(pl);
    cado_poly_clear (cpoly);

    timingstats_dict_add_mythread(stats, "main");
    timingstats_dict_disp(stats);
    timingstats_dict_clear(stats);

    return 0;
}
