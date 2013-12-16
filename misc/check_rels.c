/* check_rels: check the factorization of the norm on both alg and rat side.
   Can also, in option, check the primality of ideal and correct uncomplete
   relation or relation with non prime ideal */

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
#include "utils_with_io.h"
#include "relation.h"

#include "implicit_mpz_poly.h"

static index_t nrels_read, nrels_ok, nrels_err, nrels_completed, nrels_noprime;
static index_t nrels_toolarge;
cado_poly cpoly;
unsigned long lpb[2] = {0, 0};
int verbose = 0;
int abhexa = 0;
int check_primality = 0; /* By default no primality check */
int complete_rels = 0; /* By default, we just check the rels */

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

static inline void
print_error_line(uint64_t num, int64_t a, uint64_t b, int err, int already)
{
  if (!already)
  {
    char *str = (err) ? "Error" : "Warning";
    fprintf (stderr, "%s with relation %" PRIu64 " (a,b) = "
                     "(%" PRId64 ",%" PRIu64 "):\n", str, num, a, b);
  }
}

static inline void
factor_nonprime_ideal(earlyparsed_relation_ptr rel, weight_t i)
{
  exponent_t e = rel->primes[i].e;
  p_r_values_t p = rel->primes[i].p;
  unsigned int side = rel->primes[i].h, pr = 2;
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
      rel->primes[i] = (prime_t) {.h= side, .p= pr, .e= e * e_pr_in_p};
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

/* return 0 if everything is ok (factorization, primality, and complete)
 * return -1 if a ideal does not divide the norm (primality and completeness
 * are not checked)
 * return r = r2*2^2 + r1*2 + r0 > 0 if all ideals divides the norm but:
 *   r0 = 1 means that an ideal was not prime
 *   r1 = 1 means that the factorization of a norm was not complete
 *   r2 = 1 means that an ideal was larger that a lpb
 * If check_primality = 0, primality checks are not done and all ideals are
 * supposed prime (so r0 is always 0).
 * If complete_rels != 0, at the end, if r > 0, the relation is correct
 * (it is completed, all ideals are below the lpbs and all non-prime are
 * factored if check_primality != 0)
 */
int
process_one_relation (earlyparsed_relation_ptr rel)
{
  mpz_t norm[2];
  mpz_init (norm[0]);
  mpz_init (norm[1]);
  unsigned int fac_error = 0, need_completion = 0, nonprime = 0, toolarge = 0;
  /* If we complete relations, do not print error msg but warning. */
  int err = (complete_rels) ? 0 : 1;

  /* compute the norm on alg and rat sides */
  for(unsigned int side = 0 ; side < 2 ; side++)
  {
    cado_poly_side_ptr ps = cpoly->pols[side];
    mp_poly_homogeneous_eval_siui(norm[side], ps->f, ps->degree, rel->a, rel->b);
  }

  /* check for correctness of the factorization of the norms */
  for(weight_t i = 0; i < rel->nb ; i++)
  {
    unsigned int side = rel->primes[i].h;
    p_r_values_t p = rel->primes[i].p;
    exponent_t e = rel->primes[i].e;
    for (int j = 0; j < e; ++j)
    {
      //if (mpz_fdiv_q_ui (norm[side], norm[side], p) != 0)
      if (!mpz_divisible_ui_p (norm[side], p))
      {
        if (verbose != 0)
        {
          print_error_line (rel->num, rel->a, rel->b, 1, fac_error);
          fprintf (stderr, "    given factor %" PRpr " with exponent %u does "
                           "not divide the norm on side %u\n", p, e, side);
        }
        fac_error = 1;
      }
      else
        mpz_divexact_ui (norm[side], norm[side], p);
    }
  }

  /* With an error in the factorization of the norm, no need to continue */
  if (fac_error)
  {
    mpz_clear(norm[0]);
    mpz_clear(norm[1]);
    return -1;
  }

  /* check primality of all ideals appearing in the relations */
  if (check_primality != 0)
  {
    weight_t len = rel->nb; // no need to check the prime that can be added
    for(weight_t i = 0; i < len ; i++)
    {
      p_r_values_t p = rel->primes[i].p;
      if (!modul_isprime((const long unsigned int *)&p))
      {
        if (verbose != 0)
        {
          unsigned int side = (unsigned int) rel->primes[i].h;
          print_error_line (rel->num, rel->a, rel->b, err, nonprime);
          fprintf (stderr, "    given factor %" PRpr " is not prime on side "
                           "%u\n", p, side);
        }
        nonprime = 1;
        if (complete_rels) // if complete_rels = 1, we factor it
          factor_nonprime_ideal(rel, i);
      }
    }
  }

  /* complete relation if it is asked and necessary */
  if (complete_rels)
  {
    for (unsigned long p = 2; mpz_cmp_ui (norm[0], 1) != 0 ||
                              mpz_cmp_ui (norm[1], 1) != 0 ; p = getprime(p))
    {
      for(unsigned int side = 0 ; side < 2 ; side++)
      {
        exponent_t e = 0;
        while (mpz_divisible_ui_p (norm[side], p))
        {
          e++;
          mpz_divexact_ui (norm[side], norm[side], p);
        }
        if (e != 0)
        {
          rel_add_prime (rel, side, p, e);
          if (verbose != 0)
          {
            print_error_line (rel->num, rel->a, rel->b, err,
                                                   need_completion | nonprime);
            gmp_fprintf (stderr, "    factorization of the norm on side %u is "
                                 "not complete, need to add %" PRpr "^%u\n",
                                 side, p, e);
          }
          need_completion = 1;
        }
      }
    }

    getprime(0);
  }
  else
  {
    for(unsigned int side = 0 ; side < 2 ; side++)
    {
      if (mpz_cmp_ui (norm[side], 1) != 0)
      {
        if (verbose != 0)
        {
          print_error_line (rel->num, rel->a, rel->b, err,
                                                 need_completion | nonprime);
          gmp_fprintf (stderr, "    factorization of the norm on side %u is not "
                               "complete (%Zu is missing)\n", side, norm[side]);
        }
        need_completion = 1;
      }
    }
  }

  /* check that ideals appearing in the relations are below the lpb */
  if (lpb[0] != 0 || lpb[1] != 0)
  {
    for(weight_t i = 0; i < rel->nb ; i++)
    {
      p_r_values_t p = rel->primes[i].p;
      unsigned int side = rel->primes[i].h;
      if (lpb[side] != 0 && p > lpb[side])
      {
        if (verbose != 0)
        {
          print_error_line (rel->num, rel->a, rel->b, err,
                            nonprime | need_completion | toolarge);
          fprintf (stderr, "    given factor %" PRpr " is greater than lpb "
                           "(=%lu) on side %u\n", p, lpb[side], side);
        }
        toolarge = 1;
      }
    }
  }


  mpz_clear(norm[0]);
  mpz_clear(norm[1]);
  return nonprime + 2*need_completion + 4*toolarge;
}

static inline void
print_relation (FILE *outfile, earlyparsed_relation_ptr rel)
{
  char buf[1 << 12], *p, *op;
  size_t t;
  unsigned int i, j;

  if (!abhexa)
  {
    p = d64toa10(buf, rel->a);
    *p++ = ',';
    p = u64toa10(p, rel->b);
    *p++ = ':';
  }
  else
  {
    p = d64toa16(buf, rel->a);
    *p++ = ',';
    p = u64toa16(p, rel->b);
    *p++ = ':';
  }

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

/* Callback function called by filter_rels */

static void *
thread_callback (void * context_data, earlyparsed_relation_ptr rel)
{
  FILE *outfile = (FILE *) context_data;

  nrels_read++;
  int ret = process_one_relation (rel);

  if (complete_rels)
  {
    /* if complete_rels != 0, the only two cases where a rel is wrong are
       when ret=-1 or ret >= 4 */
    if (ret == -1)
      nrels_err++;
    else if (ret >= 4)
    {
      nrels_err++;
      nrels_toolarge++;
    }
    else
    {
      nrels_ok++;
      if (ret % 2) // rel contained a non prime ideal
        nrels_noprime++;
      if (ret >= 2) // factorization of a norm of rel was no complete
        nrels_completed++;

      print_relation (outfile, rel);
    }
  }
  else
  {
    /* if complete_rels = 0, the only case where a rel is ok is ret = 0 */
    if (ret == 0)
      nrels_ok++;
    else
    {
      nrels_err++;
      if (ret > 0)
      {
        if (ret % 2) // rel contain a non prime ideal
          nrels_noprime++;
        ret /= 2;
        if (ret % 2) // factorization of a norm of rel is no complete
          nrels_completed++;
        ret /= 2;
        if (ret % 2) // an ideal was larger than a lpb
          nrels_toolarge++;
      }
    }
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
    fprintf (stderr, "    -abhexa          - read and write a and b as hexa "
                                             "(instead of decimal)\n");
    fprintf (stderr, "    -check_primality - check primality of ideal "
                                             "(by default no checking)\n");
    fprintf (stderr, "    -complete <file> - write rels in file. If possible "
                                             "incorrect rels are corrected\n");
    fprintf (stderr, "    -lpb0 <l>        - chech that ideals on side 0 are "
                                             "below 2^l\n");
    fprintf (stderr, "    -lpb1 <l>        - chech that ideals on side 1 are "
                                             "below 2^l\n");
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
    }

    /* Update parameter list at least once to register argc/argv pointers. */
    param_list_update_cmdline (pl, &argc, &argv);
    /* print command-line arguments */
    param_list_print_command_line (stdout, pl);
    fflush(stdout);

    param_list_parse_ulong(pl, "lpb0", &lpb[0]);
    param_list_parse_ulong(pl, "lpb1", &lpb[1]);
    lpb[0] = (lpb[0] == 0) ? 0 : 1UL << lpb[0];
    lpb[1] = (lpb[1] == 0) ? 0 : 1UL << lpb[1];

    const char * polyfilename = param_list_lookup_string(pl, "poly");
    const char *outfilename = param_list_lookup_string(pl, "complete");
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

    if (outfilename)
      complete_rels = 1;

    if ((filelist != NULL) + (argc != 0) != 1)
    {
      fprintf(stderr, "Provide either -filelist or freeform file names\n");
      usage(argv0);
    }

    if (!polyfilename)
      usage(argv0);

    cado_poly_init (cpoly);
    if (!cado_poly_read (cpoly, polyfilename))
    {
      fprintf (stderr, "Error reading polynomial file\n");
      exit (EXIT_FAILURE);
    }

    char ** files = filelist ? filelist_from_file(basepath, filelist, 0) : argv;

    nrels_read = nrels_err = nrels_ok = nrels_completed = nrels_noprime = 0;
    nrels_toolarge = 0;
    if (complete_rels)
    {
      outfile = fopen_maybe_compressed (outfilename, "w");
      printf("# Correct relations will be written in %s\n", outfilename);
    }

    printf ("# Verbose output: %s\n", verbose ? "yes" : "no");
    printf ("# Correct wrong relations if possible: %s\n",
                                              complete_rels ? "yes" : "no");
    printf ("# Check primality of ideal: %s\n", check_primality ? "yes" : "no");

    timingstats_dict_init(stats);
    filter_rels(files, (filter_rels_callback_t) &thread_callback, (void*)outfile,
                EARLYPARSE_NEED_PRIMES |
                (abhexa ? EARLYPARSE_NEED_AB_HEXA : EARLYPARSE_NEED_AB_DECIMAL),
                NULL, stats);

    printf("Number of read relations: %" PRid "\n", nrels_read);
    if (complete_rels)
    {
      printf("Number of deleted relations: %" PRid "\n", nrels_err);
      printf("    among which %" PRid " contained an ideal larger than a lpb\n",
             nrels_toolarge);
      printf("Number of keeped relations: %" PRid "\n", nrels_ok);
      printf("    among which %" PRid " were completed\n"
             "            and %" PRid " contained at least one "
             "non-primes ideal\n", nrels_completed, nrels_noprime);
    }
    else
    {
      printf("Number of correct relations: %" PRid "\n", nrels_ok);
      printf("Number of wrong relations: %" PRid "\n", nrels_err);
      printf("    among which %" PRid " were not complete\n"
             "            and %" PRid " contained an ideal larger than a lpb\n"
             "            and %" PRid " contained at least one "
             "non-primes ideal\n", nrels_completed, nrels_toolarge,
             nrels_noprime);
    }

    if (filelist)
      filelist_clear(files);

    if (complete_rels)
      fclose_maybe_compressed (outfile, outfilename);

    param_list_clear(pl);
    cado_poly_clear (cpoly);

    timingstats_dict_add_mythread(stats, "main");
    timingstats_dict_disp(stats);
    timingstats_dict_clear(stats);

    return (nrels_err) ? EXIT_FAILURE : EXIT_SUCCESS ;
}
