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
#include "relation.h"
#include "params.h"
#include "cado_poly.h"
#include "gzip.h"
#include "las-duplicate.h"
#include "las-coordinates.h"
#include "las-norms.h"
#include "verbose.h"

static void *
dupsup (FILE *output, relation & rel, const mpz_t sq, const mpz_t rho, const int side, const int is_dupe)
{
  if (0)
    gmp_fprintf (output, "# sq = %Zd, rho = %Zd, side = %d\n", sq, rho, side);
  rel.print(output, is_dupe ? "# DUPE " : "");
  return NULL;
}

/* If the line is a special-q comment, sets sq and rho and returns 1.
   Otherwise returns 0. */
static int
read_sq_comment(mpz_t sq, mpz_t rho, int *side, const char *line)
{
  if (gmp_sscanf(line, "# Sieving algebraic q=%Zd; rho=%Zd;", sq, rho) == 2) {
    *side = ALGEBRAIC_SIDE;
    return 1;
  }
  if (gmp_sscanf(line, "# Sieving rational q=%Zd; rho=%Zd;", sq, rho) == 2) {
    *side = RATIONAL_SIDE;
    return 1;
  }
  return 0;
}

static void
read_poly(cado_poly_ptr cpoly, param_list pl)
{
    cado_poly_init(cpoly);
    const char *cpoly_filename;
    if ((cpoly_filename = param_list_lookup_string(pl, "poly")) == NULL) {
        fprintf(stderr, "Error: -poly is missing\n");
        param_list_print_usage(pl, NULL, stderr);
	cado_poly_clear(cpoly);
	param_list_clear(pl);
        exit(EXIT_FAILURE);
    }

    if (!cado_poly_read(cpoly, cpoly_filename)) {
	fprintf(stderr, "Error reading polynomial file %s\n", cpoly_filename);
	cado_poly_clear(cpoly);
	param_list_clear(pl);
	exit(EXIT_FAILURE);
    }
}

static int
parse_config(siever_config_ptr sc, param_list pl)
{
    sc->side = param_list_parse_switch(pl, "-ratq") ? RATIONAL_SIDE : ALGEBRAIC_SIDE;
    param_list_parse_double(pl, "lambda0", &(sc->sides[RATIONAL_SIDE]->lambda));
    param_list_parse_double(pl, "lambda1", &(sc->sides[ALGEBRAIC_SIDE]->lambda));
    int seen = 1;
    seen  = param_list_parse_int   (pl, "I",       &(sc->logI));
    seen &= param_list_parse_ulong (pl, "lim0",    &(sc->sides[RATIONAL_SIDE]->lim));
    seen &= param_list_parse_int   (pl, "lpb0",    &(sc->sides[RATIONAL_SIDE]->lpb));
    seen &= param_list_parse_int   (pl, "mfb0",    &(sc->sides[RATIONAL_SIDE]->mfb));
    seen &= param_list_parse_ulong (pl, "lim1",    &(sc->sides[ALGEBRAIC_SIDE]->lim));
    seen &= param_list_parse_int   (pl, "lpb1",    &(sc->sides[ALGEBRAIC_SIDE]->lpb));
    seen &= param_list_parse_int   (pl, "mfb1",    &(sc->sides[ALGEBRAIC_SIDE]->mfb));
    return seen;
}


static void declare_usage(param_list pl)
{
  param_list_decl_usage(pl, "path_antebuffer", "path to antebuffer program");
  param_list_decl_usage(pl, "poly", "polynomial file");
  param_list_decl_usage(pl, "I",    "set sieving region to 2^I");
  param_list_decl_usage(pl, "skew", "(alias S) skewness");
  param_list_decl_usage(pl, "lim0", "(alias rlim) rational factor base bound");
  param_list_decl_usage(pl, "lim1", "(alias alim) algebraic factor base bound");
  param_list_decl_usage(pl, "lpb0", "(alias lpbr) set rational large prime bound to 2^lpb0");
  param_list_decl_usage(pl, "lpb1", "(alias lpba) set algebraic large prime bound to 2^lpb1");
  param_list_decl_usage(pl, "mfb0", "(alias mfbr) set rational cofactor bound 2^mfb0");
  param_list_decl_usage(pl, "mfb1", "(alias mfba) set algebraic cofactor bound 2^mfb1");
  param_list_decl_usage(pl, "lambda0", "(alias rlambda) rational lambda value");
  param_list_decl_usage(pl, "lambda1", "(alias alambda) algebraic lambda value");
  param_list_decl_usage(pl, "powlim0", "(alias rpowlim) limit on powers on rat side");
  param_list_decl_usage(pl, "powlim1", "(alias apowlim) limit on powers on alg side");
  param_list_decl_usage(pl, "ratq", "(switch) use rational special-q");
  param_list_decl_usage(pl, "mt",   "number of threads to use");
  verbose_decl_usage(pl);
}

static void
usage (param_list pl, char *argv0)
{
    param_list_print_usage(pl, argv0, stderr);
    exit(EXIT_FAILURE);
}


int
main (int argc, char * argv[])
{
    char * argv0 = argv[0];
    facul_strategy_t *strategy[2];
    siever_config conf;
    int nb_threads = 1;

    param_list pl;
    param_list_init(pl);
    declare_usage(pl);
    argv++,argc--;

#ifdef HAVE_MINGW
    _fmode = _O_BINARY;     /* Binary open for all files */
#endif

    if (argc == 0)
      usage (pl, argv0);

    for( ; argc ; ) {
        if (param_list_update_cmdline(pl, &argc, &argv)) { continue; }
        /* Since we accept file names freeform, we decide to never abort
         * on unrecognized options */
        break;
        // fprintf (stderr, "Unknown option: %s\n", argv[0]);
        // abort();
    }
    /* print command-line arguments */
    verbose_interpret_parameters(pl);
    param_list_print_command_line (stdout, pl);
    fflush(stdout);

    if (argc == 0) {
      fprintf(stderr, "Error, provide freeform file names\n");
      usage(pl, argv0);
    }

    cado_poly cpoly;
    read_poly(cpoly, pl);
    int ok = parse_config(conf, pl);
    ok = ok && param_list_parse_double(pl, "skew",    &(cpoly->skew));
    if (!ok) {
        fprintf(stderr, "Error: mandatory parameter missing.\n");
	param_list_clear(pl);
        exit(EXIT_FAILURE);
    }
    param_list_parse_int(pl, "mt", &nb_threads);

    if (param_list_warn_unused(pl))
    {
      fprintf(stderr, "Error, unused parameters are given\n");
      usage(pl, argv0);
    }

    tune_las_memset();

    for (int side = 0; side < 2; side++)
      strategy[side] = facul_make_strategy(conf->sides[side]->lim,
                                           conf->sides[side]->lpb, 0, 0);

    mpz_t sq, rho;
    mpz_init(sq);
    mpz_init(rho);

    for (int argi = 0; argi < argc; argi++) {
      FILE *f = fopen_maybe_compressed(argv[argi], "rb");
      sieve_info_ptr si = NULL;
      while (!feof(f)) {
        char line[1024];
        int side = 0;
        if (fgets(line, sizeof(line), f) == NULL)
          break;
        if (read_sq_comment(sq, rho, &side, line)) {
          unsigned long limits[2] = {conf->sides[0]->lim, conf->sides[1]->lim};
          if (si != NULL)
            clear_sieve_info(si);
          si = fill_in_sieve_info(sq, rho, side, 1U << conf->logI, 1U << (conf->logI - 1),
                                  limits, strategy, cpoly, conf);
        } else {
            relation rel;
            if (rel.parse(line)) {
                ASSERT_ALWAYS(si != NULL);
                int is_dupe = relation_is_duplicate(rel, nb_threads, si);
                dupsup(stdout, rel, sq, rho, side, is_dupe);
            }
        }
      }
      fclose_maybe_compressed(f, argv[argi]);
      clear_sieve_info(si);
    }
    
    for (int side = 0; side < 2; side++)
      facul_clear_strategy(strategy[side]);
    cado_poly_clear(cpoly);
    mpz_clear(sq);
    mpz_clear(rho);
    param_list_clear(pl);

    return 0;
}
