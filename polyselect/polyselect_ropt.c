#include "cado.h"
#include <stdio.h>
#include <pthread.h>
#include "portability.h"
#include "utils.h"
#include "area.h"
#include "ropt.h"

/* thread structure for ropt */
typedef struct
{
  cado_poly_ptr poly;
  unsigned int id;
  unsigned int poly_id;
  double ropt_time;
} __ropt_thread_struct;
typedef __ropt_thread_struct ropt_thread_t[1];
typedef __ropt_thread_struct * ropt_thread_ptr;
typedef const __ropt_thread_struct * ropt_thread_srcptr;

pthread_mutex_t lock=PTHREAD_MUTEX_INITIALIZER; /* used as mutual exclusion
                                                   lock for those variables */
unsigned int nthreads = 1;
int tot_found = 0; /* total number of polynomials */
static int verbose = 0;
int ropteffort = DEFAULT_RSEFFORT; /* sieving effort, among 1-5 */
cado_poly best_poly;
double best_MurphyE = 0.0; /* Murphy's E (the larger the better) */

static inline void
ropt_wrapper (cado_poly_ptr input_poly, unsigned int poly_id, double *ropt_time)
{
  double curr_MurphyE, st;
  mpz_t t;
  cado_poly ropt_poly;
  cado_poly_init (ropt_poly);

  if (nthreads > 1)
    pthread_mutex_lock (&lock);
  printf ("\n### input polynomial %u ###\n", poly_id);
  cado_poly_fprintf_with_info (stdout, input_poly, "# ");
  if (nthreads > 1)
    pthread_mutex_unlock (&lock);

  /* If the content of the algebraic polynomial has content <> 1, then print a
     warning (this should not be frequent) and divides all coefficients of the
     polynomial by the content. */
  mpz_init (t);
  mpz_poly_content (t, input_poly->alg);
  if (mpz_cmp_ui (t, 1) != 0)
  {
    gmp_printf ("# WARNING: the content of the algebraic side of polynomial %u "
                "is not 1 (%Zd). The input polynomial will be divided by its "
                "content.\n", poly_id, t);
    mpz_poly_divexact_mpz (input_poly->alg, input_poly->alg, t);
  }
  mpz_clear (t);

  st = seconds_thread ();
  ropt_polyselect (ropt_poly, input_poly, ropteffort, verbose);
  *ropt_time += seconds_thread () - st;

  /* MurphyE */
  ropt_poly->skew = L2_skewness (ropt_poly->alg, SKEWNESS_DEFAULT_PREC);
  curr_MurphyE = MurphyE (ropt_poly, bound_f, bound_g, area, MURPHY_K);
  if (nthreads > 1)
    pthread_mutex_lock (&lock);
  if (curr_MurphyE > best_MurphyE)
  {
    best_MurphyE = curr_MurphyE;
    cado_poly_set (best_poly, ropt_poly);
  }
  printf ("\n### root-optimized polynomial %u ###\n", poly_id);
    cado_poly_fprintf_with_info_and_MurphyE (stdout, ropt_poly, curr_MurphyE,
                                             bound_f, bound_g, area, "# ");
  printf ("### Best MurphyE so far is %.2e\n", best_MurphyE);
  fflush (stdout);
  if (nthreads > 1)
    pthread_mutex_unlock (&lock);

  cado_poly_clear (ropt_poly);
}

void *
thread_ropt (void *args)
{
  ropt_thread_ptr data = (ropt_thread_ptr) args;
  ropt_wrapper (data->poly, data->poly_id, &(data->ropt_time));
  return NULL;
}


static void
declare_usage(param_list pl)
{
  char str[200];
  param_list_decl_usage(pl, "inputpolys", "root-sieve the size-optimized "
                                          "polynomials given in this file");
  snprintf (str, 200, "root-sieve effort ranging from 1 to 10 (default %d)",
            DEFAULT_RSEFFORT);
  param_list_decl_usage(pl, "ropteffort", str);
  param_list_decl_usage(pl, "t", "number of threads to use (default 1)");
  snprintf (str, 200, "sieving area (default %.2e)", AREA);
  param_list_decl_usage(pl, "area", str);
  snprintf (str, 200, "algebraic smoothness bound (default %.2e)", BOUND_F);
  param_list_decl_usage(pl, "Bf", str);
  snprintf (str, 200, "rational smoothness bound (default %.2e)", BOUND_G);
  param_list_decl_usage(pl, "Bg", str);
  param_list_decl_usage(pl, "v", "(switch) verbose mode");
  verbose_decl_usage(pl);
}

static void
usage (const char *argv, const char * missing, param_list pl)
{
  if (missing)
  {
    fprintf(stderr, "\nError: missing or invalid parameter \"-%s\"\n",
                    missing);
  }
  param_list_print_usage (pl, argv, stderr);
  param_list_clear (pl);
  exit (EXIT_FAILURE);
}

int
main (int argc, char *argv[])
{
  char **argv0 = argv;
  const char *polys_filename = NULL;
  double st0 = seconds ();
  FILE *polys_file = NULL;
  cado_poly *input_polys = NULL;
  unsigned int nb_input_polys = 0; /* number of input polynomials */
  unsigned int size_input_polys = 16; /* Size of input_polys tab. */
  double rootsieve_time = 0.0;

  cado_poly_init (best_poly);
  input_polys = (cado_poly *) malloc (size_input_polys * sizeof (cado_poly));
  ASSERT_ALWAYS (input_polys != NULL);
  for (unsigned int i = 0; i < size_input_polys; i++)
    cado_poly_init (input_polys[i]);

  /* read params */
  param_list pl;
  param_list_init (pl);

  declare_usage(pl);

  param_list_configure_switch (pl, "-v", &verbose);

  if (argc == 1)
    usage (argv0[0], NULL, pl);

  argv++, argc--;
  for ( ; argc; ) {
    if (param_list_update_cmdline (pl, &argc, &argv)) continue;
    fprintf (stderr, "Unhandled parameter %s\n", argv[0]);
    usage (argv0[0], NULL, pl);
  }

  param_list_parse_uint (pl, "t", &nthreads);
  if (param_list_parse_double (pl, "area", &area) == 0) /* no -area */
    area = AREA;
  if (param_list_parse_double (pl, "Bf", &bound_f) == 0) /* no -Bf */
    bound_f = BOUND_F;
  if (param_list_parse_double (pl, "Bg", &bound_g) == 0) /* no -Bg */
    bound_g = BOUND_G;
  /* sieving effort that passed to ropt */
  param_list_parse_int (pl, "ropteffort", &ropteffort);
  /* filename of the file with the list of polynomials to root-sieve */
  polys_filename = param_list_lookup_string (pl, "inputpolys");

  if (param_list_warn_unused(pl))
    usage (argv0[0], NULL, pl);

  /* print command line */
  verbose_interpret_parameters(pl);
  param_list_print_command_line (stdout, pl);

  /* Check that ropteffort is in [1,10] */
  if (ropteffort < 1 || ropteffort > 10)
  {
    fprintf (stderr, "Error, -ropteffort should be in [1,10]\n");
    usage (argv0[0], NULL, pl);
  }

  /* Check that nthreads is >= 1 */
  if (nthreads == 0)
  {
    fprintf (stderr, "Error, -t must be non-zero.\n");
    usage (argv0[0], NULL, pl);
  }

  if (polys_filename == NULL)
    usage (argv0[0], "inputpolys", pl);

  printf ("# Info: Will use %u thread%s\n# Info: ropteffort = %d\n", nthreads,
          (nthreads > 1) ? "s": "", ropteffort);

  /* detect L1 cache size */
  ropt_L1_cachesize ();

  printf ("# Info: L1_cachesize/2 = %u, size_tune_sievearray = %u\n",
          L1_cachesize, size_tune_sievearray);

  /* Open file containing polynomials. */
  printf ("# Reading polynomials from %s\n", polys_filename);
  polys_file = fopen(polys_filename, "r");
  if (polys_file == NULL)
  {
    perror("Could not open file");
    exit(EXIT_FAILURE);
  }

  /* Read all polynomials from file. Store then in input_polys. */
  while (cado_poly_read_next_poly_from_stream (input_polys[nb_input_polys],
                                                                    polys_file))
  {
    nb_input_polys++;
    if (nb_input_polys == size_input_polys) /* Realloc if needed */
    {
      if (verbose > 0)
        fprintf (stderr, "# Reallocating input_polys\n");
      unsigned int new_size = 2 * size_input_polys;
      input_polys = (cado_poly *) realloc (input_polys, new_size * sizeof (cado_poly));
      ASSERT_ALWAYS (input_polys != NULL);
      for (unsigned int i = size_input_polys; i < new_size; i++)
        cado_poly_init (input_polys[i]);
      size_input_polys = new_size;
    }
  }
  /* Did we stop at the end of the file or was there an error. */
  if (!feof (polys_file))
  {
    fprintf (stderr, "Error while reading file %s.\n", polys_filename);
    abort ();
  }
  fclose (polys_file);
  printf ("# %u polynomials read.\n", nb_input_polys);

  /* Main loop: do root-optimization on input_polys. */
  if (nthreads > 1) /* multi thread version */
  {
    /* Allocated memory for threads and threads_data */
    pthread_t *threads = NULL;
    ropt_thread_t *threads_data = NULL;
    threads = (pthread_t *) malloc (nthreads * sizeof (pthread_t));
    ASSERT_ALWAYS (threads != NULL);
    threads_data = (ropt_thread_t *) malloc (nthreads * sizeof (ropt_thread_t));
    ASSERT_ALWAYS (threads_data != NULL);
    for (unsigned int i = 0; i < nthreads; i++)
    {
      threads_data[i]->ropt_time = 0.0;
      threads_data[i]->id = i;
    }

    for (unsigned int i = 0; i < nb_input_polys; )
    {
      unsigned int j;
      for (j = 0; j < nthreads && i < nb_input_polys; j++, i++)
      {
        threads_data[j]->poly = input_polys[i];
        threads_data[j]->poly_id = i;
        pthread_create (&threads[j], NULL, thread_ropt,
                                                    (void *) (threads_data[j]));
      }

      /* we have created j threads, with j = nthreads usually, except at the
         end of we might have j < nthreads */
      while (j > 0)
        pthread_join (threads[--j], NULL);
    }

    for (unsigned int i = 0; i < nthreads; i++)
    {
      rootsieve_time += threads_data[i]->ropt_time;
      if (verbose > 0)
        printf ("# Stat: rootsieve on thread %u took %.2fs\n", i,
                threads_data[i]->ropt_time);
    }

    free (threads);
    free (threads_data);
  }
  else /* mono thread version */
  {
    for (unsigned int i = 0; i < nb_input_polys; i++)
      ropt_wrapper (input_polys[i], i, &rootsieve_time);
  }

  /* print total time and rootsieve time.
     These two lines gets parsed by the script. */
  printf ("\n\n# Stat: total phase took %.2fs\n", seconds () - st0);

#ifndef HAVE_RUSAGE_THREAD /* rootsieve_time is correct only if RUSAGE_THREAD
                              works or in mono-thread mode */
  if (nthreads == 1)
#endif
    printf ("# Stat: rootsieve took %.2fs\n", rootsieve_time);

  if (best_MurphyE == 0.0)
  {
    if (nb_input_polys > 0)
    {
      /* At least one poly was parsed and the best MurphyE is still 0 => Error */
      fprintf (stderr, "Error, the best MurphyE value is still 0 at the end.\n");
      abort ();
    }
    else
    {
      /* No polynomials were parsed in the input files but the whole file was
         parsed without errors (see !feof above) => Warning */
      printf ("# WARNING: No polynomials were found in the input file %s\n",
              polys_filename);
    }
  }
  else
  {
    printf ("# Best polynomial found:\n");
    cado_poly_fprintf_with_info_and_MurphyE (stdout, best_poly, best_MurphyE,
                                             bound_f, bound_g, area, NULL);
  }

  cado_poly_clear (best_poly);
  for (unsigned int i = 0; i < size_input_polys; i++)
    cado_poly_clear (input_polys[i]);
  free (input_polys);
  param_list_clear (pl);

  return EXIT_SUCCESS;
}
