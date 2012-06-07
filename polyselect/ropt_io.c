/**
 * @file ropt_io.c
 * Ease input and output effort for ropt_main.c.
 */


#include "ropt_io.h"
#include "ropt.h"

void usage (char **argv);


/**
 * Get L1 cache size in the beginning.
 */
void
ropt_L1_cachesize ()
{
  /* L1 cache size */
  unsigned int ret = cachesize_cpuid (0);
  if ((2048 <= ret) && (ret <= (1U << 20))) {
    L1_cachesize = ret / 2;
    size_tune_sievearray = L1_cachesize / 2;
  }
  else {
    ret = cachesize_guess (0);
    if ((2048 <= ret) && (ret <= (1 << 20))) {
      L1_cachesize = ret / 2;
      size_tune_sievearray = L1_cachesize / 2;
    }
  }
}


/**
 * Ropt on a single polynomials, called by ropt_on_*() functions.
 */
static void
ropt_common ( ropt_poly_t poly,
              ropt_param_t param )
{
  /* bestpoly, info */
  ropt_info_t info;
  ropt_info_init (info);

  ropt_bestpoly_t bestpoly;
  ropt_poly_setup (poly);
  ropt_bestpoly_init (bestpoly, poly->d);
  ropt_bestpoly_setup (bestpoly, poly->f, poly->g, poly->d);

  /* print f, g */
  print_poly_fg (poly->f, poly->g, poly->d, poly->n, 1);

  /* call ropt */
  ropt (poly, bestpoly, param, info);

  fprintf (stderr, "\n# Info: Best E is:\n");
  print_poly_fg (bestpoly->f, bestpoly->g, poly->d, poly->n, 1);

  ropt_info_free (info);
  ropt_bestpoly_free (bestpoly, poly->d);
}


/**
 * Read a line.
 */
static inline int
ropt_on_stdin_readpoly_readline ( char *s )
{
  int c;

  while ((c = getchar ()) != '\n' && c != EOF)
    *s++ = c;
  *s = '\0';
  fprintf (stderr, "%s", s);
  return c;
}


/**
 * Read a polynomial from stdin.
 */
static void
ropt_on_stdin_readpoly ( mpz_t N,
                         mpz_t *f,
                         mpz_t *g,
                         mpz_t M )
{
  fflush(stdin);
  int i, ret, d;

  /* input buffer */
  char s[MAX_LINE_LENGTH];

  while (feof (stdin) == 0) {
    ret = ropt_on_stdin_readpoly_readline(s);
    if (ret == EOF)
      break;
    if (strlen (s) + 1 >= MAX_LINE_LENGTH) {
      fprintf (stderr, "Error, buffer overflow\n");
      exit (1);
    }
    /* n: input number */
    if (strncmp (s, "n:", 2) == 0) {
      if (mpz_set_str (N, s + 2, 0) != 0) {
        fprintf (stderr, "Error while reading N\n`");
        exit (1);
      }
    }
    /* ci: coeff of degree i */
    else if (sscanf (s, "c%d:", &i) == 1) {
      if (i > MAX_DEGREE) {
        fprintf (stderr, "Error, too large degree %d\n", i);
        exit (1);
      }
      if (mpz_set_str (f[i], s + 3, 0) != 0) {
        fprintf (stderr, "Error while reading f[%d]", i);
        exit (1);
      }
    }
    else if (strncmp (s, "Y1:", 3) == 0) {
      if (mpz_set_str (g[1], s + 3, 0) != 0) {
        fprintf (stderr, "Error while reading Y1");
        exit (1);
      }
    }
    else if (strncmp (s, "Y0:", 3) == 0) {
      if (mpz_set_str (g[0], s + 3, 0) != 0)
      {
        fprintf (stderr, "Error while reading Y0");
        exit (1);
      }
    }
    else if (strncmp (s, "# M", 3) == 0 || strncmp (s, "m:", 2) == 0) {
      if (mpz_set_str (M, s + 2 + (s[0] == '#'), 0) != 0) {
        fprintf (stderr, "Error while reading M or m");
        exit (1);
      }
    }
  }

  for (d = MAX_DEGREE; d > 0 && mpz_cmp_ui (f[d], 0) == 0; d --);
  if (mpz_cmp_ui (M, 0) == 0) {
    mpz_t t;
    /* M = -Y0/Y1 mod N */
    mpz_invert (M, g[1], N);
    mpz_neg (M, M);
    mpz_mul (M, M, g[0]);
    mpz_mod (M, M, N);
    /* check M is also a root of the algebraic polynomial mod N */
    mpz_init_set (t, f[d]);
    for (i = d - 1; i >= 0; i --) {
      mpz_mul (t, t, M);
      mpz_add (t, t, f[i]);
      mpz_mod (t, t, N);
    }
    if (mpz_cmp_ui (t, 0) != 0) {
      fprintf (stderr, "Polynomials have no common root mod N\n");
      exit (1);
    }
    mpz_clear (t);
  }
}


/**
 * Ropt on one polynomial from stdin.
 */
void
ropt_on_stdin ( ropt_param_t param )
{

  /* poly struct */
  ropt_poly_t poly;
  ropt_poly_init (poly);

  /* read poly from stdin */
  ropt_on_stdin_readpoly (poly->n, poly->f, poly->g, poly->m);

  /* run ropt */
  fprintf (stderr, "\n# Polynomial (# 0).\n");
  ropt_common (poly, param);

  ropt_poly_free (poly);
}


/**
 * Ropt on all CADO-format polynomials in the file.
 */
void
ropt_on_cadopoly ( FILE *file,
                   ropt_param_t param )
{
  unsigned int flag = 0U, count = 0U;
  char str[MAX_LINE_LENGTH];

  /* poly struct */
  ropt_poly_t poly;
  ropt_poly_init (poly);

  /* parse each line */
  while (1) {

    if (fgets(str, MAX_LINE_LENGTH, file) == NULL)
      break;

    if ( str[0] == 'Y' ) {
      if ( str[1] == '1' ) {
        gmp_sscanf (str, "Y1: %Zd\n", poly->g[1]);
        (flag) ^= (1<<8);
      }
      else if ( str[1] == '0' ) {
        gmp_sscanf (str, "Y0: %Zd\n", poly->g[0]);
        (flag) ^= (1<<7);
      }
      else {
        fprintf (stderr, "Error in parsing line %s.\n", str);
        exit (1);
      }
    }
    else if ( str[0] == 'c') {
      if ( str[1] == '0' ) {
        gmp_sscanf (str, "c0: %Zd\n", poly->f[0]);
        (flag) ^= 1;
      }
      else if ( str[1] == '1' ) {
        gmp_sscanf (str, "c1: %Zd\n", poly->f[1]);
        (flag) ^= (1<<1);
      }
      else if ( str[1] == '2' ) {
        gmp_sscanf (str, "c2: %Zd\n", poly->f[2]);
        (flag) ^= (1<<2);
      }
      else if ( str[1] == '3' ) {
        gmp_sscanf (str, "c3: %Zd\n", poly->f[3]);
        (flag) ^= (1<<3);
      }
      else if ( str[1] == '4' ) {
        gmp_sscanf (str, "c4: %Zd\n", poly->f[4]);
        (flag) ^= (1<<4);
      }
      else if ( str[1] == '5' ) {
        gmp_sscanf (str, "c5: %Zd\n", poly->f[5]);
        (flag) ^= (1<<5);
      }
      else if ( str[1] == '6' ) {
        gmp_sscanf (str, "c6: %Zd\n", poly->f[6]);
        (flag) ^= (1<<6);
      }
      else {
        fprintf (stderr, "Error in parsing line %s.\n", str);
        exit (1);
      }
    }
    else if ( str[0] == 'n') {
      gmp_sscanf (str, "n: %Zd\n", poly->n);
      (flag) ^= (1<<9);
    }
    else
      continue;

    if ( flag == 1023 || // deg 6
         flag == 959 || // deg 5
         flag == 927 ) // deg 4
    {
      fprintf (stderr, "\n# Polynomial (# %5d).\n", count);
      ropt_common (poly, param);
      count ++;
      flag = 0U;
    }

  } // next line

  ropt_poly_free (poly);
}


/**
 * Ropt on all Msieve-format polynomials in the file.
 */
void
ropt_on_msievepoly ( FILE *file,
                     ropt_param_t param )
{
  unsigned int count = 0U;
  char str[MAX_LINE_LENGTH];

  /* poly struct */
  ropt_poly_t poly;
  ropt_poly_init (poly);

  mpz_t ad, l, m;
  mpz_init (ad);
  mpz_init (l);
  mpz_init (m);

  /* parse each line */
  while (1) {

    /* read ad, l, m */   
    if (fgets(str, MAX_LINE_LENGTH, file) == NULL)
      break;
    int ret = gmp_sscanf (str, "%Zd %Zd %Zd", ad, l, m);
    if (ret != 3)
      continue;

    /* set degree, n, ad, l, m */
    poly->d = param->d;
    mpz_set (poly->n, param->n);
    mpz_set (poly->f[poly->d], ad);
    Lemma21 (poly->f, poly->n, poly->d, l, m);
    mpz_set (poly->g[1], l);
    mpz_neg (poly->g[0], m);
    optimize (poly->f, poly->d, poly->g, 0, 1);

    fprintf (stderr, "\n# Polynomial (# %5d).\n", count);

#if SKIP_ROPT
    /* may be used for polyselect2l results (skipped polynomials) */
    print_poly_info_short (poly->f, poly->g, poly->d, poly->n);
#else
    ropt_common (poly, param);
#endif
    count ++;
  }

  /* free */
  mpz_clear (ad);
  mpz_clear (l);
  mpz_clear (m);
  ropt_poly_free (poly);
}



/**
 * parse manually input parameters to param.
 */
void
ropt_parse_param ( int argc,
                   char **argv,
                   ropt_param_t param )
{

  /* filename only */
  if (argc > 1) {
    /* stage 2 (root sieve) parameters only */
    if (strcmp (argv[1], "--s2") == 0) {

      param->stage_flag = 2;
      argv += 1;
      argc -= 1;
      while (argc >= 2 && argv[1][0] == '-') {
        if (strcmp (argv[1], "-v") == 0)
        {
          param->verbose ++;
          argv += 1;
          argc -= 1;
        }
        else if (argc >= 3 && strcmp (argv[1], "-a") == 0)
        {
          param->s2_w = atoi (argv[2]);
          argv += 2;
          argc -= 2;
        }
        else if (argc >= 3 && strcmp (argv[1], "-b") == 0)
        {
          mpz_set_str (param->s2_u, argv[2], 10);
          argv += 2;
          argc -= 2;
        }
        else if (argc >= 3 && strcmp (argv[1], "-c") == 0)
        {
          mpz_set_str (param->s2_v, argv[2], 10);
          argv += 2;
          argc -= 2;
        }
        else if (argc >= 3 && strcmp (argv[1], "-mod") == 0)
        {
          mpz_set_str (param->s2_mod, argv[2], 10);
          argv += 2;
          argc -= 2;
        }
        else if (argc >= 3 && strcmp (argv[1], "-bmax") == 0)
        {
          param->s2_Amax = atol (argv[2]);
          argv += 2;
          argc -= 2;
        }
        else if (argc >= 3 && strcmp (argv[1], "-cmax") == 0)
        {
          param->s2_Bmax = atol (argv[2]);
          argv += 2;
          argc -= 2;
        }
        else if (argc >= 3 && strcmp (argv[1], "-n") == 0)
        {
          mpz_set_str (param->n, argv[2], 10);
          argv += 2;
          argc -= 2;
        }
        else if (argc >= 3 && strcmp (argv[1], "-norm") == 0)
        {
          param->bound_lognorm = atof (argv[2]);
          argv += 2;
          argc -= 2;
        }
        else if (argc >= 3 && strcmp (argv[1], "-d") == 0)
        {
          param->d = atoi (argv[2]);
          argv += 2;
          argc -= 2;
        }
        else {
          usage (argv);
        }
      }
    }
    /* stage 1 and stage 2 parameters */
    else {
      while (argc >= 2 && argv[1][0] == '-') {
        if (strcmp (argv[1], "-v") == 0)
        {
          param->verbose ++;
          argv += 1;
          argc -= 1;
        }
        else if (argc >= 3 && strcmp (argv[1], "-amin") == 0)
        {
          param->w_left_bound = atoi (argv[2]);
          argv += 2;
          argc -= 2;
        }
        else if (argc >= 3 && strcmp (argv[1], "-amax") == 0)
        {
          param->w_length = atoi (argv[2]) - param->w_left_bound + 1;
          if (param->w_length < 0) {
            fprintf (stderr, "Error in options -amin and/or -amax.\n");
            exit(1);
          }
          argv += 2;
          argc -= 2;
        }
        else if (argc >= 3 && strcmp (argv[1], "-bmax") == 0)
        {
          param->s2_Amax = atol (argv[2]);
          argv += 2;
          argc -= 2;
        }
        else if (argc >= 3 && strcmp (argv[1], "-cmax") == 0)
        {
          param->s2_Bmax = atol (argv[2]);
          argv += 2;
          argc -= 2;
        }
        else if (argc >= 3 && strcmp (argv[1], "-n") == 0)
        {
          mpz_set_str (param->n, argv[2], 10);
          argv += 2;
          argc -= 2;
        }
        else if (argc >= 3 && strcmp (argv[1], "-norm") == 0)
        {
          param->bound_lognorm = atof (argv[2]);
          argv += 2;
          argc -= 2;
        }
        else if (argc >= 3 && strcmp (argv[1], "-d") == 0)
        {
          param->d = atoi (argv[2]);
          argv += 2;
          argc -= 2;
        }
        else if (argc >= 3 && strcmp (argv[1], "-e") == 0)
        {
          param->stage_flag = 1;
          param->s1_num_e_sl = atoi (argv[2]);
          argv += 2;
          argc -= 2;

          int j = 0;
          for (j = 0; j < param->s1_num_e_sl; j ++) {
            param->s1_e_sl[j] = (unsigned short) atoi (argv[1]);
            argv += 1;
            argc -= 1;
          }
        }
        else {
          usage (argv);
        }
      }
    }
  }
}


#if SKIP_ROPT
/**
 * Print ad, l, m and polynomial information.
 */
double
print_poly_info_short ( mpz_t *f,
                        mpz_t *g,
                        int d,
                        mpz_t N )
{
  /* print info about the polynomial */
  unsigned int nroots = 0;
  double skew, logmu, alpha, e, alpha_proj;
  int i;
  double exp_rot[] = {0, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 0};

  /* initlize cado_poly for Murphy E */
  cado_poly cpoly;
  cado_poly_init(cpoly);

  for (i = 0; i < (d + 1); i++) {
    mpz_set(cpoly->alg->f[i], f[i]);
  }
  for (i = 0; i < 2; i++) {
    mpz_set(cpoly->rat->f[i], g[i]);
  }

  /* output original poly */
  gmp_printf ("%Zd ", f[d]);
  mpz_neg (g[0], g[0]);
  for (i = 1; i >= 0; i --) {
    gmp_printf ("%Zd ", g[i]);
  }
  printf ("\n");
  mpz_neg (g[0], g[0]);

  /* compute skew, logmu, nroots */
  nroots = numberOfRealRoots (f, d, 0, 0);
  skew = L2_skewness (f, d, SKEWNESS_DEFAULT_PREC, DEFAULT_L2_METHOD);
  logmu = L2_lognorm (f, d, skew, DEFAULT_L2_METHOD);
  alpha = get_alpha (f, d, ALPHA_BOUND);
  alpha_proj = get_biased_alpha_projective (f, d, ALPHA_BOUND);

  mpz_set (cpoly->n, N);
  cpoly->alg->degree = d;
  cpoly->rat->degree = 1;
  cpoly->skew = skew;
  e = MurphyE (cpoly, BOUND_F, BOUND_G, AREA, MURPHY_K);

  printf ("# lognorm: %.2f, alpha: %.2f, (alpha_proj: %.2f) E: %.2f, nr: %u, exp_E: %1.2f, MurphyE: %1.2e\n",
          logmu,
          alpha,
          alpha_proj,
          logmu + alpha,
          nroots,
          logmu - sqrt (2.0 * exp_rot[d] * log (skew) ),
          e );

  fflush( stdout );
  cado_poly_clear (cpoly);

  return e;
}
#endif
