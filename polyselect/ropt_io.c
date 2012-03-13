/**
 * @file ropt_io.c
 * Input-output for ropt_main.c.
 */

#include "ropt_io.h"

/*
  Some information-related routines
*/
static int
readline ( char *s )
{
  int c;

  while ((c = getchar ()) != '\n' && c != EOF)
    *s++ = c;
  *s = '\0';
  fprintf (stderr, "%s", s);
  return c;
}


void
read_ggnfs ( mpz_t N,
             mpz_t *f,
             mpz_t *g,
             mpz_t M )
{
  fflush(stdin);
  int i, ret, d;
  char s[MAX_LINE_LENGTH]; /* input buffer */

  while (feof (stdin) == 0) {
    ret = readline (s);
    if (ret == EOF)
      break;
    if (strlen (s) + 1 >= MAX_LINE_LENGTH) {
      fprintf (stderr, "Error, buffer overflow\n");
      exit (1);
    }
    if (strncmp (s, "n:", 2) == 0) { /* n: input number */
      if (mpz_set_str (N, s + 2, 0) != 0) {
        fprintf (stderr, "Error while reading N\n`");
        exit (1);
      }
    }
    else if (sscanf (s, "c%d:", &i) == 1) {/* ci: coeff of degree i */
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


/*
  Call print_poly, given f, g.
*/
double
print_poly_fg ( mpz_t *f,
                mpz_t *g,
                int d,
                mpz_t N,
                int verbose )
{
  double e;
  int i;

  cado_poly cpoly;
  cado_poly_init(cpoly);
  for (i = 0; i < (d + 1); i++)
    mpz_set(cpoly->alg->f[i], f[i]);
  for (i = 0; i < 2; i++)
    mpz_set(cpoly->rat->f[i], g[i]);
  mpz_set(cpoly->n, N);
  cpoly->skew = L2_skewness (f, d, SKEWNESS_DEFAULT_PREC, DEFAULT_L2_METHOD);
  cpoly->alg->degree = d;
  cpoly->rat->degree = 1;

  if (verbose == 2) {
    e = print_cadopoly (stdout, cpoly, 1);
    fflush(stdout);
  }
  else {
    e = MurphyE (cpoly, BOUND_F, BOUND_G, AREA, MURPHY_K);
  }

  cado_poly_clear (cpoly);
  return e;
}


#if SKIP_ROOTSIEVE_M
/*
  Print ad, l, m and information as above
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
