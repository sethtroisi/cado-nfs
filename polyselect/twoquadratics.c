#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <assert.h>
#include <math.h>
#include "portability.h"
#include "utils.h"
#include "auxiliary.h"

typedef struct {
  unsigned int dim; /* dimension of the vector */
  mpz_t *c;         /* its coordinates */
} mpz_vector_struct_t;

typedef mpz_vector_struct_t mpz_vector_t[1];
typedef mpz_vector_struct_t * mpz_vector_ptr;

/* Management of the structure, set and print coefficients. */
void mpz_vector_init(mpz_vector_t v, unsigned int d)
{
  ASSERT_ALWAYS (d > 0);
  v->dim = d;
  v->c = (mpz_t *) malloc (d * sizeof(mpz_t));
  ASSERT_ALWAYS (v->c != NULL);
  for (unsigned int i = 0; i < d; i++)
    mpz_init (v->c[i]);
}

void mpz_vector_clear(mpz_vector_t v)
{
  for (unsigned int i = 0; i < v->dim; i++)
    mpz_clear (v->c[i]);
  free (v->c);
}

void mpz_vector_swap (mpz_vector_t v1, mpz_vector_t v2)
{
  ASSERT_ALWAYS (v1->dim == v2->dim);
  mpz_t *tmp = v1->c;
  v1->c = v2->c;
  v2->c = tmp;
}

void mpz_vector_set (mpz_vector_t v, mpz_vector_t s)
{
  ASSERT_ALWAYS (v->dim == s->dim);
  for (unsigned int i = 0; i < v->dim; i++)
    mpz_set (v->c[i], s->c[i]);
}

void mpz_vector_get_mpz_poly (mpz_poly_t p, mpz_vector_t v)
{
  for (unsigned int i = 0; i < v->dim; i++)
    mpz_poly_setcoeff (p, i, v->c[i]);
}

void mpz_vector_setcoordinate (mpz_vector_t v, unsigned int i, mpz_t z)
{
  ASSERT_ALWAYS (i < v->dim);
  mpz_set(v->c[i], z);
}

void mpz_vector_setcoordinate_ui (mpz_vector_t v, unsigned int i,
                                   unsigned int z)
{
  ASSERT_ALWAYS (i < v->dim);
  mpz_set_ui (v->c[i], z);
}

int
mpz_vector_is_coordinate_zero (mpz_vector_t v, unsigned int i)
{
  ASSERT_ALWAYS (i < v->dim);
  return (mpz_sgn (v->c[i]) == 0);
}

void mpz_vector_dot_product (mpz_t res, mpz_vector_t a, mpz_vector_t b)
{
  ASSERT_ALWAYS (a->dim == b->dim);
  mpz_mul(res, a->c[0], b->c[0]);
  for (unsigned int i = 1; i < a->dim; i++)
    mpz_addmul(res, a->c[i], b->c[i]);
}

void mpz_vector_norm (mpz_t res, mpz_vector_t a)
{
  mpz_vector_dot_product (res, a, a);
}

void mpz_vector_skew_dot_product (mpz_t res, mpz_vector_t a, mpz_vector_t b,
                                  mpz_t skewness)
{
  ASSERT_ALWAYS (a->dim == b->dim);
  mpz_t tmp, s, s2;
  mpz_init (tmp);
  mpz_init (s2);
  mpz_init_set_ui (s, 1);
  mpz_mul (s2, skewness, skewness);

  mpz_mul(res, a->c[0], b->c[0]);
  for (unsigned int i = 1; i < a->dim; i++)
  {
    mpz_mul(tmp, a->c[i], b->c[i]);
    mpz_mul(s, s, s2);
    mpz_addmul(res, tmp, s);
  }

  mpz_clear(tmp);
  mpz_clear(s2);
  mpz_clear(s);
}

void mpz_vector_skew_norm (mpz_t res, mpz_vector_t a, mpz_t skewness)
{
  mpz_vector_skew_dot_product (res, a, a, skewness);
}


/* compute r <- r - q*v */
void
mpz_vector_submul (mpz_vector_t r, mpz_t q, mpz_vector_t v)
{
  ASSERT_ALWAYS (r->dim == v->dim);
  for (unsigned int i = 0; i < r->dim; i++)
    mpz_submul(r->c[i], q, v->c[i]);
}



typedef struct {
  cado_poly poly;
  mpz_t m;
  mpz_t p;  // common root is m/p mod N
  mpz_t skew;
  double note;
} cado_poly_extended_s;

typedef cado_poly_extended_s cado_poly_extended[1];
typedef cado_poly_extended_s *cado_poly_extended_ptr;

void
cado_poly_extended_init (cado_poly_extended poly)
{
  cado_poly_init (poly->poly);
  mpz_init (poly->m);
  mpz_init (poly->p);
  mpz_init (poly->skew);
  poly->note = 0.0;
}

void
cado_poly_extended_clear (cado_poly_extended poly)
{
  cado_poly_clear (poly->poly);
  mpz_clear (poly->m);
  mpz_clear (poly->p);
  mpz_clear (poly->skew);
  poly->note = 0.0;
}

void
cado_poly_extended_set (cado_poly_extended poly, mpz_poly_t f, mpz_poly_t g,
                        mpz_t N, mpz_t m, mpz_t p, mpz_t skew, double note)
{
  mpz_poly_set (poly->poly->pols[0], f->coeff, f->deg);
  mpz_poly_set (poly->poly->pols[1], g->coeff, g->deg);
  mpz_set (poly->poly->n, N);
  poly->poly->skew = mpz_get_d (skew);
  mpz_invert (poly->poly->m, p, N);
  mpz_mul (poly->poly->m, poly->poly->m, m);
  mpz_mod (poly->poly->m, poly->poly->m, N);

  mpz_set (poly->m, m);
  mpz_set (poly->p, p);
  mpz_set (poly->skew, skew);
  poly->note = note;
}

void
cado_poly_extended_set_note (cado_poly_extended poly, double note)
{
  poly->note = note;
}

double
cado_poly_extended_get_note (cado_poly_extended poly)
{
  return poly->note;
}

void
cado_poly_extended_print (FILE *out, cado_poly_extended poly, char *pre)
{
  mpz_t tmp;
  mpz_init (tmp);
  double s = poly->poly->skew;

  mpz_poly_ptr f0 = poly->poly->pols[0];
  mpz_poly_ptr f1 = poly->poly->pols[1];

  gmp_fprintf (out, "%sN = %Zd\n", pre, poly->poly->n);

  fprintf(out, "%sf0 = ", pre);
  mpz_poly_fprintf (out, f0);
  fprintf(out, "%sf1 = ", pre);
  mpz_poly_fprintf (out, f1);

  fprintf (out, "%snote = %e\n", pre, poly->note);
  fprintf (out, "%salpha_f0 = %.2f\n", pre, get_alpha (f0, ALPHA_BOUND));
  fprintf (out, "%salpha_f1 = %.2f\n", pre, get_alpha (f1, ALPHA_BOUND));
  gmp_fprintf (out, "%sskewness = %Zd\n", pre, poly->skew);
  fprintf (out, "%sL2_skew_norm_f0 = %.2f\n", pre, L2_lognorm (f0, s));
  fprintf (out, "%sL2_skew_norm_f1 = %.2f\n", pre, L2_lognorm (f1, s));
  gmp_fprintf (out, "%sp = %Zd\n", pre, poly->p);
  gmp_fprintf (out, "%sm = %Zd\n", pre, poly->m);
  mpz_mod (tmp, poly->m, poly->p);
  gmp_fprintf (out, "%sr = %Zd\n", pre, tmp);
  gmp_fprintf (out, "%scommon_root = %Zd\n", pre, poly->poly->m);

  mpz_clear (tmp);
}



/* return 0 if gcd(a, b) != 1, return non-zero if a and b are coprime */
int
mpz_coprime_p (mpz_t a, mpz_t b)
{
  mpz_t g;
  mpz_init(g);
  mpz_gcd (g, a, b);
  int ret = mpz_cmp_ui (g, 1);
  mpz_clear(g);
  return (ret == 0) ? 1 : 0;
}

/* Set m to the nearest integer to m0, such that
   m = r mod P, i.e. find the integer in [m-p/2,m+p/2] which is congruent to
   r modulo P.
 */
void
compute_m (mpz_t m, mpz_t m0, mpz_t r, mpz_t P)
{
  mpz_t t;
  mpz_init (t);
  mpz_tdiv_q_2exp (m, P, 1); // m = floor(p/2)
  mpz_add (m, m ,m0);        // m = m0 + floor(p/2)
  mpz_sub (t, m ,r);         // t = m0 + floor(p/2) - r
  mpz_mod (t, t, P);         // t = (m0 + floor(p/2) -r) mod P
  mpz_sub (m, m, t);         // m = m0 + floor(p/2) - t
                             // By definition, (m0 + floor(p/2) - t) = r mod P
  mpz_clear (t);
}

/* Compute maximun skewness, which in this case is (1/6)^(1/4) * m^(1/2).
   We approximate (1/6)^(1/4) by 3124081/4889451 (the difference is less than
   10^-14).
 */
void
compute_default_max_skew (mpz_t skew, mpz_t m)
{
  mpz_sqrt (skew, m);
  mpz_mul_ui (skew, skew, 3124081);
  mpz_tdiv_q_ui (skew, skew, 4889451);
}

/* Returns the square root of N modulo P, using Tonelli-Shanks' algorithm.
   Assume P is prime.
   Returns 0 in case of error (e.g. N is not a square), else non-zero.
*/
int
mpz_msqrt (mpz_t r, mpz_t N, mpz_t P)
{
  if (mpz_cmp_ui (P, 2) == 0)
  {
    mpz_mod (r, N, P);
    return mpz_sgn(r); // Return 0 if r == 0
  }
  else if (mpz_legendre (N, P) != 1)
    return 0;
  else
  {
    mpz_t Q, n, z, c, R, t, t0, b;

    mpz_init(n);
    mpz_init_set_ui(z, 2);
    mpz_init(c);
    mpz_init(R);
    mpz_init(t);
    mpz_init(t0);
    mpz_init(b);
    mpz_init(Q);

    mpz_mod(n, N, P);

    mpz_sub_ui (Q, P, 1); // Q = P-1
    int M, S;
    for (S = 0; mpz_even_p(Q); S++)
      mpz_tdiv_q_2exp (Q, Q, 1);
    M = S;

    while (mpz_legendre(z, P) != -1)
        mpz_add_ui(z, z, 1);

    // We want c = z^Q mod P , t = n^Q mod P and R = n^((Q+1)/2) mod P
    mpz_powm (c, z, Q, P); // c = z^Q mod P
    mpz_sub_ui (Q, Q, 1);
    mpz_tdiv_q_2exp (Q, Q, 1);
    mpz_powm (R, n, Q, P); // R = n ^ ((Q-1)/2) mod P
    mpz_mul (t, R, R);
    mpz_mod (t, t, P);       // t = n^(Q-1) mod P
    mpz_mul (R, R, n);
    mpz_mod (R, R, P);       // R = n ^ ((Q+1)/2) mod P
    mpz_mul (t, t, n);
    mpz_mod (t, t, P);       // t = n^Q mod P

    while (mpz_cmp_ui (t, 1) != 0)
    {
      int i = 0;
      mpz_set(t0, t);
      while(mpz_cmp_ui(t, 1) != 0)
      {
        mpz_pow_ui(t, t, 2);
        mpz_mod(t, t, P);
        i++;
      }
      mpz_set(b, c);
      int j;
      for (j = 0; j < M - i - 1; j++)
      {
        mpz_pow_ui(b, b, 2);
        mpz_mod(b, b, P);
      }
      mpz_mul(R, R, b);
      mpz_mod(R, R, P);
      mpz_mul(c, b, b);
      mpz_mod(c, c, P);
      mpz_mul(t, t0, c);
      mpz_mod(t, t, P);
      M = i;
    }

  mpz_set (r, R);

  mpz_clear(Q);
  mpz_clear(n);
  mpz_clear(z);
  mpz_clear(c);
  mpz_clear(t);
  mpz_clear(t0);
  mpz_clear(b);
  mpz_clear(R);
  return 1;
  }
}

/* Lagrange Algorithm.
   Reduce lattice of rank 2.
   The obtained basis is always the two shortest vector.
   See, Nguyen And Vallee, "The LLL algorithm: Survey and applications" p.41
   */
void LagrangeAlgo (mpz_vector_t a, mpz_vector_t b, mpz_t skewness)
{
  mpz_t norm_a, norm_b, q, r, tmp;
  mpz_init (norm_a);
  mpz_init (norm_b);
  mpz_init (q);
  mpz_init (r);
  mpz_init (tmp);

  int skew = (mpz_cmp_ui (skewness, 1) > 0);

  if (skew)
  {
    mpz_vector_skew_norm (norm_a, a, skewness);
    mpz_vector_skew_norm (norm_b, b, skewness);
  }
  else
  {
    mpz_vector_norm (norm_a, a);
    mpz_vector_norm (norm_b, b);
  }

  if (mpz_cmp (norm_a, norm_b) < 0)
  {
    mpz_vector_swap (a, b);
    mpz_swap (norm_a, norm_b);
  }
  while (mpz_cmp (norm_a, norm_b) > 0)
  {
    if (skew)
      mpz_vector_skew_dot_product (tmp, a, b, skewness);
    else
      mpz_vector_dot_product (tmp, a, b);
    mpz_ndiv_q (q, tmp, norm_b);
    mpz_vector_submul (a, q, b);
    mpz_vector_swap (a, b);
    mpz_set (norm_a, norm_b);
    if (skew)
      mpz_vector_skew_norm (norm_b, b, skewness);
    else
      mpz_vector_norm (norm_b, b);
  }

  mpz_clear (norm_a);
  mpz_clear (norm_b);
  mpz_clear (q);
  mpz_clear (r);
  mpz_clear (tmp);
}


/* Compute f and g the two polynomials found via Montgomery's Two
   Quadratics Method, given N, P and m.
   P and N should be coprime and m^2-N should be divisible by P.
   The geometric progression is [c0, c1, c2] = [ p, m, (m^2-N)/p ]
   We start with a = [m, -p, 0] and b = [(m*t-c2)/P, -t, 1]
      where t = c2/m mod P
   We compute a reduced basis of the lattice spanned by a and b with the
   Lagrange algorithm.
 */
void MontgomeryTwoQuadratics (mpz_poly_t f, mpz_poly_t g, mpz_t skew, mpz_t N,
                              mpz_t P, mpz_t m, mpz_t max_skewness)
{
  ASSERT_ALWAYS (mpz_coprime_p(P, N));

  mpz_vector_t a, b, reduced_a, reduced_b;
  mpz_vector_init (a, 3);
  mpz_vector_init (b, 3);
  mpz_vector_init (reduced_a, 3);
  mpz_vector_init (reduced_b, 3);

  mpz_t tmp, c2, t, max_skew, min_skew;
  mpz_init(tmp);
  mpz_init(c2);
  mpz_init(t);
  mpz_init(max_skew);
  mpz_init(min_skew);

  // compute c2
  mpz_mul (tmp, m, m);
  mpz_sub (tmp, tmp, N);
  ASSERT_ALWAYS (mpz_divisible_p(tmp, P));
  mpz_divexact (c2, tmp, P);    // c2 = (m^2 - N) / P

  //compute t
  int ret = mpz_invert (tmp, m, P);
  ASSERT_ALWAYS (ret != 0);
  mpz_mul (t, tmp, c2);          // t = c2 / m mod P

  // Set vector a
  mpz_vector_setcoordinate (a, 0, m); // a[0] = m
  mpz_neg (tmp, P);
  mpz_vector_setcoordinate (a, 1, tmp); // a[1] = -P
  mpz_vector_setcoordinate_ui (a, 2, 0); // a[2] = 0

  // Set vector b
  mpz_mul (tmp, m, t);
  mpz_sub (tmp, tmp, c2);
  ASSERT_ALWAYS (mpz_divisible_p(tmp, P));
  mpz_divexact (tmp, tmp, P);
  mpz_vector_setcoordinate (b, 0, tmp); // b[0] = (m*t - c2) / P
  mpz_neg (tmp, t);
  mpz_vector_setcoordinate (b, 1, tmp); // b[1] = -t
  mpz_vector_setcoordinate_ui (b, 2, 1); // b[2] = 1

  /* Find the maximun value of skew for which LagrangeAlgo returns 2 polynomials
   * of degree 2. By taking skew <= default_max_skew, we already know that the
   * first polynomial (the one corresponding to a) will have degree 2, so we
   * just check that b[2] != 0.
   */
  mpz_set_ui (min_skew, 1);
  compute_default_max_skew (max_skew, m);
  // Overwrite the maximum default value if max_skewness argument is greater
  // than 0 and lesser than default value
  if (mpz_cmp_ui(max_skewness, 0) > 0 && mpz_cmp(max_skewness, max_skew) < 0)
    mpz_set (max_skew, max_skewness);

  mpz_set (skew, max_skew);

  do
  {
    mpz_vector_set (reduced_a, a);
    mpz_vector_set (reduced_b, b);

    // Perform Lagrange algorithm a and b
    LagrangeAlgo (reduced_a, reduced_b, skew);

    if (mpz_vector_is_coordinate_zero (reduced_b, 2))
      mpz_set (max_skew, skew);
    else
    {
      mpz_sub (tmp, max_skew, skew);
      if (mpz_cmp_ui (tmp, 1) == 0)
      {
        mpz_set (min_skew, skew);
        mpz_set (max_skew, skew);
      }
      else
        mpz_set (min_skew, skew);
    }
    mpz_add (skew, min_skew, max_skew);
    mpz_tdiv_q_2exp (skew, skew, 1);
  } while (mpz_cmp(min_skew, max_skew) < 0); // min_skew < max_skew

  mpz_vector_get_mpz_poly(f, reduced_a);
  mpz_vector_get_mpz_poly(g, reduced_b);

  mpz_clear(tmp);
  mpz_clear(c2);
  mpz_clear(t);
  mpz_clear(max_skew);
  mpz_clear(min_skew);
  mpz_vector_clear (a);
  mpz_vector_clear (b);
  mpz_vector_clear (reduced_a);
  mpz_vector_clear (reduced_b);
}


// note := exp(alpha(f)) * exp(alpha(g)) * sqrt(exp(norm(f)) * exp(norm(g)))
static inline double
notation (mpz_poly_t f, mpz_poly_t g)
{
  double alphaf = get_alpha (f, ALPHA_BOUND);
  double alphag = get_alpha (g, ALPHA_BOUND);
  double skewness = L2_skewness (f, SKEWNESS_DEFAULT_PREC);
  double normf = L2_lognorm (f, skewness);
  double normg = L2_lognorm (g, skewness);
  return exp(alphaf + alphag + (normf + normg)/2);
}



static void declare_usage(param_list pl)
{
  param_list_decl_usage(pl, "N", "input number (default c59)");
  param_list_decl_usage(pl, "minP", "Use P > minP (default 2)");
  param_list_decl_usage(pl, "maxP", "Use P <= maxP (default nextprime(minP))");
  param_list_decl_usage(pl, "skewness", "maximun skewness possible "
                                        "(default (1/6)^(1/4) * m^(1/2))");
  param_list_decl_usage(pl, "v", "verbose output (print all polynomials)");
  param_list_decl_usage(pl, "q", "quiet output (print only best polynomials)");
}

static void usage (const char *argv, param_list pl)
{
  param_list_print_usage(pl, argv, stderr);
  exit (EXIT_FAILURE);
}

int
main (int argc, char *argv[])
{
  char *argv0 = argv[0];
  int verbose = 0, quiet = 0;

  mpz_t N, minP, maxP, max_skewness;
  mpz_init (N);
  mpz_init (minP);
  mpz_init (maxP);
  mpz_init (max_skewness);

  param_list pl;

  /* read params */
  param_list_init(pl);
  declare_usage(pl);

  param_list_configure_switch (pl, "-v", &verbose);
  param_list_configure_switch (pl, "-q", &quiet);

  if (argc == 1)
    usage (argv[0], pl);

  argc--,argv++;
  for ( ; argc ; )
  {
    if (param_list_update_cmdline (pl, &argc, &argv)) { continue; }
    fprintf (stderr, "Unhandled parameter %s\n", argv[0]);
    usage (argv0, pl);
  }

  if (!param_list_parse_mpz(pl, "skewness", max_skewness))
    mpz_set_ui (max_skewness, 0);
  else if (mpz_cmp_ui (max_skewness, 1) < 0)
  {
    gmp_fprintf(stderr, "Error, skewness (%Zd) should be greater or equal "
                        "to 1\n", max_skewness);
    abort();
  }

  if (!param_list_parse_mpz(pl, "minP", minP))
    mpz_set_ui (minP, 2);
  if (!param_list_parse_mpz(pl, "maxP", maxP))
    mpz_nextprime(maxP, minP);
  if (!param_list_parse_mpz(pl, "N", N))
    mpz_set_str (N,
        "71641520761751435455133616475667090434063332228247871795429", 10);

  if (mpz_cmp_ui (minP, 0) < 0)
  {
    gmp_fprintf(stderr, "Error, minP (%Zd) should be greater or equal to 0\n",
                        minP);
    abort();
  }
  if (mpz_cmp (maxP, minP) < 0)
  {
    gmp_fprintf(stderr, "Error, maxP (%Zd) should be greater or equal to minP "
                        "(%Zd)\n", maxP, minP);
    abort();
  }
  if (mpz_cmp_ui (N, 3) < 0)
  {
    gmp_fprintf(stderr, "Error, N (%Zd) should be greater or equal to 3\n", N);
    abort();
  }
  if (mpz_divisible_ui_p (N, 2))
  {
    gmp_fprintf(stderr, "Error, N (%Zd) should not be divisible by 2\n", N);
    abort();
  }

  if (quiet)
    verbose = -1;

  if (param_list_warn_unused(pl))
    usage (argv0, pl);
  param_list_print_command_line (stdout, pl);

  gmp_printf("### N = %Zd\n", N);

/* Given an integer N and bounds on P, tests all pairs of polynomials
   with P such that minP < P <= maxP and print the pair having the best
   rating, along with its skewness, norm, alpha(f), alpha(g), and the
   prime number P which gave it
 */
  cado_poly_extended best_poly, poly;
  mpz_t P, sqrtN, r, m, skew_used;
  mpz_poly_t f, g;

  cado_poly_extended_init (best_poly);
  cado_poly_extended_init (poly);
  cado_poly_extended_set_note (best_poly, DBL_MAX);
  mpz_init(P);
  mpz_init(sqrtN);
  mpz_init(r);
  mpz_init(m);
  mpz_init(skew_used);
  mpz_poly_init (f, 2);
  mpz_poly_init (g, 2);

  // sqrtN = floor(sqrt(N))
  mpz_sqrt(sqrtN, N);

  mpz_nextprime (P, minP);
  while(mpz_cmp (P, maxP) <= 0)
  {
    if (mpz_kronecker(N, P) == 1)
    {
      for (unsigned int k = 0; k < 2; k++)
      {
        if (k == 0) // For the first time compute r = sqrt of N modulo P
        {
          int ret = mpz_msqrt (r, N, P);
          ASSERT_ALWAYS (ret != 0);
        }
        else // The second root mod P is P-r
          mpz_sub(r, P, r);

        if (verbose >= 0)
          gmp_printf("### P = %Zd: compute polynomials for r = %Zd\n", P, r);
        // m is the first integer >= sqrtN such that m = r (mod P)
        compute_m (m, sqrtN, r, P);
        MontgomeryTwoQuadratics (f, g, skew_used, N, P, m, max_skewness);
        double note = notation (f, g);
        if(note < cado_poly_extended_get_note(best_poly))
          cado_poly_extended_set (best_poly, f, g, N, m, P, skew_used, note);
        if (verbose >= 1)
        {
          cado_poly_extended_set (poly, f, g, N, m, P, skew_used, note);
          cado_poly_extended_print (stdout, poly, "# ");
        }
      }
    }
    else if (verbose >= 0)
      gmp_printf("### P = %Zd: N is not a square modulo P, skipping it.\n", P);
    mpz_nextprime(P, P);
  }

  printf("### Best polynomials found:\n");
  cado_poly_extended_print (stdout, best_poly, "");

  mpz_clear(P);
  mpz_clear(sqrtN);
  mpz_clear(r);
  mpz_clear(m);

  mpz_poly_clear (f);
  mpz_poly_clear (g);

  cado_poly_extended_clear (poly);
  cado_poly_extended_clear (best_poly);
  mpz_clear (N);
  mpz_clear (skew_used);
  mpz_clear (max_skewness);
  mpz_clear (minP);
  mpz_clear (maxP);
  param_list_clear(pl);
  return EXIT_SUCCESS;
}

