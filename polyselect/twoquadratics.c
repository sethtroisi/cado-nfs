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
typedef const mpz_vector_struct_t * mpz_vector_srcptr;

/* Management of the structure, set and print coefficients. */
void mpz_vector_init(mpz_vector_t v, unsigned int l)
{
  v->dim = l;
  ASSERT_ALWAYS (l > 0);
  v->c = (mpz_t *) malloc (l * sizeof(mpz_t));
  ASSERT_ALWAYS (v->c != NULL);
  for (unsigned int i = 0; i < l; i++)
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

/* compute r <- r - q*v */
void
mpz_vector_submul (mpz_vector_t r, mpz_t q, mpz_vector_t v)
{
  ASSERT_ALWAYS (r->dim == v->dim);
  for (unsigned int i = 0; i < r->dim; i++)
    mpz_submul(r->c[i], q, v->c[i]);
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

/* Set res to the smaller integer greater or equal to b, such that
   res = r mod M
 */
void
compute_m (mpz_t res, mpz_t b, mpz_t r, mpz_t M)
{
  mpz_mod (res, b, M);   // res = b mod M
  mpz_sub (res, r, res); // res = r - (b mod M)
  if (mpz_sgn (res) < 0) // add M if res < 0
    mpz_add (res, res, M);
  mpz_add (res, b, res); // res = b + r - (b mod M) (and maybe + M)
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
void LagrangeAlgo (mpz_vector_t a, mpz_vector_t b)
{
  mpz_t norm_a, norm_b, q, r, tmp;
  mpz_init (norm_a);
  mpz_init (norm_b);
  mpz_init (q);
  mpz_init (r);
  mpz_init (tmp);

  mpz_vector_norm (norm_a, a);
  mpz_vector_norm (norm_b, b);
  if (mpz_cmp (norm_a, norm_b) < 0)
  {
    mpz_vector_swap (a, b);
    mpz_swap (norm_a, norm_b);
  }
  while (mpz_cmp (norm_a, norm_b) > 0)
  {
    mpz_vector_dot_product (tmp, a, b);
    mpz_ndiv_q (q, tmp, norm_b);
    mpz_vector_submul (a, q, b); 
    mpz_vector_swap (a, b);
    mpz_set (norm_a, norm_b);
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
void MontgomeryTwoQuadratics (mpz_poly_t f, mpz_poly_t g, mpz_t N,
                              mpz_t P, mpz_t m)
{
  ASSERT_ALWAYS (mpz_coprime_p(P, N));

  mpz_vector_t a, b;
  mpz_vector_init (a, 3);
  mpz_vector_init (b, 3);

  mpz_t tmp, c2, t;
  mpz_init(tmp);
  mpz_init(c2);
  mpz_init(t);
  
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

  // Perform Lagrange algorithm a and b
  LagrangeAlgo (a, b);

  mpz_clear(tmp);
  mpz_clear(c2);
  mpz_clear(t);

  mpz_vector_get_mpz_poly(f, a);
  mpz_vector_get_mpz_poly(g, b);

  mpz_vector_clear (a);
  mpz_vector_clear (b);
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
  param_list_decl_usage(pl, "maxP", "Use P <= maxP (default 500)");
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

  mpz_t N, minP, maxP;
  mpz_init (N);
  mpz_init (minP);
  mpz_init (maxP);

  param_list pl;
  
  /* read params */
  param_list_init(pl);
  declare_usage(pl);

  if (argc == 1)
    usage (argv[0], pl);

  argc--,argv++;
  for ( ; argc ; )
  {
    if (param_list_update_cmdline (pl, &argc, &argv)) { continue; }
    fprintf (stderr, "Unhandled parameter %s\n", argv[0]);
    usage (argv0, pl);
  }

  if (!param_list_parse_mpz(pl, "minP", minP))
    mpz_set_ui (minP, 2);
  if (!param_list_parse_mpz(pl, "maxP", maxP))
    mpz_set_ui (maxP, 500);
  if (!param_list_parse_mpz(pl, "N", N))
    mpz_set_str (N,
        "71641520761751435455133616475667090434063332228247871795429", 10);

  if (param_list_warn_unused(pl))
    usage (argv0, pl);
  param_list_print_command_line (stdout, pl);

  gmp_printf("N = %Zd\n", N);

/* Given an integer N and bounds on P, tests all pairs of polynomials
   with P such that minP < P <= maxP and print the pair having the best
   rating, along with its skewness, norm, alpha(f), alpha(g), and the
   prime number P which gave it
 */
  double best_note = INFINITY, cur_note;
  mpz_t P, best_P, sqrtN, r, m, best_m;
  mpz_poly_t f, best_f, g, best_g;
  
  mpz_init(P);
  mpz_init(best_P);
  mpz_init(sqrtN);
  mpz_init(r);
  mpz_init(m);
  mpz_init(best_m);
  mpz_poly_init (f, 2);
  mpz_poly_init (best_f, 2);
  mpz_poly_init (g, 2);
  mpz_poly_init (best_g, 2);

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

        // m is the first integer >= sqrtN such that m = r (mod P)
        compute_m (m, sqrtN, r, P);
        MontgomeryTwoQuadratics (f, g, N, P, m);
        cur_note = notation (f, g);
        if(cur_note < best_note)
        {
          mpz_poly_copy (best_f, f);
          mpz_poly_copy (best_g, g);
          best_note = cur_note;
          mpz_set (best_P, P);
          mpz_set (best_m, m);
        }
      }
    }
    mpz_nextprime(P, P);
  }

  printf("# Best polynomials found:\nf = ");
  mpz_poly_fprintf (stdout, best_f);
  printf("g = ");
  mpz_poly_fprintf (stdout, best_g);

  printf ("note = %e\n", best_note);
  printf ("alpha_f = %.2f\n", get_alpha (best_f, ALPHA_BOUND));
  printf ("alpha_g = %.2f\n", get_alpha (best_g, ALPHA_BOUND));
  double skewness = L2_skewness (best_f, SKEWNESS_DEFAULT_PREC);
  printf ("Skewness = %f\n", skewness);
  printf ("L2_skew_norm_f = %.2f\n", L2_lognorm (best_f, skewness));
  printf ("L2_skew_norm_g = %.2f\n", L2_lognorm (best_g, skewness));
  gmp_printf ("P = %Zd\n", best_P);
  gmp_printf ("m = %Zd\n", best_m);
  mpz_mod (r, best_m, best_P);
  gmp_printf ("# r = m (mod P)\nr = %Zd\n", r);
  mpz_invert (r, best_P, N);
  mpz_mul (r, r, best_m);
  mpz_mod (r, r, N);
  gmp_printf ("# common_root = m / P (mod  N)\ncommon_root = %Zd\n", r);

  mpz_clear(best_P);
  mpz_clear(P);
  mpz_clear(sqrtN);
  mpz_clear(r);
  mpz_clear(m);
  mpz_clear(best_m);

  mpz_poly_clear (f);
  mpz_poly_clear (best_f);
  mpz_poly_clear (g);
  mpz_poly_clear (best_g);

  mpz_clear (N);
  mpz_clear (minP);
  mpz_clear (maxP);
  param_list_clear(pl);
  return EXIT_SUCCESS;
}

