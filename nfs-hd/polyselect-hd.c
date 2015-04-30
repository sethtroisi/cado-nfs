#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>
#include "cado.h"
#include "utils.h"
#include "functions.h"
#include "sieving_bound.h"
#include "parameters.h"
#include "polyselect/auxiliary.h"

  /*function_special_q(f0, f1, g, a, b, c, p, h, coeff0, coeff1, q, t);*/
void declare_usage(param_list pl)
{
  param_list_decl_usage(pl, "p", "prime number");
  param_list_decl_usage(pl, "h", "polynomial to build f0 and f1");
  param_list_decl_usage(pl, "coeff0", "lowest coefficient of g");
  param_list_decl_usage(pl, "coeff1", "greatest coefficient of g");
  param_list_decl_usage(pl, "q", "a typical q of a special-q");
  param_list_decl_usage(pl, "t", "dimension of the lattice");
  param_list_decl_usage(pl, "nb_times", "number of tested polynomials");
  param_list_decl_usage(pl, "weight_0",
      "weight for alpha of f0 when sum alphas");
  param_list_decl_usage(pl, "weight_1",
      "weight for alpha of f1 when sum alphas");
}

void initialise_parameters(int argc, char * argv[], mpz_ptr p, mpz_poly_ptr h,
    int * coeff0, int * coeff1, unsigned int * q, unsigned int * t,
    unsigned int * nb_times, double * weight_0, double * weight_1)
{
  param_list pl;
  param_list_init(pl);
  declare_usage(pl);
  FILE * fpl;
  char * argv0 = argv[0];

  argv++, argc--;
  for( ; argc ; ) {
    if (param_list_update_cmdline(pl, &argc, &argv)) { continue; }

    /* Could also be a file */
    if ((fpl = fopen(argv[0], "r")) != NULL) {
      param_list_read_stream(pl, fpl, 0);
      fclose(fpl);
      argv++,argc--;
      continue;
    }

    fprintf(stderr, "Unhandled parameter %s\n", argv[0]);
    param_list_print_usage(pl, argv0, stderr);
    exit (EXIT_FAILURE);
  }

  param_list_parse_mpz(pl, "p", p);

  param_list_parse_mpz_poly(pl, "h", h, ","); 

  param_list_parse_int(pl, "coeff0", coeff0);
  param_list_parse_int(pl, "coeff1", coeff1);
  ASSERT(coeff0 < coeff1);

  param_list_parse_uint(pl, "q", q);

  param_list_parse_uint(pl, "t", t);
  ASSERT(* t > 2);
  
  param_list_parse_uint(pl, "nb_times", nb_times);
  ASSERT(* nb_times > 0);

  param_list_parse_double(pl, "weight_0", weight_0);
  param_list_parse_double(pl, "weight_1", weight_1);
  ASSERT(* weight_0 >= 0.0);
  ASSERT(* weight_1 >= 0.0);
  ASSERT(* weight_0 + * weight_1 == 1.0);

  param_list_clear(pl);
}

int main(int argc, char * argv[])
{
  mpz_t p;
  mpz_init(p);
  mpz_poly_t h;
  mpz_poly_init(h, -1);
  unsigned int q;
  int coeff0;
  int coeff1;
  unsigned int t;
  unsigned int nb_times;
  double weight_0;
  double weight_1;

  initialise_parameters(argc, argv, p, h, &coeff0, &coeff1, &q, &t, &nb_times,
      &weight_0, &weight_1);

  mpz_poly_t f0;
  mpz_poly_init(f0, h->deg);
  mpz_poly_t f1;
  mpz_poly_init(f1, h->deg);
  mpz_poly_t g;
  mpz_poly_init(g, h->deg);

  mpz_t a;
  mpz_init(a);
  mpz_t b;
  mpz_init(b);
  mpz_t c;
  mpz_init(c);

  function_special_q(f0, f1, g, a, b, c, p, h, coeff0, coeff1, q, t, nb_times,
      weight_0, weight_1);
  printf("f0: ");
  mpz_poly_fprintf(stdout, f0);
  printf("alpha0: %f\n", get_alpha(f0, ALPHA_BOUND));
  printf("f1: ");
  mpz_poly_fprintf(stdout, f1);
  printf("alpha1: %f\n", get_alpha(f1, ALPHA_BOUND));
  gmp_printf("a: %Zd, b: %Zd, c: %Zd\n", a, b, c);
  printf("g: ");
  mpz_poly_fprintf(stdout, g);
  printf("h: ");
  mpz_poly_fprintf(stdout, h);

  mpz_poly_clear(g);
  mpz_poly_clear(h);

  mpz_clear(a);
  mpz_clear(b);
  mpz_clear(c);
  mpz_clear(p);

  return 0;
}
