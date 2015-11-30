#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>
#include "cado.h"
#include "utils.h"
#include "polynomials.h"
#include "sieving_bound.h"
#include "polyselect/auxiliary.h"

static void mpz_poly_coeff_fprintf(FILE * file, mpz_poly_srcptr f)
{
  for (int i = 0; i < f->deg; i++) {
    gmp_fprintf(file, "%Zd,", f->coeff[i]);
  }
  gmp_fprintf(file, "%Zd\n", f->coeff[f->deg]);
}

  /*function_special_q(f0, f1, g, a, b, c, p, h, coeff0, coeff1, q, t);*/
void declare_usage(param_list pl)
{
  param_list_decl_usage(pl, "p", "prime number");
  param_list_decl_usage(pl, "n", "degree of the extension");
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
  param_list_decl_usage(pl, "type", "0 if classical, 1 if special-q");
  param_list_decl_usage(pl, "out", "path to the output file");
}

void initialise_parameters(int argc, char * argv[], mpz_ptr p, mpz_poly_ptr h,
    int * coeff0, int * coeff1, unsigned int * q, unsigned int * t,
    unsigned int * nb_times, double * weight_0, double * weight_1,
    unsigned int * n, unsigned int * type, FILE ** outstd)
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

  * type = 1;
  param_list_parse_uint(pl, "type", type);
  param_list_parse_mpz(pl, "p", p);
  param_list_parse_uint(pl, "n", n);
  * nb_times = 100;
  param_list_parse_uint(pl, "nb_times", nb_times);
  ASSERT(* nb_times > 0);

  * weight_0 = 0.5;
  * weight_1 = 0.5;
  param_list_parse_double(pl, "weight_0", weight_0);
  param_list_parse_double(pl, "weight_1", weight_1);
  ASSERT(* weight_0 >= 0.0);
  ASSERT(* weight_1 >= 0.0);
  ASSERT(* weight_0 + * weight_1 == 1.0);

  param_list_parse_int(pl, "coeff0", coeff0);
  param_list_parse_int(pl, "coeff1", coeff1);
  ASSERT(coeff0 < coeff1);

  if (* type == 1) {
    param_list_parse_mpz_poly(pl, "h", h, ",");
    ASSERT(* n == (unsigned int )h->deg);

    param_list_parse_uint(pl, "q", q);

    param_list_parse_uint(pl, "t", t);
    ASSERT(* t > 2);
  }

  unsigned int size_path = 1024;
  char path [size_path];
  path[0] = '\0';
  param_list_parse_string(pl, "out", path, size_path);
  if (path[0] == '\0') {
    * outstd = stdout;
  } else {
    * outstd = fopen(path, "w");
  }

  param_list_clear(pl);
}

int main(int argc, char * argv[])
{
  mpz_t p;
  mpz_init(p);
  unsigned int n;
  mpz_poly_t h;
  mpz_poly_init(h, -1);
  unsigned int q;
  int coeff0;
  int coeff1;
  unsigned int t;
  unsigned int nb_times;
  double weight_0;
  double weight_1;
  unsigned int type;
  FILE * outstd;

  initialise_parameters(argc, argv, p, h, &coeff0, &coeff1, &q, &t, &nb_times,
      &weight_0, &weight_1, &n, &type, &outstd);

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

  gmp_fprintf(outstd, "# p^n: %Zd^%u\n", p, n);
  gmp_fprintf(outstd, "n: %Zd\n", p);
  //TODO: check if it is needed or not for cado_poly.
  fprintf(outstd, "skew: 1.000\n");
  fprintf(outstd, "# t = %u\n", t);

  if (type == 0) {
    function_classical(f0, f1, p, n, coeff0, coeff1, nb_times, weight_0,
        weight_1);
  } else if (type == 1) {
    function_special_q(f0, f1, g, a, b, c, p, h, coeff0, coeff1, q, t, nb_times,
        weight_0, weight_1);
  }

#ifndef NDEBUG
  mpz_t res;
  mpz_init(res);
  mpz_t tmp;
  mpz_init(tmp);
  mpz_poly_resultant(res, f0, f1);
  mpz_mod(tmp, res, p);
  ASSERT(mpz_cmp_ui(tmp, 0) == 0);
  mpz_pow_ui(tmp, p, n);
  mpz_mod(tmp, res, tmp);
  ASSERT(mpz_cmp_ui(tmp, 0) == 0);
  gmp_fprintf(outstd, "# Resultant of poly0 and poly1: %Zd\n", res);
  mpz_clear(res);
  mpz_clear(tmp);
#endif // NDEBUG

  fprintf(outstd, "# poly0: ");
  mpz_poly_fprintf(outstd, f0);
  fprintf(outstd, "poly0: ");
  mpz_poly_coeff_fprintf(outstd, f0);
  fprintf(outstd, "# alpha0: %f\n", get_alpha(f0, ALPHA_BOUND));
  fprintf(outstd, "# poly1: ");
  mpz_poly_fprintf(outstd, f1);
  fprintf(outstd, "poly1: ");
  mpz_poly_coeff_fprintf(outstd, f1);
  fprintf(outstd, "# alpha1: %f\n", get_alpha(f1, ALPHA_BOUND));

  if (type == 1) {
    fprintf(outstd, "# log2(q) = %u\n", q);
    gmp_fprintf(outstd, "# a: %Zd, b: %Zd, c: %Zd\n", a, b, c);
    fprintf(outstd, "# g: ");
    mpz_poly_fprintf(outstd, g);
    fprintf(outstd, "# h: ");
    mpz_poly_fprintf(outstd, h);
  }

  mpz_poly_clear(g);
  mpz_poly_clear(h);
  mpz_poly_clear(f0);
  mpz_poly_clear(f1);

  mpz_clear(a);
  mpz_clear(b);
  mpz_clear(c);
  mpz_clear(p);

  fclose(outstd);

  return 0;
}
