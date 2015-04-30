/* Test */

#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>
#include "cado.h"
#include "utils.h"
#include "functions.h"
#include "sieving_bound.h"
#include "parameters.h"

  /*function_special_q(f0, f1, g, a, b, c, p, h, coeff0, coeff1, q, t);*/
void declare_usage(param_list pl)
{
  param_list_decl_usage(pl, "p", "prime number");
  param_list_decl_usage(pl, "h", "polynomial to build f0 and f1");
  param_list_decl_usage(pl, "coeff0", "lowest coefficient of g");
  param_list_decl_usage(pl, "coeff1", "greatest coefficient of g");
  param_list_decl_usage(pl, "q", "a typical q of a special-q");
  param_list_decl_usage(pl, "t", "dimension of the lattice");

}

void initialise_parameters(int argc, char * argv[], mpz_ptr p, mpz_poly_ptr h,
    int * coeff0, int * coeff1, mpz_ptr q, unsigned int * t)
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

  param_list_parse_mpz(pl, "q", q);

  param_list_parse_uint(pl, "t", t);
  ASSERT(* t > 2);

  param_list_clear(pl);
}

int main(int argc, char * argv[])
{
  mpz_t p;
  mpz_init(p);
  mpz_poly_t h;
  mpz_poly_init(h, -1);
  mpz_t q;
  mpz_init(q);
  int coeff0;
  int coeff1;
  unsigned int t;

  initialise_parameters(argc, argv, p, h, &coeff0, &coeff1, q, &t);

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

  function_special_q(f0, f1, g, a, b, c, p, h, coeff0, coeff1, q, t);
  printf("f0: ");
  mpz_poly_fprintf(stdout, f0);
  printf("f1: ");
  mpz_poly_fprintf(stdout, f1);
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
  mpz_clear(q);

  return 0;
}

/*int main(int argc, char * argv[])*/
/*{*/
  /*mpz_poly_t f0;*/
  /*mpz_poly_init(f0, 6);*/
  /*mpz_poly_t f1;*/
  /*mpz_poly_init(f1, 6);*/
  /*mpz_poly_t g;*/
  /*mpz_poly_init(g, 6);*/
  /*mpz_poly_t h;*/
  /*mpz_poly_init(h, 6);*/
  /*for (int i = 0; i < h->alloc; i++) {*/
    /*mpz_poly_setcoeff_si(h, i, 1);*/
  /*}*/

  /*mpz_t a;*/
  /*mpz_init(a);*/
  /*mpz_t b;*/
  /*mpz_init(b);*/
  /*mpz_t c;*/
  /*mpz_init(c);*/
  /*mpz_t p;*/
  /*mpz_init(p);*/
  /*[>mpz_set_str(p, "10822639589", 10);<]*/
  /*mpz_set_str(p, "122663", 10);*/
  /*mpz_t q;*/
  /*mpz_init(q);*/
  /*[>mpz_set_str(q, "33554432", 10);<]*/
  /*mpz_set_str(q, "1", 10);*/

  /*int coeff0 = -2;*/
  /*int coeff1 = 2;*/
  /*unsigned int t = 7;*/

  /*function_special_q(f0, f1, g, a, b, c, p, h, coeff0, coeff1, q, t);*/
  /*printf("f0: ");*/
  /*mpz_poly_fprintf(stdout, f0);*/
  /*printf("f1: ");*/
  /*mpz_poly_fprintf(stdout, f1);*/
  /*gmp_printf("a: %Zd, b: %Zd, c: %Zd\n", a, b, c);*/
  /*printf("g: ");*/
  /*mpz_poly_fprintf(stdout, g);*/
  /*printf("h: ");*/
  /*mpz_poly_fprintf(stdout, h);*/

  /*mpz_poly_clear(g);*/
  /*mpz_poly_clear(h);*/

  /*mpz_clear(a);*/
  /*mpz_clear(b);*/
  /*mpz_clear(c);*/
  /*mpz_clear(p);*/
  /*mpz_clear(q);*/

  /*Other tests*/

  /*uint64_t number_element = 8388608;*/
  /*sieving_bound_t H;*/
  /*sieving_bound_init(H, 3);*/
  /*sieving_region_adapted(H, number_element);*/
  /*sieving_bound_fprintf(stdout, H);*/

  /*unsigned int number_fields = 2;*/
  /*mpz_t * mean = malloc(sizeof(mpz_t) * number_fields);*/
  /*mpz_poly_t * f = malloc(sizeof(mpz_poly_t) * number_fields);*/
  /*for (unsigned int i = 0; i < number_fields; i++) {*/
    /*mpz_init(mean[i]);*/
    /*mpz_set_ui(mean[i], 0);*/
    /*mpz_poly_init(f[i], -1);*/
  /*}*/
  /*mpz_poly_set(f[0], f0);*/
  /*mpz_poly_set(f[1], f1);*/
  /*uint64_t number_a = 100000;*/
  /*mean_approx_number(mean, f, number_fields, number_a, H);*/
  /*for (unsigned int i = 0; i < number_fields; i++) {*/
    /*printf("%zu\n", mpz_sizeinbase(mean[i], 2));*/
  /*}*/
  /*for (unsigned int i = 0; i < number_fields; i++) {*/
    /*mpz_clear(mean[i]);*/
    /*mpz_poly_clear(f[i]);*/
  /*}*/
  /*sieving_bound_clear(H)*/
  /*mpz_poly_clear(f0);*/
  /*mpz_poly_clear(f1);*/
  /*free(f); */
  /*free(mean); */

  /*return 0;*/
/*}*/
