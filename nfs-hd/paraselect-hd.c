#include <stdlib.h>
#include "cado.h"
#include "utils.h"
#include "parameters.h"

int main()
{
  mpz_t p;
  mpz_init(p);
  mpz_set_str(p, "10822639589", 10);
  unsigned int n = 6;

  uint64_t number_a = 100;

  unsigned int lpb_min = 25;
  unsigned int lpb_max = 26;

  unsigned int t_min = 3;
  unsigned int t_max = 4;

  uint64_t size_start = 34359738368;

  mpz_poly_t h;
  mpz_poly_init(h, n);
  for (int k = 0; k <= (int)n; k++) {
    mpz_poly_setcoeff_si(h, k, 1);
  }

  int coeff0 = -2;
  int coeff1 = 2;

  unsigned int nb_times = 2;
  
  double weight_0 = 0.5;
  double weight_1 = 0.5;

  mpz_poly_t f0;
  mpz_poly_init(f0, n);
  mpz_poly_t f1;
  mpz_poly_init(f1, n);

  find_parameters_adapted(p, n, number_a, lpb_min, lpb_max, t_min, t_max,
      size_start, h, coeff0, coeff1, nb_times, weight_0, weight_1);

  mpz_poly_clear(h);
  mpz_clear(p);
}
