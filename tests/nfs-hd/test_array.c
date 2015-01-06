#include "sieving_interval.h"
#include "array.h"
#include "cado.h"
#include "utils.h"
#include "macros.h"

int main()
{
  sieving_interval_t H;
  unsigned int t = 4;

  sieving_interval_init(H, t);
  for (unsigned int i = 0; i < t; i++) {
    sieving_interval_set_hi(H, i, 2);
  }

  uint64_t nb;
  sieving_interval_number_element(&nb, H);

  mpz_vector_t vector;
  mpz_vector_init(vector, t);

  for (unsigned int i = 0; i < t - 1; i++) {
    mpz_vector_setcoordinate_si(vector, i, -H->h[i]);
  }

  uint64_t index;
  uint64_t index1;

  mpz_vector_t vector1;
  mpz_vector_init(vector1, t);

  mpz_poly_t poly;
  mpz_poly_init(poly, 0);

  mpz_poly_t poly1;
  mpz_poly_init(poly1, 0);

  array_mpz_vector_index(&index, vector, H, nb);
  ASSERT_ALWAYS(index == 0);
  for (unsigned int i = 1; i < nb; i++) {
    mpz_vector_add_one(vector, H);
    array_mpz_vector_index(&index, vector, H, nb);
    ASSERT_ALWAYS(index == i);
    array_index_mpz_vector(vector1, index, H, nb);
    array_index_mpz_poly(poly, index, H, nb);
    ASSERT_ALWAYS(mpz_vector_cmp(vector, vector1) == 0);
    mpz_vector_get_mpz_poly(poly1, vector1);
    ASSERT_ALWAYS(mpz_poly_cmp(poly1, poly) == 0);
    array_mpz_poly_index(&index1, poly1, H, nb);
    ASSERT_ALWAYS(index1 == index);
  }

  mpz_vector_clear(vector1);
  mpz_vector_clear(vector);
  mpz_poly_clear(poly);
  mpz_poly_clear(poly1);
  sieving_interval_clear(H);

  return 0;
}
