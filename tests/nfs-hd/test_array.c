#include <stdint.h>
#include <stdio.h>
#include <inttypes.h>
#include "sieving_interval.h"
#include "array.h"
#include "cado.h"
#include "utils.h"
#include "macros.h"
#include "int64_vector.h"

int main()
{
  sieving_interval_t H;
  unsigned int t = 4;

  sieving_interval_init(H, t);
  for (unsigned int i = 0; i < t; i++) {
    sieving_interval_set_hi(H, i, 2);
  }

  uint64_t nb = sieving_interval_number_element(H);

  int64_vector_t vector;
  int64_vector_init(vector, t);

  mpz_vector_t v;
  mpz_vector_init(v, t);

  for (unsigned int i = 0; i < t - 1; i++) {
    int64_vector_setcoordinate(vector, i, -(int64_t)H->h[i]);
  }

  uint64_t index;
  uint64_t index1;

  array_int64_vector_index(&index, vector, H, nb);
  ASSERT_ALWAYS(index == 0);
  for (unsigned int i = 1; i < nb; i++) {
    int64_vector_add_one(vector, H);
    array_int64_vector_index(&index, vector, H, nb);
    ASSERT_ALWAYS(index == i);
    array_index_mpz_vector(v, index, H, nb);
    for (unsigned int j = 0; j < t; j++) {
      ASSERT_ALWAYS(mpz_cmp_si(v->c[j], vector->c[j] == 0));
    }
    array_mpz_vector_index(&index1, v, H, nb);
    ASSERT_ALWAYS(index == index1);
  }

  mpz_vector_clear(v);
  int64_vector_clear(vector);
  sieving_interval_clear(H);

  return 0;
}
