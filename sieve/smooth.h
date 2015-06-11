/* for batch cofactorization */

#include <stdint.h>
#include <gmp.h>

typedef struct {
  int64_t *a;
  uint64_t *b;
  mpz_t *R;    /* cofactors on side 0 (the "rational" side) */
  mpz_t *A;    /* cofactors on side 1 (the "algebraic" side) */
  size_t alloc;
  size_t size;
} cofac_list_t;
typedef cofac_list_t cofac_list[1];

typedef struct {
  mpz_t *l;
  size_t alloc;
  size_t size;
} mpz_list_t;
typedef mpz_list_t mpz_list[1];
