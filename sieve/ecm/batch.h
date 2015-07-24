#ifndef COFAC_LIST_H
#define COFAC_LIST_H

#include <stdlib.h>
#include <inttypes.h>
#include <gmp.h>
#include "facul.h"
#include "facul_doit.h"
#include "utils.h"

#define STATUS_SMOOTH  0
#define STATUS_UNKNOWN 1
#define STATUS_USELESS 2

/* TODO: should replace by mpz_array_t from utils/mpz_array.h */
typedef struct {
  mpz_t *l;
  size_t alloc;
  size_t size;
} mpz_list_t;
typedef mpz_list_t mpz_list[1];

typedef struct {
  int64_t *a;
  uint64_t *b;
  mpz_t *R;    /* cofactors on side 0 */
  mpz_t *A;    /* cofactors on side 1 */
  size_t alloc;
  size_t size;
} cofac_list_t;
typedef cofac_list_t cofac_list[1];

#ifdef __cplusplus
extern "C" {
#endif

void cofac_list_init (cofac_list);
void cofac_list_realloc (cofac_list, size_t);
void cofac_list_clear (cofac_list);
void cofac_list_add (cofac_list, long, unsigned long, mpz_t, mpz_t);
unsigned long prime_product (mpz_t, prime_info, unsigned long, unsigned long);
void find_smooth (cofac_list, int, int, unsigned long, unsigned long,
                  FILE*, FILE*, int);
void factor (cofac_list, cado_poly, int, int, int);
void create_batch_file (const char*, unsigned long, unsigned long,
                        cado_poly, int);

#ifdef __cplusplus
}
#endif

#endif /* COFAC_LIST_H */
