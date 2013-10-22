#ifndef FILTER_UTILS_H_
#define FILTER_UTILS_H_

#define NB_PRIMES_OPT 31

typedef struct {
  index_t h;
  p_r_values_t p;
  exponent_t e;
} prime_t;

typedef struct {
  index_t id;
  int32_t e; //or exponent_t ??
} ideal_merge_t;

typedef struct {
  uint64_t nrels;
  uint64_t nprimes;
  double W; //weight of the active part of the matrix
} info_mat_t;

#define CA_DUP1 314159265358979323UL
#define CB_DUP1 271828182845904523UL

#define CA_DUP2 271828182845904523UL
#define CB_DUP2 577215664901532889UL

#include <time.h>
#include <pthread.h>

#include "filter_io.h"
#include "filter_memalloc.h"


unsigned int insert_rel_in_table_no_e (earlyparsed_relation_ptr, index_t, index_t **, weight_t *);
unsigned int insert_rel_in_table_with_e (earlyparsed_relation_ptr, index_t,
                                         ideal_merge_t **, int32_t *);

#endif /* FILTER_UTILS_H_ */
