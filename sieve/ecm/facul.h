#include <stdio.h>
#include <gmp.h>

typedef struct {
  long method; /* Which method to use (P-1, P+1 or ECM) */
  void *plan;  /* Parameters for that method */
} facul_strategy_t;

facul_strategy_t * facul_make_strategy (const int n);
void facul_clear_strategy (facul_strategy_t *);
void facul_print_stats (FILE *);
int facul (unsigned long *, const mpz_t, facul_strategy_t *);
