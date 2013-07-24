#ifndef FACUL_H
#define FACUL_H

#include <stdio.h>
#include <stdint.h>
#include <gmp.h>

#define PM1_METHOD 1
#define PP1_27_METHOD 2
#define PP1_65_METHOD 3
#define EC_METHOD 4

#define FACUL_NOT_SMOOTH (-1)

#define STATS_LEN 128

typedef struct {
  long method; /* Which method to use (P-1, P+1 or ECM) */
  void *plan;  /* Parameters for that method */
} facul_method_t;


/* All prime factors in the input number must be > fb. A factor of the 
   input number is assumed to be prime if it is < fb^2.
   The input number is taken to be not smooth if it has a 
   prime factor > 2^lpb. */

typedef struct {
  unsigned long lpb;        /* Large prime bound 2^lpb */
  uint64_t assume_prime_thresh; /* The factor base bound squared. If the 
                               square exceeds UINT64_MAX, store UINT64_MAX. 
                               We assume that primes <= fbb have already been 
                               removed, thus any factor <= assume_prime_thresh 
                               is assumed prime without further test. */
  facul_method_t *methods;  /* List of methods to try */
} facul_strategy_t;


int nb_curves (const unsigned int);
facul_strategy_t * facul_make_strategy (const unsigned long,
					const unsigned int);
void facul_clear_strategy (facul_strategy_t *);
void facul_print_stats (FILE *);
int facul (unsigned long *, const mpz_t, const facul_strategy_t *);

#endif /* FACUL_H */
