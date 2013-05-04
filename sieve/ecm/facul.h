#ifndef FACUL_H
#define FACUL_H

#include <stdio.h>
#include <stdint.h>
#include <gmp.h>

#define PM1_METHOD 1
#define PP1_METHOD 2
#define EC_METHOD 3

#define FACUL_NOT_SMOOTH (-1)

#define STATS_LEN 128

typedef struct {
  long method; /* Which method to use (P-1, P+1 or ECM) */
  void *plan;  /* Parameters for that method */
} facul_method_t;


/* All prime factors in the input number must be > fb. A factor of the 
   input number is assumed to be prime if it is < fb^2.
   The input number is taken to be not smooth if it has a 
   prime factor > lpb. */

typedef struct {
  unsigned long lpb;        /* Large prime bound as an unsigned long 
			       (the integer value, not the bit size!) */
  uint64_t assume_prime_thresh; /* The factor base bound squared. If the 
                               square exceeds UINT64_MAX, store UINT64_MAX. 
                               We assume that primes <= fbb have already been 
                               removed, thus any factor <= assume_prime_thresh 
                               is assumed prime without further test. */
  facul_method_t *methods;  /* List of methods to try */
} facul_strategy_t;


/* Here we take lpb in bits */
facul_strategy_t * facul_make_strategy (const int n, const unsigned long fbb, 
					const unsigned int lpb);
void facul_clear_strategy (facul_strategy_t *);
void facul_print_stats (FILE *);
int facul (unsigned long *, const mpz_t, const facul_strategy_t *);

#endif /* FACUL_H */
