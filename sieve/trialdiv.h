#include <gmp.h>
#include "fb.h"

/* The maximum number of words in the numbers to be trial-divided.
   The $l$ value from the thesis text is equal to TRIALDIV_MAXLEN - 1 */
#define TRIALDIV_MAXLEN 6

typedef struct {
  unsigned long p;
  unsigned long w[TRIALDIV_MAXLEN - 1]; /* w[i] = w^{i+1} mod p */
  unsigned long pinv; /* pinv == 1/p (mod w) */
  unsigned long plim; /* plim = (w-1)/p */
} trialdiv_divisor_t;

/* Divides primes in d out of N and stores them (with multiplicity) in f */
size_t trialdiv (unsigned long *, mpz_t, const trialdiv_divisor_t *,
		 const size_t);

/* Initialise a trialdiv_divisor_t array with the primes stored in an 
   fbprime_t array. Allocates memory. */
trialdiv_divisor_t *trialdiv_init (const fbprime_t *, const unsigned int);

/* Frees the memory allocated by trialdiv_init() */
void trialdiv_clear (trialdiv_divisor_t *d);
