#ifndef RELATION_H_
#define RELATION_H_

// #define STR_LEN_MAX 2048 /* maximal number of characters of a relation line */
static const int STR_LEN_MAX __attribute__((deprecated)) = 2048;
static const int RELATION_MAX_BYTES = 4096;

#include <stdint.h>
#include <stdio.h>

#include "typedefs.h"

typedef struct {
  unsigned long p;      /* rational prime */
  int e;                /* exponent (may want negative exponent in sqrt) */
} rat_prime_t;

typedef struct {
  unsigned long p;      /* algebraic prime */
  unsigned long r;      /* corresponding root: r = a/b mod p */
  int e;                /* exponent (may want negative exponent in sqrt) */
} alg_prime_t;

typedef struct {
  int64_t a;		/* only a is allowed to be negative */
  uint64_t b;
  rat_prime_t *rp;	/* array of rational primes */
  alg_prime_t *ap;	/* array of algebraic primes */
  uint8_t nb_rp;	/* number of rational primes */
  uint8_t nb_ap;        /* number of algebraic primes */
  uint8_t nb_rp_alloc;	/* allocated space for rp */
  uint8_t nb_ap_alloc;	/* allocated space for ap */
} relation_t;

#ifdef __cplusplus
extern "C" {
#endif

// Relation I/O
extern void relation_init(relation_t *rel);
extern void relation_clear(relation_t *rel);
p_r_values_t relation_compute_r (int64_t a, uint64_t b, p_r_values_t p);
extern void relation_compute_all_r (relation_t * rel);
extern char * u64toa16 (char *p, uint64_t m);
extern char * u64toa10 (char *p, uint64_t m);
extern char * d64toa10 (char *p, int64_t m);
extern char * d64toa16 (char *p, int64_t m);


extern void fprint_relation(FILE *file, relation_t * rel, const char *);

/* FIXME: The following interface still strongly relies on the fact that
 * the rational side is [0] and the algebraic side is [1] */
extern void relation_add_prime (relation_t *rel, int side, unsigned long p);



extern void relation_copy (relation_t *s, relation_t * r);


/* Not clear whether we're willing to insist on _not_ exposing these
 * three functions or not */
extern void relation_provision_for_primes(relation_t * rel, int nr, int na);
extern void relation_compress_rat_primes(relation_t * rel);
extern void relation_compress_alg_primes(relation_t * rel);

#ifdef __cplusplus
}
#endif

#endif	/* RELATION_H_ */
