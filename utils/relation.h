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

struct relation_stream_s {
    FILE * source;
    relation_t rel;     // This object is owned by the relation_stream
                        // object. Callers may manipulate the relation_t
                        // member, but not clear it.

    // parameters. may be set by the caller after calling
    // relation_stream_init.
    
    int parse_only_ab;  // used when only the (a,b) pair information
                        // is interesting, e.g. for duplicates removal.
                        
    // note that if parse_only_ab is unset, the relation returned by
    // relation_stream_get is sorted.
    
    // various stats stuff. May be used by the caller.
    size_t pos;
    index_t nrels;
    index_t lnum;
    // only valid after relation_stream_trigger_disp_progress
    // of relation_stream_disp_progress_now_p
    double dt, mb_s, rels_s;

    // temporaries, + various stuff for internal use.
    int pipe;   /* whether source was popen()ed */
    double t0, t1;
};

typedef struct relation_stream_s relation_stream[1];
typedef struct relation_stream_s * relation_stream_ptr;
typedef const struct relation_stream_s * relation_stream_srcptr;

#ifdef __cplusplus
extern "C" {
#endif

// Relation I/O
extern void relation_init(relation_t *rel);
extern void relation_clear(relation_t *rel);
// extern int read_relation(relation_t *rel, const char *str);
// extern int fread_relation(FILE *file, relation_t *rel);
extern unsigned long findroot(long a, unsigned long b, unsigned long p);
extern void computeroots(relation_t * rel);
extern char * u64toa16 (char *p, uint64_t m);
extern char * u64toa10 (char *p, uint64_t m);
extern char * d64toa10 (char *p, int64_t m);
extern char * d64toa16 (char *p, int64_t m);


extern void fprint_relation(FILE *file, relation_t * rel);
extern void fprint_relation_raw (FILE *file, relation_t * rel);
extern void reduce_exponents_mod2 (relation_t *rel);

/* FIXME: The following interface still strongly relies on the fact that
 * the rational side is [0] and the algebraic side is [1] */
extern void relation_add_prime (relation_t *rel, int side, unsigned long p);



extern void relation_copy (relation_t *s, relation_t * r);


/* reads over relations in a file, just discarding them */
extern void skip_relations_in_file(FILE * file, int n) __attribute__((deprecated));

/*
extern int read_relation_quick_and_dirty(int64_t * pa, uint64_t * pb, char * line, FILE * f);
extern int read_relation_from_stream(relation_t * rel, char * line, FILE * f);
*/

/* Not clear whether we're willing to insist on _not_ exposing these
 * three functions or not */
extern void relation_provision_for_primes(relation_t * rel, int nr, int na);
extern void relation_compress_rat_primes(relation_t * rel);
extern void relation_compress_alg_primes(relation_t * rel);

extern void relation_stream_init(relation_stream_ptr rs);
extern void relation_stream_closefile(relation_stream_ptr rs);
extern void relation_stream_clear(relation_stream_ptr rs);
extern void relation_stream_openfile(relation_stream_ptr rs, const char * name);
extern void relation_stream_bind(relation_stream_ptr rs, FILE * f);
extern void relation_stream_unbind(relation_stream_ptr rs);
extern int relation_stream_disp_progress_now_p(relation_stream_ptr rs);
extern int relation_stream_get(relation_stream_ptr rs, char * line, size_t lsize, int forced_read, unsigned int ab_base, int allow_comment);
extern int relation_stream_get_skip(relation_stream_ptr rs);
extern void relation_stream_trigger_disp_progress(relation_stream_ptr rs);

#ifdef __cplusplus
}
#endif

#endif	/* RELATION_H_ */
