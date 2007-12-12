#ifndef	CADO_UTILS_H_
#define	CADO_UTILS_H_

#include "cado.h"
#include <stdint.h>

#ifdef	__cplusplus
extern "C" {
#endif

// This reads a file created by polyselect and fill in the structure
// accordingly. Return 1 if success, 0 if failure (and diagnostic on
// stderr)
extern int read_polynomial (cado_poly, char *filename);
extern void fprint_polynomial (FILE *, mpz_t *, const int);

// Relation I/O
extern void copy_rel(relation_t *Rel, relation_t rel);
extern int read_relation(relation_t *rel, const char *str);
extern int fread_relation(FILE *file, relation_t *rel);
extern unsigned long findroot(long a, unsigned long b, unsigned long p);
extern void computeroots(relation_t * rel);
extern void fprint_relation(FILE *file, relation_t rel);


extern uint64_t microseconds();

/* cputime */
static inline int cputime(void) { return (int) microseconds() / 1000; }
static inline double seconds(void) { return (double) microseconds() /1.0e6; }

/* long_poly arithmetic */
typedef struct {
  int alloc;    /* number of allocated coefficients */
  int degree;   /* degree < alloc */
  long *coeff; /* coefficient list */
} __long_poly_struct;
typedef __long_poly_struct long_poly_t[1];
extern int roots_mod_long (long*, mpz_t*, int, const long);
extern int nbits (unsigned long);

/* getprime */
extern unsigned long getprime (unsigned long);

#ifdef	__cplusplus
}	// extern "C"
#endif

#endif	/* CADO_UTILS_H_ */
