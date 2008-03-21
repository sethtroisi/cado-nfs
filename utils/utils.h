#ifndef	CADO_UTILS_H_
#define	CADO_UTILS_H_

#include "cado.h"
#include "mod_ul.h"
#include <stdint.h>

#define LONG int64_t
#define ULONG uint64_t

#include <limits.h>

/* It's awful, but I confess that this ULONG_BITS is not ``portable'',
 * norm-wise.  GMP_LIMB_BITS is at hand, but could differ. An #ifdef
 * switch depending on macros like __x86_64 is considerably more fragile.
 */
#define	ULONG_BITS	((int) (sizeof(unsigned long) * CHAR_BIT))

/* the following function is missing in GMP */
#ifndef mpz_addmul_si
#define mpz_addmul_si(a, b, c)                  \
  do {                                          \
    if (c >= 0)                                 \
      mpz_addmul_ui (a, b, c);                  \
    else                                        \
      mpz_submul_ui (a, b, -c);                 \
  }                                             \
  while (0)
#endif

#ifdef	__cplusplus
extern "C" {
#endif

/* file polyfile.c */
// This reads a file created by polyselect and fill in the structure
// accordingly. Return 1 if success, 0 if failure (and diagnostic on
// stderr)
extern int read_polynomial (cado_poly, char *filename);
extern void fprint_polynomial (FILE *, mpz_t *, const int);
void clear_polynomial (cado_poly);
void cado_poly_init (cado_poly);
void cado_poly_clear (cado_poly);
void check_polynomials (cado_poly);

  /* fpoly.c */
double fpoly_eval (const double *, const int, const double);
double fpoly_dichotomy (double *, int, double, double, double, unsigned int);
void   fpoly_print (FILE *, const double *f, const int deg, char *name);

// Relation I/O
extern void copy_rel(relation_t *Rel, relation_t rel);
extern void clear_relation(relation_t *rel);
extern int read_relation(relation_t *rel, const char *str);
extern int fread_relation(FILE *file, relation_t *rel);
extern unsigned long findroot(long a, unsigned long b, unsigned long p);
extern void computeroots(relation_t * rel);
extern void fprint_relation(FILE *file, relation_t rel);
extern void reduce_exponents_mod2 (relation_t *rel);


extern uint64_t microseconds();

/* cputime */
static inline int cputime(void) { return (int) microseconds() / 1000; }
static inline double seconds(void) { return (double) microseconds() /1.0e6; }

/* long_poly: long_poly arithmetic */
typedef struct {
  int alloc;    /* number of allocated coefficients */
  int degree;   /* degree < alloc */
  LONG *coeff; /* coefficient list */
} __long_poly_struct;
typedef __long_poly_struct long_poly_t[1];
extern int roots_mod_long (LONG*, mpz_t*, int, const LONG);
extern int nbits (ULONG);
void long_poly_init (long_poly_t f, int d);
void long_poly_clear (long_poly_t f);
int long_poly_set_mod (long_poly_t fp, mpz_t *f, int d, LONG p);
int isirreducible_mod_long(long_poly_t fp, const LONG p);
int long_poly_fits (unsigned int, LONG);

/* modul_poly */
int modul_roots_mod_long  (unsigned long*, mpz_t*, int, modulus_t);
int modul_roots_mod_int64 (int64_t*,       mpz_t*, int, modulus_t);

/* getprime */
extern unsigned long getprime (unsigned long);

/* gmp_aux */
void mpz_set_uint64 (mpz_t, uint64_t);
uint64_t mpz_get_uint64 (mpz_t);
uint64_t uint64_nextprime (uint64_t);

#ifdef	__cplusplus
}	// extern "C"
#endif

#endif	/* CADO_UTILS_H_ */
