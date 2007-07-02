#include "cado.h"

// This reads a file created by polyselect and fill in the structure
// accordingly. Return 1 if success, 0 if failure (and diagnostic on
// stderr)
int read_polynomial (cado_poly, char *filename);

// Relation I/O
void copy_rel(relation_t *Rel, relation_t rel);
int read_relation(relation_t *rel, const char *str);
int fread_relation(FILE *file, relation_t *rel);
unsigned long findroot(long a, unsigned long b, unsigned long p);
void computeroots(relation_t * rel);
void fprint_relation(FILE *file, relation_t rel);

/* cputime */
int cputime (void);

/* long_poly arithmetic */
typedef struct {
  int alloc;    /* number of allocated coefficients */
  int degree;   /* degree < alloc */
  long *coeff; /* coefficient list */
} __long_poly_struct;
typedef __long_poly_struct long_poly_t[1];

