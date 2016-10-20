#ifndef _ECM_PRAC_BC_H_
#define _ECM_PRAC_BC_H_

#include <gmp.h>
#include "bytecoder.h"

struct prac_cost_s
{
  double dbl; /* cost of a doubling */
  double dadd; /* cost of a differential addition */
};

typedef struct prac_cost_s prac_cost_t;

unsigned int prac_bytecode (char **, unsigned int, unsigned int,
                            unsigned int, const prac_cost_t *, int, int);
void prac_bytecode_fprintf (FILE *, const char *, unsigned int);
int prac_bytecode_check (const char *, unsigned int, mpz_srcptr, int);
void prac_cache_free ();
#endif
