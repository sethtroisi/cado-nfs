#ifndef _ECM_ADDCHAIN_BC_H_
#define _ECM_ADDCHAIN_BC_H_

#include "bytecoder.h"

struct addchain_cost_s
{
  double dbl;
  double add;
  double dbladd;
  double dbl_precomp;
  double add_precomp;
};

typedef struct addchain_cost_s addchain_cost_t;

/* ADDCHAIN_Q_MAX must be < 255.
 * Set it to 101, because for B1 < 1024, it is not necessary to try bigger q,
 * they never produce better additions chains.
 */
#define ADDCHAIN_Q_MAX 101
/* The code for ADDCHAIN_[n]DBL must be different from 0 because if we use a
 * dict to combine multiple ADDCHAIN_DBL into ADDCHAIN_[n]DBL, 0 has a special
 * meaning (which is remove from the bytecode).
 */
#define ADDCHAIN_DBL ((literal_t) 0x7f)
#define ADDCHAIN_DBL_STR "\x7f"
#define ADDCHAIN_2DBL ((literal_t) 0xff)
/* max number of consecutive DBL that we can encode */
#define ADDCHAIN_MAX_CONSECUTIVE_DBL 2

unsigned int addchain_bytecode (char **, unsigned int, unsigned int,
                                unsigned int, const addchain_cost_t*, int, int);
void addchain_bytecode_fprintf (FILE *, const char *, unsigned int);
int addchain_bytecode_check (const char *, unsigned int, mpz_srcptr, int);
#endif
