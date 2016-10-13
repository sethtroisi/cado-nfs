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

#define ADDCHAIN_Q_MAX 253 /* Must be < 255 */
#define ADDCHAIN_DBL ((literal_t) 0x7f)
#define ADDCHAIN_DBL_STR "\x7f"
#define ADDCHAIN_2DBL ((literal_t) 0xff)

unsigned int addchain_bytecode (char **, const unsigned int, const unsigned int,
                                const addchain_cost_t *, bc_dict_t *, int);
void addchain_bytecode_fprintf (FILE *, const char *, unsigned int);
int addchain_bytecode_check (const char *, unsigned int, mpz_srcptr, int);
#endif
