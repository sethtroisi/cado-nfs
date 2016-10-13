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

#define ADDCHAIN_Q_MAX 253 /* Must always fit in a char and be < 255 */
#define ADDCHAIN_DBL ((literal_t) 0x7f)

double addchain_bytecode (const unsigned int, const addchain_cost_t *,
                          bc_state_t *, int);
#endif
