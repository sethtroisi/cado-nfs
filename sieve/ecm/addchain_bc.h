#ifndef _ECM_ADDCHAIN_BC_H_
#define _ECM_ADDCAHIN_BC_H_

struct addchain_cost_s
{
  double dbl;
  double add;
  double dbladd;
  double dbl_precomp;
  double add_precomp;
};

typedef struct addchain_cost_s addchain_cost_t[1];
typedef struct addchain_cost_s * addchain_cost_ptr;
typedef const struct addchain_cost_s * addchain_cost_srcptr;

#define ADDCHAIN_Q_MAX 253 /* Must always fit in a char */
#define ADDCHAIN_DBL ((literal_t) 0x7f)

/* Op cost for Twisted Edwards Curves with a=-1
 * TODO explain different models and why dbladd > add+dbl*/
addchain_cost_t TwEdwards_minus1_opcost = [ {.dbl = 7.} ];
#endif
