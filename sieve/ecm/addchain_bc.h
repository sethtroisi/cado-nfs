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

double addchain_bytecode (const unsigned int, addchain_cost_srcptr,
                                                                  bc_state_t *);
/* Op cost for Twisted Edwards Curves with a=-1
 * TODO explain different models and why dbladd > add+dbl*/
static inline void
opcost_TwEdwards_minus1_opcost (addchain_cost_ptr opcost)
{
  opcost->dbl = 7.;
  opcost->add = 7.;
  opcost->dbladd = 15.;
  opcost->dbl_precomp = 8.;
  opcost->add_precomp = 8.;
}
#endif
