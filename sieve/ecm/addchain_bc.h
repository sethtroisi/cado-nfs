#ifndef _ECM_ADDCHAIN_BC_H_
#define _ECM_ADDCHAIN_BC_H_

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

/* Costs of operations for Twisted Edwards Curves with a=-1
 *  For those curves, we use 3 different models:
 *    projective, extended and (only internally) completed
 *  The costs corresponds to operations on different models:
 *    - dbl corresponds to a doubling projective -> projective
 *    - add corresponds to an addition extended,extended -> projective
 *    - dbladd corresponds to a doubling and an addition
 *        projective, extended -> projective
 *      It is more costly than add + dbl but it is due to the fact that if we
 *      wanted to do the same operation with 1 add and 1 dbl, we would need to
 *      convert the output of the dbl from projective to extended in order to be
 *      able to use the add, which is costly.
 *    - dbl_precomp corresponds to a doubling extended -> extended
 *    - add_precomp corresponds to an addition extended, extended, -> extended
 *  We count 1 for a multiplication and 1 for a squaring on the base field.
 */
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
