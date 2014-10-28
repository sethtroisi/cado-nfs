#ifndef GENERATE_STRATEGIES
#define GENERATE_STRATEGIES

#include "tab_decomp.h"
#include "tab_strategy.h"
#include "tab_point.h"
#include "convex_hull.h"

/************************************************************************/
/*                 COLLECT DATA FOR SEVERAL (R1,R2)                     */
/*     WITH R1 ET R2 THE LENGHT OF COFACTORS (1: rat, 2: alg)           */
/************************************************************************/

double
compute_proba_strategy_r(tabular_decomp_t * init_tab, strategy_t * strat,
			 int len_p_min, int len_p_max, int side);

double
compute_time_strategy_r(tabular_decomp_t * init_tab, strategy_t * strat,
			int side, double proba_failed_facto_other_side);

void
collect_data_all_strategy_r1_r2(tabular_strategy_t * tab_strat,
				tabular_decomp_t * init_tab1,
				tabular_decomp_t * init_tab2, int lb, int ub);

/************************************************************************/
/*                   GENERATE MATRIX                                    */
/************************************************************************/

tabular_strategy_t *generate_strategies_oneside(fm_t * zero, tabular_fm_t * pm1,
						tabular_fm_t * pp1,
						tabular_fm_t * ecm_m16,
						tabular_fm_t * ecm_rc, int lb,
						int ub, int r);

tabular_strategy_t ***generate_matrix(tabular_fm_t * pm1, tabular_fm_t * pp1,
				      tabular_fm_t * ecm_m16,
				      tabular_fm_t * ecm_rc, int rlb, int rub,
				      int rmfb, int alb, int aub, int amfb);

/************************************************************************/
/*                      CONVEX_HULL_ST                                  */
/************************************************************************/

tabular_point_t *convert_tab_point_to_tab_strategy(tabular_strategy_t * t);

tabular_strategy_t *convert_tab_strategy_to_tab_point(tabular_point_t * t,
						      tabular_strategy_t *
						      init);

tabular_strategy_t *convex_hull_strategy(tabular_strategy_t * t);

#endif				/* GENERATE_STRATEGIES */
