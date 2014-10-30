#ifndef GENERATE_STRATEGIES
#define GENERATE_STRATEGIES

#include "tab_decomp.h"
#include "tab_strategy.h"
#include "tab_point.h"
#include "convex_hull.h"

/************************************************************************/
/*                      COLLECT DATA FOR ONLY ONE COFACTOR              */
/************************************************************************/

double compute_proba_strategy(tabular_decomp_t * init_tab, strategy_t * strat,
			      int len_p_min, int len_p_max);

double compute_time_strategy(tabular_decomp_t * init_tab, strategy_t * strat);

/************************************************************************/
/*                   GENERATE MATRIX                                    */
/************************************************************************/

tabular_strategy_t *generate_strategies_oneside(tabular_decomp_t * tab_decomp,
						fm_t * zero,
						tabular_fm_t * pm1,
						tabular_fm_t * pp1,
						tabular_fm_t * ecm_m16,
						tabular_fm_t * ecm_rc, int lb,
						int ub, int r);

tabular_strategy_t *generate_strategy_r0_r1(tabular_strategy_t * strat_r0,
					    tabular_strategy_t * strat_r1);

tabular_strategy_t ***generate_matrix(const char *name_directory_decomp,
				      tabular_fm_t * pm1, tabular_fm_t * pp1,
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
