#include "cado.h"
#include "modredc_ul_default.h"

#define ecm ecm_ul

#define ell_pointorder ell_pointorder_ul
#define ellM_curveorder_jacobi ellM_curveorder_jacobi_ul
#define ell_curveorder ell_curveorder_ul

/* The next ones are static so there's no need to rename them, but it's 
   nice to have these functions distinguishable e.g. in profiler output */

#define ecm_stage2 ecm_stage2_ul

#define ec_arith ec_arith_ul
#define ec_point_init ec_point_init_ul
#define ec_point_clear ec_point_clear_ul
#define ec_point_set ec_point_set_ul
#define ec_point_swap ec_point_swap_ul
#define ec_point_print ec_point_print_ul

#define edwards_neg edwards_neg_ul
#define edwards_add edwards_add_ul
#define edwards_sub edwards_sub_ul
#define edwards_dbl edwards_dbl_ul
#define edwards_tpl edwards_tpl_ul
#define edwards_smul_ui edwards_smul_ui_ul

#define montgomery_dbl montgomery_dbl_ul
#define montgomery_dadd montgomery_dadd_ul
#define montgomery_smul_ui montgomery_smul_ui_ul
#define bytecode_prac_interpret_ec_montgomery bytecode_prac_interpret_ec_montgomery_ul

#define weierstrass_add weierstrass_add_ul
#define weierstrass_dbl weierstrass_dbl_ul
#define weierstrass_smul_ui weierstrass_smul_ui_ul


#include "ecm.c"

