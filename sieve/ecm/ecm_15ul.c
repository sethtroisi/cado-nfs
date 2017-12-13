#include "cado.h"
#include "modredc_15ul_default.h"

#define ecm ecm_15ul
#define ell_pointorder ell_pointorder_15ul
#define ellM_curveorder_jacobi ellM_curveorder_jacobi_15ul
#define ell_curveorder ell_curveorder_15ul


#define ecm_stage2 ecm_stage2_15ul

#define ec_arith ec_arith_15ul
#define ec_point_init ec_point_init_15ul
#define ec_point_clear ec_point_clear_15ul
#define ec_point_set ec_point_set_15ul
#define ec_point_swap ec_point_swap_15ul
#define ec_point_print ec_point_print_15ul

#define edwards_neg edwards_neg_15ul
#define edwards_add edwards_add_15ul
#define edwards_sub edwards_sub_15ul
#define edwards_dbl edwards_dbl_15ul
#define edwards_tpl edwards_tpl_15ul
#define edwards_smul_ui edwards_smul_ui_15ul

#define montgomery_dbl montgomery_dbl_15ul
#define montgomery_dadd montgomery_dadd_15ul
#define montgomery_smul_ui montgomery_smul_ui_15ul
#define bytecode_prac_interpret_ec_montgomery bytecode_prac_interpret_ec_montgomery_15ul

#define weierstrass_add weierstrass_add_15ul
#define weierstrass_dbl weierstrass_dbl_15ul
#define weierstrass_smul_ui weierstrass_smul_ui_15ul


#include "ecm.c"
