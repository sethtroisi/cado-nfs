#include "cado.h"
#include "modredc_2ul2_default.h"

#define ecm ecm_2ul2
#define ell_pointorder ell_pointorder_2ul2
#define ellM_curveorder_jacobi ellM_curveorder_jacobi_2ul2
#define ell_curveorder ell_curveorder_2ul2
#define ellM_double ellM_double_2ul2
#define ellM_add ellM_add_2ul2
#define ellM_interpret_bytecode ellM_interpret_bytecode_2ul2
#define ecm_stage2 ecm_stage2_2ul2

#define ecmE ecmE_2ul2
#define ellE_double ellE_double_2ul2
#define ellE_add ellE_add_2ul2
#define ecmE_stage2 ecmE_stage2_2ul2

#define ec_arith ec_arith_2ul2
#define ec_point_init ec_point_init_2ul2
#define ec_point_clear ec_point_clear_2ul2
#define ec_point_set ec_point_set_2ul2
#define ec_point_swap ec_point_swap_2ul2
#define ec_point_print ec_point_print_2ul2

#define edwards_neg edwards_neg_2ul2
#define edwards_add edwards_add_2ul2
#define edwards_sub edwards_sub_2ul2
#define edwards_dbl edwards_dbl_2ul2
#define edwards_tpl edwards_tpl_2ul2
#define edwards_smul_ui edwards_smul_ui_2ul2

#define montgomery_dbl montgomery_dbl_2ul2
#define montgomery_dadd montgomery_dadd_2ul2
#define montgomery_smul_ui montgomery_smul_ui_2ul2

#define weierstrass_add weierstrass_add_2ul2
#define weierstrass_dbl weierstrass_dbl_2ul2
#define weierstrass_smul_ui weierstrass_smul_ui_2ul2

#include "ecm.c"
