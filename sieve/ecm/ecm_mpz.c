#include "cado.h"
#include "mod_mpz_default.h"

#define ecm ecm_mpz
#define ell_pointorder ell_pointorder_mpz
#define ellM_curveorder_jacobi ellM_curveorder_jacobi_mpz
#define ell_curveorder ell_curveorder_mpz
#define ellM_double ellM_double_mpz
#define ellM_add ellM_add_mpz
#define ellM_interpret_bytecode ellM_interpret_bytecode_mpz
#define ecm_stage2 ecm_stage2_mpz

#define ecmE ecmE_mpz
#define ellE_double ellE_double_mpz
#define ellE_add ellE_add_mpz
#define ecmE_stage2 ecmE_stage2_mpz

#define ec_arith ec_arith_mpz
#define ec_point_init ec_point_init_mpz
#define ec_point_clear ec_point_clear_mpz
#define ec_point_set ec_point_set_mpz
#define ec_point_swap ec_point_swap_mpz
#define ec_point_print ec_point_print_mpz

#define edwards_neg edwards_neg_mpz
#define edwards_add edwards_add_mpz
#define edwards_sub edwards_sub_mpz
#define edwards_dbl edwards_dbl_mpz
#define edwards_tpl edwards_tpl_mpz
#define edwards_smul_ui edwards_smul_ui_mpz

#define montgomery_dbl montgomery_dbl_mpz
#define montgomery_dadd montgomery_dadd_mpz
#define montgomery_smul_ui montgomery_smul_ui_mpz

#define weierstrass_add weierstrass_add_mpz
#define weierstrass_dbl weierstrass_dbl_mpz
#define weierstrass_smul_ui weierstrass_smul_ui_mpz

#include "ecm.c"
